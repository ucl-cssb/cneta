/*
This program aims to simulate SVs along a phylogenetic tree.

Tree representation:
if we have n samples these trees always have the germline variation as the n+1 samples
the final edge length is always zero
edges are always directed from -> to going from the node/MRCA
therefore leaf nodes are always in the second column of edges and the root is only in the first
Internally nodes and edges have ids starting from 0
to match up with ape in R we always print id+1
*/

#include "model.hpp"
#include "parse_cn.hpp"
#include "evo_tree.hpp"
#include "tree_op.hpp"
#include "genome.hpp"


#include <boost/program_options.hpp>

// using namespace std;


// The number of bins for each chromosome. Each bin corresponds to a window of size 500,000 bp. 4401 bins in total
const int BIN_SIZE[] = {367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};
vector<int> CHR_BIN_SIZE(BIN_SIZE, BIN_SIZE + sizeof(BIN_SIZE) / sizeof(BIN_SIZE[0]));
const int NODE_MIN_TIME = 1000;

enum SIM_METHOD {SIM_TIME, SIM_SEQ};
enum MUT_TYPE {DUP, DEL, GAIN, LOSS, WGD};

// key: chr, seg, copy_number
typedef map<int, map<int,int>> copy_number;

int debug = 0;
gsl_rng* r;
unsigned seed;

// mean of exp size distributions
double mean_dup_size;
double mean_del_size;


// unary function and pointer to unary function
// allows use of gsl rng for standard template algorithms
inline long unsigned myrng(long unsigned n){
  return gsl_rng_uniform_int(r, n);
}

long unsigned (*fp_myrng)(long unsigned);


// rates: The rate of all events (segment duplication, segment deletion) for finite sites model
double get_site_rates(const vector<double>& rate_consts, int model, int stateA, int stateB, int cn_max, double& site_dup_rate, double& site_del_rate){
    double rate = 0.0;

    if(model == MK){ // Mk
        // assume duplication rate equals deletion rate
        assert(rate_consts[0] == rate_consts[1]);
        rate = 4 * rate_consts[0] / 5;
    }else if(model == BOUNDT){  // bounded total
        double dup_rate = rate_consts[0];
        double del_rate = rate_consts[1];
        int state = stateA + stateB;

        if(state == 0){
            site_dup_rate = 0.0;
            site_del_rate = 0.0;
        }else if(state == cn_max){
            site_dup_rate = 0.0;
            site_del_rate = 2 * cn_max * del_rate;
        }else{
            site_dup_rate = 2 * state * dup_rate;
            site_del_rate = 2 * state * del_rate;
        }
        rate = site_dup_rate + site_del_rate;
    }else if(model == BOUNDA){   // bounded allele specific
        double dup_rate = rate_consts[0];
        double del_rate = rate_consts[1];

        if(stateA == 0 && stateB == 0){
            site_dup_rate = 0.0;
            site_del_rate = 0.0;
        }else if(stateA + stateB == cn_max){   // only deletion can occur
            site_dup_rate = 0.0;
            if(stateA == 0 || stateB == 0){   // deletion can only occur at one allele
                site_del_rate = del_rate;
            }else{
                site_del_rate = 2 * del_rate;
            }
        }else{
            if(stateA == 0 || stateB == 0){
                site_dup_rate = dup_rate;
                site_del_rate = del_rate;
            }else{
                site_dup_rate = 2 * dup_rate;
                site_del_rate = 2 * del_rate;
                // cout << "assign rate " << site_dup_rate << "\t" << site_del_rate << endl;
            }
        }
        rate = site_dup_rate + site_del_rate;
    }else{

    }

    // cout << stateA << "\t" << stateB << "\t" << site_dup_rate << "\t" << site_del_rate << endl;

    return rate;
}

// The WGD rate of the genome changes as the copy numbers change, due to the upper limit on total copy number
// WGD is impossible when it is approaching cn_max
void update_wgd_rate(genome& g, double& wgd_rate, int cn_max){
    int max_cn = get_max_cn_genome(g);
    if(max_cn <= 0 || 2 * max_cn > cn_max){
        // cout << "Allowed MAX CN is " << cn_max << endl;
        // cout << "Maximum copy number on genome " << g.node_id << " is " << max_cn << endl;
        wgd_rate = 0.0;
    }
}

// The chromosome gain rate of the genome changes as the copy numbers change, due to the upper limit on total copy number
void update_chr_gain_rate(genome& g, int i, double& chr_gain_rate, int cn_max){
    double max_cn_c = get_max_cn_chr(g, i);
    double max_snum = get_max_seg_num_chr(g, i);
    if(max_cn_c <= 0 || max_snum >= cn_max){
        chr_gain_rate = 0.0;
    }
}


void update_chr_loss_rate(genome& g, int i, double& chr_loss_rate){
    double max_cn_c = get_max_cn_chr(g, i);
    if(max_cn_c <= 0){
        chr_loss_rate = 0.0;
    }
}


// rates: The total mutation rate at each site in the genome, used for randomly selecting a site
// The site is pre-specified at the beginning, denoted by each segment ID
// site_*_rates: duplication/deletion probabilities for all sites on the genome, used for picking a site for duplication/deletion
// chr_*_rates: chromosome gain/loss probabilities for all sites on the genome, used for picking a site for chromosome gain/loss
// type_rates: used for selecting different mutational types
double get_total_rates_allele_specific(genome& g, vector<double>& site_dup_rates, vector<double>& site_del_rates, vector<double>& chr_gain_rates, vector<double>& chr_loss_rates, vector<double>& type_rates, const vector<double>& rate_consts, int model, int cn_max){
    assert(rate_consts.size() == 5);
    double chr_gain_rate = rate_consts[2];
    double chr_loss_rate = rate_consts[3];
    double wgd_rate = rate_consts[4];
    double dup_rate_all = 0.0;
    double del_rate_all = 0.0;
    double gain_rate_all = 0.0;
    double loss_rate_all = 0.0;

    g.calculate_cn(); // Find the current state of the genome before computing rates (affected by copy number)
    g.calculate_allele_cn();

    vector<double> rates_chrs;  // mutation rates on each chromosome
    int s = 0;
    // Tranverse all the sites in the genome to get their rates
    for(int i = 0; i < g.cn_profile.size(); ++i){   // For each chromosome
        vector<double> rates_chr_sites; // mutation rates on all sites of one chromosome
        for(int j = 0; j < g.cn_profile[i].size(); ++j){    // For each segment in the chromosome
            // Get the copy number of this segment
            int stateA = g.allele_cn_profile[i][j];
            int stateB = g.allele_cn_profile[i + NUM_CHR][j];
            double site_dup_rate = 0.0;
            double site_del_rate = 0.0;
            double rate = get_site_rates(rate_consts, model, stateA, stateB, cn_max, site_dup_rate, site_del_rate);
            // if(debug){
            //     cout << "Rate on chr " << i+1 << " seg " << j+1 << " site " << s++ << " stateA " << stateA << " stateB " << stateB << " is " << rate << endl;
            //     cout << "Site duplication rate " << site_dup_rate << "\tSite deletion rate " << site_del_rate << endl;
            // }

            dup_rate_all += site_dup_rate;
            del_rate_all += site_del_rate;

            site_dup_rates.push_back(site_dup_rate);     // mutation rates on all sites of all chromosomes
            site_del_rates.push_back(site_del_rate);     // mutation rates on all sites of all chromosomes

            rates_chr_sites.push_back(rate);
        }
        // Get the rate on one chromosome
        double rate_c = accumulate(rates_chr_sites.begin(), rates_chr_sites.end(), 0.0);    // rates for duplication/deletion

        update_chr_gain_rate(g, i, chr_gain_rate, cn_max);
        update_chr_loss_rate(g, i, chr_loss_rate);
        rate_c += chr_loss_rate;
        rate_c += chr_gain_rate;

        if(debug){
            cout << "Rate on chr " << i+1 << " is " << rate_c << endl;
            cout << "Chromosome gain rate " << chr_gain_rate << "\tChromosome loss rate " << chr_loss_rate << endl;
        }
        rates_chrs.push_back(rate_c);

        chr_gain_rates.push_back(chr_gain_rate);
        chr_loss_rates.push_back(chr_loss_rate);

        gain_rate_all += chr_gain_rate;
        loss_rate_all += chr_loss_rate;
    }
    // Get the total rate of mutation for the current genome (the sum of rates accross sites, chromosomes)
    double rate = accumulate(rates_chrs.begin(), rates_chrs.end(), 0.0);
    // cout << "WGD rate before " << wgd_rate << endl;
    update_wgd_rate(g,  wgd_rate, cn_max);
    // cout << "WGD rate after " << wgd_rate << endl;
    rate += wgd_rate;

    type_rates.push_back(dup_rate_all);
    type_rates.push_back(del_rate_all);
    type_rates.push_back(gain_rate_all);
    type_rates.push_back(loss_rate_all);
    type_rates.push_back(wgd_rate);

    if(debug){
        cout << "The total rate at all sites in the genome: " << rate << endl;
        // cout << site_dup_rates.size() << endl;
        // cout << site_del_rates.size() << endl;
        // cout << chr_gain_rates.size() << endl;
        // cout << chr_loss_rates.size() << endl;
        // cout << type_rates.size() << endl;

        cout << "Site rates:";
        for(int i = 0; i <  site_dup_rates.size(); i++){
            cout << "\t" << site_dup_rates[i] << "\t" << site_del_rates[i];
        }
        cout << endl;

        cout << "Chromosome rates:";
        for(int i = 0; i <  chr_gain_rates.size(); i++){
            cout << "\t" << chr_gain_rates[i] << "\t" << chr_loss_rates[i];
        }
        cout << endl;

        cout << "Rates of different mutational types on genome " << g.node_id;
        for(int i = 0; i < type_rates.size(); i++){
            cout << "\t" << type_rates[i];
        }
        cout << endl;
    }

    return rate;
}


// Note: the sites are counted for all the haplotypes. The positions on different haplotypes are counted as one site.
// The starting point of sites is 0
void site2chr(int site, int& chr, int& seg, const vector<int>& chr_lengths){
    if(debug) cout << "Site to convert is " << site << endl;
    if(site >= 0 && site < chr_lengths[0]){
        chr = 0;
        seg = site;
    }else{
        if(debug) cout << "There are " << chr_lengths.size() << " chromosomes" << endl;
        for(int i = 0; i < chr_lengths.size(); i++){
            int sum1 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 1, 0);
            int sum2 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 2, 0);
            // if(debug){
            //     cout << "Size for chr " << i+1 << " is " << chr_lengths[i] << endl;
            //     cout << "sum until chromosome " << i + 1 << " is " << sum1 << endl;
            //     cout << "sum until chromosome " << i + 2 << " is " << sum2 << endl;
            // }
            if(site >= sum1 && site < sum2){
                chr = i + 1;
                seg = site - sum1;
                break;
            }
        }
    }
    assert(seg < chr_lengths[chr]);
}




// input: c -- chromosome number (for each haplotype); seg_id -- location on a chromosome (should be intepreted as segment ID)
// if return value > 0,  mutations are generated successfully
int generate_duplication(genome& g, int c, int seg_id, int model, int cn_max, int debug){
    if(g.chrs[c].size() <= 0){
        cout << "All segments on chr " << c+1 << " has lost" << endl;
        return 0;
    }
    if(debug){
      cout << "GENOME BEFORE duplication " << endl;
      // g.print();
      g.print_cn();
    }
    // mean variant size in bins
    double mdup = mean_dup_size;
    int len = 0;
    if(mdup > 1) len = gsl_ran_exponential(r, mdup);

    if(debug){
        cout << "before duplication, Chr " << c+1 << " has " << g.chrs[c].size() << " segments" << endl;
        // for(int i = 0; i < g.chrs[c].size(); i++){
        //     cout << "\t" << g.chrs[c][i].seg_id;
        // }
        // cout << endl;
        cout << "dup duplication:" << len << endl;

        //cout << "SV: insert: " << loc + len + 1 << "\t" << loc << "\t" << loc + len + 1 << endl;
        cout << "\tSV: duplicating segment, chr, seg_id, len: " << c+1 << "\t" << seg_id << "\t" << len + 1 << endl;
        cout << "Previous state: ";
        // for(int i = 0; i < len + 1; i++){
            cout << "\t " << g.cn_profile[c%NUM_CHR][seg_id];
        // }
        cout << endl;
    }
    // Find the number of copies
    // int state = g.cn_profile[c%NUM_CHR][g.chrs[c][loc].seg_id];
    int state = g.cn_profile[c%NUM_CHR][seg_id];
    if(state >= cn_max) {
        // cout << "cn_max is " << cn_max << endl;
        return 0;
    }
    int ins_start = 0;
    double u = runiform(r, 0, 1);
    if(u < 0.5){
        // randomly select a position
        if(debug) cout << "nontandem duplication" << endl;
        ins_start = gsl_rng_uniform_int(r, g.chrs[c].size());
    }else{ // tandem duplication
        // Find the location of seg_id in the genome
        for(int i = 0; i < g.chrs[c].size(); i++){
            if(g.chrs[c][i].seg_id == seg_id){
                ins_start = i;
                break;
            }
        }
    }
    int ins_end = ins_start + len + 1; //assume tandem duplication
    if(debug) cout << "   start " << ins_start << ", end " << ins_end << endl;

    if(model == MK){
        vector<int> possible_states;
        for (int i = state + 1; i <= cn_max; i++){
            possible_states.push_back(i);
        }
        // Randomly pick a possible state
        int ssize = possible_states.size();
        int pstate;
        if(debug){
            cout << "There are " << ssize << " possible state(s) for duplication" << endl;
        }
        if(ssize > 1){
            int k = gsl_rng_uniform_int(r, ssize);
            pstate = possible_states[k];
        }
        else{
            pstate = state + 1;
        }
        // Ensure there are "state" copy of the region
        // if(debug){
        //     cout << "copy start " << g.chrs[c][loc].seg_id << endl;
        //     cout << "copy end " << g.chrs[c][loc + len].seg_id << endl;
        // }
        if(pstate - state <= 0){
            cout << "no possible states" << '\n';
            return 0;
        }
        // assume tandem duplication
        for(int j = 0; j < pstate-state; j++){
            // g.chrs[c].insert(g.chrs[c].begin() + loc+(len + 1)*(j+1), g.chrs[c].begin() + loc, g.chrs[c].begin() + loc + len + 1);
            g.chrs[c].insert(g.chrs[c].begin()+ins_start+(len + 1)*j, len + 1, segment(c%NUM_CHR, seg_id));

            if(debug){
                cout << j+1 << "th dup" << endl;
                int it;
                it = g.chrs[c][ins_start+(len + 1)*j].seg_id;
                cout << "dup start " << it << endl;
                it = g.chrs[c][ins_start+(len + 1)*j+len].seg_id;
                cout << "dup end " << it << endl;
            }
        }
    }else if(model == BOUNDT || model == BOUNDA){
        // int max_state = get_max_cn_chr_seg(g, c, seg_id, len);
        // int max_snum = get_max_seg_num_seg(g, c, seg_id, len);
        // // cout << "cn_max is " << cn_max << endl;
        // if(max_snum > cn_max){
        //     return 0;
        // }
        // g.chrs[c].insert(g.chrs[c].begin()+ins_end, g.chrs[c].begin() + loc, g.chrs[c].begin() + loc + len + 1);
        g.chrs[c].insert(g.chrs[c].begin()+ins_start, len + 1, segment(c%NUM_CHR, seg_id));

        if(debug){
            cout << "End position " << ins_start+len + 1 << endl;
            cout << "Insertion position " << ins_start << endl;
            cout << "Current state: ";
            for(int i = 0; i < len + 1; i++){
                cout << "\t " << g.cn_profile[c%NUM_CHR][g.chrs[c][ins_start + i].seg_id];
            }
            cout << endl;
            set<double> segs;
            for(int i = 0; i < g.chrs[c].size(); i++){
                segs.insert(g.chrs[c][i].seg_id);
            }
            cout << "After duplication, Chr " << c+1 << " has " << g.chrs[c].size() << " segments (unique: " << segs.size() << ")" << endl;
        }
    }else{
        cerr << "Model not supported!" << endl;
        exit(EXIT_FAILURE);
    }

    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
        cout << "GENOME AFTER duplication " << endl;
        // g.print();
        g.print_cn();
    }
    return 1;
}


// Assume the chosen segment is in the specified chromosome
int generate_deletion(genome& g, int c, int seg_id, int model, int cn_max, int debug){
    if(g.chrs[c].size() <= 0){
        cout << "All segments on chr " << c+1 << " has lost" << endl;
        return 0;
    }
    if(debug){
      cout << "GENOME BEFORE deletion " << endl;
      // g.print();
      g.print_cn();
    }

    double mdel = mean_del_size;
    int len = 0;
    if(mdel > 1)   len = (int) gsl_ran_exponential(r, mdel);

    if(debug){
      cout << "before deletion, Chr " << c+1 << " has " << g.chrs[c].size() << " segments" << endl;
      // for(int i = 0; i < g.chrs[c].size(); i++){
      //     cout << "\t" << g.chrs[c][i].seg_id;
      // }
      // cout << endl;
      cout << "deletion length:" << len << endl;
    }

    // Find the regions to delete
    int del_start = 0;
    int del_msize = 0;
    for(int i = 0; i < g.chrs[c].size(); i++){
        if(g.chrs[c][i].seg_id == seg_id){
            if(g.chrs[c][i-1].seg_id != seg_id){
                del_start = i;
                // cout << i << "\t" << g.chrs[c][i].seg_id << "\t" << seg_id << endl;
                del_msize = 0;
            }
            del_msize++;
        }
    }
    // if(len + del_start >= g.chrs[c].size()){
    if(del_msize < len){
        cout << "Not enough size to delete " << len << ", delete " << del_msize << endl;
        len = del_msize;
        // return 0;
    }
    int del_end = del_start + len + 1;
    if(debug) cout << "deletion start " << del_start << ", end " << del_end << endl;

    if(model == MK){
        // int state = g.cn_profile[c%NUM_CHR][g.chrs[c][loc].seg_id];
        int state = g.cn_profile[c % NUM_CHR][seg_id];
        // Find the number of segments in other haplotypes
        int hap_state = 0;
        int ploidy_at_c = get_ploidy_at_chr(g, c);
        // cout << "Ploidy at " << c << " is " << ploidy_at_c << endl;
        for(int j = 0; j < ploidy_at_c; j++){
            if(j == c/NUM_CHR) continue;
            int hap = c % NUM_CHR + j * NUM_CHR;
            // cout << "Alternative haplotype " << hap << endl;
            for(int i = 0; i < g.chrs[hap].size();i++){
                if(g.chrs[hap][i].seg_id == seg_id)
                {
                    hap_state++;
                }
            }
        }
        if(debug){
            cout << "Copies in other haplotypes: " << hap_state << endl;
            cout << "\tSV: deleting segment, chrs, seg_id, len: " << c << "\t" << seg_id << "\t" << len + 1 << endl;
        }

        vector<int> possible_states;
        int pstate;
        for (int i = 0; i <= state - hap_state - 1; i++){
            possible_states.push_back(i);
        }
        // Randomly pick a possible state
        int ssize = possible_states.size();
        if(debug){
            cout << "There are " << ssize << " possible state(s) for deletion" << endl;
        }
        if(ssize <= 0){
            if(debug) cout << "Impossible to do deletion on Chr " << c << endl;
            return 0;
        }else{
            if(ssize > 1){
                int k = gsl_rng_uniform_int(r, ssize);
                pstate = possible_states[k];
            }else{
                pstate = possible_states[0];
            }
            // Ensure there are "state" copy of the region
            int start = g.chrs[c][del_start].seg_id;
            int end = g.chrs[c][del_start + len].seg_id;
            if(debug){
                cout << "del start " << start << endl;
                cout << "del end " << end << endl;
            }

            if(pstate < 0){
                return 0;
            }

            for(int j = 0; j <= pstate; j++){
                // erase all the segment with seg_id in the specified region
                for(int i = 0; i < g.chrs[c].size();i++){
                    int seg = g.chrs[c][i].seg_id;
                    if(seg >= start && seg <= end){
                        g.chrs[c].erase(g.chrs[c].begin() + i);
                        if(debug){
                            cout << "deleting " << seg << " at position " << i << endl;
                        }
                        break;
                    }
                }
            }
        }
    }else if(model == BOUNDT || model == BOUNDA){
        // Find the minimum copy number of the region to delete when the interval to delete has different segments
        // int min_state = g.cn_profile[c%NUM_CHR][g.chrs[c][loc].seg_id];
        // for (int i = loc; i < loc + len + 1; i++){
        //     int state = g.cn_profile[c%NUM_CHR][g.chrs[c][i].seg_id];
        //     if(state < min_state){
        //         min_state = state;
        //     }
        // }
        // if(min_state <= 0){
        //     return 0;
        // }
        g.chrs[c].erase(g.chrs[c].begin() + del_start, g.chrs[c].begin() + del_end);
    }else{
        cerr << "Model not supported!" << endl;
        exit(EXIT_FAILURE);
    }

    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
        cout << "Current state: " << g.cn_profile[c%NUM_CHR][seg_id] << endl;
        set<double> segs;
        for(int i = 0; i < g.chrs[c].size(); i++){
            segs.insert(g.chrs[c][i].seg_id);
        }
        cout << "After deletion, Chr " << c+1 << " has " << g.chrs[c].size() << " segments (unique: " << segs.size() << ")" << endl;

        cout << "GENOME AFTER deletion " << endl;
        // g.print();
        g.print_cn();
    }
    return 1;
}


int generate_chr_gain(genome& g, int c, int debug){
    // vector<int> available_chrs;
    // get_available_chr(g, available_chrs);
    // if(available_chrs.size()<=0){
    //     // cout << "No available chromosomes!" << endl;
    //     return 0;
    // }
    // int ci = gsl_rng_uniform_int(r, available_chrs.size());
    // int c = available_chrs[ci];
    // int max_cn_c = get_max_cn_chr(g, c);
    // Some segments with maximum copy number may have several copies and hence have copy number larger than specified value when having a chromosome gain (with overlapping events)
    // int max_snum = get_max_seg_num_chr(g, c);
    // if(g.chrs[c].size() <= 0 || max_snum > cn_max){
    //     return 0;
    // }

    // Make sure that the maximum copy number is smaller than the specified value
    int orig_num_chr = get_num_chr(g);
    int orig_size = g.chrs.size();
    assert(orig_size % NUM_CHR == 0);
    // cout << "Size of chrs " << orig_size << endl;
    // assert(orig_size % NUM_CHR == 0);
    int is_in = 0;
    int new_chr = (orig_size / NUM_CHR) * NUM_CHR + c%NUM_CHR; // If c is not in the genome
    // To make sure the chromosome IDs match in different set of copies, add all the haploid chromosomes with others being empty
    for(int i = 0; i < orig_size%NUM_CHR; i++){  // If c is in the genome, but has lost all segments
      if(g.chrs[c%NUM_CHR + i*NUM_CHR].size() <= 0){
          new_chr = c%NUM_CHR + i*NUM_CHR;
          is_in = 1;
          break;
      }
    }
    if(!is_in){
        for(int i = orig_size; i < orig_size+NUM_CHR; i++){ // Fill the position for other chromosomes before the one to duplicate so that the chromosome ID follows the pattern
          vector<segment> e;
          g.chrs.insert(g.chrs.end(), e);
          // cout << "Size of chr " <<  i << " is " << g.chrs[i].size() << endl;
        }
    }

    g.chrs[new_chr].insert(g.chrs[new_chr].begin(), g.chrs[c].begin(), g.chrs[c].end());
    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
      cout << "Chromosome gain in Chromosome ID " << c << endl;
      cout << "ID for gained Chromosome " << new_chr << endl;
      cout << "Copies to make in chr " << c << endl;
      for(int i = 0; i < g.chrs[c].size(); i++){
          cout << "\t" << g.chrs[c][i].chr << "," << g.chrs[c][i].seg_id;
      }
      cout << endl;
      cout << "There are " << orig_num_chr << " non-empty chromosome IDs before chr gain" << endl;
      cout << "There are " << get_num_chr(g) << " non-empty chromosome IDs now" << endl;
      g.print_cn();
    }
    return 1;
}


int generate_chr_loss(genome& g, int c, int debug) {
    // vector<int> available_chrs;
    // get_available_chr(g, available_chrs);
    // if(available_chrs.size()<=0){
    //     // cout << "No available chromosomes!" << endl;
    //     return 0;
    // }
    // int ci = gsl_rng_uniform_int(r, available_chrs.size());
    // int c = available_chrs[ci];
    // int c = gsl_rng_uniform_int(r, g.chrs.size());
    int orig_num_chr = get_num_chr(g);
    // Delete the segments in the chromosome, but keep the chromosome ID so that it can be remapped by NUM_CHR
    g.chrs[c].erase(g.chrs[c].begin(), g.chrs[c].end());
    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
        cout << "Chromosome loss in " << c+1 << endl;
        cout << "There are " << orig_num_chr << " non-empty chromosomes before chr loss" << endl;
        cout << "There are " << get_num_chr(g) << " non-empty chromosomes now" << endl;
        g.print_cn();
    }
    return 1;
}


int generate_wgd(genome& g, int cn_max, int debug) {
    int max_cn_g = get_max_cn_genome(g);
    if(g.chrs.size() <= 0 || 2 * max_cn_g > cn_max){
        return 0;
    }

    int orig_num_chr = get_num_chr(g);
    int orig_size = g.chrs.size();
    // g.chrs.insert(g.chrs.end(), g.chrs.begin(), g.chrs.end());   // cause duplicated insertions on Mac
    vector<segment> s;
    g.chrs.insert(g.chrs.end(), orig_size, s);
    // replicate the segments for each chromosome
    for(int i = 0; i < orig_size; i++){
         // cout << "Insert chr " << i << " for chr " << i + orig_size << endl;
         g.chrs[i + orig_size].insert(g.chrs[i + orig_size].begin(), g.chrs[i].begin(), g.chrs[i].end());
    }
    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
        cout << "Whole genome doubling" << endl;
        cout << "There are " << orig_num_chr << " non-empty chromosomes before WGD" << endl;
        cout << "There are " << get_num_chr(g) << " non-empty chromosomes now" << endl;
        g.print_cn();
        g.print();
    }
    return 1;
}


// Set the initial copy number matrix
copy_number initialize_cn(const vector<int>& chr_lengths, int num_seg, int model){
    copy_number cn;

    for(int i = 0; i < chr_lengths.size(); ++i){
      int num_seg = chr_lengths[i];

      for(int j = 0; j < num_seg; j++){
          if(model == BOUNDT){
             cn[i][j] = 2;     // index of normal state
          }else{
             cn[i][j] = 4;
          }
       }
    }

    return cn;
}


// double*  get_rate_matrix(const vector<double>& rate_consts, int cn_max, int nstate, int model) {
//     if(debug){
//         cout << "Getting rate matrix" << endl;
//     }
//     return qmat;
// }


void print_cn_state(const copy_number& curr_cn){
    for(auto cn_profile : curr_cn){
        for(auto seg: cn_profile.second){
            // chr, seg, cn
            cout <<  cn_profile.first << "\t"  << seg.first << "\t" << seg.second << endl;
        }
    }
}



// Envolving sequences along the tree (available for bounded model, implemented according to approach in book Yang, 2004, P437)
// only support segment duplication and deletion
void evolve_sequences(map<int, copy_number>& cn_matrix, const int& node_id, const evo_tree& tree, double* qmat, int nstate, int num_seg, const vector<double>& rate_consts, map<int, int>& num_muts){
    if(!tree.nodes[node_id].isRoot){
        double rate = rate_consts[0] + rate_consts[1];
        int edge_id = tree.nodes[node_id].e_in;
        double blen = tree.edges[edge_id].length;
        int parent = tree.edges[edge_id].start;
        copy_number curr_cn = cn_matrix[parent];
        copy_number next_cn(curr_cn);

        if(debug) cout << "\nEvolve for node " << node_id << " , parent " << parent << " , branch length " << blen << ", mutation rate " << rate << endl;

        if(blen > 0){
            double *pmat = new double[(nstate) * (nstate)];
            memset(pmat, 0.0, (nstate) * (nstate) * sizeof(double));

            get_transition_matrix_bounded(qmat, pmat, blen, nstate);

            if(debug){
                r8mat_print(nstate, nstate, pmat, "  Transition matrix:");
            }
            // compute gsl_ran_discrete_preproc for all possible states at one time to save computation
            vector<gsl_ran_discrete_t*> dis_all;
            for(int si = 0; si < nstate; ++si){
                double *probs = new double[nstate];
                memset(probs, 0.0, nstate);
                // get the discrete distribution for each state si
                for(int sk = 0; sk < nstate; ++sk){
                    double prob = pmat[si + sk * nstate];
                    probs[sk] = prob;
                    // cout << si << "\t" << sk << "\t" << probs[sk] << endl;
                }
                // Randomly select the next state
                gsl_ran_discrete_t*  dis = gsl_ran_discrete_preproc(nstate, probs);
                dis_all.push_back(dis);

                delete [] probs;
            }

            delete [] pmat;

            if(debug) cout << "size of current CNP " << curr_cn.size() << endl;

            int i = 0;
            for(auto cn_profile : curr_cn){
                int chr = cn_profile.first;
                for(auto seg: cn_profile.second){
                    i++;
                    int seg_id = seg.first;
                    int si = seg.second;
                    // sampling at once based on the transition matrix
                    int sel = gsl_ran_discrete(r, dis_all[si]);
                    if(debug) cout << "Selected state for site " << i << " is " << sel << endl;
                    next_cn[chr][seg_id] = sel;
                }
            }
            // number of mutations is not available when evolving sequencing along the tree. If counting changed sites, multiple hits may be overlooked
            // Compute the number of mutations by Poisson distribution, rate: per allele per site per year
            int num_mut_mean = num_seg * rate * blen * NORM_PLOIDY;
            int num_mut = gsl_ran_poisson(r, num_mut_mean);
            num_muts[edge_id] = num_mut;

            if(debug) cout << "Number of mutations: " << edge_id << "\t" << num_mut << "\t" << num_seg << "\t" << rate << "\t" << blen << endl;
        }

        if(debug){
            cout << "Current copy number state is: ";
            print_cn_state(curr_cn);

            cout << "Next copy number state after time " << blen << " is: ";
            print_cn_state(next_cn);
        }

        cn_matrix[node_id] = next_cn;
    }

    if(!tree.nodes[node_id].daughters.empty()){
      evolve_sequences(cn_matrix, tree.nodes[node_id].daughters[0], tree, qmat, nstate, num_seg, rate_consts, num_muts);
      evolve_sequences(cn_matrix, tree.nodes[node_id].daughters[1], tree, qmat, nstate, num_seg, rate_consts, num_muts);
    }
}


void write_cn(map<int, copy_number>& cn_matrix, int node_id, ogzstream& out, int cn_max, int model){
    copy_number cn_profile = cn_matrix[node_id];

    for(auto c : cn_profile){
        // cout << c.first << endl;
        for(auto s : c.second){
            int cn = s.second;
            if(model == BOUNDA) cn = state_to_total_cn(s.second, cn_max);
            // node_id, chr, seg_id, copy number
            if(debug) cout << node_id + 1 << "\t" << c.first << "\t" << s.first << "\t" << s.second << "\t" << cn << endl;
            out << node_id + 1 << "\t" << c.first << "\t" << s.first << "\t" << cn << endl;
        }
    }
}


void write_allele_cn(map<int, copy_number>& cn_matrix, int node_id, ogzstream& out, int cn_max){
    copy_number cn_profile = cn_matrix[node_id];

    for(auto c : cn_profile){
        // cout << c.first << endl;
        for(auto s : c.second){
            int state = s.second;
            int cnA = 0;
            int cnB = 0;
            state_to_allele_cn(state, cn_max, cnA, cnB);
            // node_id, chr, seg_id, copy number
            if(debug) cout << node_id + 1 << "\t" << c.first << "\t" << s.first << "\t" << s.second << "\t" << cnA  << "\t" << cnB << endl;
            out << node_id + 1 << "\t" << c.first << "\t" << s.first << "\t" << cnA  << "\t" << cnB << endl;
        }
    }
}

// Print the simulated copy numbers
void print_sequences(map<int, copy_number>& cn_matrix, int cn_max, int model, int num_seg, map<int, int> num_muts, evo_tree& test_tree, string dir, string prefix, int age, int print_nex, int print_allele = 1, int cons = 1){
    stringstream sstm;

    sstm << dir << prefix << "-cn.txt.gz";
    ogzstream out_cn(sstm.str().c_str());
    for(int j = 0; j < test_tree.nleaf; ++j){
        write_cn(cn_matrix, j, out_cn, cn_max, model);
    }
    out_cn.close();

    // to get the number of segments
    int num_total_bins = 0;
    int num_invar_bins = 0;
    cout << "reading data and calculating CNA regions" << endl;
    vector<vector<vector<int>>> s_info = read_cn(sstm.str(), test_tree.nleaf-1, num_total_bins, cn_max);
    vector<vector<int>> segs;
    segs = get_all_segs(s_info, test_tree.nleaf - 1, num_total_bins, num_invar_bins, 1);
    int seg_size = segs.size();
    sstm.str("");

    if(print_allele && model == BOUNDA){
        // cout << "Writing allele specific copy number " << endl;
        sstm << dir << prefix << "-allele-cn.txt.gz";
        ogzstream out_allele_cn(sstm.str().c_str());
        for(int j = 0; j < test_tree.nleaf; ++j){
          write_allele_cn(cn_matrix, j, out_allele_cn, cn_max);
        }
        out_allele_cn.close();
        sstm.str("");
    }

    // Output the copy numbers for internal nodes
    sstm << dir << prefix << "-inodes-cn.txt.gz";
    //ofstream out_cn(sstm.str());
    ogzstream out_cn_inodes(sstm.str().c_str());
    for(int j=test_tree.nleaf; j < (2 * test_tree.nleaf - 1); ++j){
        write_cn(cn_matrix, j, out_cn_inodes, cn_max, model);
    }
    out_cn_inodes.close();
    sstm.str("");


    sstm << dir << prefix << "-info.txt";
    ofstream out_info(sstm.str());
    out_info << "NODE TIMES:";
    out_info << "\tnid\ttime" << endl;
    for(int j = 0; j <  test_tree.nodes.size(); ++j){
        out_info << "\t" << j+1 << "\t" << test_tree.nodes[j].time << endl;
    }
    out_info << endl;
    out_info << "SEGMENTS: " << seg_size << endl;
    out_info << endl;
    for(int i = 0; i < test_tree.nleaf; ++i){
        test_tree.print_ancestral_edges(i, out_info);
    }
    out_info << endl;

    out_info.close();
    sstm.str("");

    sstm << dir << prefix << "-tree.txt";
    ofstream out_tree(sstm.str());
    //test_tree.write(out_tree);
    out_tree << "start\tend\tlength\teid\tnmut" << endl;
    for(int i = 0; i < test_tree.edges.size(); ++i){
        out_tree << test_tree.edges[i].start+1 << "\t" << test_tree.edges[i].end+1 << "\t" << test_tree.edges[i].length << "\t" << test_tree.edges[i].id+1 << "\t" << num_muts[test_tree.edges[i].id] << endl;
    }
    out_tree.close();
    sstm.str("");

    if(print_nex){
        sstm << dir << prefix << "-tree.nex";
        ofstream nex_tree(sstm.str());
        string newick = test_tree.make_newick();
        test_tree.write_nexus(newick, nex_tree);
        nex_tree.close();
        sstm.str("");
    }

    if(cons){
        sstm << dir << prefix << "-rel-times.txt";
        double node_min = NODE_MIN_TIME;
        for(int j = 0; j < test_tree.nleaf-1; ++j){
            if(test_tree.nodes[j].time < node_min) node_min = test_tree.nodes[j].time;
        }
        //cout << "leaf minimum: " << node_min << endl;
        ofstream out_rel(sstm.str());
        for(int j = 0; j < test_tree.nleaf-1; ++j){
            double delta = test_tree.nodes[j].time - node_min;
            if(delta < SMALL_VAL) delta = 0;
            out_rel << j+1 << "\t" << delta << "\t" << ceil(age + delta) << endl;
        }
        out_rel.close();
        sstm.str("");
    }

    // Print total number of mutations and total branch length
    int total_mut = 0;
    double total_blen = 0;
    for(int i = 0; i < test_tree.edges.size(); ++i){
        total_blen += test_tree.edges[i].length;
        total_mut += num_muts[test_tree.edges[i].id];
    }
    double mu_est = total_mut / total_blen / num_seg / NORM_PLOIDY;
    cout << "The total number of mutations along branches is: " << total_mut << endl;
    cout << "The total branch length is: " << total_blen << endl;
    cout << "The estimated mutation rate per segment per year is: " << mu_est << endl;
}


// Randomly choose a haplotype of a chromosome
void select_haplotype(genome& g, int& c, long unsigned (*fp_myrng)(long unsigned)){
    // Find the number of possible haplotype of this chromosome
    int ploidy_at_c = get_ploidy_at_chr(g, c);
    assert(ploidy_at_c > 0);

    int hap_c = c;
    vector<int> haps;   // possible haplotypes
    for(int h = 0; h < ploidy_at_c; h++){
        haps.push_back(h);
    }

    // shuffle(haps.begin(), haps.end(), default_random_engine(seed));
    random_shuffle(haps.begin(), haps.end(), fp_myrng);

    for(int h = 0; h < ploidy_at_c; h++){   // randomly choose a haplotype
        hap_c = c + NUM_CHR * haps[h];
        // cout << "Size of chr " << hap_c << " is " << g.chrs[hap_c].size() << endl;
        if(g.chrs[hap_c].size() > 0)    break;
    }
    c = hap_c;
}


// Randomly choose a haplotype of a chromosome based on available segments
void select_haplotype_by_seg(genome& g, int& c, int seg_id, long unsigned (*fp_myrng)(long unsigned)){
    // Find the number of possible haplotype of this chromosome
    int ploidy_at_c = get_ploidy_at_chr(g, c);
    assert(ploidy_at_c > 0);

    int hap_c = c;
    vector<int> haps;   // possible haplotypes
    for(int h = 0; h < ploidy_at_c; h++){
        haps.push_back(h);
    }

    // shuffle(haps.begin(), haps.end(), default_random_engine(seed));
    random_shuffle(haps.begin(), haps.end(), fp_myrng);

    for(int h = 0; h < ploidy_at_c; h++){   // randomly choose a haplotype
        hap_c = c + NUM_CHR * haps[h];
        // cout << "Size of chr " << hap_c << " is " << g.chrs[hap_c].size() << endl;
        // check if haplotype has the chosen segment
        int has_seg = 0;
        if(g.chrs[hap_c].size() > 0){
            for(int i = 0; i < g.chrs[hap_c].size(); i++){
                if(g.chrs[hap_c][i].seg_id == seg_id){
                    has_seg = 1;
                    break;
                }
            }
            if(has_seg) break;
        }
    }
    c = hap_c;
}


// Simulate mutations under different models by simulating waiting times of a Markov chain
// Each mutation occur at some sites on one chromosome (allele-specific)
vector<mutation> generate_mutation_by_model(genome& g, const int& edge_id, const double& blength, const double& node_time, const vector<int>& chr_lengths, const vector<double>& rate_consts, int model, int cn_max, unsigned seed, int& num_fail){
    vector<mutation> ret;
    if(debug){
        cout << "\tGenerate_mutations, blength:" << "\t" << blength << endl;
        cout << "node ID " << g.node_id << endl;
        cout << "edge ID " << edge_id << endl;
        cout << "node time " << node_time << endl;
    }
    double avg_rate = 0.0;  // The mean of genome mutation rates along an edge;

    // Simulate the waiting times of a Markov chain
    double time = 0.0;
    int count = 0;  // Count the number of mutations
    while(time < blength){
        vector<double> site_dup_rates, site_del_rates; // The rate is dynamically changing according to the status of each site, used for selecting duplication/deletion
        vector<double> chr_gain_rates, chr_loss_rates; // The updated rates of chromosome gain/loss
        vector<double> type_rates;  // The updated rates of each mutation type
        if(debug){
          cout << "\tComputing total mutation rate on the genome " << endl;
        }
        double rate = get_total_rates_allele_specific(g, site_dup_rates, site_del_rates, chr_gain_rates, chr_loss_rates, type_rates, rate_consts, model, cn_max);
        if(rate <= 0){
            if(debug) cout << "This genome cannot be mutated any more on edge " << edge_id << endl;
            break;
        }
        avg_rate = avg_rate + rate;
        count++;

        double tevent = gsl_ran_exponential(r, 1 / rate);
        time += tevent;

        // choose type of event
        int e = rchoose(r, type_rates);     // assign the event to a site with probabilities proportional to the rates at the site
        if(debug){
            cout << "There are " << g.chrs.size() << " chromosome IDs (including empty ones) now" << endl;
            cout << "current time " << time << endl;
            cout << "Total mutation rate on the genome " << rate << endl;
            cout << "chosen event " << e << endl;
            // cout << "Sizes of rates vector " << site_dup_rates.size() << "\t" << site_del_rates.size() << "\t" << chr_gain_rates.size() << "\t" << chr_loss_rates.size() << endl;
        }

        int res = 0;
        int site = 0;
        int c = 0;
        int seg_id = 0;

        // Randomly select a site for duplication or deletion. Convert the site position back to chromosome ID and segment ID
        if(e == DUP){   //duplications
            site = rchoose(r, site_dup_rates);
            site2chr(site, c, seg_id, chr_lengths);
            select_haplotype_by_seg(g, c, seg_id, fp_myrng);
            res = generate_duplication(g, c, seg_id, model, cn_max, debug);
            if(debug){
                cout << "chosen site " << site << " is at chr " << c+1 << " seg " << seg_id << endl;
            }
        }else if(e == DEL){   //deletions
            // make sure a deletion somewhere is possible
            site = rchoose(r, site_del_rates);
            site2chr(site, c, seg_id, chr_lengths);
            select_haplotype_by_seg(g, c, seg_id, fp_myrng);
            res = generate_deletion(g, c, seg_id, model, cn_max, debug);
            if(debug){
                // cout << "Deletion rates " << endl;
                // for(int i = 0; i < site_del_rates.size(); i++){
                //        cout << i << "\t" << site_del_rates[i] << endl;
                // }
                cout << "chosen site " << site << " is at chr " << c+1 << " seg " << seg_id << endl;
            }
        }else if(e == GAIN){   //chr gain
            c = rchoose(r, chr_gain_rates);
            seg_id = -1;
            select_haplotype(g, c, fp_myrng);
            res = generate_chr_gain(g, c, debug);
            if(debug){
                cout << "chosen chr " << c+1 << " for gain event" << endl;
            }
        }else if(e == LOSS){   //chr loss
            c = rchoose(r, chr_loss_rates);
            seg_id = -1;
            select_haplotype(g, c, fp_myrng);
            res = generate_chr_loss(g, c, debug);
            if(debug){
                cout << "chosen chr " << c+1 << " for loss event" << endl;
            }
        }else if(e == WGD){   // whole genome duplication
            c = -1;
            seg_id = -1;
            res = generate_wgd(g, cn_max, debug);
        }
        else{
            cerr << "Unknown mutation type " << e << endl;
            exit(EXIT_FAILURE);
        }

        if(res > 0){    // mutations are generated successfully
            // time/blength: relative time in the branch; node_time+time: real time when the mutation generates
            mutation mut(edge_id, e, time / blength, node_time + time, c, seg_id);
            ret.push_back(mut);
            g.nmuts[e]++;
            g.mutations.push_back(mut);
            // time += tevent;
        }else{
             // if(debug)
             // cout << "\tMutation failed:" << endl;
             cout << "failed: tevent, total time, time/branch len, event, chr, site, loc/seg_id, copy_number\n" << g.node_id+1 << "\t" << edge_id+1 << "\t" << tevent << "\t" << time << "\t" << blength << "\t" << e << "\t" << c+1 << "\t" << site << "\t" << seg_id << "\t" << g.cn_profile[c%NUM_CHR][seg_id] << endl;
             cout << "all segments on chr " << c << ":";
             for(int i = 0; i < g.chrs[c].size(); i++){
                 cout << "\t" << g.chrs[c][i].seg_id;
             }
             cout << endl;

             g.print_cn();
             // time -= tevent;
             num_fail++;
             continue;
        }

        if(debug){
            cout << "node ID, edge ID, mut times, tevent, total time, time/branch len, event, chr, loc\t" << g.node_id+1 << "\t" << edge_id+1 << "\t" << tevent << "\t" << time << "\t" << blength << "\t" << e << "\t" << c+1 << "\t" << seg_id << endl;
        }
    }

    if(count > 0){
        avg_rate = avg_rate / count;
    }
    cout << "Average genome mutation rate on edge (per year) " << edge_id << " is: " << avg_rate << endl;

    return ret;
}


// simulate mutations under infinite sites model
vector<mutation> generate_mutation_times(const int& edge_id, const double& blength, const double& node_time, const vector<double>& rate_consts){
  if(debug) cout << "\tgenerate_mutations, blength:" << "\t" << blength << endl;

  vector<mutation> ret;

  // model this as a poisson process
  double time = 0.0;
  while(time < blength){
    vector<double> rates;
    for(int i = 0; i < rate_consts.size(); ++i){
      rates.push_back(rate_consts[i] / 2.0);        // in coalescent Nmut ~ Pois( theta*l/2)
      //cout << "##:\t" << rates[i] << endl;
    }

    double rate = accumulate(rates.begin(), rates.end(), 0.0);
    double tevent = gsl_ran_exponential (r, 1 / rate);
    time += tevent;

    // choose type of event
    int e = rchoose(r, rates);
    if(debug) cout << "mut times, tevent, total time, time/branch len, event\t" << tevent << "\t" << time << "\t" << blength << "\t" << e << endl;

    ret.push_back(mutation( edge_id, e, time / blength, node_time + time, 0, 0));
  }
  ret.pop_back();

  return ret;
}


void apply_mutations(const vector<mutation>& muts, genome& g){
  // mean variant size in bins
  double mdup = mean_dup_size;
  double mdel = mean_del_size;

  for(int i = 0; i < muts.size(); ++i){
    g.nmuts[ muts[i].type]++;
    g.mutations.push_back(muts[i]);

    if(muts[i].type == DUP){   //duplications
        if(debug){
            cout << "GENOME BEFORE duplication " << endl;
            //g.print();
            g.print_cn();
        }

        // make sure a duplication somewhere is possible
        bool possible = false;
        for(int j = 0; j < g.chrs.size(); ++j){
            if(g.chrs[j].size() > 1) possible = true;
        }

        if(possible){
            bool done = false;
            while(!done){
                int c = gsl_rng_uniform_int(r, g.chrs.size());
                int len = gsl_ran_exponential(r, mdup);
                //int len =  1 + gsl_ran_poisson (r, ldup);
                //cout << "dup len:" << len << endl;

                if(len < g.chrs[c].size()){
                    // choose a random location up to size-len
                    int loc = gsl_rng_uniform_int(r, g.chrs[c].size() - len);

                    if(debug) cout << "\tSV: duplicating segment, chr, start, len: " << c << "\t" << loc << "\t" << len + 1 << endl;
                    //cout << "SV: insert: " << loc + len + 1 << "\t" << loc << "\t" << loc + len + 1 << endl;
                    g.chrs[c].insert(g.chrs[c].begin() + loc + len + 1, g.chrs[c].begin() + loc, g.chrs[c].begin() + loc + len + 1);
                    done = true;
                }
            }
        }else{
            if(debug) cout << "\tSV: duplication failed:" << endl;
        }

        if(debug){
            cout << "GENOME AFTER duplication " << endl;
            //g.print();
            g.print_cn();
        }
    }else if(muts[i].type == DEL){   //deletions
      if(debug){
        cout << "GENOME BEFORE deletion " << endl;
        //g.print();
        g.print_cn();
      }

      // make sure a deletion somewhere is possible
      bool possible = false;
      for(int j = 0; j < g.chrs.size(); ++j){
          if(g.chrs[j].size() > 1) possible = true;
      }

      if(possible){
        bool done = false;
        while(!done){
          int c = gsl_rng_uniform_int(r, g.chrs.size());
          int len = (int) gsl_ran_exponential(r, mdel);
          //int len =  1 + gsl_ran_poisson (r, ldel);
          //cout << "del len:" << len << endl;

          if(len < g.chrs[c].size()){
            int loc = gsl_rng_uniform_int(r,g.chrs[c].size()-len);
            if(debug) cout << "\tSV: deleting segment, chrs, start, len: " << c << "\t" << loc << "\t" << len + 1 << endl;
            g.chrs[c].erase(g.chrs[c].begin() + loc, g.chrs[c].begin() + loc + len + 1);
            done = true;
          }
        }
      }else{
          if(debug) cout << "\tSV: deletion failed:" << endl;
      }

      if(debug){
        cout << "GENOME AFTER deletion " << endl;
        //g.print();
        g.print_cn();
      }
    }else if(muts[i].type == LOSS){   //chr loss
      if(g.chrs.size() > 0){
        int c = gsl_rng_uniform_int(r, g.chrs.size());
        // g.chrs.erase(g.chrs.begin()+c);
        g.chrs[c].erase(g.chrs[c].begin(), g.chrs[c].end());
      }

    }else if(muts[i].type == GAIN){   //chr gain
      if(g.chrs.size() > 0){
        int c = gsl_rng_uniform_int(r, g.chrs.size());
        int orig_num_chr = g.chrs.size();
        g.chrs.insert(g.chrs.end(), g.chrs.begin()+c, g.chrs.begin()+c+1);
        g.chrs[orig_num_chr].insert(g.chrs[orig_num_chr].begin(), g.chrs[c].begin(), g.chrs[c].end());
      }
    }else if(muts[i].type == WGD){   // whole genome duplication
      if(g.chrs.size() > 0){
          int orig_num_chr = g.chrs.size();
          g.chrs.insert(g.chrs.end(), g.chrs.begin(), g.chrs.end());
          // replicate the segments for each chromosome
          for(int i = 0; i < orig_num_chr; i++){
               g.chrs[i + orig_num_chr].insert(g.chrs[i + orig_num_chr].begin(), g.chrs[i].begin(), g.chrs[i].end());
          }
      }
    }else{
      cerr << "Unknown mutation type!" << endl;
      exit(EXIT_FAILURE);
    }
  }
}



void traverse_tree_mutating(const int& node_id, const evo_tree& tree, const vector<int>& chr_lengths, const vector<double>& rate_consts, map<int, vector<mutation>>& all_muts, map<int, int>& failed_muts, vector<genome>& genomes, int model, int cn_max, unsigned seed){
  //cout << "\ttraverse_tree: " << node_id+1 << endl;
  if(!tree.nodes[node_id].isRoot){
    // copy the parent
    genomes[node_id] = genomes[tree.nodes[node_id].parent];
    genomes[node_id].node_id = node_id;

    // apply the mutations from parent -> daughter
    //cout << "MUTATING genome: " << tree.nodes[node_id].parent+1 << " -> " << node_id+1 << "\t edge id: " << tree.nodes[node_id].e_in+1 << endl;
    int edge_id = tree.nodes[node_id].e_in;
    vector<mutation> muts;
    int num_fail = 0;
    if(tree.edges[edge_id].length > 0){
        if(debug){
            cout << "Generating mutations for node " << node_id << endl;
        }
        if(model == INFINITE){  
            muts = generate_mutation_times(edge_id, tree.edges[edge_id].length, tree.nodes[tree.nodes[node_id].parent].time, rate_consts);
            apply_mutations( muts, genomes[node_id]);
        }else{
            muts = generate_mutation_by_model(genomes[node_id], edge_id, tree.edges[edge_id].length, tree.nodes[tree.nodes[node_id].parent].time, chr_lengths, rate_consts, model, cn_max, seed, num_fail);
        }
        // Print out copy number for checking
        if(debug){
            cout << "Generated mutations for node " << node_id << endl;
            map<int, map<int,int> >::const_iterator it1;
            map<int,int>::const_iterator it2;
            for(it1 = genomes[node_id].cn_profile.begin(); it1 != genomes[node_id].cn_profile.end(); ++it1){
              for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
                cout << node_id+1 << "\t" << it1->first+1 << "\t" << it2->first << "\t" << it2->second << endl;
              }
            }
        }
    }
    all_muts[edge_id] = muts;
    failed_muts[edge_id] = num_fail;
  }

  if(!tree.nodes[node_id].daughters.empty()){
    traverse_tree_mutating(tree.nodes[node_id].daughters[0], tree, chr_lengths, rate_consts, all_muts, failed_muts, genomes, model, cn_max, seed);
    if(debug) cout << "Finished for node " << tree.nodes[node_id].daughters[0] << endl;
    traverse_tree_mutating(tree.nodes[node_id].daughters[1], tree, chr_lengths, rate_consts, all_muts, failed_muts, genomes, model, cn_max, seed);
    if(debug) cout << "Finished for node " << tree.nodes[node_id].daughters[1] << endl;
  }

}



void simulate_samples(vector<genome>& genomes, map<int,vector<mutation>>& muts, map<int,int>& failed_muts, const evo_tree& tree, genome& germline, const vector<int>& chr_lengths, const vector<double>& rate_consts, int model, int cn_max, unsigned seed, int verbose = 0){
  // assign the germline to the root of the tree
  germline.node_id = tree.root_node_id;

  for(int i = 0; i < (2 * tree.nleaf - 1); ++i){
    if(i == tree.root_node_id){
      genomes.push_back(germline);
    }else{
      genomes.push_back(genome());
    }
  }

  // move through the evolutionary tree mutating genomes
  traverse_tree_mutating(tree.root_node_id, tree, chr_lengths, rate_consts, muts, failed_muts, genomes, model, cn_max, seed);

  // final samples returned, print out leaf nodes
  if(verbose){
    cout << "MUTATIONS:" << endl;
    for(map<int, vector<mutation> >::iterator it = muts.begin(); it != muts.end(); it++){
      cout << "EDGE, id: " << it->first+1
       << "\t" << tree.edges[it->first].start+1 << " -> " << tree.edges[it->first].end+1
       << "\t" << it->second.size() << endl;
      vector<mutation> v = it->second;
      for(int i = 0; i < v.size(); ++i){
          v[i].print();
      }
    }
    cout << endl;
    cout << "LEAF GENOMES:" << endl;
    for(int i = 0; i < tree.nleaf; ++i){
      //genomes[i].print();
      genomes[i].print_muts();
      //genomes[i].print_cn();
      tree.print_ancestral_edges( genomes[i].node_id);
      cout << endl;
    }
  }
}


void run_sample_set(int Ns, gsl_rng* r, unsigned seed, int Ne, double beta, const vector<double>& rate_consts, double* pvs, int* ret){
  // static const int arr[] = {367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};
  // vector<int> chr_lengths(arr, arr + sizeof(arr) / sizeof(arr[0]));
  vector<int> chr_lengths{367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};

  int model = BOUNDA;
  int cn_max = 4;

  // vector<double> rate_consts (prc, prc + sizeof(prc) / sizeof(prc[0]));
  genome germline(chr_lengths, 2);

  vector<int> edges;

  vector<double> lengths;
  vector<double> epoch_times;
  vector<double> node_times;
  stringstream sstm;
  vector<genome> results;
  map<int, vector<mutation>> muts;
  map<int, int> failed_muts;
  //cout << "\n\n###### New sample collection ######" << endl;
  //cout << "###### Ns+1= " << Ns+1 << endl;

  generate_coal_tree(Ns, r, fp_myrng, edges, lengths, epoch_times, node_times, Ne, beta);
  evo_tree test_tree(Ns+1, edges, lengths);

  //for(int i = 0; i < 6; ++i) epars.push_back( ptree[i]);
  //for(int i = 0; i < 8; ++i) lengths.push_back( pl[i]);
  //evo_tree test_tree = construct_tree(Ns, epars, lengths, node_times);

  simulate_samples(results, muts, failed_muts, test_tree, germline, chr_lengths, rate_consts, model, cn_max, seed, false);

  int nbins = 4401;
  for(int i = 0; i < (test_tree.nleaf-1); ++i){
    vector<int> cns = results[i].get_cn_vector();
    for(int j = 0; j < nbins; ++j){
      ret[nbins*i + j] = cns[j];
    }
  }

  if(0){
    sstm << "test-data-cn.txt.gz";
    ogzstream out_cn(sstm.str().c_str());
    for(int j = 0; j < test_tree.nleaf; ++j){
      results[j].write(out_cn);
    }
    out_cn.close();
    sstm.str("");
  }

}


int get_num_seg(const vector<int>& chr_lengths){
    int total_seg = 0;
    for(int i = 0; i < NUM_CHR; i++){
        total_seg += chr_lengths[i];
        // cout << alpha[i] << "\t" << theta[i] << "\t" << chr_lengths[i] << endl;
    }
    return total_seg;
}


void run_test(int mode, string dir, unsigned seed){
    //create test tree
    if(mode == 2){
      //static const int arr1[] = {7,8, 6,7, 8,1, 8,2, 7,9, 9,3, 9,4, 6,5 };
      static const int arr1[] = {6,7, 5,6, 7,0, 7,1, 6,8, 8,2, 8,3, 5,4};
      vector<int> e(arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]));
      static const double arr2[] = {5.1,6.3,10.2,9.5,5.2,3.2,5.4,0};
      vector<double> l(arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]));
      evo_tree test_tree(5, e, l);
      test_tree.print();
    }

    // test out the genome classes and functionality
    if(mode == 3){
      if(0){
        vector<mutation> dels;
        dels.push_back( mutation( 0, 1, 0.5, 0, 1, 2));
        dels.push_back( mutation( 0, 1, 0.5, 0, 2, 3));
        dels.push_back( mutation( 0, 1, 0.5, 0, 3, 4));
        dels.push_back( mutation( 0, 1, 0.5, 0, 4, 5));

        vector<mutation> dups;
        dups.push_back( mutation( 0, 0, 0.5, 0, 11, 2));
        dups.push_back( mutation( 0, 0, 0.5, 0, 12, 3));
        dups.push_back( mutation( 0, 0, 0.5, 0, 13, 4));
        dups.push_back( mutation( 0, 0, 0.5, 0, 14, 5));

        vector<mutation> wgds;
        wgds.push_back( mutation( 0, 4, 0.5, 0, -1, -1));

        genome g1(2,10);
        g1.node_id = 0;
        g1.print();
        g1.print_cn();
        apply_mutations(dups, g1);
        apply_mutations(dels, g1);
        apply_mutations(wgds, g1);
        g1.print();
        g1.print_cn();
      }

      if(1){
        static const int arr3[] = {10,4,3};
        vector<int> chr_lengths (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]));

        genome g2(chr_lengths);
        g2.node_id = 0;
        g2.print();
        g2.print_cn();

        genome g2d(chr_lengths,2);
        g2d.node_id = 0;
        g2d.print();
        g2d.print_cn();

      }

      //genome g3 = g2;
      //g3.print();
      //g3.print_cn();
    }
    if(mode == 4){
      stringstream sstm;
      int Ns = 4;
      int Ne = 0;
      double beta = 0;
      // double prates[] = {0.1, 0.1, 0.1, 0.1, 0.05};
      vector<double> rate_consts{0.1, 0.1, 0.1, 0.1, 0.05};
      double pvs[] = {30.0, 30.0};

      for(int i = 0; i < 10; ++i){
        int* ret = new int[Ns*4401];
        run_sample_set(Ns, r, seed, Ne, beta, rate_consts, pvs, &(ret[0]));

        sstm << dir << "sim-data-" << i+1 << "-cn.txt";
        ofstream out_cn(sstm.str());
        for(int i = 0; i < (Ns*4401); ++i){
            out_cn << ret[i] << endl;
        }
        out_cn.close();
        sstm.str("");
        delete[] ret;
      }
    }
}


void print_simulations(int mode, int num_seg, vector<genome>& results, const vector<int>& chr_lengths, map<int, vector<mutation>>& muts, map<int, int>& failed_muts, evo_tree& test_tree, string dir, string prefix, int age, int cn_max, int print_allele, int print_mut, int print_nex, int cons = 1){
    stringstream sstm;

    sstm << dir << prefix << "-cn.txt.gz";
    ogzstream out_cn(sstm.str().c_str());
    for(int j = 0; j < test_tree.nleaf; ++j){
        results[j].write(out_cn);
    }
    out_cn.close();
    int num_total_bins = 0;
    int num_invar_bins = 0;
    cout << "reading data and calculating CNA regions" << endl;
    vector<vector<vector<int>>> s_info = read_cn(sstm.str(), test_tree.nleaf-1, num_total_bins, cn_max);
    vector<vector<int>> segs;
    if(mode == 0){  // merge consecutive bins with the same CN
        segs = get_invar_segs(s_info, test_tree.nleaf-1, num_total_bins, num_invar_bins);
    }else{
        segs = get_all_segs(s_info, test_tree.nleaf-1, num_total_bins, num_invar_bins, 1);
    }
    int seg_size = segs.size();
    sstm.str("");

    if(print_allele){
        // cout << "Writing allele specific copy number " << endl;
        sstm << dir << prefix << "-allele-cn.txt.gz";
        ogzstream out_allele_cn(sstm.str().c_str());
        for(int j = 0; j < test_tree.nleaf; ++j){
          // cout << "Chr size " << results[j].chrs.size() << endl;
          results[j].write_allele_cn(out_allele_cn, chr_lengths);
        }
        out_allele_cn.close();
        sstm.str("");
    }

    int ntotn = 2 * test_tree.nleaf - 1;
    // Output the copy numbers for internal nodes
    sstm << dir << prefix << "-inodes-cn.txt.gz";
    //ofstream out_cn(sstm.str());
    ogzstream out_cn_inodes(sstm.str().c_str());
    for(int j=test_tree.nleaf; j < ntotn; ++j){
      results[j].write(out_cn_inodes);
    }
    out_cn_inodes.close();
    sstm.str("");

    if(print_allele){
        sstm << dir << prefix << "-inodes-allele-cn.txt.gz";
        ogzstream out_allele_cn_inodes(sstm.str().c_str());
        for(int j=test_tree.nleaf; j < ntotn; ++j){
          // cout << "Chr size " << results[j].chrs.size() << endl;
          results[j].write_allele_cn(out_allele_cn_inodes, chr_lengths);
        }
        out_allele_cn_inodes.close();
        sstm.str("");
    }

    sstm << dir << prefix << "-info.txt";
    ofstream out_info(sstm.str());
    out_info << "NODE TIMES:";
    out_info << "\tnid\ttime" << endl;
    for(int j = 0; j <  test_tree.nodes.size(); ++j){
        out_info << "\t" << j+1 << "\t" << test_tree.nodes[j].time << endl;
    }
    out_info << endl;
    out_info << "SEGMENTS: " << seg_size << endl;
    out_info << endl;
    for(int i = 0; i < test_tree.nleaf; ++i){
        test_tree.print_ancestral_edges(results[i].node_id, out_info);
    }
    out_info << endl;
    // number of mutations is known when simulating waiting time
    for(int i = 0; i < test_tree.nleaf; ++i){
        results[i].print_muts(out_info);
    }
    out_info.close();
    sstm.str("");

    sstm << dir << prefix << "-tree.txt";
    vector<int> nmuts;
    ofstream out_tree(sstm.str());
    //test_tree.write(out_tree);
    out_tree << "start\tend\tlength\teid\tnmut\tnfailed" << endl;
    for(int i = 0; i < test_tree.edges.size(); ++i){
        out_tree << test_tree.edges[i].start+1 << "\t" << test_tree.edges[i].end+1 << "\t" << test_tree.edges[i].length << "\t" << test_tree.edges[i].id+1 << "\t" << muts[test_tree.edges[i].id].size() << "\t" << failed_muts[test_tree.edges[i].id] << endl;
      nmuts.push_back(muts[test_tree.edges[i].id].size());
    }
    out_tree.close();
    sstm.str("");

    if(print_nex){
        sstm << dir << prefix << "-tree.nex";
        ofstream nex_tree(sstm.str());
        string newick = test_tree.make_newick();
        test_tree.write_nexus(newick, nex_tree);
        nex_tree.close();
        sstm.str("");

        sstm << dir << prefix << "-tree-nmut.nex";
        ofstream nex_tree2(sstm.str());
        newick = test_tree.make_newick_nmut(0, nmuts);
        test_tree.write_nexus(newick, nex_tree2);
        nex_tree2.close();
        sstm.str("");
    }

    if(print_mut){
        sstm << dir << prefix << "-mut.txt";
        ofstream out_mut(sstm.str());
        for(int j = 0; j < test_tree.nleaf; ++j){
            for(int i = 0; i < results[j].mutations.size(); ++i){
                int chr = results[j].mutations[i].chr;
                out_mut << j+1 << "\t" << results[j].mutations[i].edge_id+1
        << "\t" << results[j].mutations[i].type << "\t" << results[j].mutations[i].btime << "\t" << results[j].mutations[i].gtime << "\t" << chr+1 << "\t" << (chr+1)%22 << "\t" << results[j].mutations[i].seg << endl;
            }
        }
        out_mut.close();
        sstm.str("");
    }

    if(cons){
        sstm << dir << prefix << "-rel-times.txt";
        double node_min = NODE_MIN_TIME;
        for(int j = 0; j < test_tree.nleaf-1; ++j){
            if(test_tree.nodes[j].time < node_min) node_min = test_tree.nodes[j].time;
        }
        //cout << "leaf minimum: " << node_min << endl;
        ofstream out_rel(sstm.str());
        for(int j = 0; j < test_tree.nleaf-1; ++j){
            double delta = test_tree.nodes[j].time - node_min;
            if(delta < SMALL_VAL) delta = 0;
            out_rel << j+1 << "\t" << delta << "\t" << ceil(age + delta) << endl;
        }
        out_rel.close();
        sstm.str("");
    }

    // Print total number of mutations and total branch length
    int total_mut = 0;
    int failed_mut = 0;
    int succ_mut = 0;
    double total_blen = 0;
    for(int i = 0; i < test_tree.edges.size(); ++i){
        total_blen += test_tree.edges[i].length;
        total_mut += muts[test_tree.edges[i].id].size();
        total_mut += failed_muts[test_tree.edges[i].id];
        succ_mut += muts[test_tree.edges[i].id].size();
        failed_mut += failed_muts[test_tree.edges[i].id];
    }
    double mu_est = total_mut / total_blen / num_seg / NORM_PLOIDY;
    cout << "The total number of mutations along branches is: " << total_mut << endl;
    cout << "The total number of successful mutations is: " << succ_mut << endl;
    cout << "The total number of failed mutations (due to limitations on copy number) is: " << failed_mut << endl;
    cout << "The total branch length is: " << total_blen << endl;
    cout << "The estimated mutation rate per segment per year is: " << mu_est << endl;
}


void run_simulations(string tree_file, int mode, int method, const vector<int>& chr_lengths, int num_seg, int Ns, int Nsims, int cn_max, int model, int cons, int Ne, double beta, double gtime, double delta_t, int age, const vector<double>& rate_consts, string dir, string prefix, int print_allele, int print_mut, int print_nex, int verbose = 0){
    // can specify the germline in a number of different ways
    //genome germline(1,10);
    //genome germline(42,1000);
    //genome germline(chr_lengths);
    genome germline(chr_lengths, NORM_PLOIDY);

    string orig_prefix = prefix;
    cout << "\nNumber of datasets to simulate " << Nsims << endl;

    for(int i = 0; i < Nsims; ++i){
        vector<genome> genomes;      // one genome for each node on the tree
        map<int, vector<mutation>> muts;   // hold all mutations by edge id
        map<int, int> failed_muts;  // the number of failed mutations on each edge

        cout << "\nSimulation " << i + 1 << endl;

        //cout << "\n\n###### New sample collection ######" << endl;
        if(orig_prefix.empty()) {
          prefix = "sim-data-" + to_string(i + 1);
        }
        cout << "Prefix of output file: " << prefix <<endl;

        evo_tree test_tree;
        if(tree_file != ""){
            test_tree = read_tree_info(tree_file, Ns);
        }else{
            test_tree = generate_random_tree(Ns, r, fp_myrng, Ne, age, beta, gtime, delta_t, cons, debug);
        }
        if(debug){
          test_tree.print();
          cout << "Simulated tree height " << get_tree_height(test_tree.get_node_times()) << endl;
        }


        if(method == SIM_TIME){     // applicable for all models 
            simulate_samples(genomes, muts, failed_muts, test_tree, germline, chr_lengths, rate_consts, model, cn_max, seed, verbose);
            int num_seg = get_num_seg(chr_lengths);
            print_simulations(mode, num_seg, genomes, chr_lengths, muts, failed_muts, test_tree, dir, prefix, age, cn_max, print_allele, print_mut, print_nex, cons);
        }else{  // model can only be BOUNDA or BOUNDT
            assert(model == BOUNDA || model == BOUNDT);
            int root = Ns + 1;
            int nstate = cn_max + 1;
            if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;
            // cout << "number of states now " << nstate << "\t" << cn_max << endl;

            map<int, copy_number> cn_matrix; // copy number matrix for each node
            map<int, int> num_muts; // number of mutations on each edge
            copy_number init_cn = initialize_cn(chr_lengths, num_seg, model);
            if(debug){
                cout << "There are " << nstate << " possible states" << endl;
                cout << "Initial copy number state at " << root << ": " << endl;
                print_cn_state(init_cn);
            }
            cn_matrix[root] = init_cn;

            int msize = (nstate) * (nstate);
            double *qmat = new double[msize];
            memset(qmat, 0.0, msize * sizeof(double));

            double dup_rate = rate_consts[0];
            double del_rate = rate_consts[1];
            if(model == BOUNDA){
                if(debug){
                    cout << "\tGetting allele specific rate matrix" << endl;
                }
                get_rate_matrix_allele_specific(qmat, dup_rate, del_rate, cn_max);
            }else{
                if(debug){
                    cout << "\tGetting total rate matrix" << endl;
                }
                get_rate_matrix_bounded(qmat, dup_rate, del_rate, cn_max);
            }

            if(debug){
                r8mat_print(nstate, nstate, qmat, "  Rate matrix:");
                check_matrix_row_sum(qmat, nstate);
            }

            evolve_sequences(cn_matrix, root, test_tree, qmat, nstate, num_seg, rate_consts, num_muts);
            print_sequences(cn_matrix, cn_max, model, num_seg, num_muts, test_tree, dir, prefix, age, print_nex, print_allele, cons);

            delete [] qmat;
        }
    }
    // cout << "finish simulations" << endl;
}

//////////////////////////////////////////////////////////
///                                                    ///
///   MAIN                                             ///
///                                                    ///
//////////////////////////////////////////////////////////
int main (int argc, char** const argv) {
    /********* major input parameters ***********/
    int Ns;   // number of samples
    int age;

    int model;
    int cons;
    int cn_max;
    /*************************************/

    string dir; // output directory
    string prefix; // prefix of output file
    int Nsims;  // number of multi-region samples
    // five event types: duplication, deletion, chromosome gain, chromosome loss, wgd
    // rates are 1/mean
    double dup_rate, del_rate, chr_gain, chr_loss, wgd;
    // parameters for mean of dup/del size distributions
    double dup_size, del_size;
    // effective population size
    int Ne;
    double beta, gtime;    // population growth rate
    double delta_t;    // relative timing difference
    int mode, method;
    int seg_max, fix_nseg;
    int print_allele, print_mut, print_nex;
    string tree_file, stime;

    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ;
    po::options_description required("Required parameters");
    required.add_options()
      ("odir,o", po::value<string>(&dir)->required(), "output directory")
       ;
    po::options_description optional("Optional parameters");
    optional.add_options()
      ("nsim,n", po::value<int>(&Nsims)->default_value(1), "number of multi-region samples")

      ("tree_file", po::value<string>(&tree_file)->default_value(""), "input tree file. Mutations will be generated along this tree if provided.")
      ("nregion,r", po::value<int>(&Ns)->default_value(5), "number of regions")
      ("age,a", po::value<int>(&age)->default_value(MAX_AGE), "age of the patient to simulate")

      ("mode", po::value<int>(&mode)->default_value(0), "running mode of the program (0: simulating genome in fix-sized bins (4401 bins of size 500 Kbp by default), 1: simulating genome in segments of unknown size, 2 to 4: test)")
      ("method", po::value<int>(&method)->default_value(0), "method of simulation (0: simulating waiting times, 1: simulating sequences directly)")
      ("model,d", po::value<int>(&model)->default_value(0), "model of evolution (0: Mk, 1: one-step bounded (total), 2: one-step bounded (allele-specific), 3: infinite sites)")

      // for generating random coalescence tree
      ("epop,e", po::value<int>(&Ne)->default_value(2), "effective population size of cell populations")
      ("gtime", po::value<double>(&gtime)->default_value(0.002739726), "generation time in year")
      ("beta,b", po::value<double>(&beta)->default_value(0.0), "population growth rate")
      ("tdiff,t", po::value<double>(&delta_t)->default_value(0), "relative timing difference to earliest sample")
      // ("stime", po::value<string>(&stime)->default_value(""), "sampling time of different samples (Format: numTipDates Date1 from to Date2 from to ... DateN from to)")
      ("constrained", po::value<int>(&cons)->default_value(1), "whether or not to constrain tree height by patient age. If yes (1), the initial branch lengths will be adjusted by specified patient age so that the tree height is smaller than patient age.")

      // segment options   
      ("cn_max", po::value<int>(&cn_max)->default_value(4), "maximum copy number of a segment")
      ("seg_max", po::value<int>(&seg_max)->default_value(100), "maximum number of segments to simulate")
      ("fix_nseg", po::value<int>(&fix_nseg)->default_value(1), "whether or not to fix the number of segments to simulate. If not, the number of segments on each chromosome is proportional to chromosome size")

      // rate and size of mutations
      ("dup_rate", po::value<double>(&dup_rate)->default_value(0.001), "duplication rate (allele/locus/year)")
      ("del_rate", po::value<double>(&del_rate)->default_value(0.001), "deletion rate (allele/locus/year)")
      ("chr_gain", po::value<double>(&chr_gain)->default_value(0.0), "chromosome gain rate (haplotype/chr/year)")
      ("chr_loss", po::value<double>(&chr_loss)->default_value(0.0), "chromosome loss rate (haplotype/chr/year)")
      ("wgd", po::value<double>(&wgd)->default_value(0.0), "WGD (whole genome doubling) rate (year)")
      ("dup_size", po::value<double>(&mean_dup_size)->default_value(1), "mean of duplication size distributions (in bins)")
      ("del_size", po::value<double>(&mean_del_size)->default_value(1), "mean of deletion size distributions (in bins)")

      // output options
      ("prefix,p", po::value<string>(&prefix)->default_value(""), "prefix of output file (sim-data-N if not specified")
      ("print_allele", po::value<int>(&print_allele)->default_value(1), "whether or not to output allele-specific copy numbers")
      ("print_mut", po::value<int>(&print_mut)->default_value(1), "whether or not to output the list of mutations")
      ("print_nex", po::value<int>(&print_nex)->default_value(1), "whether or not to output the tree in NEXUS format")

      ("seed", po::value<unsigned>(&seed)->default_value(0), "seed used for generating random numbers")
      ("verbose", po::value<int>(&debug)->default_value(0), "verbose level (0: default, 1: debug)")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(required).add(optional);
    po::variables_map vm;

    try{
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
        if(vm.count("help")){
            cout << cmdline_options << endl;
            return 1;
        }
        if(vm.count("version")){
            cout << "sveta [version 0.1], a program to simulate structural variations along a phylogenetic tree" << endl;
            return 1;
        }
        po::notify(vm);
    }catch(const std::exception& e){
        std::cerr << e.what() << std::endl;
        return 1;
    }

    if(delta_t > 0){
        cout << "Simulating time differences among samples with baseline time difference " << delta_t << endl;
    }

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    setup_rng(r, seed);

    fp_myrng = &myrng;

    // output directory
    if(dir.back() != '/'){
        dir = dir + "/";
    }

    vector<int> chr_lengths = CHR_BIN_SIZE;
    int num_seg = 0;

    if(method == SIM_SEQ) {   // when simulating sequences directly, each site is a final segment
        cout << "Simulating sequences directly.\n\tEach site is a final segment in this way of simulation" << endl;
        mode = 1;
    }else{
        cout << "\nSimulating waiting times" << endl;
    }

    if(mode == 1){
        cout << "Under this mode, each site is treated as a final segment. The mean duplication/deletion size is currently fixed to be 1" << endl;
        mean_dup_size = 1;
        mean_del_size = 1;
        // Randomly generate the number of segments
        if(fix_nseg){
             num_seg = seg_max;
        }else{
             num_seg = runiform(r, 0.2, 0.8) * seg_max;      // use 0.2 to avoid too few segments; use 0.8 to avoid more than seg_max segments
             if(num_seg > seg_max)  num_seg = seg_max;
        }
        // cout << "Approximate number of segments to simulate is " << num_seg << endl;
        // Distribute the segments according to the size of chromosomes
        double theta[NUM_CHR] = {};
        double alpha[NUM_CHR] = {};
        // int bin_size = accumulate(CHR_BIN_SIZE.begin(), CHR_BIN_SIZE.end(), 0);
        // cout << "Total number of bins is " << bin_size << endl;
        int total_seg = 0;
        for(int i = 0; i < NUM_CHR; i++){
            alpha[i] = CHR_BIN_SIZE[i];
        }
        gsl_ran_dirichlet(r, NUM_CHR, alpha, theta);
        for(int i = 0; i < NUM_CHR; i++){
            chr_lengths[i] = ceil(theta[i] * num_seg);
            assert(chr_lengths[i] > 0);
            total_seg += chr_lengths[i];
            // cout << alpha[i] << "\t" << theta[i] << "\t" << chr_lengths[i] << endl;
        }
        if(total_seg > num_seg){
            int diff = total_seg - num_seg;
            while(diff != 0){
                for(int i = 0; i < NUM_CHR; i++){
                    if(chr_lengths[i] > 1){
                        chr_lengths[i] -= 1;
                        diff--;
                    }
                    if(diff == 0) break;
                }
            }
        }
        cout << "Number of segments simulated is " << num_seg << endl;
    }else{    // When simulating waiting times, only simulate allele-specific model
      if(model == BOUNDT){
          cout << "When simulating waiting times, only allele-specific model is allowed" << endl;
          model = BOUNDA;
      }
    }

    vector<double> rate_consts = {dup_rate, del_rate, chr_gain, chr_loss, wgd};

    cout << "\nEvolution model simulated: " << model << endl;
    cout << "\tMaximum copy number of a segment is " << cn_max << endl;
    cout << "\tMutation rates:\t" << rate_consts[0] << "\t" << rate_consts[1]  << "\t" << rate_consts[2]  << "\t" << rate_consts[3]  << "\t" << rate_consts[4] << endl;
    cout << "\tSizes of segment duplication/deletion:\t" << mean_dup_size << "\t" << mean_del_size << endl;

    // simulate coalescent tree and apply SVs
    if(mode < 2){
        if(beta > 0){
            cout << "\nSimulating exponential growth" << endl;
        }

        run_simulations(tree_file, mode, method, chr_lengths, num_seg, Ns, Nsims, cn_max, model, cons, Ne, beta, gtime, delta_t, age, rate_consts, dir, prefix, print_allele, print_mut, print_nex);
        // cout << "FINISHED" << endl;
    }else{
        cout << "\nRunning test" << endl;
        run_test(mode, dir, seed);
    }

    gsl_rng_free(r);
}
