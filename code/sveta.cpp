#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <map>
#include <random>
#include <sstream>
#include <ctime>
#include <map>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include "gzstream.h"

#include "sveta.hpp"
#include "evo_tree.hpp"
#include "genome.hpp"
#include "stats.hpp"
#include "utilities.hpp"

#include <boost/program_options.hpp>

using namespace std;

int debug = 0;
// The number of bins for each chromosme. Each bin corresponds to a window of size 500,000 bp.
static const int BIN_SIZE[] = {367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};
vector<int> CHR_BIN_SIZE(BIN_SIZE, BIN_SIZE + sizeof(BIN_SIZE) / sizeof(BIN_SIZE[0]));
const int NODE_MIN_TIME = 1000;


// Find the number of a specific chromosme
int get_ploidy_at_chr(genome& g, int c){
    int ploidy_at_c = 0;
    for(int i=0; i<g.chrs.size(); i++){
        if(i%22 == c && g.chrs[i].size()>0){
            ploidy_at_c += 1;
        }
    }
    return ploidy_at_c;
}

// Find the current maximum copy number for one chromosme
int get_max_cn_chr(genome& g, int c){
    g.calculate_cn();
    vector<segment> segs = g.chrs[c];
    int max_cn = 0;
    for (int i=0; i<segs.size(); i++){
        int cn = g.cn_profile[c%22][g.chrs[c][i].seg_id];
        if(cn > max_cn){
            max_cn = cn;
        }
    }
    if(debug){
        cout << "maximum copy number for chr " << c << " is " << max_cn << endl;
    }
    return max_cn;
}

int get_max_seg_num_chr(genome& g, int c){
    g.calculate_cn();
    map<int, int> seg_cn;     // record seg_id and current copy number
    vector<int> segs;     // record seg_id and current multiplicity
    for (int i=0; i<g.chrs[c].size(); i++){
        int sid = g.chrs[c][i].seg_id;
        int cn = g.cn_profile[c%22][sid];
        seg_cn.insert(pair<int, int>(sid, cn));
        segs.push_back(sid);
    }
    map<int,int> dup;
    for_each(segs.begin(), segs.end(), [&dup]( int val ){ dup[val]++; } );

    int num_mseg = 0;
    for( auto p : dup ) {
        int total_cn = p.second + seg_cn[p.first];
        // if(debug){
        //     cout << "Number of segments: " << p.first << ' ' << p.second << ' ' << total_cn << endl;
        // }
        if(total_cn > num_mseg){
            num_mseg = total_cn;
        }
    }
    return num_mseg;
}

// Find the current maximum copy number for one segment on a chromosme
int get_max_cn_chr_seg(genome& g, int c, int loc, int len){
    if(g.chrs[c].size() <= 0){
        return 0;
    }
    g.calculate_cn();
    int max_state = g.cn_profile[c%22][g.chrs[c][loc].seg_id];
    // Check the state of all the regions
    for (int i = loc + 1; i < loc + len + 1; i++){
        int state = g.cn_profile[c%22][g.chrs[c][i].seg_id];
        // cout << "copy number for " << i << " is " << state << endl;
        if(state > max_state){
            max_state = state;
        }
    }
    if(debug){
        cout << "maximum copy number for chr " << c << ", loc " << loc << ", length " << len << " is " << max_state << endl;
    }
    return max_state;
}


int get_max_seg_num_seg(genome& g, int c, int loc, int len){
    g.calculate_cn();
    // Some segments may not have the largest copy number for now. But there current copy number plus the multiplicity will make the new copy number larger than specified.
    map<int, int> seg_cn;     // record seg_id and current copy number
    vector<int> segs;     // record seg_id and current multiplicity
    for (int i=loc; i < loc + len + 1; i++){
        int sid = g.chrs[c][i].seg_id;
        int cn = g.cn_profile[c%22][sid];
        seg_cn.insert(pair<int, int>(sid, cn));
        segs.push_back(sid);
    }
    map<int,int> dup;
    for_each(segs.begin(), segs.end(), [&dup]( int val ){ dup[val]++; } );

    int num_mseg = 0;
    for( auto p : dup ) {
        int total_cn = p.second + seg_cn[p.first];
        // if(debug){
        //     cout << "Number of segments: " << p.first << ' ' << p.second << ' ' << total_cn << endl;
        // }
        if(total_cn > num_mseg){
            num_mseg = total_cn;
        }
    }
    return num_mseg;
}


// Find the current maximum copy number for the whole genome
int get_max_cn_genome(genome& g){
    g.calculate_cn();
    int max_cn = 0;
    for(int i=0; i<g.cn_profile.size(); i++){
        int cn = get_max_cn_chr(g, i);
        if(cn > max_cn){
            max_cn = cn;
        }
    }
    if(debug){
        cout << "maximum copy number for the whole genome is " << max_cn << endl;
    }
    return max_cn;
}


// rates: The rate of all events (segment duplication, segment deletion)
double get_site_rates(vector<double>& rates, const vector<double>& rate_constants, int model, int state, int cn_max){
    double rate;

    if(model == 0){ // JC69
        // assume duplication rate equals deletion rate
        assert(rate_constants[0] == rate_constants[1]);
        for(int i=0; i<2; ++i){
          rates.push_back(rate_constants[i]);
        }
        rate = 4 * rate_constants[0] / 5;
    }

    if(model == 1){ // bounded
        double dup_rate = rate_constants[0];
        double del_rate = rate_constants[1];

        if(state == 0){
            rate =  0;
        }
        else if(state == cn_max){
            rate =  2 * cn_max * del_rate;
        }
        else{
            rate =  2 * state * dup_rate + 2 * state * del_rate;
        }
    }
    return rate;
}


// rates: The total rate at each site in the genome, used for randomly selecting a site
double get_total_rates(genome& g, vector<double>& rates, const vector<double>& rate_constants, int model, int cn_max){
    assert(rate_constants.size() == 5);
    double chr_gain_rate = rate_constants[2];
    double chr_loss_rate = rate_constants[3];
    double wgd_rate = rate_constants[4];

    g.calculate_cn(); // Find the current state of the genome before computing rates (affected by copy number)
    vector<double> rates_chrs;
    // Tranverse all the sites in the genome to get their rates
    for(int i=0; i < g.cn_profile.size(); ++i){   // For each chromosome
        vector<double> rates_chr_sites;
        for(int j=0; j < g.cn_profile[i].size(); ++j){    // For each segment in the chromosome
            // Get the copy number of this segment
            int state = g.cn_profile[i][j];
            vector<double> site_rates;
            double rate = get_site_rates(site_rates, rate_constants, model, state, cn_max);
            // cout << "Rate on chr " << i << " seg " << j << " state " << state << " is " << rate << endl;
            rates.push_back(rate);
            rates_chr_sites.push_back(rate);
        }
        double rate_c = accumulate(rates_chr_sites.begin(), rates_chr_sites.end(), 0.0);
        double max_cn_c = get_max_cn_chr(g, i);
        double max_snum = get_max_seg_num_chr(g, i);
        if(max_cn_c > 0){
            rate_c += chr_loss_rate;

            if(max_snum < cn_max){
                rate_c += chr_gain_rate;
            }
        }
        if(debug){
            cout << "Rate on chr " << i << " is " << rate_c << endl;
        }
        rates_chrs.push_back(rate_c);
    }
    // Get the total rate of substitution (the sum of rates accross sites)
    double rate = accumulate(rates_chrs.begin(), rates_chrs.end(), 0.0);
    int max_cn = get_max_cn_genome(g);
    if(max_cn > 0 && 2 * max_cn <= cn_max){
        rate += wgd_rate;
    }
    if(debug){
        cout << "The total rate at all sites in the genome: " << rate << endl;
    }
    return rate;
}


// Note: the sites are counted for all the haplotypes. The positions on different haplotypes are counted as one site.
void site2chr(int site, int& chr, int& seg, const vector<int>& chr_lengths){
    if(debug) cout << "Site to convert is " << site << endl;
    if(site >= 0 && site < chr_lengths[0]){
        chr = 0;
        seg = site;
    }
    else{
        if(debug) cout << "There are " << chr_lengths.size() << " chromosomes" << endl;
        for(int i=0; i<chr_lengths.size(); i++){
            int sum1 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 1, 0);
            int sum2 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 2, 0);
            // cout << "Size for chr " << i+1 << " is " << chr_lengths[i] << endl;
            // cout << "sum until chromosome " << i + 1 << " is " << sum1 << endl;
            // cout << "sum until chromosome " << i + 2 << " is " << sum2 << endl;
            if(site > sum1 && site < sum2){
                chr = i + 1;
                seg = site - sum1;
                break;
            }
        }
    }
}

// Find the number of chromosmes with segments
int get_num_chr(genome& g){
    int n = 0;
    for(int i=0; i<g.chrs.size(); i++){
        if(g.chrs[i].size()>0){
            n++;
        }
    }
    return n;
}

// Find the ID of chromosmes with segments
void get_available_chr(genome& g, vector<int>& available_chrs){
    available_chrs.clear();
    for(int i=0; i<g.chrs.size(); i++){
        // cout << "size of chr " << i << " is " << g.chrs[i].size() << endl;
        if(g.chrs[i].size()>0){
            available_chrs.push_back(i);
        }
    }
}

int generate_duplication(genome& g, int c, int loc, int model, int cn_max, int debug){
    if( g.chrs[c].size() <= 0 || loc >= g.chrs[c].size()){
        return 0;
    }
    if(debug){
      cout << "GENOME BEFORE duplication " << endl;
      // g.print();
      g.print_cn();
    }
    // mean variant size in bins
    double mdup = g.mean_dup_size;
    int len = 0;
    if(mdup > 1)    len = gsl_ran_exponential(r, mdup);

    if(debug){
      cout << "before dup, Chr " << c << " has " << g.chrs[c].size() << " segments" << endl;
      // for(int i=0; i<g.chrs[c].size(); i++){
      //     cout << "\t" << g.chrs[c][i].seg_id;
      // }
      // cout << endl;
      cout << "dup len:" << len << endl;
    }

    if( loc + len >= g.chrs[c].size() ){
      return 0;
    }

    if(debug) cout << "\tSV: duplicating segment, chr, start, len: " << c << "\t" << loc << "\t" << len+1 << endl;
    //cout << "SV: insert: " << loc+len+1 << "\t" << loc << "\t" << loc+len+1 << endl;
    // Find the number of copies
    int state = g.cn_profile[c%22][g.chrs[c][loc].seg_id];
    if(debug){
        cout << "Previous state: ";
        for(int i=0; i<len+1; i++){
            cout << "\t " << g.cn_profile[c%22][g.chrs[c][loc + i].seg_id];
        }
        cout << endl;
    }

    int li = loc+len+1; //assume tandom duplication

    if(model == 0){
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
        if(debug){
            cout << "copy start " << g.chrs[c][loc].seg_id << endl;
            cout << "copy end " << g.chrs[c][loc+len].seg_id << endl;
        }

        if(pstate-state <= 0){
            return 0;
        }
        // assume tandom duplication
        for(int j=0; j < pstate-state; j++){
            g.chrs[c].insert(g.chrs[c].begin()+loc+(len+1)*(j+1), g.chrs[c].begin()+loc, g.chrs[c].begin()+loc+len+1);

            if(debug){
                cout << j+1 << "th dup" << endl;
                int it;
                it = g.chrs[c][loc+(len+1)*(j+1)].seg_id;
                cout << "dup start " << it << endl;
                it = g.chrs[c][loc+(len+1)*(j+1)+len].seg_id;
                cout << "dup end " << it << endl;
            }
        }
    }

    if(model == 1){
        // int max_state = get_max_cn_chr_seg(g, c, loc, len);
        int max_snum = get_max_seg_num_seg(g, c, loc, len);
        // cout << "cn_max is " << cn_max << endl;
        if(max_snum > cn_max){
            return 0;
        }
        // tandom duplication
        // double u = runiform(r, 0, 1);
        // if(u<0.5){
        //     // randomly select a position
        //     li = gsl_rng_uniform_int(r, g.chrs[c].size());
        // }
        g.chrs[c].insert(g.chrs[c].begin()+li, g.chrs[c].begin()+loc, g.chrs[c].begin()+loc+len+1);
        g.calculate_cn();   // compute copy number after each mutation event

        if(debug){
            cout << "End position " << loc+len+1 << endl;
            cout << "Insertion position " << li << endl;
            cout << "Current state: ";
            for(int i=0; i<len+1; i++){
                cout << "\t " << g.cn_profile[c%22][g.chrs[c][li + i].seg_id];
            }
            cout << endl;
            set<double> segs;
            for(int i=0; i<g.chrs[c].size(); i++){
                segs.insert(g.chrs[c][i].seg_id);
            }
            cout << "After duplication, Chr " << c << " has " << g.chrs[c].size() << " segments (unique: " << segs.size() << ")" << endl;
        }
    }
    // }
    if(debug){
        cout << "GENOME AFTER duplication " << endl;
        // g.print();
        g.print_cn();
    }
    return 1;
}


int generate_deletion(genome& g, int c, int loc, int model, int cn_max, int debug){
    if( g.chrs[c].size() <= 0 || loc >= g.chrs[c].size()){
        return 0;
    }
    if(debug){
      cout << "GENOME BEFORE deletion " << endl;
      // g.print();
      g.print_cn();
    }

    double mdel = g.mean_del_size;
    int len = 0;
    if(mdel > 1)   len = (int) gsl_ran_exponential(r, mdel);

    if(debug){
      cout << "before del, Chr " << c << " has " << g.chrs[c].size() << " segments" << endl;
      // for(int i=0; i<g.chrs[c].size(); i++){
      //     cout << "\t" << g.chrs[c][i].seg_id;
      // }
      // cout << endl;
      cout << "del len:" << len << endl;
    }

    if( len + loc >= g.chrs[c].size() ){
        return 0;
    }

    if(model == 0){
        int state = g.cn_profile[c%22][g.chrs[c][loc].seg_id];
        // Find the number of segments in other haplotypes
        int hap_state = 0;
        int ploidy_at_c = get_ploidy_at_chr(g, c);
        // cout << "Ploidy at " << c << " is " << ploidy_at_c << endl;
        for(int j=0; j<ploidy_at_c; j++){
            if(j == c/22) continue;
            int hap = c % 22 + j * 22;
            // cout << "Alternative haplotype " << hap << endl;
            for(int i=0; i<g.chrs[hap].size();i++){
                if(g.chrs[hap][i].seg_id == loc)
                {
                    hap_state++;
                }
            }
        }
        if(debug){
            cout << "Copies in other haplotypes: " << hap_state << endl;
            cout << "\tSV: deleting segment, chrs, start, len: " << c << "\t" << loc << "\t" << len+1 << endl;
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
            // cout << "Impossible to do deletion on Chr " << c << endl;
            return 0;
        }
        else{
            if(ssize > 1){
                int k = gsl_rng_uniform_int(r, ssize);
                pstate = possible_states[k];
            }
            else{
                pstate = possible_states[0];
            }
            // Ensure there are "state" copy of the region
            int start = g.chrs[c][loc].seg_id;
            int end = g.chrs[c][loc+len].seg_id;
            if(debug){
                cout << "del start " << start << endl;
                cout << "del end " << end << endl;
            }

            if(pstate < 0){
                return 0;
            }

            for(int j=0; j <= pstate; j++){
                // erase all the segment with seg_id in the specified region
                for(int i=0; i<g.chrs[c].size();i++){
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
    } // model == 0

    if(model == 1){
        // Find the minimum copy number of the region to delete
        int min_state = g.cn_profile[c%22][g.chrs[c][loc].seg_id];
        for (int i = loc; i < loc + len + 1; i++){
            int state = g.cn_profile[c%22][g.chrs[c][i].seg_id];
            if(state < min_state){
                min_state = state;
            }
        }
        if(min_state <= 0){
            return 0;
        }
        g.chrs[c].erase(g.chrs[c].begin()+loc, g.chrs[c].begin()+loc+len+1);
    }

    g.calculate_cn();   // compute copy number after each mutation event

    if(debug){
        cout << "Current state: " << g.cn_profile[c%22][g.chrs[c][loc].seg_id] << endl;
        set<double> segs;
        for(int i=0; i<g.chrs[c].size(); i++){
            segs.insert(g.chrs[c][i].seg_id);
        }
        cout << "After deletion, Chr " << c << " has " << g.chrs[c].size() << " segments (unique: " << segs.size() << ")" << endl;

        cout << "GENOME AFTER deletion " << endl;
        // g.print();
        g.print_cn();
    }
    return 1;
}


int generate_chr_gain(genome& g, int cn_max, int debug){
    vector<int> available_chrs;
    get_available_chr(g, available_chrs);
    if(available_chrs.size()<=0){
        // cout << "No available chromosmes!" << endl;
        return 0;
    }
    int ci = gsl_rng_uniform_int(r, available_chrs.size());
    int c = available_chrs[ci];
    // int max_cn_c = get_max_cn_chr(g, c);
    // Some segments with maximum copy number may have several copies and hence have copy number larger than specified value when having a chromosme gain
    int max_snum = get_max_seg_num_chr(g, c);
    if( g.chrs[c].size() <= 0 || max_snum > cn_max){
        return 0;
    }

    // Make sure that the maximum copy number is smaller than the specified value
    int orig_num_chr = get_num_chr(g);
    int orig_size = g.chrs.size();
    // cout << "Size of chrs " << orig_size << endl;
    // assert(orig_size % 22 == 0);
    int new_chr = (orig_size / 22 ) * 22 + c%22; // If c is not in chrs
    for(int i=0; i<orig_size%22; i++){
      if(g.chrs[c%22 + i*22].size() <= 0){
          new_chr = c%22 + i*22;
          break;
      }
    }
    // cout << "ID of new chr " << new_chr << endl;
    for(int i=orig_size; i<=new_chr; i++){ // Fill the position for other chromosmes before the one to duplicate so that the chromosme ID follows the pattern
      vector<segment> e;
      g.chrs.insert(g.chrs.end(), e);
      // cout << "Size of chr " <<  i << " is " << g.chrs[i].size() << endl;
    }
    g.chrs[new_chr].insert(g.chrs[new_chr].begin(), g.chrs[c].begin(), g.chrs[c].end());
    g.calculate_cn();   // compute copy number after each mutation event

    if(debug){
      cout << "Chromosome gain in " << c << endl;
      cout << "ID for gained Chromosome " << new_chr << endl;
      cout << "There are " << orig_num_chr << " chromosomes before chr gain" << endl;
      cout << "There are " << get_num_chr(g) << " chromosomes now" << endl;
      g.print_cn();
    }
    return 1;
}


int generate_chr_loss(genome& g, int debug) {
    vector<int> available_chrs;
    get_available_chr(g, available_chrs);
    if(available_chrs.size()<=0){
        // cout << "No available chromosmes!" << endl;
        return 0;
    }
    int ci = gsl_rng_uniform_int(r, available_chrs.size());
    int c = available_chrs[ci];
    // int c = gsl_rng_uniform_int(r, g.chrs.size() );
    int orig_num_chr = get_num_chr(g);
    // Delete the segments in the chromosme, but keep the chromosme ID so that it can be remapped by 22
    g.chrs[c].erase(g.chrs[c].begin(), g.chrs[c].end());
    g.calculate_cn();   // compute copy number after each mutation event

    if(debug){
        cout << "Chromosome loss in " << c << endl;
        cout << "There are " << orig_num_chr << " chromosomes before chr loss" << endl;
        cout << "There are " << get_num_chr(g) << " chromosomes now" << endl;
        g.print_cn();
    }
    return 1;
}


int generate_wgd(genome& g, int cn_max, int debug) {
    int max_cn_g = get_max_cn_genome(g);
    if( g.chrs.size() <= 0 || 2 * max_cn_g > cn_max){
        return 0;
    }

    int orig_num_chr = get_num_chr(g);
    int orig_size = g.chrs.size();
    g.chrs.insert(g.chrs.end(), g.chrs.begin(), g.chrs.end());
    // replicate the segments for each chromosme
    for(int i=0; i<orig_size; i++){
         g.chrs[i + orig_size].insert(g.chrs[i + orig_size].begin(), g.chrs[i].begin(), g.chrs[i].end());
    }
    g.calculate_cn();   // compute copy number after each mutation event

    if(debug){
        cout << "Whole genome doubling" << endl;
        cout << "There are " << orig_num_chr << " chromosomes before WGD" << endl;
        cout << "There are " << get_num_chr(g) << " chromosomes now" << endl;
        g.print_cn();
        // g.print();
    }
    return 1;
}


// Simulate mutations under different models
vector<mutation> generate_mutation_by_model(genome& g, const int& edge_id, const double& blength, const double& node_time, const vector<int>& chr_lengths, const vector<double>& rate_constants, int model, int cn_max){
    vector<mutation> ret;
    if(debug) cout << "\tGenerate_mutations, blength:" << "\t" << blength << endl;

    vector<double> site_rates;
    if(debug){
      cout << "\tComputing total mutation rate on the genome " << endl;
    }
    double rate = get_total_rates(g, site_rates, rate_constants, model, cn_max);
    // cout << "Total mutation rate on the genome " << rate << endl;

    // Simulate the waiting times of a Markov chain
    double time = 0.0;
    while( time < blength ){
        double tevent = gsl_ran_exponential(r, 1/rate);
        time += tevent;
        // choose type of event
        int e = rchoose(r, rate_constants);     // assign the event to a site with probabilities proportional to the rates at the site

        // Randomly select a site for duplication or deletion. Convert the site position back to chromosome ID and segment ID
        int site = 0;
        int c = 0;
        int loc = 0;
        int ploidy_at_c = 0;
        if(e == 0 || e == 1){
            site = rchoose(r, site_rates);
            site2chr(site, c, loc, chr_lengths);
            // Randomly choose a haplotype
            double x = runiform(r, 0, 1);
            // Find the number of possible haplotype of this chromosme
            ploidy_at_c = get_ploidy_at_chr(g, c);
            if(ploidy_at_c == 0) {  // The chromosme has been completely lost
                continue;
            }
            int haplotype = myrng(ploidy_at_c);
            c += 22 * haplotype;
            if(debug){
                cout << "There are " << g.chrs.size() << " chromosome IDs now" << endl;
                cout << "site " << site << " is at chr " << c << " pos " << loc << " haplotype " << haplotype << " ploidy " << ploidy_at_c << endl;
            }
        }

        int res = 0;
        if(e == 0){   //duplications
            res = generate_duplication(g, c, loc, model, cn_max, debug);
        }
        else if( e == 1 ){   //deletions
            // make sure a deletion somewhere is possible
            res = generate_deletion(g, c, loc, model, cn_max, debug);
        }
        else if( e == 2 ){   //chr gain
            res = generate_chr_gain(g, cn_max, debug);
        }
        else if( e == 3 ){   //chr loss
            res = generate_chr_loss(g, debug);
        }
        else if( e == 4 ){   // whole genome duplication
            res = generate_wgd(g, cn_max, debug);
        }
        else{
            cerr << "Unknown mutation type " << e << endl;
            exit(1);
        }

        if(res > 0){
            mutation mut(edge_id, e, time/blength, node_time+time);
            ret.push_back(mut);
            g.nmuts[e]++;
            g.mutations.push_back(mut);
            // time += tevent;
        }else{
             if(debug) cout << "\tMutation failed:" << endl;
             continue;
        }

        if(debug){
            cout << "mut times, tevent, total time, time/branch len, event, chr, loc\t" << tevent << "\t" << time << "\t" << blength << "\t" << e << "\t" << c << "\t" << loc << endl;
        }
    }

    return ret;
}


vector<mutation> generate_mutation_times(const int& edge_id, const double& blength, const double& node_time, const vector<double>& rate_constants ){
  bool print = false;

  if(print) cout << "\tgenerate_mutations, blength:" << "\t" << blength << endl;

  vector<mutation> ret;

  // model this as a poisson process
  double time = 0.0;
  while( time < blength ){
    vector<double> rates;
    for(int i=0; i<rate_constants.size(); ++i){
      rates.push_back(rate_constants[i]/2.0);        // in coalescent Nmut ~ Pois( theta*l/2 )
      //cout << "##:\t" << rates[i] << endl;
    }

    double rate = accumulate(rates.begin(), rates.end(), 0.0);
    double tevent = gsl_ran_exponential (r, 1/rate);
    time += tevent;

    // choose type of event
    int e = rchoose(r,rates);
    if(print) cout << "mut times, tevent, total time, time/branch len, event\t" << tevent << "\t" << time << "\t" << blength << "\t" << e << endl;

    ret.push_back( mutation( edge_id, e, time/blength, node_time+time ) );
  }
  ret.pop_back();

  return ret;
}

void apply_mutations( const vector<mutation>& muts, genome& g, bool print = false ){
  //print = true;

  // mean variant size in bins
  double mdup = g.mean_dup_size;
  double mdel = g.mean_del_size;

  for(int i=0; i<muts.size(); ++i){

    g.nmuts[ muts[i].type ]++;
    g.mutations.push_back( muts[i] );

    if(muts[i].type == 0){   //duplications
      if(print == true){
        cout << "GENOME BEFORE duplication " << endl;
        //g.print();
        g.print_cn();
      }

      // make sure a duplication somewhere is possible
      bool possible = false;
      for(int j=0; j < g.chrs.size(); ++j){
          if( g.chrs[j].size() > 1 ) possible = true;
    }

    if( possible == true ){
        bool done = false;
        while( done == false ){
          int c = gsl_rng_uniform_int(r, g.chrs.size() );
          int len = gsl_ran_exponential(r, mdup);
          //int len =  1 + gsl_ran_poisson (r, ldup);
          //cout << "dup len:" << len << endl;

          if( len < g.chrs[c].size() ){
            // choose a random location up to size-len
            int loc = gsl_rng_uniform_int(r,g.chrs[c].size()-len);

            if(print == true) cout << "\tSV: duplicating segment, chr, start, len: " << c << "\t" << loc << "\t" << len+1 << endl;
            //cout << "SV: insert: " << loc+len+1 << "\t" << loc << "\t" << loc+len+1 << endl;
            g.chrs[c].insert(g.chrs[c].begin()+loc+len+1, g.chrs[c].begin()+loc, g.chrs[c].begin()+loc+len+1);
            done = true;
          }
        }
     }else{
         if(print == true) cout << "\tSV: duplication failed:" << endl;
     }

      if(print == true){
        cout << "GENOME AFTER duplication " << endl;
        //g.print();
        g.print_cn();
      }
    }
    else if( muts[i].type == 1){   //deletions
      if(print == true){
        cout << "GENOME BEFORE deletion " << endl;
        //g.print();
        g.print_cn();
      }

      // make sure a deletion somewhere is possible
      bool possible = false;
      for(int j=0; j < g.chrs.size(); ++j){
          if( g.chrs[j].size() > 1 ) possible = true;
      }

      if( possible == true ){
        bool done = false;
        while( done == false ){
          int c = gsl_rng_uniform_int(r, g.chrs.size() );
          int len = (int) gsl_ran_exponential(r, mdel);
          //int len =  1 + gsl_ran_poisson (r, ldel);
          //cout << "del len:" << len << endl;

          if( len < g.chrs[c].size() ){
            int loc = gsl_rng_uniform_int(r,g.chrs[c].size()-len);
            if(print == true) cout << "\tSV: deleting segment, chrs, start, len: " << c << "\t" << loc << "\t" << len+1 << endl;
            g.chrs[c].erase(g.chrs[c].begin()+loc, g.chrs[c].begin()+loc+len+1);
            done = true;
          }
        }
      }else{
          if(print == true) cout << "\tSV: deletion failed:" << endl;
      }

      if(print == true){
        cout << "GENOME AFTER deletion " << endl;
        //g.print();
        g.print_cn();
      }
    }
    else if( muts[i].type == 2){   //chr loss
      if( g.chrs.size() > 0 ){
        int c = gsl_rng_uniform_int(r, g.chrs.size() );
        // g.chrs.erase(g.chrs.begin()+c);
        g.chrs[c].erase(g.chrs[c].begin(), g.chrs[c].end());
      }

    }
    else if( muts[i].type == 3){   //chr gain
      if( g.chrs.size() > 0 ){
        int c = gsl_rng_uniform_int(r, g.chrs.size() );
        int orig_num_chr = g.chrs.size();
        g.chrs.insert(g.chrs.end(), g.chrs.begin()+c, g.chrs.begin()+c+1);
        g.chrs[orig_num_chr].insert(g.chrs[orig_num_chr].begin(), g.chrs[c].begin(), g.chrs[c].end());
      }
    }
    else if( muts[i].type == 4){   // whole genome duplication
      if( g.chrs.size() > 0 ){
          int orig_num_chr = g.chrs.size();
          g.chrs.insert(g.chrs.end(), g.chrs.begin(), g.chrs.end());
          // replicate the segments for each chromosme
          for(int i=0; i<orig_num_chr; i++){
               g.chrs[i + orig_num_chr].insert(g.chrs[i + orig_num_chr].begin(), g.chrs[i].begin(), g.chrs[i].end());
          }
      }
    }
    else{
      cerr << "Unknown mutation type" << endl;
      exit(1);
    }
  }
}

void traverse_tree_mutating(const int& node_id, const evo_tree& tree, const vector<int>& chr_lengths, const vector<double>& rate_constants,
        map<int, vector<mutation>>& all_muts, vector<genome>& genomes, int model, int cn_max){
  //cout << "\ttraverse_tree: " << node_id+1 << endl;
  if( !tree.nodes[node_id].isRoot ){
    // copy the parent
    genomes[ node_id ] = genomes[ tree.nodes[node_id].parent ];
    genomes[ node_id ].node_id = node_id;

    // apply the mutations from parent -> daughter
    //cout << "MUTATING genome: " << tree.nodes[node_id].parent+1 << " -> " << node_id+1 << "\t edge id: " << tree.nodes[node_id].e_in+1 << endl;
    int edge_id = tree.nodes[node_id].e_in;
    vector<mutation> muts;
    if( tree.edges[edge_id].length > 0 ){
        if(debug){
            cout << "Generating mutations for node " << node_id << endl;
        }
        if(model == 2){
            muts = generate_mutation_times(edge_id, tree.edges[edge_id].length, tree.node_times[ tree.nodes[node_id].parent ], rate_constants);
            apply_mutations( muts, genomes[ node_id ] );
        }
        else{
            muts = generate_mutation_by_model(genomes[node_id], edge_id, tree.edges[edge_id].length, tree.node_times[ tree.nodes[node_id].parent], chr_lengths, rate_constants, model, cn_max);
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
  }

  if( tree.nodes[ node_id ].daughters.size() == 0 ){
    // we are done
    return;
  }else{
    traverse_tree_mutating(tree.nodes[ node_id ].daughters[0], tree, chr_lengths, rate_constants, all_muts, genomes, model, cn_max);
    traverse_tree_mutating(tree.nodes[ node_id ].daughters[1], tree, chr_lengths, rate_constants, all_muts, genomes, model, cn_max);
  }
  return;
}

void simulate_samples(vector<genome>& genomes, map<int,vector<mutation> >& muts, const evo_tree& tree, genome& germline, const vector<int>& chr_lengths, const vector<double>& rate_constants, int model, int cn_max, bool print = true){
  // assign the germline to the root of the tree
  germline.node_id = tree.root_node_id;

  for(int i=0; i<(tree.nnode + tree.nleaf); ++i){
    if(i == tree.root_node_id){
      genomes.push_back( germline );
    }
    else{
      genomes.push_back( genome() );
    }
  }

  // move through the evolutionary tree mutating genomes
  traverse_tree_mutating(tree.root_node_id, tree, chr_lengths, rate_constants, muts, genomes, model, cn_max);

  // final samples returned, print out leaf nodes
  if(print == true){
    cout << "MUTATIONS:" << endl;
    for(map<int, vector<mutation> >::iterator it = muts.begin(); it != muts.end(); it++){
      cout << "EDGE, id: " << it->first+1
       << "\t" << tree.edges[it->first].start+1 << " -> " << tree.edges[it->first].end+1
       << "\t" << it->second.size() << endl;
      vector<mutation> v = it->second;
      for(int i=0; i<v.size(); ++i){
          v[i].print();
      }
    }
    cout << endl;
    cout << "LEAF GENOMES:" << endl;
    for(int i=0; i<tree.nleaf; ++i){
      //genomes[i].print();
      genomes[i].print_muts();
      //genomes[i].print_cn();
      tree.print_ancestral_edges( genomes[i].node_id );
      cout << endl;
    }
  }
}


//void run_sample_set(int Ns, double* prc, double* pvs, double* ptree, double* pl, int* ret){
void run_sample_set(int Ns, double* prc, double* pvs, int* ret){
  static const int arr[] = {367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};
  vector<int> chr_lengths (arr, arr + sizeof(arr) / sizeof(arr[0]) );

  int model = 2;
  int cn_max = 4;

  vector<double> rate_consts (prc, prc + sizeof(prc) / sizeof(prc[0]) );
  genome germline(chr_lengths,2);
  germline.mean_dup_size = pvs[0];
  germline.mean_del_size = pvs[1];

  vector<int> edges;

  vector<double> lengths;
  vector<double> epoch_times;
  vector<double> node_times;
  stringstream sstm;
  vector<genome> results;
  map<int, vector<mutation> > muts;
  //cout << "\n\n###### New sample collection ######" << endl;
  //cout << "###### Ns+1= " << Ns+1 << endl;

  generate_coal_tree(Ns, edges, lengths, epoch_times, node_times);
  evo_tree test_tree(Ns+1, edges, lengths);

  //for(int i=0; i<6; ++i) epars.push_back( ptree[i] );
  //for(int i=0; i<8; ++i) lengths.push_back( pl[i] );
  //evo_tree test_tree = construct_tree(Ns, epars, lengths, node_times );

  test_tree.node_times = node_times;
  simulate_samples(results, muts, test_tree, germline, chr_lengths, rate_consts, model, cn_max, false);

  int nbins = 4401;
  for(int i=0; i<(test_tree.nleaf-1); ++i){
    vector<int> cns = results[i].get_cn_vector();
    for(int j=0; j<nbins; ++j){
      ret[nbins*i + j] = cns[j];
    }
  }

  if(0){
    sstm << "test-data-cn.txt.gz";
    ogzstream out_cn(sstm.str().c_str());
    for(int j=0; j<test_tree.nleaf; ++j){
      results[j].write(out_cn);
    }
    out_cn.close();
    sstm.str("");
  }
  return;
}



void run_test(int mode, string dir){
    //create test tree
    if(mode == 2){
      //static const int arr1[] = {7,8, 6,7, 8,1, 8,2, 7,9, 9,3, 9,4, 6,5 };
      static const int arr1[] = {6,7, 5,6, 7,0, 7,1, 6,8, 8,2, 8,3, 5,4 };
      vector<int> e (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
      static const double arr2[] = {5.1,6.3,10.2,9.5,5.2,3.2,5.4,0};
      vector<double> l (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
      evo_tree test_tree(5, e, l);
      test_tree.print();
    }

    // test out the genome classes and functionality
    if(mode == 3){
      if(0){
        vector<mutation> dels;
        dels.push_back( mutation( 0, 1, 0.5, 0 ) );
        dels.push_back( mutation( 0, 1, 0.5, 0 ) );
        dels.push_back( mutation( 0, 1, 0.5, 0 ) );
        dels.push_back( mutation( 0, 1, 0.5, 0 ) );

        vector<mutation> dups;
        dups.push_back( mutation( 0, 0, 0.5, 0 ) );
        dups.push_back( mutation( 0, 0, 0.5, 0 ) );
        dups.push_back( mutation( 0, 0, 0.5, 0 ) );
        dups.push_back( mutation( 0, 0, 0.5, 0 ) );

        vector<mutation> wgds;
        wgds.push_back( mutation( 0, 4, 0.5, 0 ) );

        genome g1(2,10);
        g1.node_id = 0;
        g1.print();
        g1.print_cn();
        apply_mutations( dups, g1 );
        apply_mutations( dels, g1 );
        apply_mutations( wgds, g1 );
        g1.print();
        g1.print_cn();
      }

      if(1){
        static const int arr3[] = {10,4,3};
        vector<int> chr_lengths (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]) );

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
      double prc[] = {0.1, 0.1, 0.1, 0.1, 0.05};
      double pvs[] = {30.0, 30.0};

      for(int i=0; i<10; ++i){
        int* ret = new int[Ns*4401];
        run_sample_set(Ns, prc, pvs, &(ret[0]) );

        sstm << dir << "sim-data-" << i+1 << "-cn.txt";
        ofstream out_cn(sstm.str());
        for(int i=0; i<(Ns*4401); ++i){
            out_cn << ret[i] << endl;
        }
        out_cn.close();
        sstm.str("");
        delete[] ret;
      }
    }
}


// randomly assign leaf edges to time points t0, t1, t2, t3, ...
void assign_tip_times(double delta_t, int Ns, vector<double>& tobs, const vector<int>& edges, vector<double>& lengths){
    tobs.clear(); // clear the sample time differences
    if(delta_t > 0){
      cout << "assigning temporal structure" << endl;
      bool assign0 = false;
      int bcount = 0;
      int num_diff = Ns; // maximum of tips with different times
      for(int l=1; l<edges.size(); l=l+2){
        if(edges[l] < Ns){
          int ind = 0;
          if(assign0 == false){
            ind = 0;
            assign0 = true;
          }
          else{
            ind = gsl_rng_uniform_int(r, num_diff);
          }
          double stime = ind * delta_t;
          cout << "\t sample / sampling time point :" << edges[l]+1 << "\t" << stime << endl;
          lengths[bcount] = lengths[bcount] + stime;
          // Store the sampling time for time scaling by the age of a patient
          tobs.push_back(stime);
        }
        bcount++;
      }
    }else{
      for(int i=0; i<Ns; i++){
          tobs.push_back(0);
      }
    }
}


void print_simulations(int mode, vector<genome>& results, map<int, vector<mutation>>& muts, evo_tree& test_tree, string dir, string prefix, int age, int print_allele, int print_mut, int print_nex){
    stringstream sstm;

    sstm << dir << prefix << "-cn.txt.gz";
    ogzstream out_cn(sstm.str().c_str());
    for(int j=0; j<test_tree.nleaf; ++j){
        results[j].write(out_cn);
    }
    out_cn.close();
    int num_total_bins = 0;
    int num_invar_bins = 0;
    cout << "reading data and calculating CNA regions" << endl;
    vector<vector<vector<int>>> s_info = read_cn(sstm.str(), test_tree.nleaf-1, num_total_bins);
    // We now need to convert runs of variable bins into segments of constant cn values, grouped by chromosme
    vector<vector<int>> segs;
    if(mode == 0){
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
        for(int j=0; j<test_tree.nleaf; ++j){
          // cout << "Chr size " << results[j].chrs.size() << endl;
          results[j].write_allele_cn(out_allele_cn);
        }
        out_allele_cn.close();
        sstm.str("");
    }

    // Output the copy numbers for internal nodes
    sstm << dir << prefix << "-inodes-cn.txt.gz";
    //ofstream out_cn(sstm.str());
    ogzstream out_cn_inodes(sstm.str().c_str());
    for(int j=test_tree.nleaf; j<test_tree.ntotn; ++j){
      results[j].write(out_cn_inodes);
    }
    out_cn_inodes.close();
    sstm.str("");

    if(print_allele){
        sstm << dir << prefix << "-inodes-allele-cn.txt.gz";
        ogzstream out_allele_cn_inodes(sstm.str().c_str());
        for(int j=test_tree.nleaf; j<test_tree.ntotn; ++j){
          // cout << "Chr size " << results[j].chrs.size() << endl;
          results[j].write_allele_cn(out_allele_cn_inodes);
        }
        out_allele_cn_inodes.close();
        sstm.str("");
    }


    sstm << dir << prefix << "-info.txt";
    ofstream out_info(sstm.str());
    out_info << "NODE TIMES:";
    out_info << "\tnid\ttime" << endl;
    for(int j=0; j< test_tree.node_times.size(); ++j){
        out_info << "\t" << j+1 << "\t" << test_tree.node_times[j] << endl;
    }
    out_info << endl;
    out_info << "SEGMENTS: " << seg_size << endl;
    out_info << endl;
    for(int i=0; i<test_tree.nleaf; ++i){
        test_tree.print_ancestral_edges( results[i].node_id, out_info );
    }
    out_info << endl;
    for(int i=0; i<test_tree.nleaf; ++i){
        results[i].print_muts(out_info);
    }
    out_info.close();
    sstm.str("");

    sstm << dir << prefix << "-tree.txt";
    vector<int> nmuts;
    ofstream out_tree(sstm.str());
    //test_tree.write(out_tree);
    out_tree << "start\tend\tlength\teid\tnmut" << endl;
    for(int i=0; i<test_tree.nedge; ++i){
        out_tree << test_tree.edges[i].start+1 << "\t" << test_tree.edges[i].end+1 << "\t" << test_tree.edges[i].length << "\t" << test_tree.edges[i].id+1 << "\t" << muts[test_tree.edges[i].id].size() << endl;
      nmuts.push_back(muts[test_tree.edges[i].id].size());
    }
    out_tree.close();
    sstm.str("");

    if(print_nex){
        sstm << dir << prefix << "-tree.nex";
        ofstream nex_tree(sstm.str());
        int precision = 5;
        string newick = test_tree.make_newick(precision);
        test_tree.write_nexus(newick, nex_tree);
        nex_tree.close();
        sstm.str("");

        sstm << dir << prefix << "-tree-nmut.nex";
        ofstream nex_tree2(sstm.str());
        precision = 0;
        newick = test_tree.make_newick_nmut(precision, nmuts);
        test_tree.write_nexus(newick, nex_tree2);
        nex_tree2.close();
        sstm.str("");
    }

    if(print_mut){
        sstm << dir << prefix << "-mut.txt";
        ofstream out_mut(sstm.str());
        for(int j=0; j<test_tree.nleaf; ++j){
            for(int i=0; i<results[j].mutations.size(); ++i){
                out_mut << j+1 << "\t" << results[j].mutations[i].edge_id+1
        << "\t" << results[j].mutations[i].type << "\t" << results[j].mutations[i].btime << "\t" << results[j].mutations[i].gtime << endl;
            }
        }
        out_mut.close();
        sstm.str("");
    }

    sstm << dir << prefix << "-rel-times.txt";
    double node_min = NODE_MIN_TIME;
    for(int j=0; j<test_tree.nleaf-1; ++j){
        if( test_tree.node_times[j] < node_min ) node_min = test_tree.node_times[j];
    }
    //cout << "leaf minimum: " << node_min << endl;
    ofstream out_rel(sstm.str());
    for(int j=0; j<test_tree.nleaf-1; ++j){
        double delta = test_tree.node_times[j] - node_min;
        out_rel << j+1 << "\t" << delta << "\t" << ceil(age + delta) << endl;
    }
    out_rel.close();
    sstm.str("");

}


void run_simulations(int mode, const vector<int>& chr_lengths, int Ns, int Nsims, int model, int cons, int cn_max, int Ne, double delta_t, int age, const vector<double>& rate_consts, const vector<double>& var_size, string dir, string prefix, int print_allele, int print_mut, int print_nex){
    // can specify the germline in a number of different ways
    //genome germline(1,10);
    //genome germline(42,1000);
    //genome germline(chr_lengths);
    genome germline(chr_lengths, NORM_PLOIDY);
    germline.mean_dup_size = var_size[0];
    germline.mean_del_size = var_size[1];

    vector<int> edges;
    vector<double> lengths;
    vector<double> epoch_times;
    vector<double> node_times;

    vector<genome> results;
    map<int, vector<mutation>> muts;   // hold all mutations by edge id
    string orig_prefix = prefix;

    for(int i=0; i<Nsims; ++i){
        cout << "Simulation " << i+1 << endl;
        //cout << "\n\n###### New sample collection ######" << endl;
        if(orig_prefix.empty()) {
          prefix = "sim-data-" + to_string(int(i+1));
        }
        cout << "Prefix of output file: " << prefix <<endl;
        generate_coal_tree(Ns, edges, lengths, epoch_times, node_times);
        if(debug){
          cout << "Initial coalescence tree: " << endl;
          evo_tree test_tree0(Ns+1, edges, lengths);
          test_tree0.print();
        }

        // Scale branch lengths by Ne
        for(int l = 0; l < lengths.size(); ++l) lengths[l] = lengths[l] * Ne;

        assign_tip_times(delta_t, Ns, tobs, edges, lengths);

        evo_tree test_tree(Ns+1, edges, lengths);
        if(debug){
          cout << "Tree after tip sampling: " << endl;
          test_tree.print();
        }
        // Ensure that the rescaled tree by age has a height smaller than the age of last sample
        // Need tobs to get the miminmal tree height
        if(cons){
          double min_height = *max_element(tobs.begin(), tobs.end());
          double max_height = age + min_height;
          double tree_height = runiform(r, min_height, max_height);
          double old_tree_height = test_tree.get_tree_height();
          double ratio = tree_height/old_tree_height;
          // test_tree.scale_time(ratio);
          // Only scaling internal nodes to allow large differences at the tip
          test_tree.scale_time_internal(ratio);
          if(debug){
              cout << "Tree height before scaling " << old_tree_height << endl;
              cout << "Scaling ratio " << ratio << endl;
              test_tree.print();
          }
          cout << "Simulated tree height " << test_tree.get_tree_height() << endl;
        }

        simulate_samples(results, muts, test_tree, germline, chr_lengths, rate_consts, model, cn_max, false);
        print_simulations(mode, results, muts, test_tree, dir, prefix, age, print_allele, print_mut, print_nex);

        edges.clear();
        lengths.clear();
        results.clear();
        muts.clear();
        epoch_times.clear();
        node_times.clear();
    }
}

//////////////////////////////////////////////////////////
///                                                    ///
///   MAIN                                             ///
///                                                    ///
//////////////////////////////////////////////////////////

// Tree representation
// if we have n samples these trees always have the germline variation as the n+1 samples
// the final edge length is always zero
// edges are always directed from -> to going from the node/MRCA
// therefore leaf nodes are always in the second column of edges and the root is only in the first

// Internally nodes and edges have ids starting from 0
// to match up with ape in R we always print id+1
int main (int argc, char ** const argv) {
    string dir; // output directory
    string prefix; // prefix of output file
    int Ns; // number of regions
    int Nsims;  // number of multi-region samples
    // five event types: duplication, deletion, chromosome gain, chromosome loss, wgd
    // rates are 1/mean
    double dup_rate, del_rate, chr_gain, chr_loss, wgd;
    // parameters for mean of dup/del size distributions
    double dup_size, del_size;
    // effective population size
    int Ne;
    // relative timing difference
    double delta_t;
    int seed, mode;
    int model, cons;
    int cn_max, seg_max;
    int print_allele, print_mut, print_nex;

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
      ("mode", po::value<int>(&mode)->default_value(0), "running mode of the program (0: Simuting genomes in fix-sized bins, 1: Simulating genomes in random segments, 2 to 4: Test)")

      ("model,d", po::value<int>(&model)->default_value(0), "model of evolution (0: JC69, 1: 1-step bounded, 2: Poisson)")
      ("age,a", po::value<int>(&age)->default_value(100), "age of the patient to simulate")
      ("epop,e", po::value<int>(&Ne)->default_value(2), "effective population size")
      ("tdiff,t", po::value<double>(&delta_t)->default_value(1), "relative timing difference")
      ("constrained", po::value<int>(&cons)->default_value(1), "constraints on branch length (0: none, 1: fixed total time)")
      ("cn_max", po::value<int>(&cn_max)->default_value(4), "maximum copy number of a segment")
      ("seg_max", po::value<int>(&seg_max)->default_value(100), "maximum number of segments to simulate")

      ("nregion,r", po::value<int>(&Ns)->default_value(5), "number of regions")
      ("nsim,n", po::value<int>(&Nsims)->default_value(1), "number of multi-region samples")
      ("prefix,p", po::value<string>(&prefix)->default_value(""), "prefix of output file (it will be sim-data-N if not specified")

      ("dup_rate", po::value<double>(&dup_rate)->default_value(0.0001), "duplication rate (allele/locus/year)")
      ("del_rate", po::value<double>(&del_rate)->default_value(0.0002), "deletion rate (allele/locus/year)")
      ("chr_gain", po::value<double>(&chr_gain)->default_value(0), "chromosome gain rate (haplotype/chr/year)")
      ("chr_loss", po::value<double>(&chr_loss)->default_value(0), "chromosome loss rate (haplotype/chr/year)")
      ("wgd", po::value<double>(&wgd)->default_value(0), "WGD (whole genome doubling) rate (year)")
      ("dup_size", po::value<double>(&dup_size)->default_value(30), "mean of duplication size distributions")
      ("del_size", po::value<double>(&del_size)->default_value(30), "mean of deletion size distributions")

      ("print_allele", po::value<int>(&print_allele)->default_value(1), "whether or not to output allele-specific copy numbers")
      ("print_mut", po::value<int>(&print_mut)->default_value(1), "whether or not to output the list of mutations")
      ("print_nex", po::value<int>(&print_nex)->default_value(1), "whether or not to output the tree in NEXUS format")

      ("seed", po::value<int>(&seed)->default_value(0), "seed used for generating random numbers")
      ("verbose", po::value<int>(&debug)->default_value(0), "verbose level (0: default, 1: debug)")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(required).add(optional);
    po::variables_map vm;

    try {
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
    } catch (const std::exception& e) {
          std::cerr << e.what() << std::endl;
          return 1;
    }

  setup_rng(seed);

  // output directory
  if(dir.back() != '/'){
      dir = dir + "/";
  }

  vector<int> chr_lengths = CHR_BIN_SIZE;
  if(mode==1){
      cout << "Under this model, each site is treated as a final segment. The mean duplication/deletion size is currently fixed to be 1" << endl;
      dup_size = 1;
      del_size = 1;
      // Randomly generate the number of segments
      int num_seg = runiform(r, 0.2, 0.8) * seg_max;
      // cout << "Approximate number of segments to simulate is " << num_seg << endl;
      // Distribute the segments according to the size of chromosmes
      double theta[NUM_CHR] = {};
      double alpha[NUM_CHR] = {};
      // int bin_size = accumulate(CHR_BIN_SIZE.begin(), CHR_BIN_SIZE.end(), 0);
      // cout << "Total number of bins is " << bin_size << endl;
      int total_seg = 0;
      for(int i = 0; i<NUM_CHR; i++){
          alpha[i] = CHR_BIN_SIZE[i];
      }
      gsl_ran_dirichlet(r, NUM_CHR, alpha, theta);
      for(int i = 0; i<NUM_CHR; i++){
          chr_lengths[i] = ceil(theta[i] * num_seg);
          assert(chr_lengths[i] > 0);
          total_seg += chr_lengths[i];
          // cout << alpha[i] << "\t" << theta[i] << "\t" << chr_lengths[i] << endl;
      }
      cout << "Number of segments simulated is " << total_seg << endl;
  }
  vector<double> var_size = {dup_size, del_size};
  vector<double> rate_consts = {dup_rate, del_rate, chr_gain, chr_loss, wgd};

  // priors on rates
  //vector<double> pars;
  //read_params(argv[4],pars);
  cout << "rates:\t" << rate_consts[0] << "\t" << rate_consts[1]  << "\t" << rate_consts[2]  << "\t" << rate_consts[3]  << "\t" << rate_consts[4] << endl;
  cout << "sizes of segment duplication/deletion:\t" << var_size[0] << "\t" << var_size[1] << endl;

  // simulate coalescent tree and apply SVs
  if(mode >= 0){
      run_simulations(mode, chr_lengths, Ns, Nsims, model, cons, cn_max, Ne, delta_t, age, rate_consts, var_size, dir, prefix, print_allele, print_mut, print_nex);
  }
  else{
      run_test(mode, dir);
  }

}
