/*
This program aims to simulate SVs across a genome along a phylogenetic tree.

Tree representation:
if we have n samples these trees always have the germline variation as the n + 1 samples
the final edge length is always zero
edges are always directed from -> to going from the node/MRCA
therefore leaf nodes are always in the second column of edges and the root is only in the first
Internally nodes and edges have ids starting from 0
to match up with ape in R we always print id + 1

Output copy number representation:
sample ID, chromosome ID, segment ID, CN.
Here, the ID for sample, chromosome and segment start from 1.
*/

#include "model.hpp"
#include "parse_cn.hpp"
#include "evo_tree.hpp"
#include "tree_op.hpp"
#include "genome.hpp"


#include <boost/program_options.hpp>

// using namespace std;

// The number of bins for each chromosome. Each bin corresponds to a window of size 500,000 bp. 4401 bins in total
vector<int> CHR_BIN_SIZE(BIN_SIZE, BIN_SIZE + sizeof(BIN_SIZE) / sizeof(BIN_SIZE[0]));

const int NODE_MIN_TIME = 1000;

enum SIM_METHOD {SIM_TIME, SIM_SEQ};
enum MUT_TYPE {DUP, DEL, GAIN, LOSS, WGD};

struct PRINT_LEVEL{
    int print_allele;
    int print_mut;
    int print_nex;
    int print_relative;
    int print_baseline;
};

struct SV_SIZE{
  double mean_dup_size;
  double mean_del_size;
};

gsl_rng* r;

// unary function and pointer to unary function
// allows use of gsl rng for standard template algorithms
inline long unsigned myrng(long unsigned n){
  return gsl_rng_uniform_int(r, n);
}

long unsigned (*fp_myrng)(long unsigned);


// Set the initial copy number matrix
copy_number initialize_cn(const vector<int>& chr_lengths, int model){
    copy_number cn;

    for(int i = 0; i < chr_lengths.size(); ++i){
      int num_seg = chr_lengths[i];

      for(int j = 0; j < num_seg; j++){
          if(model == BOUNDT){
             cn[i][j] = NORM_PLOIDY;     // index of normal state
          }else{
             cn[i][j] = NORM_ALLElE_STATE;
          }
       }
    }

    return cn;
}


// used for debugging
void print_cn_state(const copy_number& curr_cn){
    for(auto cn_profile : curr_cn){
        for(auto seg: cn_profile.second){
            // chr, seg, cn
            cout << cn_profile.first + 1 << "\t" << seg.first + 1 << "\t" << seg.second << endl;
        }
    }
}


// Envolving sequences along the tree (available for bounded model, implemented according to approach in book Yang, 2004, P437)
// only support site duplication and deletion
void evolve_sequences(map<int, copy_number>& cn_matrix, const int& node_id, const evo_tree& tree, double* qmat, int nstate, int num_seg, const vector<double>& rate_consts, map<int, int>& num_muts, int debug){
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
                gsl_ran_discrete_t* dis = gsl_ran_discrete_preproc(nstate, probs);
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
      evolve_sequences(cn_matrix, tree.nodes[node_id].daughters[0], tree, qmat, nstate, num_seg, rate_consts, num_muts, debug);
      evolve_sequences(cn_matrix, tree.nodes[node_id].daughters[1], tree, qmat, nstate, num_seg, rate_consts, num_muts, debug);
    }
}


void write_cn(map<int, copy_number>& cn_matrix, int node_id, ogzstream& out, int cn_max, int model){
    copy_number cn_profile = cn_matrix[node_id];

    for(auto c : cn_profile){
        for(auto s : c.second){
            int cn = s.second;
            if(model == BOUNDA) cn = state_to_total_cn(s.second, cn_max);
            // node_id, chr, seg_id, copy number
            // if(debug) cout << node_id + 1 << "\t" << c.first + 1 << "\t" << s.first + 1 << "\t" << s.second << "\t" << cn << endl;
            out << node_id + 1 << "\t" << c.first + 1 << "\t" << s.first + 1 << "\t" << cn << endl;
        }
    }
}


// Only segment-level mutations are allowed when simluating sequences directly
void write_rcn(map<int, copy_number>& cn_matrix, int node_id, ogzstream& out, int cn_max, int model){
    copy_number cn_profile = cn_matrix[node_id];

    for(auto c : cn_profile){
        for(auto s : c.second){
            int cn = s.second;
            if(model == BOUNDA) cn = state_to_total_cn(s.second, cn_max);
            // node_id, chr, seg_id, copy number
            int rcn = cn - NORM_PLOIDY;
            if(rcn < -2) rcn = -2;
            if(rcn > 2) rcn = 2;
            // if(debug) cout << node_id + 1 << "\t" << c.first + 1 << "\t" << s.first + 1 << "\t" << s.second << "\t" << cn << "\t" << rcn << endl;
            out << node_id + 1 << "\t" << c.first + 1 << "\t" << s.first + 1 << "\t" << rcn << endl;
        }
    }
}


// Write relative CNs based on baseline strategy in PISCA
void write_allele_rcn(map<int, copy_number>& cn_matrix, int node_id, ogzstream& out, int cn_max, int model, gsl_rng* r){
    copy_number cn_profile = cn_matrix[node_id];
    int ploidy = NORM_PLOIDY / 2;

    for(auto c : cn_profile){
        for(auto s : c.second){
            int state = s.second;
            int cnA = 0;
            int cnB = 0;
            state_to_allele_cn(state, cn_max, cnA, cnB);
            // node_id, chr, seg_id, copy number
            int rcnA = get_rcn_baseline(cnA, ploidy, r);
            int rcnB = get_rcn_baseline(cnB, ploidy, r);
            // int rcn = rcnA + rcnB;
            // if(debug) cout << node_id + 1 << "\t" << c.first + 1 << "\t" << s.first + 1 << "\t" << s.second << "\t" << rcnA << "\t" << rcnB << endl;
            out << node_id + 1 << "\t" << c.first + 1 << "\t" << s.first + 1 << "\t" << rcnA << "\t" << rcnB << endl;
        }
    }
}


void write_allele_cn(map<int, copy_number>& cn_matrix, int node_id, ogzstream& out, int cn_max){
    copy_number cn_profile = cn_matrix[node_id];

    for(auto c : cn_profile){
        for(auto s : c.second){
            int state = s.second;
            int cnA = 0;
            int cnB = 0;
            state_to_allele_cn(state, cn_max, cnA, cnB);
            // node_id, chr, seg_id, copy number
            // if(debug) cout << node_id + 1 << "\t" << c.first + 1 << "\t" << s.first + 1 << "\t" << s.second << "\t" << cnA << "\t" << cnB << endl;
            out << node_id + 1 << "\t" << c.first + 1 << "\t" << s.first + 1 << "\t" << cnA << "\t" << cnB << endl;
        }
    }
}

void write_sample_times(stringstream& sstm, evo_tree& test_tree, string dir, string prefix, int age){
    cout << "Writing sampling times" << endl;
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
        out_rel << j + 1 << "\t" << delta << "\t" << ceil(age + delta) << endl;
    }
    out_rel.close();
    sstm.str("");
}

// Print the simulated copy numbers
void print_sequences(map<int, copy_number>& cn_matrix, int cn_max, int model, int num_seg, map<int, int> num_muts, evo_tree& test_tree, string dir, string prefix, int age, PRINT_LEVEL print_level, gsl_rng* r, int cons = 1){
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
    cout << "Reading data and calculating CNA regions" << endl;
    vector<vector<vector<int>>> s_info = read_cn(sstm.str(), test_tree.nleaf - 1, num_total_bins, cn_max);

    vector<vector<int>> segs;
    segs = get_all_segs(s_info, test_tree.nleaf - 1, num_total_bins, num_invar_bins, 1);
    int seg_size = segs.size();
    sstm.str("");

    if(print_level.print_allele && model == BOUNDA){
        cout << "Writing allele specific copy number " << endl;
        sstm << dir << prefix << "-allele-cn.txt.gz";
        ogzstream out_allele_cn(sstm.str().c_str());
        for(int j = 0; j < test_tree.nleaf; ++j){
          write_allele_cn(cn_matrix, j, out_allele_cn, cn_max);
        }
        out_allele_cn.close();
        sstm.str("");
    }

    if(print_level.print_relative){
        cout << "Writing relative total copy number " << endl;
        sstm << dir << prefix << "-rcn.txt.gz";
        ogzstream out_rcn(sstm.str().c_str());
        for(int j = 0; j < test_tree.nleaf; ++j){
            write_rcn(cn_matrix, j, out_rcn, cn_max, model);
        }
        out_rcn.close();
        sstm.str("");
    }

    if(print_level.print_baseline){
        cout << "Writing allele specific total copy number " << endl;
        sstm << dir << prefix << "-allele-rcn.txt.gz";
        ogzstream out_rcn_baseline(sstm.str().c_str());
        for(int j = 0; j < test_tree.nleaf; ++j){
            write_allele_rcn(cn_matrix, j, out_rcn_baseline, cn_max, model, r);
        }
        out_rcn_baseline.close();
        sstm.str("");
    }


    cout << "Writing total copy number of internal nodes" << endl;
    sstm << dir << prefix << "-inodes-cn.txt.gz";
    int ntotn = (2 * test_tree.nleaf - 1);
    ogzstream out_cn_inodes(sstm.str().c_str());
    for(int j = test_tree.nleaf; j < ntotn; ++j){
        write_cn(cn_matrix, j, out_cn_inodes, cn_max, model);
    }
    out_cn_inodes.close();
    sstm.str("");

    if(print_level.print_allele){
        cout << "Writing allele specific copy number of internal nodes" << endl;
        sstm << dir << prefix << "-inodes-allele-cn.txt.gz";
        ogzstream out_allele_cn_inodes(sstm.str().c_str());
        for(int j = test_tree.nleaf; j < ntotn; ++j){
          write_allele_cn(cn_matrix, j, out_allele_cn_inodes, cn_max);
        }
        out_allele_cn_inodes.close();
        sstm.str("");
    }

    cout << "Writing all the mutational events" << endl;
    sstm << dir << prefix << "-info.txt";
    ofstream out_info(sstm.str());
    out_info << "NODE TIMES:";
    out_info << "\tnid\ttime" << endl;
    for(int j = 0; j <  test_tree.nodes.size(); ++j){
        out_info << "\t" << j + 1 << "\t" << test_tree.nodes[j].time << endl;
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

    cout << "Writing the tree in TXT format" << endl;
    sstm << dir << prefix << "-tree.txt";
    ofstream out_tree(sstm.str());
    //test_tree.write(out_tree);
    out_tree << "start\tend\tlength\teid\tnmut" << endl;
    for(int i = 0; i < test_tree.edges.size(); ++i){
        out_tree << test_tree.edges[i].start + 1 << "\t" << test_tree.edges[i].end + 1 << "\t" << test_tree.edges[i].length << "\t" << test_tree.edges[i].id + 1 << "\t" << num_muts[test_tree.edges[i].id] << endl;
    }
    out_tree.close();
    sstm.str("");

    if(print_level.print_nex){
        cout << "Writing the tree in NEXUS format" << endl;
        sstm << dir << prefix << "-tree.nex";
        ofstream nex_tree(sstm.str());
        string newick = test_tree.make_newick();
        test_tree.write_nexus(newick, nex_tree);
        nex_tree.close();
        sstm.str("");
    }

    if(cons){
        write_sample_times(sstm, test_tree, dir, prefix, age);
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


// Simulate mutations under different models by simulating waiting times of a Markov chain
// Each mutation occur at some sites on one chromosome (haplotype-specific)
vector<mutation> generate_mutation_by_model(genome& g, const int& edge_id, const double& blength, const double& node_time, const vector<int>& chr_lengths, const vector<double>& rate_consts, const SV_SIZE& sv_size, int model, int cn_max, int& num_fail, int debug){
    if(debug){
        cout << "\tGenerate_mutations at node ID " << g.node_id + 1 << " (time " << node_time << "), edge ID " << edge_id + 1 << " with branch length:" << "\t" << blength << endl;
    }

    vector<mutation> ret;
    double avg_rate = 0.0;  // The mean of genome mutation rates along an edge;
    int nsite = accumulate(chr_lengths.begin(), chr_lengths.end(), 0);
    double time = 0.0;
    int count = 0;  // Count the number of mutations

    while(time < blength){
        // The rate is dynamically changing according to the status of each site
        vector<double> site_dup_rates(nsite, 0.0), site_del_rates(nsite, 0.0);
        vector<double> chr_gain_rates(NUM_CHR, 0.0), chr_loss_rates(NUM_CHR, 0.0);
        vector<double> type_rates;  // The rates of each mutation type

        if(debug){
          cout << "\tComputing total mutation rate on the genome with " << nsite << " sites " << endl;
          cout << "GENOME BEFORE " << endl;
          g.print();
          g.print_cn();
          cout << "There are " << g.chrs.size() << " chromosome IDs (including empty ones) before" << endl;
          cout << "There are " << get_num_available_chr(g) << " non-empty chromosomes before" << endl;
        }

        double rate = get_total_rates_allele_specific(g, site_dup_rates, site_del_rates, chr_gain_rates, chr_loss_rates, type_rates, rate_consts, model, cn_max, debug);
        if(rate <= 0){
            if(debug) cout << "This genome cannot be mutated any more on edge " << edge_id + 1 << endl;
            break;
        }
        avg_rate = avg_rate + rate;

        double tevent = gsl_ran_exponential(r, 1 / rate);
        time += tevent;
        count++;

        int e = rchoose(r, type_rates);     // assign the event to a site with probabilities proportional to the rates at the site

        if(debug){
            cout << "Current time " << time << endl;
            cout << "Total mutation rate on the genome " << rate << endl;
            cout << "Chosen event " << MUT_TYPES[e] << endl;
        }

        int res = 0;
        int site = 0;  // site across the whole reference genome
        int c = 0;
        int seg_id = 0;

        // Randomly select a site for duplication or deletion. Convert the site position back to chromosome ID and segment ID
        if(e == DUP){   //duplications
            site = rchoose(r, site_dup_rates);
            seg_id = site;
            int loc = site2chr(site, c, chr_lengths, debug);    // loc is the relative location of a site on c, previously used as seg_id
            select_haplotype_by_seg(g, c, seg_id, fp_myrng);
            res = generate_duplication(g, c, seg_id, sv_size.mean_dup_size, model, cn_max, r, debug);
            if(debug){
                cout << "duplication of site " << site << ", chr " << c + 1 << " seg " << seg_id + 1 << " location " << loc + 1 << endl;
            }
        }else if(e == DEL){   //deletions
            site = rchoose(r, site_del_rates);
            seg_id = site;
            int loc = site2chr(site, c, chr_lengths, debug);
            select_haplotype_by_seg(g, c, seg_id, fp_myrng);
            res = generate_deletion(g, c, seg_id, sv_size.mean_del_size, model, cn_max, r, debug);
            if(debug){
                cout << "deletion of site " << site << ", chr " << c + 1 << " seg " << seg_id + 1 << " location " << loc + 1 << endl;
            }
        }else if(e == GAIN){   //chr gain
            c = rchoose(r, chr_gain_rates);
            seg_id = -1;
            select_haplotype(g, c, fp_myrng);
            res = generate_chr_gain(g, c, cn_max, debug);
            if(debug){
                cout << "chosen chr " << c + 1 << " for gain event" << endl;
            }
        }else if(e == LOSS){   //chr loss
            c = rchoose(r, chr_loss_rates);
            seg_id = -1;
            select_haplotype(g, c, fp_myrng);
            res = generate_chr_loss(g, c, debug);
            if(debug){
                cout << "chosen chr " << c + 1 << " for loss event" << endl;
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
            // time/blength: relative time in the branch; node_time + time: real time when the mutation generates
            mutation mut(edge_id, e, time / blength, node_time + time, c, seg_id);
            ret.push_back(mut);
            g.nmuts[e]++;
            g.mutations.push_back(mut);
            if(debug) cout << "Update copy number after each mutation event" << endl;
            g.calculate_cn();
            g.calculate_allele_cn();
        }else{
             if(debug){
                cout << "\tMutation failed on chr " << c + 1 << " with segments below:";
                for(int i = 0; i < g.chrs[c].size(); i++){
                    cout << "\t" << g.chrs[c][i].seg_id;
                }
                cout << endl;
             }
             // time -= tevent;
             num_fail++;
             continue;
        }

        if(debug){
            cout << "node\tedge\ttime2event\ttotalTime\tbranchLength\tevent\tchr\tsite\tcopyNumber\n\t" << g.node_id + 1 << "\t" << edge_id + 1 << "\t" << tevent << "\t" << time << "\t" << blength << "\t  " << MUT_TYPES[e] << "\t" << c + 1 << "\t" << seg_id + 1 << "\t" << g.cn_profile[c % NUM_CHR][seg_id] << endl;

            cout << "GENOME AFTER " << MUT_TYPES[e] << endl;
            g.print_cn();
            g.print();

            cout << "There are " << get_num_available_chr(g) << " non-empty chromosomes now" << endl;
        }
    }

    if(count > 0){
        avg_rate = avg_rate / count;
    }
    cout << "Average genome mutation rate (per year) on edge " << edge_id << " is: " << avg_rate << endl;

    return ret;
}


// simulate mutations under infinite sites model
vector<mutation> generate_mutation_times(const int& edge_id, const double& blength, const double& node_time, const vector<double>& rate_consts, int debug){
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
    double tevent = gsl_ran_exponential(r, 1 / rate);
    time += tevent;

    // choose type of event
    int e = rchoose(r, rates);
    if(debug) cout << "mut times, tevent, total time, time/branch len, event\t" << tevent << "\t" << time << "\t" << blength << "\t" << e << endl;

    ret.push_back(mutation( edge_id, e, time / blength, node_time + time, 0, 0));
  }

  ret.pop_back();

  return ret;
}


void apply_mutations(const vector<mutation>& muts, genome& g, const SV_SIZE& sv_size, int debug){
  // mean variant size in bins
  double mdup = sv_size.mean_dup_size;
  double mdel = sv_size.mean_del_size;

  for(int i = 0; i < muts.size(); ++i){
    g.nmuts[ muts[i].type]++;
    g.mutations.push_back(muts[i]);

    if(muts[i].type == DUP){   //duplications
        if(debug){
            cout << "GENOME BEFORE duplication " << endl;
            g.print();
            // g.print_cn();
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
            int loc = gsl_rng_uniform_int(r, g.chrs[c].size() - len);
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
        g.chrs.insert(g.chrs.end(), g.chrs.begin()+c, g.chrs.begin()+c + 1);
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


void traverse_tree_mutating(const int& node_id, const evo_tree& tree, const vector<int>& chr_lengths, const vector<double>& rate_consts, const SV_SIZE& sv_size, map<int, vector<mutation>>& all_muts, map<int, int>& failed_muts, vector<genome>& genomes, int model, int cn_max, int debug){
  //cout << "\ttraverse_tree: " << node_id + 1 << endl;
  if(!tree.nodes[node_id].isRoot){
    // copy the parent
    genomes[node_id] = genomes[tree.nodes[node_id].parent];
    genomes[node_id].node_id = node_id;

    int edge_id = tree.nodes[node_id].e_in;
    vector<mutation> muts;
    int num_fail = 0;

    // apply the mutations from parent -> daughter
    if(tree.edges[edge_id].length > 0){
        if(debug){
            cout << "MUTATING genome: " << tree.nodes[node_id].parent + 1 << " -> " << node_id + 1 << "\t edge id: " << tree.nodes[node_id].e_in + 1 << endl;
            cout << "Generating mutations for node " << node_id << endl;
        }

        if(model == INFINITE){
            muts = generate_mutation_times(edge_id, tree.edges[edge_id].length, tree.nodes[tree.nodes[node_id].parent].time, rate_consts, debug);
            apply_mutations(muts, genomes[node_id], sv_size, debug);
        }else{
            muts = generate_mutation_by_model(genomes[node_id], edge_id, tree.edges[edge_id].length, tree.nodes[tree.nodes[node_id].parent].time, chr_lengths, rate_consts, sv_size, model, cn_max, num_fail, debug);
        }

        // Print out copy number for checking
        if(debug){
            cout << "Generated mutations for node " << node_id << endl;
            print_cnp(genomes[node_id].cn_profile);
        }
    }

    all_muts[edge_id] = muts;
    failed_muts[edge_id] = num_fail;
  }

  if(!tree.nodes[node_id].daughters.empty()){
    traverse_tree_mutating(tree.nodes[node_id].daughters[0], tree, chr_lengths, rate_consts, sv_size, all_muts, failed_muts, genomes, model, cn_max, debug);
    if(debug) cout << "Finished for node " << tree.nodes[node_id].daughters[0] << endl;

    traverse_tree_mutating(tree.nodes[node_id].daughters[1], tree, chr_lengths, rate_consts, sv_size, all_muts, failed_muts, genomes, model, cn_max, debug);
    if(debug) cout << "Finished for node " << tree.nodes[node_id].daughters[1] << endl;
  }

}



void simulate_samples(vector<genome>& genomes, map<int,vector<mutation>>& muts, map<int,int>& failed_muts, const evo_tree& tree, genome& germline, const vector<int>& chr_lengths, const vector<double>& rate_consts, const SV_SIZE& sv_size, int model, int cn_max, int debug = 0){
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
  traverse_tree_mutating(tree.root_node_id, tree, chr_lengths, rate_consts, sv_size, muts, failed_muts, genomes, model, cn_max, debug);

  // final samples returned, print out leaf nodes
  if(debug){
    cout << "MUTATIONS:" << endl;
    for(map<int, vector<mutation> >::iterator it = muts.begin(); it != muts.end(); it++){
      cout << "EDGE, id: " << it->first + 1
       << "\t" << tree.edges[it->first].start + 1 << " -> " << tree.edges[it->first].end + 1
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
      genomes[i].print_muts(cout);
      //genomes[i].print_cn();
      tree.print_ancestral_edges(genomes[i].node_id);
      cout << endl;
    }
  }
}


void run_sample_set(int Ns, gsl_rng* r, unsigned seed, const ITREE_PARAM& itree_param, const vector<double>& rate_consts, const SV_SIZE& sv_size, double* pvs, int* ret){
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
  //cout << "###### Ns + 1= " << Ns + 1 << endl;

  generate_coal_tree(Ns, r, fp_myrng, edges, lengths, epoch_times, node_times, itree_param);
  evo_tree test_tree(Ns + 1, edges, lengths);

  //for(int i = 0; i < 6; ++i) epars.push_back( ptree[i]);
  //for(int i = 0; i < 8; ++i) lengths.push_back( pl[i]);
  //evo_tree test_tree = construct_tree(Ns, epars, lengths, node_times);

  simulate_samples(results, muts, failed_muts, test_tree, germline, chr_lengths, rate_consts, sv_size, model, cn_max, false);

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


void run_test(int mode, string dir, unsigned seed, const ITREE_PARAM& itree_param, const SV_SIZE& sv_size, int debug){
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
        dels.push_back(mutation(0, 1, 0.5, 0, 1, 2));
        dels.push_back(mutation(0, 1, 0.5, 0, 2, 3));
        dels.push_back(mutation(0, 1, 0.5, 0, 3, 4));
        dels.push_back(mutation(0, 1, 0.5, 0, 4, 5));

        vector<mutation> dups;
        dups.push_back(mutation(0, 0, 0.5, 0, 11, 2));
        dups.push_back(mutation(0, 0, 0.5, 0, 12, 3));
        dups.push_back(mutation(0, 0, 0.5, 0, 13, 4));
        dups.push_back(mutation(0, 0, 0.5, 0, 14, 5));

        vector<mutation> wgds;
        wgds.push_back(mutation(0, 4, 0.5, 0, -1, -1));

        genome g1(2,10);
        g1.node_id = 0;
        g1.print();
        g1.print_cn();
        apply_mutations(dups, g1, sv_size, debug);
        apply_mutations(dels, g1, sv_size, debug);
        apply_mutations(wgds, g1, sv_size, debug);
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

        genome g2d(chr_lengths, 2);
        g2d.node_id = 0;
        g2d.print();
        g2d.print_cn();
      }
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
        run_sample_set(Ns, r, seed, itree_param, rate_consts, sv_size, pvs, &(ret[0]));

        sstm << dir << "sim-data-" << i + 1 << "-cn.txt";
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


void print_simulations(int mode, int num_seg, vector<genome>& results, const vector<int>& chr_lengths, map<int, vector<mutation>>& muts, map<int, int>& failed_muts, evo_tree& test_tree, string dir, string prefix, int age, int cn_max, PRINT_LEVEL print_level, gsl_rng* r, int cons, int debug){
    stringstream sstm;

    sstm << dir << prefix << "-cn.txt.gz";
    ogzstream out_cn(sstm.str().c_str());
    for(int j = 0; j < test_tree.nleaf; ++j){
        results[j].write(out_cn);
    }
    out_cn.close();

    int num_total_bins = 0;
    int num_invar_bins = 0;
    int is_total = 1;
    int is_rcn = 0;
    cout << "Reading data and calculating CNA regions" << endl;
    vector<vector<vector<int>>> s_info = read_cn(sstm.str(), test_tree.nleaf - 1, num_total_bins, cn_max, is_total, is_rcn, debug);

    vector<vector<int>> segs;
    if(mode == 0){
        cout << "Merging consecutive bins with the same CN" << endl;
        segs = get_invar_segs(s_info, test_tree.nleaf - 1, num_total_bins, num_invar_bins, is_total, debug);
    }else{
        segs = get_all_segs(s_info, test_tree.nleaf - 1, num_total_bins, num_invar_bins, 1, is_total, debug);
    }
    int seg_size = segs.size();
    sstm.str("");

    if(print_level.print_allele){
        cout << "Writing allele specific copy number " << endl;
        sstm << dir << prefix << "-allele-cn.txt.gz";
        ogzstream out_allele_cn(sstm.str().c_str());
        for(int j = 0; j < test_tree.nleaf; ++j){
          // cout << "Chr size " << results[j].chrs.size() << endl;
          results[j].write_allele_cn(out_allele_cn);
        }
        out_allele_cn.close();
        sstm.str("");
    }

    if(print_level.print_relative){
        cout << "Writing relative total copy number " << endl;
        sstm << dir << prefix << "-rcn.txt.gz";
        ogzstream out_rcn(sstm.str().c_str());
        for(int j = 0; j < test_tree.nleaf; ++j){
          results[j].write_rcn(out_rcn, r);
        }
        out_rcn.close();
        sstm.str("");
    }

    if(print_level.print_baseline){
        cout << "Writing allele specific total copy number " << endl;
        sstm << dir << prefix << "-allele-rcn.txt.gz";
        ogzstream out_rcn_baseline(sstm.str().c_str());
        for(int j = 0; j < test_tree.nleaf; ++j){
            // cout << "writing baseline RCN for node " << j + 1 << endl;
            results[j].write_allele_rcn(out_rcn_baseline, r);
        }
        out_rcn_baseline.close();
        sstm.str("");
    }

    int ntotn = 2 * test_tree.nleaf - 1;
    cout << "Writing total copy number of internal nodes" << endl;
    sstm << dir << prefix << "-inodes-cn.txt.gz";
    ogzstream out_cn_inodes(sstm.str().c_str());
    for(int j = test_tree.nleaf; j < ntotn; ++j){
      results[j].write(out_cn_inodes);
    }
    out_cn_inodes.close();
    sstm.str("");

    if(print_level.print_allele){
        cout << "Writing allele specific copy number of internal nodes" << endl;
        sstm << dir << prefix << "-inodes-allele-cn.txt.gz";
        ogzstream out_allele_cn_inodes(sstm.str().c_str());
        for(int j = test_tree.nleaf; j < ntotn; ++j){
          results[j].write_allele_cn(out_allele_cn_inodes);
        }
        out_allele_cn_inodes.close();
        sstm.str("");
    }

    cout << "Writing all the mutational events" << endl;
    sstm << dir << prefix << "-info.txt";
    ofstream out_info(sstm.str());
    out_info << "NODE TIMES:\n";
    out_info << "\tnid\ttime" << endl;
    for(int j = 0; j <  test_tree.nodes.size(); ++j){
        out_info << "\t" << j + 1 << "\t" << test_tree.nodes[j].time << endl;
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

    cout << "Writing the tree in TXT format" << endl;
    sstm << dir << prefix << "-tree.txt";
    vector<int> nmuts;
    ofstream out_tree(sstm.str());
    out_tree << "start\tend\tlength\teid\tnmut\tnfailed" << endl;
    for(int i = 0; i < test_tree.edges.size(); ++i){
        out_tree << test_tree.edges[i].start + 1 << "\t" << test_tree.edges[i].end + 1 << "\t" << test_tree.edges[i].length << "\t" << test_tree.edges[i].id + 1 << "\t" << muts[test_tree.edges[i].id].size() << "\t" << failed_muts[test_tree.edges[i].id] << endl;
        nmuts.push_back(muts[test_tree.edges[i].id].size());
    }
    out_tree.close();
    sstm.str("");

    if(print_level.print_nex){
        cout << "Writing the tree in NEXUS format" << endl;
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

    if(print_level.print_mut){
        sstm << dir << prefix << "-mut.txt";
        ofstream out_mut(sstm.str());
        for(int j = 0; j < test_tree.nleaf; ++j){
            for(int i = 0; i < results[j].mutations.size(); ++i){
                int chr = results[j].mutations[i].chr;
                out_mut << j + 1 << "\t" << results[j].mutations[i].edge_id + 1
        << "\t" << results[j].mutations[i].type << "\t" << results[j].mutations[i].btime << "\t" << results[j].mutations[i].gtime << "\t" << chr + 1 << "\t" << (chr + 1)%22 << "\t" << results[j].mutations[i].seg << endl;
            }
        }
        out_mut.close();
        sstm.str("");
    }

    if(cons){
        write_sample_times(sstm, test_tree, dir, prefix, age);
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
    if(failed_mut > 0){
        assert(total_mut == succ_mut + failed_mut);
        cout << "The total number of successful mutations is: " << succ_mut << endl;
        cout << "The total number of failed mutations (due to limitations on copy number) is: " << failed_mut << endl;
    }
    cout << "The total branch length is: " << total_blen << endl;
    cout << "The estimated mutation rate per segment per year is: " << mu_est << endl;
}


void run_simulations(string tree_file, int mode, int method, const vector<int>& chr_lengths, int num_seg, int Ns, int Nsims, int cn_max, int model, int cons, const ITREE_PARAM& itree_param, double delta_t, int age, const vector<double>& rate_consts, const SV_SIZE& sv_size, string dir, string prefix, PRINT_LEVEL print_level, gsl_rng* r, int debug = 0){
    genome germline(chr_lengths, NORM_PLOIDY);

    string orig_prefix = prefix;
    cout << "\nNumber of datasets to simulate " << Nsims << endl;

    for(int i = 0; i < Nsims; ++i){
        vector<genome> genomes;      // one genome for each node on the tree
        map<int, vector<mutation>> muts;   // hold all mutations by edge id
        map<int, int> failed_muts;  // the number of failed mutations on each edge

        cout << "\nSimulation " << i + 1 << endl;

        //cout << "\n\n###### New sample collection ######" << endl;
        if(orig_prefix.empty()){
          prefix = "sim-data-" + to_string(i + 1);
        }
        cout << "Prefix of output file: " << prefix <<endl;

        evo_tree test_tree;
        if(tree_file != ""){
            test_tree = read_tree_info(tree_file, Ns);
        }else{
            test_tree = generate_random_tree(Ns, r, fp_myrng, age, itree_param, delta_t, cons, debug);
        }
        if(debug){
          test_tree.print();
          cout << "Simulated tree height " << get_tree_height(test_tree.get_node_times()) << endl;
        }

        if(method == SIM_TIME){     // applicable for all models
            simulate_samples(genomes, muts, failed_muts, test_tree, germline, chr_lengths, rate_consts, sv_size, model, cn_max, debug);
            int num_seg = get_num_seg(chr_lengths);
            print_simulations(mode, num_seg, genomes, chr_lengths, muts, failed_muts, test_tree, dir, prefix, age, cn_max, print_level, r, cons, debug);
        }else{  // model can only be BOUNDA or BOUNDT
            assert(model == BOUNDA || model == BOUNDT);
            int root = Ns + 1;
            int nstate = cn_max + 1;
            if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;
            // cout << "number of states now " << nstate << "\t" << cn_max << endl;

            map<int, copy_number> cn_matrix; // copy number matrix for each node
            map<int, int> num_muts; // number of mutations on each edge
            copy_number init_cn = initialize_cn(chr_lengths, model);
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

            evolve_sequences(cn_matrix, root, test_tree, qmat, nstate, num_seg, rate_consts, num_muts, debug);
            print_sequences(cn_matrix, cn_max, model, num_seg, num_muts, test_tree, dir, prefix, age, print_level, r, cons);

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
    int debug;
    unsigned seed;

    int Ns;   // number of multi-region samples
    int Nsims;  // number of simulations

    int age;
    int model;
    int mode, method;
    int cons;

    int seg_max, fix_nseg;
    int cn_max;

    /*************************************/
    // five event types: duplication, deletion, chromosome gain, chromosome loss, wgd
    // rates are 1/mean
    double dup_rate, del_rate, chr_gain, chr_loss, wgd;
    // mean of exp size distributions for site duplication/deletion
    double mean_dup_size;
    double mean_del_size;

    int Ne;     // effective population size
    double beta, gtime;    // population growth rate
    double delta_t;    // relative timing difference
    // string stime;

    int print_allele, print_mut, print_nex, print_relative, print_baseline;
    string tree_file;
    string dir; // output directory
    string prefix; // prefix of output file

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

      ("mode", po::value<int>(&mode)->default_value(0), "running mode of the program (0: simulating genome in fix-sized bins (4401 bins of size 500 Kbp by default), 1: simulating genome in segments of variable size, 2 to 4: test)")
      ("method", po::value<int>(&method)->default_value(0), "method of simulation (0: simulating waiting times, 1: simulating sequences directly)")
      ("model,d", po::value<int>(&model)->default_value(0), "model of evolution (0: Mk, 1: one-step bounded (total), 2: one-step bounded (haplotype-specific), 3: infinite sites)")

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
      ("fix_nseg", po::value<int>(&fix_nseg)->default_value(1), "whether or not to fix the number of segments to simulate. If not, the number of segments (<= seg_max) is randomly chosen")

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
      ("print_allele", po::value<int>(&print_allele)->default_value(1), "whether or not to output haplotype-specific copy numbers")
      ("print_mut", po::value<int>(&print_mut)->default_value(1), "whether or not to output the list of mutations")
      ("print_nex", po::value<int>(&print_nex)->default_value(1), "whether or not to output the tree in NEXUS format")
      ("print_relative", po::value<int>(&print_relative)->default_value(0), "whether or not to print relative copy numbers, similar to output of CGHcall: -2 (double deletion), -1 (single deletion), 0 (normal), 1 (gain), 2 (double gain)")
      ("print_baseline", po::value<int>(&print_baseline)->default_value(0), "whether or not to print relative copy numbers obtained from baseline strategy as used in PISCA")

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
            cout << "cnets [version 0.1], a program to simulate copy number alterations along a phylogenetic tree" << endl;
            return 1;
        }
        po::notify(vm);
    }catch(const std::exception& e){
        std::cerr << e.what() << std::endl;
        return 1;
    }

    if(delta_t > SMALL_VAL){
        cout << "Simulating time differences among samples with baseline time difference " << delta_t << endl;
    }

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    setup_rng(r, seed);

    fp_myrng = &myrng;

    PRINT_LEVEL print_level{print_allele, print_mut, print_nex, print_relative, print_baseline};

    // output directory
    if(dir.back() != '/'){
        dir = dir + "/";
    }

    vector<int> chr_lengths = CHR_BIN_SIZE;  // use 4401 bins as in example data by default (mode 0)
    int num_seg = 0;

    if(method == SIM_SEQ) {   // when simulating sequences directly, each site is a final segment
        cout << "Simulating sequences directly according to the substitution model of copy number changes at segment level.\nEach site is a final segment in this mode of simulation" << endl;
        mode = 1;
    }else{
        cout << "\nSimulating waiting times" << endl;
    }

    if(mode == 1){
        cout << "Under this mode, each site is treated as a final segment.\nThe mean duplication/deletion size is currently fixed to be 1" << endl;
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
        // due to rounding issues, actual number of segments may be larger than specified number of segments
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
    }else{    // When simulating waiting times, only simulate haplotype-specific model
      if(model == BOUNDT){
          cout << "When simulating waiting times, only haplotype-specific model is allowed" << endl;
          model = BOUNDA;
      }
    }

    vector<double> rate_consts = {dup_rate, del_rate, chr_gain, chr_loss, wgd};

    cout << "\nEvolution model simulated: " << model << endl;
    cout << "\tMaximum copy number of a segment is " << cn_max << endl;
    cout << "\tMutation rates:\t" << rate_consts[0] << "\t" << rate_consts[1]  << "\t" << rate_consts[2]  << "\t" << rate_consts[3]  << "\t" << rate_consts[4] << endl;
    cout << "\tSizes of site duplication/deletion:\t" << mean_dup_size << "\t" << mean_del_size << endl;

    ITREE_PARAM itree_param{Ne, beta, gtime};
    SV_SIZE sv_size{mean_dup_size, mean_del_size};

    // simulate coalescent tree and apply SVs
    if(mode < 2){
        if(beta > 0){
            cout << "\nSimulating exponential growth" << endl;
        }
        run_simulations(tree_file, mode, method, chr_lengths, num_seg, Ns, Nsims, cn_max, model, cons, itree_param, delta_t, age, rate_consts, sv_size, dir, prefix, print_level, r, debug);
    }else{
        cout << "\nRunning test" << endl;
        run_test(mode, dir, seed, itree_param, sv_size, debug);
    }

    gsl_rng_free(r);
}
