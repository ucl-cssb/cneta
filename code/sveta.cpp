// to build the program:
//   g++ sveta.cpp -o sveta -L/usr/local/lib/ -lgsl -I/usr/local/include
//
// to run the model:
//   ./sveta


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
static const int arr[] = {367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};
vector<int> chr_lengths (arr, arr + sizeof(arr) / sizeof(arr[0]) );

// rates: The rate of all events (segment duplication, segment deletion)
double get_site_rates(vector<double>& rates, const vector<double>& rate_constants, int model, int state){
    double rate;

    if(model == 0){ // JC69
        // assume duplication rate equals deletion rate
        assert(rate_constants[0] == rate_constants[1]);
        for(int i=0; i<2; ++i){
          rates.push_back(rate_constants[i]);
        }
        rate = 4 * rate_constants[0] / 5;
    }
    // if(model == 1){ // bounded
    //     // assume duplication rate equals deletion rate
    //     assert(rate_constants[0] == rate_constants[1]);
    //
    //     if(state == 0 || state == CN_MAX ){
    //         rate =  rate_constants[0] / 5;
    //         for(int i=0; i<rate_constants.size(); ++i){
    //           rates.push_back(rate_constants[0]);
    //         }
    //     }
    //     else{
    //         rate =  2 * rate_constants[0] / 5;
    //         for(int i=0; i<2; ++i){
    //           rates.push_back(rate_constants[i]);
    //         }
    //     }
    // }
    if(model == 1){ // bounded
        double dup_rate = rate_constants[0];
        double del_rate = rate_constants[1];
        if(state == 0){
            rate =  0;
        }
        else if(state == CN_MAX){
            rate =  2 * CN_MAX * del_rate;
        }
        else{
            rate =  2 * state * dup_rate + 2 * state * del_rate;
        }
    }
    return rate;
}


// rates: The total rate at each site in the genome
double get_total_rates(genome& g, vector<double>& rates, const vector<double>& rate_constants, int model){
    // Tranverse all the sites in the genome to get their rates
    for(int i=0; i < g.cn_profile.size(); ++i){   // For each chromosome
        for(int j=0; j < g.cn_profile[i].size(); ++j){    // For each segment in the chromosome
            // Get the copy number of this segment
            int state = g.cn_profile[i][j];
            vector<double> site_rates;
            double rate = get_site_rates(site_rates, rate_constants, model, state);
            // cout << "Rate on chr " << i << " seg " << j << " state " << state << " is " << rate << endl;
            rates.push_back(rate);
        }
    }
    // Get the total rate of substitution (the sum of rates accross sites)
    double rate = accumulate(rates.begin(), rates.end(), 0.0);
    if(debug){
        cout << "The total rate at all sites in the genome: " << rate << endl;
    }
    return rate;
}


// Note: the sites is counted for all the haplotypes. The positions on different haplotypes are counted as one site.
void site2chr(int site, int& chr, int& seg){
    if(site >= 0 && site < chr_lengths[0]){
        chr = 0;
        seg = site;
    }
    else{
        // cout << "There are " << chr_lengths.size() << " chromosomes" << endl;
        for(int i=0; i<chr_lengths.size(); i++){
            int sum1 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 1, 0);
            // cout << "sum until " << i + 1 << " chromosome is " << sum1 << endl;
            int sum2 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 2, 0);
            // cout << "sum until " << i + 2 << " chromosome is " << sum2 << endl;
            if(site > sum1 && site < sum2){
                chr = i + 1;
                seg = site - sum1;
                break;
            }
        }
    }
}

// Find the number of a specific chromosme
int get_ploidy_at_chr(genome& g, int c){
    int ploidy_at_c = 0;
    for(int i=0; i<g.chrs.size(); i++){
        if(i%22 == c){
            ploidy_at_c += 1;
        }
    }
    return ploidy_at_c;
}

// Simulate mutations under different models
vector<mutation> generate_mutation_by_model(genome& g, const int& edge_id, const double& blength, const double& node_time, const vector<double>& rate_constants, int model){
  vector<mutation> ret;
  if(debug) cout << "\tgenerate_mutations, blength:" << "\t" << blength << endl;
  vector<double> site_rates;
  vector<double> state_rates;

  // mean variant size in bins
  double mdup = g.mean_dup_size;
  double mdel = g.mean_del_size;

  // Find the current state of the genome
  g.calculate_cn();

  double rate = get_total_rates(g, site_rates, rate_constants, model);
  // cout << "Total mutation rate on the genome " << rate << endl;

  // Simulate the waiting times of a Markov chain
  double time = 0.0;
  while( time < blength ){
    double tevent = gsl_ran_exponential (r, 1/rate);
    time += tevent;

    // choose type of event
    int e = rchoose(r, rate_constants);

    // assign the event to a site with probabilities proportional to the rates at the site
    int site = rchoose(r,  site_rates);
    // Convert the site position back to chromosome ID and segment ID
    int c = 0;
    int loc = 0;
    site2chr(site, c, loc);
    // Randomly choose a haplotype
    double x = runiform(r, 0, 1);
    // Find the number of possible haplotype of this chromosme
    int ploidy_at_c = get_ploidy_at_chr(g, c);
    int haplotype = myrng(ploidy_at_c);
    if(debug){
        cout << "There are " << g.chrs.size() << " chromosomes now" << endl;
        cout << "site " << site << " is at chr " << c << " pos " << loc << " haplotype " << haplotype << " ploidy " << ploidy_at_c << endl;
    }
    c += 22 * haplotype;

    if(debug) cout << "mut times, tevent, total time, time/branch len, event, chr, loc\t" << tevent << "\t" << time << "\t" << blength << "\t" << e << "\t" << c << "\t" << loc << endl;

    // Generate the mutation
    mutation mut(edge_id, e, time/blength, node_time+time);
    ret.push_back(mut);
    g.nmuts[e]++;
    g.mutations.push_back(mut);

    if(e == 0){   //duplications
      if(debug){
        cout << "GENOME BEFORE duplication " << endl;
        g.print();
        g.print_cn();
      }

      if( g.chrs[c].size() > 0 && loc < g.chrs[c].size()){
        bool done = false;
        while( done == false ){
          int len = gsl_ran_exponential(r, mdup);
          // int len = 0;
          if(debug){
              cout << "before dup, Chr " << c << " has " << g.chrs[c].size() << " segments" << endl;
              for(int i=0; i<g.chrs[c].size(); i++){
                  cout << "\t" << g.chrs[c][i].seg_id;
              }
              cout << endl;
              cout << "dup len:" << len << endl;
          }

          if( loc + len < g.chrs[c].size() ){
            if(debug) cout << "\tSV: duplicating segment, chr, start, len: " << c << "\t" << loc << "\t" << len+1 << endl;
            //cout << "SV: insert: " << loc+len+1 << "\t" << loc << "\t" << loc+len+1 << endl;
            // Find the number of copies
            int state = g.cn_profile[c%22][g.chrs[c][loc].seg_id];
            vector<int> possible_states;
            if(model == 0){
                for (int i = state + 1; i <= CN_MAX; i++){
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

                // assume tandom duplication
                for(int j=0; j< pstate-state; j++){
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
                int state = g.cn_profile[c%22][g.chrs[c][loc].seg_id];
                if(state < CN_MAX ){
                    // tandom duplication
                    g.chrs[c].insert(g.chrs[c].begin()+loc+len+1, g.chrs[c].begin()+loc, g.chrs[c].begin()+loc+len+1);
                }
            }

            if(debug){
                g.calculate_cn();
                cout << "Previous state: " << state << endl;
                cout << "Current state: " << g.cn_profile[c%22][g.chrs[c][loc].seg_id] << endl;
                set<double> segs;
                for(int i=0; i<g.chrs[c].size(); i++){
                    segs.insert(g.chrs[c][i].seg_id);
                }
                cout << "After duplication, Chr " << c << " has " << g.chrs[c].size() << " segments (unique: " << segs.size() << ")" << endl;
            }
            done = true;
          }
        }
     }else{
         if(debug) cout << "\tSV: duplication failed:" << endl;
     }

      if(debug){
        cout << "GENOME AFTER duplication " << endl;
        g.print();
        g.print_cn();
      }
    }
    else if( e == 1 ){   //deletions
      if(debug){
        cout << "GENOME BEFORE deletion " << endl;
        g.print();
        g.print_cn();
      }
      // make sure a deletion somewhere is possible
      if( g.chrs[c].size() > 0 && loc < g.chrs[c].size()){
        bool done = false;
        while( done == false ){
          int len = (int) gsl_ran_exponential(r, mdel);
          // int len = 0;
          if(debug){
              cout << "before del, Chr " << c << " has " << g.chrs[c].size() << " segments" << endl;
              for(int i=0; i<g.chrs[c].size(); i++){
                  cout << "\t" << g.chrs[c][i].seg_id;
              }
              cout << endl;
              cout << "del len:" << len << endl;
          }

          if( len + loc < g.chrs[c].size() ){
            // Find the segment ID
            int state = g.cn_profile[c%22][g.chrs[c][loc].seg_id];
            // Find the number of segments in other haplotypes
            int hap_state = 0;
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
            if(model == 0){
                for (int i = 0; i <= state - hap_state - 1; i++){
                    possible_states.push_back(i);
                }
                // Randomly pick a possible state
                int ssize = possible_states.size();
                if(debug){
                    cout << "There are " << ssize << " possible state(s) for deletion" << endl;
                }
                if(ssize <= 0){
                    cout << "Impossible to do deletion on Chr " << c << endl;
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
                int state = g.cn_profile[c%22][g.chrs[c][loc].seg_id];
                if(state > 0){
                    g.chrs[c].erase(g.chrs[c].begin()+loc, g.chrs[c].begin()+loc+len+1);
                }
            }

            if(debug){
                g.calculate_cn();   // Update copy number profile
                cout << "Previous state: " << state << endl;
                cout << "Current state: " << g.cn_profile[c%22][g.chrs[c][loc].seg_id] << endl;
                set<double> segs;
                for(int i=0; i<g.chrs[c].size(); i++){
                    segs.insert(g.chrs[c][i].seg_id);
                }
                cout << "After deletion, Chr " << c << " has " << g.chrs[c].size() << " segments (unique: " << segs.size() << ")" << endl;
            }
            done = true;
          }
        }
      }else{
          if(debug) cout << "\tSV: deletion failed:" << endl;
      }

      if(debug){
        cout << "GENOME AFTER deletion " << endl;
        g.print();
        g.print_cn();
      }
    }
    else if( e == 2 ){   //chr loss
      if( g.chrs.size() > 0 ){
        int orig_num_chr = g.chrs.size();
        // Delete the segments in the chromosme, but keep the chromosme ID so that it can be remapped by 22
        g.chrs[c].erase(g.chrs[c].begin(), g.chrs[c].end());
        if(debug){
            cout << "Chromosome loss in " << c << endl;
            cout << "There are " << orig_num_chr << " chromosomes before chr loss" << endl;
            cout << "There are " << g.chrs.size() << " chromosomes now" << endl;
            g.print_cn();
        }
      }
    }
    else if( e == 3 ){   //chr gain
      if( g.chrs.size() > 0 && g.chrs[c].size() > 0){
        int orig_num_chr = g.chrs.size();
        g.chrs.insert(g.chrs.end(), g.chrs.begin()+c, g.chrs.begin()+c+1);
        g.chrs[orig_num_chr].insert(g.chrs[orig_num_chr].begin(), g.chrs[c].begin(), g.chrs[c].end());
        if(debug){
            cout << "Chromosome gain in " << c << endl;
            cout << "There are " << orig_num_chr << " chromosomes before chr gain" << endl;
            cout << "There are " << g.chrs.size() << " chromosomes now" << endl;
            g.print_cn();
            g.print_cn();
        }
      }
    }
    else if( e == 4 ){   // whole genome duplication
      if( g.chrs.size() > 0 ){
          int orig_num_chr = g.chrs.size();
          g.chrs.insert(g.chrs.end(), g.chrs.begin(), g.chrs.end());
          // replicate the segments for each chromosme
          for(int i=0; i<orig_num_chr; i++){
               g.chrs[i + orig_num_chr].insert(g.chrs[i + orig_num_chr].begin(), g.chrs[i].begin(), g.chrs[i].end());
          }
          if(debug){
              cout << "Whole genome doubling" << endl;
              cout << "There are " << orig_num_chr << " chromosomes before WGD" << endl;
              cout << "There are " << g.chrs.size() << " chromosomes now" << endl;
              g.print_cn();
              // g.print();
          }
      }
    }
    else{
      cerr << "Unknown mutation type " << e << endl;
      exit(1);
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

void traverse_tree_mutating(const int& node_id, const evo_tree& tree, const vector<double>& rate_constants,
        map<int, vector<mutation>>& all_muts, vector<genome>& genomes, int model){
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
        if(model == 2){
            muts = generate_mutation_times(edge_id, tree.edges[edge_id].length, tree.node_times[ tree.nodes[node_id].parent ], rate_constants);
            apply_mutations( muts, genomes[ node_id ] );
        }
        else{
            muts = generate_mutation_by_model(genomes[ node_id ], edge_id, tree.edges[edge_id].length, tree.node_times[ tree.nodes[node_id].parent], rate_constants, model);
        }

      if(debug){
          cout << "Generating mutations for node " << node_id << endl;
      }
    }
    all_muts[edge_id] = muts;
  }

  if( tree.nodes[ node_id ].daughters.size() == 0 ){
    // we are done
    return;
  }else{
    traverse_tree_mutating(tree.nodes[ node_id ].daughters[0], tree, rate_constants, all_muts, genomes, model );
    traverse_tree_mutating(tree.nodes[ node_id ].daughters[1], tree, rate_constants, all_muts, genomes, model );
  }
  return;
}

void simulate_samples(vector<genome>& genomes, map<int,vector<mutation> >& muts, const evo_tree& tree, genome& germline, const vector<double>& rate_constants, int model, bool print = true ){
  // assign the germline to the root of the tree
  germline.node_id = tree.root_node_id;

  for(int i=0; i<(tree.nnode+tree.nleaf); ++i){
    if(i == tree.root_node_id){
      genomes.push_back( germline );
    }
    else{
      genomes.push_back( genome() );
    }
  }

  // move through the evolutionary tree mutating genomes
  traverse_tree_mutating( tree.root_node_id, tree, rate_constants, muts, genomes, model );

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

/*
//void run_sample_set(int Ns, double* prc, double* pvs, double* ptree, double* pl, int* ret){
void run_sample_set(int Ns, double* prc, double* pvs, int* ret){

  static const int arr[] = {367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};
  vector<int> chr_lengths (arr, arr + sizeof(arr) / sizeof(arr[0]) );

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
  simulate_samples(results, muts, test_tree, germline, rate_consts, false);

  int nbins = 4401;
  for(int i=0; i<(test_tree.nleaf-1); ++i){
    vector<int> cns = results[i].get_cn_vector();
    for(int j=0; j<nbins; ++j){
      ret[nbins*i + j] = cns[j];
    }
  }

  if(0){
    sstm << "test-data-cn.txt";
    ofstream out_cn(sstm.str());
    for(int j=0; j<test_tree.nleaf; ++j){
      results[j].write(out_cn);
    }
    out_cn.close();
    sstm.str("");
  }
  return;
}
*/

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
    // four event types: duplication, deletion, chromosome gain, chromosome loss, wgd
    // rates are 1/mean
    double dup_rate, del_rate, chr_gain, chr_loss, wgd;
    // parameters for mean of dup/del size distributions
    double dup_size, del_size;
    // effective population size
    double Ne;
    // relative timing difference
    double delta_t;
    int seed;
    int model;
    int cons;

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
      ("model,d", po::value<int>(&model)->default_value(0), "model of evolution (0: JC69, 1: 1-step bounded, 2: Poisson)")
      ("age,a", po::value<int>(&age)->default_value(100), "age of the patient to simulate")
      ("epop,e", po::value<double>(&Ne)->default_value(2), "effective population size")
      ("tdiff,t", po::value<double>(&delta_t)->default_value(1), "relative timing difference")
      ("constrained", po::value<int>(&cons)->default_value(1), "constraints on branch length (0: none, 1: fixed total time)")

      ("nregion,r", po::value<int>(&Ns)->default_value(5), "number of regions")
      ("nsim,n", po::value<int>(&Nsims)->default_value(1), "number of multi-region samples")
      ("prefix,p", po::value<string>(&prefix)->default_value(""), "prefix of output file (it will be sim-data-N if not specified")

      ("dup_rate", po::value<double>(&dup_rate)->default_value(0.0005), "duplication rate (allele/locus/year)")
      ("del_rate", po::value<double>(&del_rate)->default_value(0.01), "deletion rate (allele/locus/year)")
      ("chr_gain", po::value<double>(&chr_gain)->default_value(0), "chromosome gain rate (haplotype/chr/year)")
      ("chr_loss", po::value<double>(&chr_loss)->default_value(0), "chromosome loss rate (haplotype/chr/year)")
      ("wgd", po::value<double>(&wgd)->default_value(0), "WGD (whole genome doubling) rate (year)")
      ("dup_size", po::value<double>(&dup_size)->default_value(30), "mean of duplication size distributions")
      ("del_size", po::value<double>(&del_size)->default_value(30), "mean of deletion size distributions")

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

  int mode = 0;

  // output directory
  if(dir.back() != '/'){
      dir = dir + "/";
  }
  vector<double> rate_consts = {dup_rate, del_rate, chr_gain, chr_loss, wgd};
  vector<double> var_size = {dup_size, del_size};

  // priors on rates
  //vector<double> pars;
  //read_params(argv[4],pars);

  cout << "rates:\t" << rate_consts[0] << "\t" << rate_consts[1]  << "\t" << rate_consts[2]  << "\t" << rate_consts[3]  << "\t" << rate_consts[4] << endl;
  cout << "sizes:\t" << var_size[0] << "\t" << var_size[1] << endl;

  // simulate coalescent tree and apply SVs
  if(mode == 0){
    // can specify the germline in a number of different ways
    //genome germline(1,10);
    //genome germline(42,1000);
    //genome germline(chr_lengths);

    genome germline(chr_lengths,2);
    germline.mean_dup_size = var_size[0];
    germline.mean_del_size = var_size[1];

    // if(model==1){
    //     cout << "Under this model, the mean duplication/deletion size is currently fixed to be 1" << endl;
    // }

    vector<int> edges;
    vector<double> lengths;
    vector<double> epoch_times;
    vector<double> node_times;
    stringstream sstm;
    vector<genome> results;
    map<int, vector<mutation> > muts;   // hold all mutations by edge id
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

      // Scale branch lengths
      for(int l=0; l<lengths.size(); ++l) lengths[l] = lengths[l]*Ne;

      // randomly assign leaf edges to time points t0, t1, t2, t3
      //
      cout << "assigning temporal structure" << endl;
      bool assign0 = false;
      int bcount = 0;
      for(int l=1; l<edges.size(); l=l+2){
        	if( edges[l] < Ns){
        	  int ind = 0;
        	  if(assign0 == false){
        	    ind = 0;
        	    assign0 = true;
        	  }
        	  else{
        	    ind = gsl_rng_uniform_int(r, 4);
        	  }
        	  cout << "\t sample / time point:" << edges[l]+1 << "\t" << ind << endl;
              double stime = ind*delta_t;
        	  lengths[bcount] = lengths[bcount] + stime;
              // Store the sampling time for time scaling by the age of a patient
              tobs.push_back(stime);
        	}
        	bcount++;
      }

      evo_tree test_tree(Ns+1, edges, lengths);
      if(cons){
          double min_height = *max_element(tobs.begin(), tobs.end());
          double max_height = age + min_height;
          double tree_height = runiform(r, min_height, max_height);
          double old_tree_height = test_tree.get_tree_height();
          double ratio = tree_height/old_tree_height;
          cout << "Tree height before scaling " << old_tree_height << endl;
          cout << "Scaling ratio " << ratio << endl;
          test_tree.scale_time(ratio);
          cout << "Rescaled tree height " << tree_height << endl;
          // scaling tobs accordingly
          cout << "Rescaled sampling time point:";
          for(int i = 0; i< tobs.size(); i++)
          {
              tobs[i] = tobs[i] * ratio;
              cout << "\t" << tobs[i];
          }
          cout << endl;
      }
      simulate_samples(results, muts, test_tree, germline, rate_consts, model, false);

      sstm << dir << prefix << "-cn.txt.gz";
      //ofstream out_cn(sstm.str());
      ogzstream out_cn(sstm.str().c_str());
      for(int j=0; j<test_tree.nleaf; ++j){
          results[j].write(out_cn);
      }
      out_cn.close();

      // re-read the cn data to get the segments
      vector<vector<int> > seg_data = read_data_var_regions(sstm.str().c_str(), Ns, 4);
      int Nchar = seg_data.size();
      cout << "#### Nchar: " << Nchar << endl;
      sstm.str("");

      sstm << dir << prefix << "-info.txt";
      ofstream out_info(sstm.str());
      out_info << "NODE TIMES:";
      out_info << "\tnid\ttime" << endl;
      for(int j=0; j< test_tree.node_times.size(); ++j){
          out_info << "\t" << j+1 << "\t" << test_tree.node_times[j] << endl;
      }
      out_info << endl;
      out_info << "SEGMENTS: " << Nchar << endl;
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
      ofstream out_tree(sstm.str());
      //test_tree.write(out_tree);
      out_tree << "start\tend\tlength\teid\tnmut" << endl;
      for(int i=0; i<test_tree.nedge; ++i){
          out_tree << test_tree.edges[i].start+1 << "\t" << test_tree.edges[i].end+1 << "\t" << test_tree.edges[i].length
     << "\t" << test_tree.edges[i].id+1 << "\t" << muts[test_tree.edges[i].id].size() << endl;
      }
      out_tree.close();
      sstm.str("");

      sstm << dir << prefix << "-tree.nex";
      ofstream nex_tree(sstm.str());
      int precision = 5;
      test_tree.write_nexus(precision, nex_tree);
      nex_tree.close();
      sstm.str("");

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

      sstm << dir << prefix << "-rel-times.txt";
      double node_min = 1000;
      for(int j=0; j<Ns; ++j){
          if( test_tree.node_times[j] < node_min ) node_min = test_tree.node_times[j];
      }
      //cout << "leaf minimum: " << node_min << endl;
      ofstream out_rel(sstm.str());
      if(cons){
          for(int j=0; j<Ns; ++j){
              double delta = test_tree.node_times[j] - node_min;
              out_rel << j+1 << "\t" << delta << "\t" << ceil(age + delta) << endl;
          }
      }
      else{
          for(int j=0; j<Ns; ++j){
              out_rel << j+1 << "\t" << test_tree.node_times[j] - node_min << endl;
          }
      }
      out_rel.close();
      sstm.str("");

      edges.clear();
      lengths.clear();
      results.clear();
      muts.clear();
      epoch_times.clear();
      node_times.clear();

    }

  }

  /*
  //create test tree
  if(mode == 1){
    //static const int arr1[] = {7,8, 6,7, 8,1, 8,2, 7,9, 9,3, 9,4, 6,5 };
    static const int arr1[] = {6,7, 5,6, 7,0, 7,1, 6,8, 8,2, 8,3, 5,4 };

    vector<int> e (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );

    static const double arr2[] = {5.1,6.3,10.2,9.5,5.2,3.2,5.4,0};
    vector<double> l (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

    evo_tree test_tree(5, e, l);

    test_tree.print();
  }

  // test out the genome classes and functionality
  if(mode == 2){

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
  if(mode == 3){
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
  */
}
