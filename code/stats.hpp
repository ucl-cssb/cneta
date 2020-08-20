#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>

#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <set>
#include <climits>

#include "lbfgsb/lbfgsb_new.h"
#include "matexp/matrix_exponential.hpp"
#include "matexp/r8lib.hpp"
#include <boost/numeric/ublas/io.hpp>

#include "utilities.hpp"

typedef std::numeric_limits< double > dbl;

// // Variables related to likelihood computation
// extern "C" {
//   gsl_rng * r;
//   // vector<vector<int>> vobs;
//   map<int, vector<vector<int>>> vobs;   // CNP for each site
//   vector<double> tobs;
//   vector<int> obs_num_wgd;  // possible number of WGD events
//   vector<vector<int>> obs_change_chr;
//   vector<int> sample_max_cn;
//
//   map<int, set<vector<int>>> decomp_table;  // possible state combinations for observed copy numbers
//   set<vector<int>> comps;
//
//   // vector<vector<vector<double>>> L_sk_k_mrca;  // Likelihood table for the MRCA of all sites
//
//   int model;
//   int correct_bias;
//   int cn_max;
//   int only_seg; // Whether or not to only consider segment-level mutations
//   int is_total; // whether or not the input is total copy number
//   int infer_wgd; // whether or not to infer WGD status of a sample
//   int infer_chr; // whether or not to infer chromosome gain/loss status of a sample
//   int use_repeat;   // whether or not to use repeated site patterns
//
//   // Set max_* as global variables to avoid adding more parameters in maximization
//   int m_max;
//   int max_wgd;
//   int max_chr_change;
//   int max_site_change;
//
//   int Ns;   // number of samples
//   int Nchar;  // number of sites
//   int num_invar_bins;   // number of invariant sites
//   int age;
//
//   int cons;
//   int maxj;
// }

const double ERROR_X = 1.0e-4;
const int MAX_AGE = 100;
const double LARGE_LNL = -1e9;
const double SMALL_LNL = -1e20;
const double SMALL_VAL = 1.0e-10;
// The maximum mutation rates allowed
const double MAX_MRATE = 1;
// The maximum age ratio allowed
const double MIN_RATIO = 1e-2;
const double MAX_RATIO = 0.99;
// The shortest branch length allowed (in year)
const double SMALL_BLEN = 1e-9;
// Scaling tree height to 1/HEIGHT_SCALE if the current height is larger than the upper bound (patient age at last sample)
const int HEIGHT_SCALE = 3;
// The difference from minmial height
const int HEIGHT_OFFSET = 10;

int opt_one_branch = 0;

// Pairing function: see https://en.wikipedia.org/wiki/Pairing_function
int pairInteger(int int1, int int2) {
    if (int1 <= int2) {
        return ((int1 + int2)*(int1 + int2 + 1)/2 + int2);
    } else {
        return ((int1 + int2)*(int1 + int2 + 1)/2 + int1);
    }
}

void setup_rng(int set_seed){
  gsl_rng_env_setup();

  const gsl_rng_type* T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  if( set_seed != 0 ){
    gsl_rng_set(r, set_seed);
  }else{
    int t = time(NULL);
    int pid = getpid();
    long s = t*pid;
    //cout << "pid:" << "\t" << getpid() << endl;
    // cout << "seed:" << "\t" << t << "\t" << pid << "\t" << abs(s) << endl;
    cout << "seed:" << "\t" << abs(s) << endl;
    gsl_rng_set (r, abs(s) );
  }
}

// unary function and pointer to unary function
// allows use of gsl rng for standard template algorithms
long unsigned myrng(long unsigned n)
{
  return gsl_rng_uniform_int(r, n);
}

long unsigned (*fp)(long unsigned) = myrng;

// factorial for choose(n,k)
int fact(int n){
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}

// wrapper for uniform
double runiform(gsl_rng* r, double a, double b){
  double myrandom = a + (b-a)*gsl_rng_uniform (r);

  while(myrandom==0){
    myrandom = a + (b-a)*gsl_rng_uniform (r);
  }
  return myrandom;
}


// sample an element according to a probability vector
// gsl_ran_multinomial?
int rchoose(gsl_rng* r, const vector<double>& rates){
  //cout << "rchoose, rates:";
  //for(int i=0; i<rates.size(); ++i) cout << "\t" << rates[i];
  //cout << endl;

  vector<double> p;
  double s = accumulate(rates.begin(), rates.end(), 0.0);

  for(int i=0; i<rates.size(); ++i){
    p.push_back(rates[i]/s);
  }

  vector<double> psum(p.size(),0.0);
  partial_sum(p.begin(),p.end(),psum.begin());

  double u = gsl_rng_uniform(r);
  int ret = -1;
  for(int i=0; i<rates.size(); ++i){
    if( u < psum[i] ){
      ret = i;
      break;
    }
  }

  //cout << "u=\t" << u << "\t" << ret << endl;
  return ret;
}


// generate neutral coalescent trees
// here nsample is the number of cancer samples (not including germline node)
// here we directly calculate the edges in the tree
void generate_coal_tree(const int& nsample, vector<int>& edges, vector<double>& lengths, vector<double>& epoch_times, vector<double>& times, int Ne, double beta = 0, double gtime=0.002739726){
  //cout << "GENERATING COAL TREE" << endl;
  int nlin = nsample;
  vector<int> nodes;

  // For compatibility with ape in R, root node must be labelled
  // 0 to nsample-1 are leaf, nsample is germline, nsample+1 is root
  // all other nodes ids start from here
  int node_count = nsample+2;

  // create leaf nodes
  for(int i=0; i<nsample; ++i) nodes.push_back(i);

  // create vector of event times
  // total nodes = 2*(nsample+1) -1
  for(int i=0; i<(2*nsample+1); ++i) times.push_back(0.0);

  double t_tot = 0;
  double curr_time = 0;
  while(nlin > 1){
    // sample a time from Exp( combinations(k,2) )
    double lambda = fact(nlin)/( 2*fact(nlin-2));
    double tn = gsl_ran_exponential(r, 1/lambda) * Ne;
    // cout << "Normal coalescence time  " << tn << endl;
    double t = tn;
    if(beta > 0){  // simulate exponential growth
        // t = (log(1 + beta * tn * exp(- beta * t_tot))) / beta;  // Formula in book "Gene Genealogies, Variation and Evolution", P99
        t = ((log(exp(beta*t_tot) + beta * tn)) / beta)- t_tot;    // Formula in CoalEvol7.3.5.c, equivalent to the one above
    }
    // cout << "Exponential coalescence time  " << t << endl;

    t_tot += t;
    // choose two random nodes from available list
    random_shuffle(nodes.begin(), nodes.end(),fp);

    // edge node_count -> node
    edges.push_back(node_count);
    edges.push_back(nodes[nodes.size()-1]);
    lengths.push_back(t_tot - times[ nodes[nodes.size()-1] ] );

    edges.push_back(node_count);
    edges.push_back(nodes[nodes.size()-2]);
    lengths.push_back(t_tot - times[ nodes[nodes.size()-2] ] );

    // update time for this node
    //cout << "\t" << node_count+1 << "\t" << t_tot << endl;
    epoch_times.push_back(t_tot);
    times[ node_count ] = t_tot;

    nodes.pop_back();
    nodes.pop_back();

    nodes.push_back(node_count);
    node_count++;
    nlin--;
  }

  cout << "TMRCA of tumour samples (scaled by N): " << t_tot << endl;
  // create the root and germline nodes and edges
  double lambda = runiform(r, 0, 1);
  double tn = gsl_ran_exponential(r, 1/lambda) * Ne;
  double t = tn;
  if(beta > 0){  // simulate exponential growth
      // t = (log(1 + beta * tn * exp(- beta * t_tot))) / beta;
      t = ((log(exp(beta*t_tot) + beta * tn)) / beta)- t_tot;    // Formula in CoalEvol7.3.5.c
  }
  t_tot += t;
  epoch_times.push_back(t_tot);

  // add in the time for the root and germline nodes
  times[nsample+1] = t_tot;
  times[nsample] = t_tot;

  edges.push_back(nsample+1);
  edges.push_back(node_count-1);
  lengths.push_back(t);

  edges.push_back(nsample+1);
  edges.push_back(nsample);
  lengths.push_back(0);

  // invert the times
  for(int i=0; i<times.size(); ++i){
    times[i] = t_tot - times[i];
  }

  for(int l = 0; l < lengths.size(); ++l){
      lengths[l] = lengths[l] * gtime;
  }
  //cout << "total time of tree: " << t_tot << " : ";
  //for(int i=0; i<epoch_times.size(); ++i) cout << "\t" << epoch_times[i];
  //cout << endl;

  // invert the times
  //cout << "times of nodes:" << endl;
  //for(int i=0; i<times.size(); ++i){
    //cout << i+1 << "\t" << times[i] << endl;
  //}
}


// Scale the total time by given time
evo_tree generate_coal_tree(const int& nsample, int Ne = 1, double beta = 0, double gtime=0.002739726){
  //cout << "GENERATING COAL TREE" << endl;
  vector<int> edges;
  vector<double> lengths;
  vector<double> epoch_times;
  vector<double> times;

  int nlin = nsample;
  vector<int> nodes;

  // For compatibility with ape in R, root node must be labelled
  // 0 to nsample-1 are leaf, nsample is germline, nsample+1 is root
  // all other nodes ids start from here
  int node_count = nsample+2;

  // create leaf nodes
  for(int i=0; i<nsample; ++i) nodes.push_back(i);

  // create vector of event times
  // total nodes = 2*(nsample+1) -1
  for(int i=0; i<(2*nsample+1); ++i) times.push_back(0.0);

  double t_tot = 0;
  while(nlin > 1){
    // sample a time from Exp( combinations(k,2) )
    double lambda = fact(nlin)/( 2*fact(nlin-2) );
    double tn = gsl_ran_exponential(r, 1/lambda) * Ne;
    // cout << "Normal coalescence time  " << tn << endl;
    double t = tn;
    if(beta > 0){  // simulate exponential growth
        // t = (log(1 + beta * tn * exp(- beta * t_tot))) / beta;  // Formula in book "Gene Genealogies, Variation and Evolution", P99
        t = ((log(exp(beta*t_tot) + beta * tn)) / beta)- t_tot;    // Formula in CoalEvol7.3.5.c, equivalent to the one above
    }

    t_tot += t;

    // choose two random nodes from available list
    random_shuffle(nodes.begin(), nodes.end(),fp);

    // edge node_count -> node
    edges.push_back(node_count);
    edges.push_back(nodes[nodes.size()-1]);
    lengths.push_back(t_tot - times[ nodes[nodes.size()-1] ] );

    edges.push_back(node_count);
    edges.push_back(nodes[nodes.size()-2]);
    lengths.push_back(t_tot - times[ nodes[nodes.size()-2] ] );

    // update time for this node
    //cout << "\t" << node_count+1 << "\t" << t_tot << endl;
    epoch_times.push_back(t_tot);
    times[ node_count ] = t_tot;

    nodes.pop_back();
    nodes.pop_back();

    nodes.push_back(node_count);
    node_count++;
    nlin--;
  }

  // create the root and germline nodes and edges
  // double lambda = 1;
  double lambda = runiform(r, 0, 1);
  double tn = gsl_ran_exponential(r, 1/lambda) * Ne;
  // cout << "Normal coalescence time  " << tn << endl;
  double t = tn;
  if(beta > 0){  // simulate exponential growth
      // t = (log(1 + beta * tn * exp(- beta * t_tot))) / beta;  // Formula in book "Gene Genealogies, Variation and Evolution", P99
      t = ((log(exp(beta*t_tot) + beta * tn)) / beta)- t_tot;    // Formula in CoalEvol7.3.5.c, equivalent to the one above
  }
  t_tot += t;
  epoch_times.push_back(t_tot);

  // add in the time for the root and germline nodes
  times[nsample+1] = t_tot;
  times[nsample] = t_tot;

  edges.push_back(nsample+1);
  edges.push_back(node_count-1);
  lengths.push_back(t);

  edges.push_back(nsample+1);
  edges.push_back(nsample);
  lengths.push_back(0);

  // invert the times
  for(int i=0; i<times.size(); ++i){
    times[i] = t_tot - times[i];
  }

  for(int l = 0; l < lengths.size(); ++l){
      lengths[l] = lengths[l] * gtime;
  }
  // cout << "total time of tree: " << t_tot << " : ";
  // for(int i=0; i<epoch_times.size(); ++i) cout << "\t" << epoch_times[i];
  // cout << endl;
  //
  // // invert the times
  // cout << "times of nodes:" << endl;
  // for(int i=0; i<times.size(); ++i){
  //   cout << i+1 << "\t" << times[i] << endl;
  // }

  evo_tree ret(nsample+1, edges, lengths);
  return ret;
}


// Create a new tree with the same topology as input tree but different branch lengths
evo_tree create_new_tree(gsl_vector* blens, evo_tree &rtree, int cons){
    // create a new tree
    vector<edge> enew;
    // cout << "copy original edges" << endl;
    for(int i=0; i<rtree.nedge; ++i){
      enew.push_back( rtree.edges[i] );
    }
    if(cons){
        // cout << "branches constrained" << endl;
        int count = 0;
        for(int i=0; i<rtree.nedge-1; ++i){
            if(enew[i].end > Ns){
                // cout << "count " << count << endl;
                enew[i].length = gsl_vector_get(blens, count);
                count++;
            }else{
                enew[i].length = 0;
            }
        }

        // Time to first sample
        // cout << "time to 1st sample: " << total_time << endl;
        evo_tree new_tree(rtree.nleaf, enew, rtree.get_total_time(), rtree.tobs);
        new_tree.mu = rtree.mu;
        new_tree.dup_rate = rtree.dup_rate;
        new_tree.del_rate = rtree.del_rate;

        return new_tree;
    }
    else{
        // cout << "branches unconstrained" << endl;
        for(int i=0; i<rtree.nedge-1; ++i){
            enew[i].length = gsl_vector_get(blens,i);
        }

        evo_tree new_tree(Ns+1, enew);
        new_tree.mu = rtree.mu;
        new_tree.dup_rate = rtree.dup_rate;
        new_tree.del_rate = rtree.del_rate;

        return new_tree;
    }
}


// Check if there is negative branch length
bool is_blen_valid(evo_tree &rtree){
    for(int i=0; i<rtree.nedge; ++i){
        if(rtree.edges[i].length < 0){
            return false;
        }
    }
    return true;
}


// update edge ending at i when its time is updated
void update_edge_len(evo_tree &rtree, int i){
    for(int j = 0; j < rtree.nedge; j++){
        if (rtree.edges[j].end == i){
            rtree.edges[j].length = rtree.node_times[i] - rtree.node_times[rtree.edges[j].start];
        }
        // assert(rtree.edges[j].length > 0);
    }
}


// Adjust the time of sample1 so that it is not after all other tips. Then the subsequent adjustment will be addition, not introducing negative values
void adjust_sample1(evo_tree &rtree, vector<double>& tobs, int sample1){
    // Find the maximum time differences to the first sample
    double max_diff = 0;
    for(int i = 0; i < tobs.size(); i++){
        // current differences to first sample
        double delta = rtree.node_times[i] - rtree.node_times[sample1];
        // If delta is negative, there is no need to adjust node time of first sample
        if(delta > max_diff){
            max_diff = delta;
        }
    }
    // cout << "Sample 1 time before " << rtree.node_times[sample1] << endl;
    if(max_diff > 0){
        // double increament = max_diff - max_tobs;
        rtree.node_times[sample1] =  rtree.node_times[sample1] + max_diff;
        update_edge_len(rtree, sample1);
    }
    // cout << "Sample 1 time after " << rtree.node_times[sample1] << endl;
}


// Change the time and branch length for one tip
void adjust_one_tip(evo_tree &rtree, vector<double>& tobs, int i, int sample1){
    int debug = 0;
    for(int j = 0; j < rtree.nedge; j++){
        if (rtree.edges[j].end == i){
            if(debug){
                cout << "Adjust sample " << i << endl;
                cout << "Sample time before " << rtree.node_times[i] << endl;
            }
            rtree.node_times[i] =  rtree.node_times[sample1] + tobs[i];
            if(debug) cout << "Sample time after " << rtree.node_times[i] << endl;
            rtree.edges[j].length = rtree.node_times[i] - rtree.node_times[rtree.edges[j].start];

            break;
        }
    }
}


// Adjust all tips recursively so that tree height is smaller than age
void adjust_all_tips(evo_tree &rtree, vector<double>& tobs, int age){
    assert(is_blen_valid(rtree));
    int debug = 0;
    double ttime = rtree.get_total_time();
    if(ttime < age) return;
    if(debug) {
        // cout << "Age of patient " << age << endl;
        cout << "Time until first sample before " << ttime << endl;
    }

    // Reduce to_reduce from all tips
    // Add SMALL_BLEN keep some distance from age to ensure minimal branch length after reduction
    double to_reduce = ttime - age + SMALL_BLEN;

    // Reduce some amount from terminal branches to keep all of them still positive
    vector<double> tip_blens;
    for(int j = 0; j < rtree.nedge; j++){
       if (rtree.edges[j].end < Ns ){
           tip_blens.push_back(rtree.edges[j].length);
       }
    }
    double min_tblen = *min_element(tip_blens.begin(), tip_blens.end());
    if(debug){
        cout << "Differences to reduce " << to_reduce << endl;
        cout << "Minimal terminal branch length " << min_tblen << endl;
    }

    if(min_tblen >= to_reduce + SMALL_BLEN){
        if(debug) cout << "Only need to adjust tips " << endl;
        for(int j = 0; j < rtree.nedge; j++){
           if (rtree.edges[j].end < Ns ){
               // Only need to update the time of tips
               rtree.node_times[rtree.edges[j].end] =  rtree.node_times[rtree.edges[j].end] - to_reduce;
               rtree.edges[j].length = rtree.node_times[rtree.edges[j].end] - rtree.node_times[rtree.edges[j].start];
           }
        }
    }
    else{   // min_tblen < to_reduce + SMALL_BLEN
        // Always need to adjust terminal branch lengths
        vector<double> internal_blens;
        for(int j = 0; j < rtree.nedge; j++){
            if (rtree.edges[j].end < Ns ){
                // Reduce SMALL_BLEN so that the shortest terminal branch has length SMALL_BLEN after reduction
                rtree.edges[j].length = rtree.edges[j].length - (min_tblen - SMALL_BLEN);
            }
            if (rtree.edges[j].end > Ns){
               internal_blens.push_back(rtree.edges[j].length);
            }
        }
        // If there are still residuals, find some internal edge with larger length and reduce it
        double to_reduce2 = to_reduce - (min_tblen - SMALL_BLEN);
        assert(to_reduce2 > 0);
        if(debug) cout << "Differences to reduce further " << to_reduce2 << endl;
        double max_iblen = *max_element(internal_blens.begin(), internal_blens.end());
        if(max_iblen >= to_reduce2 + SMALL_BLEN){
            if(debug){
                cout << "Adjust one internal branch length " << endl;
                cout << "Edge lengths before: ";
                for(int j = 0; j < rtree.nedge; j++){
                    cout << "\t" << rtree.edges[j].length;
                }
                cout << endl;
                cout << "Node times before: ";
                for(int i=0; i<rtree.nodes.size(); ++i){
                    cout << "\t" << rtree.node_times[i];
                }
                cout << endl;
            }
            for(int j = 0; j < rtree.nedge; j++){
                if (rtree.edges[j].end > Ns && rtree.edges[j].length >= to_reduce2 + SMALL_BLEN ){
                    if(debug) cout << "Adjusting edge " << j << endl;
                    rtree.edges[j].length = rtree.edges[j].length - to_reduce2;
                    break;
                }
            }
            if(debug){
                cout << "Edge lengths after: ";
                for(int j = 0; j < rtree.nedge; j++){
                    cout << "\t" << rtree.edges[j].length;
                }
                cout << endl;
            }
        }else{  // max_iblen < to_reduce2 + SMALL_BLEN, rarely needed
            if(debug) cout << "Adjust several internal branch lengths " << endl;
            // All internal branch lengths are smaller than to_reduce2 + SMALL_BLEN, so reduction has to be done multiple times
            double delta = to_reduce2;
            for(int j = 0; j < rtree.nedge; j++){
                if (rtree.edges[j].end > Ns && rtree.edges[j].length > 2 * SMALL_BLEN){
                    if(delta >= (rtree.edges[j].length - SMALL_BLEN)){   // need another reduction
                        rtree.edges[j].length = SMALL_BLEN;
                        delta = delta - (rtree.edges[j].length - SMALL_BLEN);
                    }
                    else{   // delta < (rtree.edges[j].length - SMALL_BLEN
                        rtree.edges[j].length = rtree.edges[j].length - delta;
                        delta = -1;
                    }
                    if(delta <= 0) break;
                }
            }
        }

        rtree.calculate_node_times();
        if(debug){
            cout << "Node times after: ";
            for(int i=0; i<rtree.nodes.size(); ++i){
                cout << "\t" << rtree.node_times[i];
            }
            cout << endl;
        }
    }
    assert(is_blen_valid(rtree));

    if(debug) {
        // cout << "Age of patient " << age << endl;
        cout << "Time until first sample after " << rtree.get_total_time() << endl;
    }
}


void adjust_blen(double& nx, double a, double b){
    // double a = BLEN_MIN;
    // double b = BLEN_MAX;
    // cout << "before adjusting " << nx << endl;
    double n, e;
    if(nx < a){
        e = a - nx;
        n = floor(e/(b-a));
        nx = a + (e - n*(b-a));
        // cout << "after adjusting lower bound " << nx << endl;
        assert(nx > a);
    }
    if(nx > b){
        e = nx - b;
        n = floor(e/(b-a));
        nx = b - (e - n*(b-a));
        // cout << "after adjusting upper bound " << nx << endl;
        assert(nx < b);
    }
    // cout << "after adjusting " << nx << endl;
    assert(nx >= a && nx <= b);
}


// Seem not work well
bool is_tip_valid(evo_tree &rtree, vector<double>& tobs, int sample1){
    int debug = 0;
    for(int i = 0; i < tobs.size(); i++){
        // current differences to first sample
        double delta = rtree.node_times[i] - rtree.node_times[sample1];
        if(fabs(delta-tobs[i]) > SMALL_VAL){
            if(debug){
                cout << delta << endl;
                cout << tobs[i] << endl;
                cout << fabs(delta-tobs[i]) << endl;
                cout << "tip " << i << " is not valid" << endl;
                rtree.print();
            }
            return false;
        }
    }
    return true;
}

// Adjust the terminal branch lengths of a tree with fixed topolgy when the constraints are violated
// The initial or perturbed tree may not have tip nodes in accord with the sampling time information
// Ensure edge length keep positive
void adjust_tree_tips(evo_tree &rtree){
    assert(is_blen_valid(rtree));
    int debug = 0;
    // find the first sample
    vector<double> tobs = rtree.tobs;
    int sample1 = 0;
    for(int i = 0; i < tobs.size(); i++){
        //
        if(tobs[i] == 0){
            sample1 = i;
            break;
        }
    }
    if(debug){
        cout << "A first sample is " << sample1 << endl;
    }

    // Check the node times for tip nodes
    adjust_sample1(rtree, tobs, sample1);
    for(int i = 0; i < tobs.size(); i++){
        if(i==sample1) continue;
        if(debug) cout << "Adjusting tip " << i << endl;
        adjust_one_tip(rtree, tobs, i, sample1);
    }
    adjust_all_tips(rtree, tobs, age);
}


// The initial branch lengths are always positive.
// But after optimization, some branch lengths may become negative.
// Increase all branch lengths with the same value to keep tip timing difference
// The result tree may break the height constraint
void adjust_tree_blens_all(evo_tree &rtree){
    int debug = 0;
    vector<double> blens;
    for(int i=0; i<rtree.nedge; ++i){
        blens.push_back(rtree.edges[i].length);
    }

    double min_blen = *min_element(blens.begin(), blens.end());
    if(min_blen < 0){
        double delta = min_blen - SMALL_BLEN;
        if(debug){
            cout << "Estimated branch lengths have negative values!" << endl;
            cout << "Increase all branch lengths (except the branch to normal node) to eliminate negative values by " << -delta << endl;
        }
        for(int i = 0; i<blens.size(); i++){
            if(rtree.edges[i].start == rtree.nleaf && rtree.edges[i].end == rtree.nleaf - 1) continue;
            blens[i] = blens[i] - delta;
            rtree.node_times[rtree.edges[i].end] = rtree.node_times[rtree.edges[i].end] - delta;
        }
    }
    assert(blens.size() == rtree.nedge);
    for(int i=0; i<rtree.nedge; ++i){
        rtree.edges[i].length = blens[i];
    }
}



// Only change negative branches. The tip differences will be adjusted later
void adjust_tree_blens(evo_tree &rtree){
    int debug = 0;
    for(int i=0; i<rtree.nedge; ++i){
        if(rtree.edges[i].length < 0){
            if(debug) cout << "Negative length for edge " << i << endl;
            rtree.edges[i].length = SMALL_BLEN;
            rtree.node_times[rtree.edges[i].end] = rtree.node_times[rtree.edges[i].start] + SMALL_BLEN;
        }
    }
}


// Scale the tree so that total tree height is in the lifetime of the patient
// This may cause tips nodes violating the given sampling time
void adjust_tree_height(evo_tree &rtree, int scale = HEIGHT_SCALE){
    int debug = 0;
    double tree_height = rtree.get_tree_height();
    double min_height = *max_element(tobs.begin(), tobs.end());
    double max_height = age + min_height;
    // cout << "ttime " << ttime << endl;
    if(tree_height > max_height){
        // Use smaller height to avoid new violations after tip adjustment
        double max_height_scaled = (max_height/scale > min_height) ? max_height/scale : min_height + scale;
        double new_tree_height = runiform(r, min_height, max_height_scaled);
        double ratio = new_tree_height/tree_height;
        if(debug){
            cout << "Estimated tree height larger than age at last sample!" << endl;
            // randomly change tree height until first sample to satifisty the constraints
            rtree.print();
            cout << "Adjusting tree height with ratio " << ratio << endl;
        }
        rtree.scale_time(ratio);
        if(debug){
            rtree.print();
        }
    }
    assert(is_blen_valid(rtree));
}





//
// Likelihood and transition probabilities
//
double get_transition_prob(const double& mu, double blength, const int& sk, const int& sj){
  //if(debug) cout << "\tget_transition_prob" << endl;

  // assume a basic Markov model here cf JC
  // Qij = u/5 when i!=j
  // we have five states 0, 1, 2, 3, 4
  // Pij(t) = exp(-ut) k_ij + (1 - exp(-ut))*pi_j
  // pi_0 = pi_1 = pi_2 = pi_3 = pi_4 = 1/5

  double prob = 0;

  if( sk == sj ){
    prob = exp(-mu*blength) + (1 - exp(-mu*blength))*0.2;
  }else{
    prob = (1 - exp(-mu*blength))*0.2;
  }

  return prob;
}


// boost::numeric::ublas::matrix<double> get_rate_matrix_bounded(double dup_rate=0.01, double del_rate=0.01) {
//     boost::numeric::ublas::matrix<double> m(cn_max+1, cn_max+1);
//     for (unsigned i = 0; i < m.size1 (); ++ i){
//         for (unsigned j = 0; j < m.size2 (); ++ j){
//             m (i, j) = 0;
//         }
//     }
//     for( unsigned i = 1; i < cn_max; i++){
//         m (i,i-1) = 2 * i * del_rate;
//         m (i,i+1) = 2 * i * dup_rate;
//         m(i, i) = 0 - m (i,i-1) - m (i,i+1);
//     }
//     m (cn_max, cn_max - 1) = 2 * cn_max * del_rate;
//     m (cn_max, cn_max) = 0 - m (cn_max, cn_max - 1);
//
//     // std::cout << m << std::endl;
//     return m;
// }

void get_rate_matrix_bounded(double* m, double dup_rate, double del_rate, const int cn_max) {
    int debug = 0;
    int ncol = cn_max + 1;

    for (unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }
    for( unsigned i = 1; i < cn_max; i++){
        m[i+(i-1)*ncol] = 2 * i * del_rate;
        m[i+(i+1)*ncol] = 2 * i * dup_rate;
        m[i+i*ncol] = 0 - m[i+(i-1)*ncol] - m[i+(i+1)*ncol];
    }
    m[cn_max + (cn_max - 1)*ncol] = 2 * cn_max * del_rate;
    m[cn_max + cn_max*ncol] = 0 - m[cn_max + (cn_max - 1)*ncol];

    if(debug){
        std::cout << m << std::endl;
        r8mat_print ( ncol, ncol, m, "  A:" );
    }
}


// A matrix with dimension (cn_max + 1) * (cn_max + 2) / 2
// Suppose copy number configuration is in specific order such as:  0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2	 2/1	3/0	0/4	 1/3	  2/2	 3/1	4/0
void get_rate_matrix_allele_specific(double* m, double dup_rate, double del_rate, int cn_max) {
    int debug = 0;
    int ncol = (cn_max + 1) * (cn_max + 2) / 2;
    if(debug){
     cout << "Total number of states for copy number " << cn_max << " is " << ncol << endl;
    }

    for (unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }

    int s = cn_max * (cn_max - 1) / 2; // start index for the last group
    int r, c;   // row index and column index
    // 1st row for the last group (only deletion is allowed)
    r = ncol - (cn_max + 1);
    if(debug) cout << "Filling row " << r << endl;
    c = s;
    m[r + c*ncol] = del_rate;
    m[r + r*ncol] = 0 - del_rate;
    // last row for the last group
    r = ncol - 1 ;
    if(debug) cout << "Filling row " << r << endl;
    c = r - (cn_max + 1);
    m[r + c*ncol] = del_rate;
    m[r + r*ncol] = 0 - del_rate;
    // middle rows with two possiblities
    for(int j = cn_max; j > 1; j--){
        r = ncol - j;
        if(debug) cout << "Filling row " << r << endl;
        c = s;
        m[r + c*ncol] = del_rate;
        m[r + (c + 1)*ncol] = del_rate;
        m[r + r*ncol] = 0 - 2 * del_rate;

        s++;
    }


    int tcn = 0;    // total copy number
    // For entries with at least one zero, filling bottom up
    for(int i = 1; i < cn_max; i++){
        tcn = cn_max - i;
        int start_r = tcn * (tcn+1)/2;

        // for zero at the front
        r = start_r;
        if(debug) cout << "Filling row " << r << endl;
        c = r - tcn;
        m[r + c*ncol] = del_rate;
        c = r + (tcn + 1);
        m[r + c*ncol] = dup_rate;
        m[r + r*ncol] = 0 - (del_rate + dup_rate);

        // for zero at the end
        r = start_r + tcn;
        if(debug) cout << "Filling row " << r << endl;
        c = r - (tcn + 1);
        m[r + c*ncol] = del_rate;
        c = r + (tcn + 2);
        m[r + c*ncol] = dup_rate;
        m[r + r*ncol] = 0 - (del_rate + dup_rate);
    }

    // For entries with none zero (from copy number 2 to cn_max-1, filling bottom up (excluding last group)
    for(int i = 1; i < cn_max - 1; i++){
        tcn = cn_max - i;
        int start_r = tcn * (tcn+1)/2 + 1;
        for(int j = 1; j < tcn; j++){   // tcn - 2 rows to fill
            r = start_r + j - 1;
            if(debug) cout << "Filling row " << r << endl;

            c = r - (tcn + 1);
            m[r + c*ncol] = del_rate;
            c = c + 1;
            m[r + c*ncol] = del_rate;

            c = r + (tcn + 1);
            m[r + c*ncol] = dup_rate;
            m[r + (c+1)*ncol] = dup_rate;

            m[r + r*ncol] = 0 - 2 * (del_rate + dup_rate);
        }
    }

    if(debug){
        r8mat_print ( ncol, ncol, m, "  A:" );
    }
}

//
void get_rate_matrix_site_change(double* m, double dup_rate, double del_rate, int site_change_max){
    int debug = 0;
    int ncol = 2 * site_change_max + 1;
    int lcol = ncol - 1;

    for (unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }
    for( unsigned i = 1; i < ncol - 1; i++){
        m[i+(i-1)*ncol] = del_rate;
        m[i+(i+1)*ncol] = dup_rate;
        m[i+i*ncol] = 0 - m[i+(i-1)*ncol] - m[i+(i+1)*ncol];
    }

    m[ncol] = dup_rate;
    m[0] = - dup_rate;

    m[lcol + (lcol - 1)*ncol] = del_rate;
    m[lcol + lcol*ncol] = - del_rate;

    if(debug){
        std::cout << m << std::endl;
        r8mat_print ( ncol, ncol, m, "  A:" );
    }
}


//
void get_rate_matrix_chr_change(double* m, double chr_gain_rate, double chr_loss_rate, int chr_change_max){
    int debug = 0;
    int ncol = 2 * chr_change_max + 1;
    int lcol = ncol - 1;

    for (unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }
    for( unsigned i = 1; i < ncol - 1; i++){
        m[i+(i-1)*ncol] = chr_loss_rate;
        m[i+(i+1)*ncol] = chr_gain_rate;
        m[i+i*ncol] = 0 - m[i+(i-1)*ncol] - m[i+(i+1)*ncol];
    }

    m[ncol] = chr_gain_rate;
    m[0] = - chr_gain_rate;

    m[lcol + (lcol - 1)*ncol] = chr_loss_rate;
    m[lcol + lcol*ncol] = - chr_loss_rate;

    if(debug){
        std::cout << m << std::endl;
        r8mat_print ( ncol, ncol, m, "  A:" );
    }
}


// Maximum allowed number of WGD events over a time interval: 2
void get_rate_matrix_wgd(double* m, double wgd_rate, int wgd_max=2){
    int debug = 0;
    int ncol = wgd_max + 1;

    for (unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }
    for(unsigned i = 0; i < wgd_max; i++){
        m[i+(i+1)*ncol] = wgd_rate;
        m[i+i*ncol] = - wgd_rate;
    }

    if(debug){
        std::cout << m << std::endl;
        r8mat_print ( ncol, ncol, m, "  A:" );
    }
}


// // http://www.guwi17.de/ublas/examples/
// // Use eigendecomposition to compute matrix exponential
// double get_transition_prob_bounded(const boost::numeric::ublas::matrix<double>& q, double t, const int& sk, const int& sj ){
//     //
//     boost::numeric::ublas::matrix<double> tmp = q * t;
//     boost::numeric::ublas::matrix<double> p = expm_pad(tmp);
//     cout << "t: " << t << endl;
//     cout << "Q matrix: " << q << endl;
//     cout << "tmp matrix: " << tmp << endl;
//     cout << "P matrix: " << p << endl;
//     return p(sk, sj);
// }
// n = cn_max + 1 for model 1 (total copy number)
void get_transition_matrix_bounded(double* q, double* p, double t, int n){
    int debug = 0;
    // int n = cn_max + 1;
    double *tmp = new double[n*n];
    memset( tmp, 0, n*n*sizeof(double) );

    for(int i = 0; i < n*n; i++){
        tmp[i] = q[i] * t;
    }
    double* res = r8mat_expm1(n, tmp);
    for(int i=0; i<n*n; i++){
        if(res[i]<0){
            p[i] = 0;
        }
        else{
            p[i] = res[i];
        }
    }

    if(debug){
        cout << "t: " << t << endl;
        r8mat_print ( n, n, q, "  Q matrix:" );
        r8mat_print ( n, n, tmp, "  TMP matrix:" );
        r8mat_print ( n, n, p, "  P matrix:" );
    }

    free(tmp);
}

double get_transition_prob_bounded(double* p, const int& sk, const int& sj, int n){
    int debug = 0;

    // int n = cn_max + 1;
    int i = sk  + sj * n;
    double v = p[i];

    if(debug){
        r8mat_print ( n, n, p, "  P matrix before access:" );
        cout << "Prob at " << i << "("<< sk << ", " << sj << ") is: " << v << endl;
    }
    // if(v<SMALL_VAL){
    //     v = 0;
    // }
    return v;
}


// Randomly swap two leaves
evo_tree perturb_tree(const int& Ns, evo_tree& tree){
    int debug = 0;
    if(debug) cout << "\tperturb one tree" << endl;

    // swap two leaf nodes: 0 to Ns-1 are the samples, Ns is germline, Ns+1 is root
    vector<int> n_0;
    for(int i=0; i<Ns; ++i) n_0.push_back( i );
    random_shuffle(n_0.begin(), n_0.end(), fp);
    int n1 = n_0[Ns-1];
    n_0.pop_back();
    random_shuffle(n_0.begin(), n_0.end(), fp);
    int n2 = n_0[0];

    if(debug){
        cout << "\toriginal tree: ";
        cout << order_tree_string_uniq(create_tree_string_uniq( tree )) << endl;
        cout << tree.make_newick(5) << endl;
        cout << "\tswapping nodes:\t" << n1+1 << "\t" << n2+1 << endl;
    }

    vector<edge> enew;
    int e1, e2;
    int end1, end2;
    double len1, len2;
    for(int i=0; i<tree.nedge; ++i){
        enew.push_back( tree.edges[i] );
        bool changed = false;
        if(enew[i].end == n1 && changed == false){
          enew[i].end = n2;
          e1 = i;
          end1 = n1;
          len1 = enew[i].length;
          changed = true;
        }
        if(enew[i].end == n2 && changed == false){
          enew[i].end = n1;
          e2 = i;
          end2 = n2;
          len2 = enew[i].length;
          changed = true;
        }
    }

    // cout << "edge 1 " << "\t" << e1  << "\t" << enew[e1].start + 1 << "\t" << end1 + 1 << "\t" <<
    // enew[e1].length << endl;
    // cout << "edge 2 " << "\t" << e2  << "\t" << enew[e2].start + 1 << "\t" << end2 + 1 << "\t" <<
    // enew[e2].length << endl;

    enew[e1].length = len2;
    enew[e2].length = len1;

    // cout << "edge 1 " << "\t" << e1  << "\t" << enew[e1].start + 1 << "\t" << enew[e1].end + 1 << "\t" <<
    // enew[e1].length << endl;
    // cout << "edge 2 " << "\t" << e2  << "\t" << enew[e2].start + 1 << "\t" << enew[e2].end + 1 << "\t" <<
    // enew[e2].length << endl;

    evo_tree new_tree(Ns+1, enew, 1);

    new_tree.mu = tree.mu;
    new_tree.dup_rate = tree.dup_rate;
    new_tree.del_rate = tree.del_rate;
    new_tree.chr_gain_rate = tree.chr_gain_rate;
    new_tree.chr_loss_rate = tree.chr_loss_rate;
    new_tree.wgd_rate = tree.wgd_rate;
    new_tree.tobs = tree.tobs;

    if(debug){
        // cout << "\tnew tree: ";
        // cout << create_tree_string_uniq( new_tree ) << endl;
        string ordered = order_tree_string_uniq(create_tree_string_uniq( new_tree ) );
        cout << "\tordered new tree: ";
        cout << ordered << endl;
        cout << new_tree.make_newick(5) << endl;
        cout << endl;
    }

    return new_tree;
}


// Use NNI for rooted (clock) trees
evo_tree do_random_NNI(evo_tree &rtree, int debug){
    if(debug){
        printf ("Before:\n");
        rtree.print();
        // Output to a file to facilitate ploting
        // ofstream out_tree("./tree-before-NNI.txt");
        // rtree.write(out_tree);
        // out_tree.close();
    }
    evo_tree ntree(rtree);
    /* pick an interior branch (excluding the top edge connecting to root Ns+1), around which it is possible to make an NNI */
    double tv;
    double tc;
    int i, u, v, c, a, b;
    edge *e_uv, *e_uc, *e_va;
    vector<int> children;
    double rnum;
    int count = 0;
    double blen_uc;
    double blen_va;
    tc = -2;
    tv = -1;
    set<double> chosen; // Stop if all the internal edges have been tested

    while (tc <= tv && chosen.size() < rtree.nintedge - 1){
        count++;
        // Exclude the edges connecting to unaltered genome
        assert(rtree.intedges[rtree.nintedge-1]->start == Ns+1);
        i = myrng(rtree.nintedge-1);
        chosen.insert(i);
        if(debug){
            cout << "Picking " << i+1 << "th internal edge: " << rtree.intedges[i]->start+1 << ", " << rtree.intedges[i]->end+1 << endl;
        }
        e_uv = rtree.intedges[i];
        // Find the children of the start node
        u = e_uv->start;
        children = rtree.nodes[u].daughters;
        // Select the child with subtree as v
        v = *max_element(children.begin(), children.end());
        c = *min_element(children.begin(), children.end());

        assert(v > Ns+1);

        tv = rtree.node_times[v];
        tc = rtree.node_times[c];

        // Find the children of v
        children = rtree.nodes[v].daughters;
        assert(children.size()==2);
        rnum = runiform(r, 0, 1);
        if(rnum > 0.5){
            a = children[0];
        }
        else{
            a = children[1];
        }
    }

    if(tc <= tv){
        if(debug) cout << "Cannot update topology!" << endl;
        return ntree;
    }
    if(debug){
        cout << "Picked an edge after " << count << " iteration(s)" << endl;
        cout << "u, v, a, c " << u+1 << "\t" << v+1 << "\t" << a+1 << "\t" << c+1 << endl;
    }

    /* record branch lengths */
    // oldALength = e_va->length;
    // oldCLength = e_uc->length;

    /* make topology change */
    if(debug){
        cout << "Making topology change" << endl;
    }
    // (u,c) --> (u,a); (v,a) --> (v,c)
    // move the subtree below
    // Find the length of edge uc and va
    for(int j=0; j < ntree.edges.size(); j++){
        edge *e = &ntree.edges[j];
        if(e->start == u && e->end == c){
            blen_uc = e->length;
            e->end = a;
            e_uc = e;
        }
        if(e->start == v && e->end == a){
            blen_va = e->length;
            e->end = c;
            e_va = e;
        }
    }

    // update branch length to keep node ages unchanged
    e_va->length = blen_uc - e_uv->length;
    e_uc->length = e_uv->length + blen_va;
    // cout << " edge v->a length " << e_va->length << endl;
    // cout << " edge u->c length " << e_uc->length << endl;

    adjust_blen(e_va->length, BLEN_MIN, BLEN_MAX);
    adjust_blen(e_uc->length, BLEN_MIN, BLEN_MAX);

    ntree.nodes.clear();
    ntree.generate_nodes();
    ntree.intedges.clear();
    ntree.generate_int_edges();

    /* adjust branch lengths */
    // e_uc->length = rtree.node_times[a] - rtree.node_times[u];
    // e_va->length = rtree.node_times[c] - rtree.node_times[v];

    if(debug){
        printf ("After:\n");
        ntree.print();
        // Output to a file to facilitate ploting
        // ofstream out_tree("./test/tree-after-NNI.txt");
        // rtree.write(out_tree);
        // out_tree.close();
    }
    return ntree;
}


// Do NNI on an internal branch i, according to P293 (Yang, 2014)
vector<evo_tree> do_NNI_on_branch(evo_tree &rtree, int i, vector<NNIMove>& moves, int debug){
    if(debug){
        printf ("Before:\n");
        rtree.print();
        // Output to a file to facilitate ploting
        // ofstream out_tree("./tree-before-NNI.txt");
        // rtree.write(out_tree);
        // out_tree.close();
    }
    vector<evo_tree> nni_trees;
    evo_tree tree1(rtree);
    evo_tree tree2(rtree);
    /* pick an interior branch (excluding the top edge connecting to root Ns+1), around which it is possible to make an NNI */
    double tv, tc;
    int u, v, c, a1, a2, b;
    edge *e_uv, *e_uc, *e_va;
    vector<int> children;
    double blen_uc, blen_va;
    tc = -2;
    tv = -1;

    if(debug){
        cout << "Picking " << i+1 << "th internal edge: " << rtree.intedges[i]->start+1 << ", " << rtree.intedges[i]->end+1 << endl;
    }
    e_uv = rtree.intedges[i];
    // Find the children of the start node
    u = e_uv->start;
    children = rtree.nodes[u].daughters;
    // Select the child with subtree as v
    v = *max_element(children.begin(), children.end());
    c = *min_element(children.begin(), children.end());
    assert(v > Ns+1);
    tv = rtree.node_times[v];
    tc = rtree.node_times[c];
    if(tc <= tv){
        cout << "Cannot update topology!" << endl;
        return nni_trees;
    }
    // Find the children of v
    children = rtree.nodes[v].daughters;
    assert(children.size()==2);
    a1 = children[0];
    a2 = children[1];

    if(debug){
        cout << "u, v, a1, a2, c " << u+1 << "\t" << v+1 << "\t" << a1+1 << "\t" << a2+1 << "\t" << c+1 << endl;
    }

    /* make topology change */
    // (u,c) --> (u,a); (v,a) --> (v,c)
    // move the subtree below
    // Find the length of edge uc and va
    for(int j=0; j < tree1.edges.size(); j++){
        edge *e = &tree1.edges[j];
        if(e->start == u && e->end == c){
            blen_uc = e->length;
            e->end = a1;
            e_uc = e;
        }
        if(e->start == v && e->end == a1){
            blen_va = e->length;
            e->end = c;
            e_va = e;
        }
    }
    // update branch length to keep node ages unchanged
    e_va->length = blen_uc - e_uv->length;
    e_uc->length = e_uv->length + blen_va;
    adjust_blen(e_va->length, BLEN_MIN, BLEN_MAX);
    adjust_blen(e_uc->length, BLEN_MIN, BLEN_MAX);
    tree1.nodes.clear();
    tree1.generate_nodes();
    tree1.intedges.clear();
    tree1.generate_int_edges();
    // NNIMove m1;
    // m1.swap_id = i;
    // m1.node1 = u;
    // m1.node2 = v;

    for(int j=0; j < tree2.edges.size(); j++){
        edge *e = &tree2.edges[j];
        if(e->start == u && e->end == c){
            blen_uc = e->length;
            e->end = a2;
            e_uc = e;
        }
        if(e->start == v && e->end == a2){
            blen_va = e->length;
            e->end = c;
            e_va = e;
        }
    }
    // update branch length to keep node ages unchanged
    e_va->length = blen_uc - e_uv->length;
    e_uc->length = e_uv->length + blen_va;
    adjust_blen(e_va->length, BLEN_MIN, BLEN_MAX);
    adjust_blen(e_uc->length, BLEN_MIN, BLEN_MAX);
    tree2.nodes.clear();
    tree2.generate_nodes();
    tree2.intedges.clear();
    tree2.generate_int_edges();
    // NNIMove m2;
    // m2.swap_id = i;
    // m2.node1 = u;
    // m2.node2 = v;

    if(debug){
        printf ("After:\n");
        tree1.print();
        tree2.print();
        // Output to a file to facilitate ploting
        // ofstream out_tree("./test/tree-after-NNI.txt");
        // rtree.write(out_tree);
        // out_tree.close();
    }

    nni_trees.push_back(tree1);
    nni_trees.push_back(tree2);

    // moves.push_back(m1);
    // moves.push_back(m2);

    return nni_trees;
}



// Find all the NNI trees for a given tree
map<int, vector<evo_tree>> get_all_NNI_trees(evo_tree &rtree, vector<NNIMove>& moves, int debug){
    map<int, vector<evo_tree>> nni_trees_all;
    // vector<evo_tree> nni_trees_all;
    for(int i = 0; i < rtree.intedges.size(); i++){
        vector<evo_tree> nni_trees = do_NNI_on_branch(rtree, i, moves, debug);
        nni_trees_all[i] = nni_trees;
        // nni_trees_all.insert(nni_trees_all.end(), nni_trees.begin(), nni_trees.end());
    }
    return nni_trees_all;
}

// Check whether the set of branch lengths is valid under the time constraints
int is_tree_valid(evo_tree &rtree, int cons){
    for(int i=0; i<rtree.edges.size(); ++i){
        if (rtree.edges[i].length < 0 ){
            cout << "Negative branch length!" << endl;
            return 0;
        }
    }
    if(cons){
        double tot = rtree.get_total_time();
        if(floor(tot) > age){
            cout << "Wrong total time (larger than age)!" << "\t" << tot << "\t" << age << endl;
            return 0;
        }
    }
    return 1;
}



// Create likelihood vectors at the tip node
// L_sk_k has one row for each tree node and one column for each possible state
vector<vector<double>> initialize_lnl_table(vector<int>& obs, evo_tree &rtree, int Ns, int model, int nstate, int is_total){
    int debug = 0;
    // construct a table for each state of each node
    vector<vector<double>> L_sk_k(rtree.ntotn, vector<double>(nstate,0));

    if(model > 1){
        for(int i=0; i<Ns; ++i){
            // For total copy number, all the possible combinations have to be considered.
            // Set related allele specific cases to be 1, with index from obs[i] * (obs[i] + 1)/2 to obs[i] * (obs[i] + 1)/2 + obs[i]. The index is computed based on pre-specified order.
            if(is_total==1){
                int si = (obs[i] * (obs[i] + 1))/2;
                int ei = si + obs[i];
                for(int k = si; k<=ei; k++){
                    L_sk_k[i][k] = 1.0;
                }
            }else{ // With allele-specific copy number, only the specific site needs to be filled
                L_sk_k[i][obs[i]] = 1.0;
            }
        }
        // set unaltered 1/1
        L_sk_k[Ns][4] = 1.0;
    }
    else{
        for(int i=0; i<Ns; ++i){
          for(int j=0; j<nstate; ++j){
    	         if( j == obs[i] ) L_sk_k[i][j] = 1.0;
          }
        }
        // set unaltered
        L_sk_k[Ns][2] = 1.0;
    }

    if(debug){
      cout << "\nCNs at tips:\n";
      for(int i=0; i<Ns; ++i){
          cout<< "\t" << obs[i];
      }
      cout << endl;
      cout << "\nLikelihood for tips:\n";
      for(int i=0; i<rtree.nleaf; ++i){
          for(int j=0; j < nstate; ++j){
            cout << "\t" << L_sk_k[i][j];
          }
          cout << endl;
      }
    }
    return L_sk_k;
}


// Create likelihood vectors and state vectors at the tip node for reconstructing joint ancestral state
// L_sk_k (S_sk_k) has one row for each tree node and one column for each possible state
void initialize_asr_table(vector<int>& obs, evo_tree &rtree, map<double, double*>& pmats, vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, int Ns, int model, int nstate, int is_total){
    int debug = 0;
    if(debug) cout << "Initializing tables for reconstructing joint ancestral state" << endl;

    for(int i=0; i<Ns; ++i){
        // cout << "node " << i + 1 << endl;
        // Find the parent node
        int parent = rtree.edges[rtree.nodes[i].e_in].start;
        double blen = rtree.edges[rtree.nodes[i].e_in].length;
        // cout << "parent " << parent + 1 << endl;
        // cout << "blen " << blen << endl;
        double* pblen = pmats[blen];
        // Find the state(s) of current node
        vector<int> tip_states;
        if(model==2){
          if(is_total==1){
              int si = (obs[i] * (obs[i] + 1))/2;
              int ei = si + obs[i];
              for(int k = si; k<=ei; k++){
                 tip_states.push_back(k);
              }
          }else{ // With allele-specific copy number, only the specific site needs to be filled
              tip_states.push_back(obs[i]);
          }
        }
        if(debug) cout << "There are " << tip_states.size() << " states for copy number " << obs[i] << endl;
        for(int j=0; j<nstate; ++j){  // For each possible parent state, find the most likely tip states
            vector<double> vec_li(nstate, SMALL_LNL);
            for(int m=0; m<tip_states.size(); ++m){
                int k = tip_states[m];
                if(debug) cout << "parent state " << j << ", child state " << k << endl;
                double li = 0;
                if (model == 0){
                    li = get_transition_prob(rtree.mu, blen, j, k);
                }
                else{
                    li =  pblen[j  + k * nstate];  // assume parent has state j
                }
                if(li > 0) li = log(li);
                else li = SMALL_LNL;
                vec_li[k] = li;
            }

            int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
            double max_li = *max_element(vec_li.begin(), vec_li.end());
            assert(max_li == vec_li[max_i]);
            // cout << "for node: i " << i + 1 << ", parent state " << j << ", max state is " << max_i << " with probability " << exp(max_li) << endl;
            S_sk_k[i][j] = max_i;
            L_sk_k[i][j] = max_li;
        }
    }

    if(debug){
      cout << "\nCNs at tips:\n";
      for(int i=0; i<Ns; ++i){
          cout<< "\t" << obs[i];
      }
      cout << endl;
      cout << "\nLikelihood for tips:\n";
      for(int i=0; i<rtree.nleaf; ++i){
          for(int j=0; j < nstate; ++j){
            cout << "\t" << L_sk_k[i][j];
          }
          cout << endl;
      }
      cout << "\nState vector for tips:\n";
      for(int i=0; i<rtree.nleaf; ++i){
          for(int j=0; j < nstate; ++j){
            cout << "\t" << S_sk_k[i][j];
          }
          cout << endl;
      }
    }

}


// Create likelihood vectors and state vectors at the tip node for reconstructing joint ancestral state, for independent chain model
// L_sk_k (S_sk_k) has one row for each tree node and one column for each possible state
void initialize_asr_table_decomp(vector<int>& obs, evo_tree &rtree, map<double, double*>& pmats_wgd, map<double, double*>& pmats_chr, map<double, double*>& pmats_seg, vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, int Ns, int model, int nstate, int is_total){
    int debug = 0;
    if(debug) cout << "Initializing tables for reconstructing joint ancestral state" << endl;
    int s_wgd, s_chr, s_seg, s_chr2, s_seg2;
    int e_wgd, e_chr, e_seg, e_chr2, e_seg2;
    double prob_wgd, prob_chr, prob_seg, prob_chr2, prob_seg2, prob_chr_all, prob_seg_all;
    int dim_wgd = max_wgd + 1;
    int dim_chr = 2 * max_chr_change + 1;
    int dim_seg = 2 * max_site_change + 1;

    for(int i=0; i<Ns; ++i){
        // cout << "node " << i + 1 << endl;
        // Find the parent node
        int parent = rtree.edges[rtree.nodes[i].e_in].start;
        double bli = rtree.edges[rtree.nodes[i].e_in].length;
        // cout << "parent " << parent + 1 << endl;
        // cout << "blen " << blen << endl;
        double *pbli_wgd;
        double *pbli_chr;
        double *pbli_seg;
        if(dim_wgd > 1){
            pbli_wgd = pmats_wgd[bli];
        }
        if(dim_chr > 1){
            pbli_chr = pmats_chr[bli];
        }
        if(dim_seg > 1){
            pbli_seg = pmats_seg[bli];
        }

        // Find the state(s) of current node
        vector<int> tip_states;
        int k = 0;
        for (auto c : comps){
            for(int m1 = 0; m1 <= m_max; m1++){
                for(int m2 = 0; m2 <= m_max; m2++){
                    int sum = pow(2, c[0] + 1) + m1 * c[1] + c[2] + 2 * m2 * c[3] + 2 * c[4];
                    if(sum==obs[i]){
                        tip_states.push_back(k);
                    }
                }
            }
            k++;
        }
        if(debug) cout << "There are " << tip_states.size() << " states for copy number " << obs[i] << endl;

        for(int j=0; j<nstate; ++j){  // For each possible parent state, find the most likely tip states
            set<vector<int>>::iterator iter = comps.begin();
            // It will move forward the passed iterator by passed value
            advance(iter, j);
            vector<int> s = *iter;
            // get_decomposition(sk, s_wgd, s_chr, s_seg, decomp_table);
            s_wgd = s[0];
            s_chr = s[1];
            s_seg = s[2];
            s_chr2 = s[3];
            s_seg2 = s[4];
            // The indices for chromosome and segment matrix have to be ajusted
            int delta_chr = (dim_chr - 1)/2;
            int delta_seg = (dim_seg - 1)/2;
            if(debug){
                cout << "Starting state " << j << "\t" << s_wgd << "\t" << s_chr << "\t" << s_seg << "\t" << s_chr2 << "\t" << s_seg2 << "\n";
                cout << "   Offset for chr and segment matrices " << delta_chr << "\t" << delta_seg << "\n";
            }

            vector<double> vec_li(nstate, SMALL_LNL);
            for(int m=0; m<tip_states.size(); ++m){
                int k = tip_states[m];
                if(debug) cout << "parent state " << j << ", child state " << k << endl;
                double li = 0;
                prob_wgd = 1;
                prob_chr = 1;
                prob_seg = 1;
                prob_chr2 = 1;
                prob_seg2 = 1;
                prob_chr_all = 1;
                prob_seg_all = 1;

                iter = comps.begin();
                advance(iter, k);
                vector<int> e = *iter;
                // get_decomposition(sk, s_wgd, s_chr, s_seg, decomp_table);
                e_wgd = e[0];
                e_chr = e[1];
                e_seg = e[2];
                e_chr2 = e[3];
                e_seg2 = e[4];
                if(dim_wgd > 1) prob_wgd = pbli_wgd[s_wgd + e_wgd * dim_wgd];
                if(dim_chr > 1){
                    prob_chr = pbli_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
                    prob_chr2 = pbli_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
                    prob_chr_all = prob_chr * prob_chr2;
                }
                if(dim_seg > 1){
                    prob_seg = pbli_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
                    prob_seg2 = pbli_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
                    prob_seg_all = prob_seg * prob_seg2;
                }
                li = prob_wgd * prob_chr_all * prob_seg_all;
                if(debug){
                    cout << "End state " << k << "\t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2  << "\n";
                    cout << "Prob for each event " << "\t" << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << li << "\n";
                }

                if(li > 0) li = log(li);
                else li = SMALL_LNL;
                vec_li[k] = li;
            }

            int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
            double max_li = *max_element(vec_li.begin(), vec_li.end());
            assert(max_li == vec_li[max_i]);
            // cout << "for node: i " << i + 1 << ", parent state " << j << ", max state is " << max_i << " with probability " << exp(max_li) << endl;
            S_sk_k[i][j] = max_i;
            L_sk_k[i][j] = max_li;
        }
    }

    if(debug){
      cout << "\nCNs at tips:\n";
      for(int i=0; i<Ns; ++i){
          cout<< "\t" << obs[i];
      }
      cout << endl;
      cout << "\nLikelihood for tips:\n";
      for(int i=0; i<rtree.nleaf; ++i){
          for(int j=0; j < nstate; ++j){
            cout << "\t" << L_sk_k[i][j];
          }
          cout << endl;
      }
      cout << "\nState vector for tips:\n";
      for(int i=0; i<rtree.nleaf; ++i){
          for(int j=0; j < nstate; ++j){
            cout << "\t" << S_sk_k[i][j];
          }
          cout << endl;
      }
    }

}

void print_comps(set<vector<int>> comps){
    cout << "\nAll possible states for copy number decomposition in order:" << endl;
    int count = 0;
    for (auto c : comps){
        count += 1;
        cout << "\tstate " << count << endl;
        for(int i = 0; i < c.size(); i++){
            cout << "\t" << c[i];
        }
        cout << endl;
    }
}


void print_decomp_table(map<int, set<vector<int>>>& decomp_table){
    int count = 0;
    for(auto item : decomp_table){
        set<vector<int>> comp = item.second;
        cout << "\nCopy number " << item.first << endl;
        for (auto c : comp){
            count += 1;
            for(int i = 0; i < c.size(); i++){
                cout << "\t" << c[i];
            }
            cout << endl;
        }
    }
    cout << "Number of possible states for copy number decomposition " << count << endl;
}

// L_sk_k has one row for each tree node and one column for each possible state; chr starting from 1
// This function is critical in obtaining correct likelihood. If one tip is not initialized, the final likelihood will be 0.
vector<vector<double>> initialize_lnl_table_decomp(vector<int>& obs, int chr, evo_tree &rtree, set<vector<int>>& comps, int Ns, int model, int cn_max, int& nstate, int is_total){
    int debug = 0;
    // construct a table for each state of each node
    vector<vector<double>> L_sk_k(rtree.ntotn, vector<double>(nstate,0));

    for(int i=0; i<Ns; ++i){
        // For total copy number, all the possible combinations have to be considered.
        // Set related allele specific cases to be 1, with index from obs[i] * (obs[i] + 1)/2 to obs[i] * (obs[i] + 1)/2 + obs[i]. The index is computed based on pre-specified order.
        int cn = obs[i];
        int num_change = 0;
        if(infer_chr == 1 && chr > 0){
            num_change = obs_change_chr[i][chr-1];
        }
        if(is_total==0){    // changing input allel-specific copy number to total
            cn = state_to_total_cn(obs[i], cn_max);
        }
        if(debug) cout << "\nStart filling likelihood table for sample "  << i + 1 << " chromosome " << chr  << " copy number " << cn << endl;
        // Fill all the possible state combinations
        int k = 0;
        for (auto c : comps){
            if(debug){
                cout << "\tcn vector:";
                for(int k = 0; k < c.size(); k++){
                    cout << "\t" << c[k];
                }
                cout << endl;
            }
            int alpha = c[0];  // wgd component
            // If a sample has one WGD event, the correponding component must be the same
            if(infer_wgd == 1 && alpha != obs_num_wgd[i]) continue;

            // WGD may occur before, at least one chromosome gain before or after WGD
            if(num_change >= 1 && (c[1] <= 0 && c[3] <= 0)){
                if(debug){
                    cout << "\t\tpotential chromosome gain in sample "  << i + 1 << " chromosome " << chr << " is " << num_change << endl;
                }
                continue;
            }
            if(num_change <= -1 && (c[1] >= 0 && c[3] >= 0)){
                if(debug){
                    cout << "\t\tpotential chromosome loss in sample "  << i + 1 << " chromosome " << chr << " is " << num_change << endl;
                }
                continue;
            }
            // int sum = pow(2, alpha + 1) + c[5] * c[1] + c[2] + 2 * c[6] * c[3] + 2 * c[4];
            // if(sum == cn){
            //     if(debug) cout << "\t\tfilling 1 here" << endl;
            //     L_sk_k[i][k] = 1.0;
            // }
            // assuming m_max >= 1. It is likely that all copies of a segment is lost before chromosome gain/loss
            for(int m1 = 0; m1 <= m_max; m1++){
                for(int m2 = 0; m2 <= m_max; m2++){
                    int sum = pow(2, alpha + 1) + m1 * c[1] + c[2] + 2 * m2 * c[3] + 2 * c[4];
                    if(sum == cn){
                        if(debug) cout << "\t\tfilling 1 here" << endl;
                        L_sk_k[i][k] = 1.0;
                    }
                }
            }
            k++;
        }
        // each row should have one entry being 1
        int sum = 0;
        for(int j = 0; j < L_sk_k[i].size(); j++){
            sum += L_sk_k[i][j];
        }
        if(sum < 1){
            cout << "Error in filling table for copy number " << cn << " in sample " << i+1 << " chromosome " << chr << endl << endl;
        }
    }
    // set likelihood for normal sample
    // Fill all the possible state combinations
    for (int j = 0; j < comps.size(); j++){
        int k = 0;
        for (auto v : comps){
            bool zeros = all_of(v.begin(), v.end(), [](int i) { return i==0; });
            if(zeros){
                L_sk_k[Ns][k] = 1.0;
                break;
            }
            k++;
        }
    }

    if(debug){
      cout << "\nCNs at tips:\n";
      for(int i=0; i<Ns; ++i){
          cout<< "\t" << obs[i];
      }
      cout << endl;
      cout << "\nLikelihood for tips:\n";
      for(int i=0; i<rtree.nleaf; ++i){
          for(int j=0; j < nstate; ++j){
            cout << "\t" << L_sk_k[i][j];
          }
          cout << endl;
      }
    }
    return L_sk_k;
}


// Get the likelihood of the tree from likelihood table
double extract_tree_lnl(const vector<vector<double>>& L_sk_k, int Ns, int model){
    int debug = 0;

    if(debug){
        for(int j=0; j<L_sk_k[Ns+1].size(); ++j){
          cout << "\t" << L_sk_k[Ns+1][j];
        }
        cout << endl;
    }

    if(model == 2){
        // The index is changed from 2 to 4 (1/1)
        if(debug) cout << "Likelihood for root is " << L_sk_k[Ns+1][4] << endl;
        if(L_sk_k[Ns+1][4] > 0) return log(L_sk_k[Ns+1][4]);
        else return LARGE_LNL;
    }else{
        if(debug) cout << "Likelihood for root is " << L_sk_k[Ns+1][2] << endl;
        if(L_sk_k[Ns+1][2] > 0) return log(L_sk_k[Ns+1][2]);
        else return LARGE_LNL;
    }
}


// Get the states of the tree from likelihood table at one site, starting from MRCA
map<int, vector<int>> extract_tree_ancestral_state(evo_tree &rtree, const vector<vector<double>>& L_sk_k, const vector<vector<int>>& S_sk_k, int Ns, int model, map<int, int> &asr_states){
    int debug = 0;
    map<int, vector<int>> asr_cns;     // copy numbers for each internal node at one site, stored in the order of decreasing IDs
    int max_id = 2 * Ns;

    if(debug){
        for(int j=0; j<L_sk_k[Ns+1].size(); ++j){
          cout << "\t" << L_sk_k[Ns+1][j];
        }
        cout << endl;
    }

    int parent_state = 2;
    if(model == 2){
        parent_state = 4;
    }
    if(model == 3) {
        int k = 0;
        for (auto v : comps){
            // if(debug) cout << k << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << "\t" << v[3] << "\t" << v[4] << endl;
            bool zeros = all_of(v.begin(), v.end(), [](int i) { return i==0; });
            if(zeros){
                break;
            }
            k++;
        }
        parent_state = k;
    }
    asr_states[Ns+1] = parent_state;
    if(is_total){
        vector<int> cns{2};
        asr_cns[Ns+1] = cns;
    }
    else{
        vector<int> cns{1,1};
        asr_cns[Ns+1] = cns;
    }

    // Traverse the tree from root to tips, need to know the parent of each node
    for(int nid = max_id; nid > Ns + 1; nid--){
        // Find the parent node of current node
        int parent = rtree.nodes[nid].parent;
        assert(asr_states.find(parent) != asr_states.end());
        parent_state = asr_states[parent];

        int state = S_sk_k[nid][parent_state];
        asr_states[nid] = state;
        if(debug){
            cout << "\t\tnode " << nid + 1 << " with state " << state  << " and parent " << parent + 1 << " whose state is " << parent_state << endl;
        }

        vector<int> cns;
        if(model==3){
            // Find the decomposition with state
            set<vector<int>>::iterator iter = comps.begin();
            advance(iter, state);
            vector<int> c = *iter;
            set<int> sums;
            for(int m1 = 0; m1 <= m_max; m1++){
                for(int m2 = 0; m2 <= m_max; m2++){
                    int cn = pow(2, c[0] + 1) + m1 * c[1] + c[2] + 2 * m2 * c[3] + 2 * c[4];
                    // only add unique copy number
                    sums.insert(cn);
                }
            }
            for(auto s : sums){
                cns.push_back(s);
            }
        }
        else{
            if(is_total == 1){
                int cn = state_to_total_cn(state, cn_max);
                cns.push_back(cn);
            }else{
                int cnA, cnB;
                state_to_allele_cn(state, cn_max, cnA, cnB);
                cns.push_back(cnA);
                cns.push_back(cnB);
            }
        }
        asr_cns[nid] = cns;
    }

    return asr_cns;
}


// Get the likelihood of the tree from likelihood table of state combinations
double extract_tree_lnl_decomp(const vector<vector<double>>& L_sk_k, set<vector<int>>& comps, int Ns){
    int debug = 0;
    if(debug) cout << "Extracting likelihood for the root" << endl;
    double likelihood = 0;
    int k = 0;
    for (auto v : comps){
        // if(debug) cout << k << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << "\t" << v[3] << "\t" << v[4] << endl;
        bool zeros = all_of(v.begin(), v.end(), [](int i) { return i==0; });
        if(zeros){
            likelihood = L_sk_k[Ns+1][k];
            break;
        }
        k++;
    }

    if(debug){
        for(int j=0; j<comps.size(); ++j){
          cout << "\t" << L_sk_k[Ns+1][j];
        }
        cout << endl;
    }

    if(likelihood > 0) return log(likelihood);
    else return LARGE_LNL;
}


void print_tree_lnl(evo_tree &rtree, vector<vector<double>>& L_sk_k, int nstate){
    cout << "\nLikelihood so far:\n";
    for(int i=0; i<rtree.ntotn; ++i){
        for(int j=0; j<nstate; ++j){
          cout << "\t" << L_sk_k[i][j];
        }
        cout << endl;
    }
}


void print_tree_state(evo_tree &rtree, vector<vector<int>>& S_sk_k, int nstate){
    cout << "\nStates so far:\n";
    for(int i=0; i<rtree.ntotn; ++i){
        for(int j=0; j<nstate; ++j){
          cout << "\t" << S_sk_k[i][j];
        }
        cout << endl;
    }
}


double get_prob_children(vector<vector<double>>& L_sk_k, evo_tree &rtree, double* pbli, double* pblj, int nsk, int nstate, int ni, int nj, int bli, int blj, int model){
    int debug = 0;

    double Li = 0;
    // loop over possible si
    for(int si=0; si<nstate; ++si){
      if(L_sk_k[ni][si] > 0){
          if (model == 0){
            Li += get_transition_prob(rtree.mu, bli, nsk, si) * L_sk_k[ni][si];
          }
          else{
            // if(debug){
            //     cout << "1st branch length " << bli << endl;
            //     cout << get_transition_prob_bounded(pbli, nsk, si, nstate) << endl;
            //     r8mat_print( nstate, nstate, pbli, " Pmatrix for 1st branch length " );
            // }
            Li += get_transition_prob_bounded(pbli, nsk, si, nstate) * L_sk_k[ni][si];
          }
       }
    }
    if(debug) cout << "\tscoring: Li\t" << Li << endl;

    double Lj = 0;
    // loop over possible sj
    for(int sj=0; sj<nstate; ++sj){
        if(L_sk_k[nj][sj] > 0){
            if (model == 0){
                 Lj += get_transition_prob(rtree.mu, blj, nsk, sj) * L_sk_k[nj][sj];
            }
            else{
                // if(debug){
                //     cout << "2nd branch length " << blj << endl;
                //     r8mat_print( nstate, nstate, pblj, " Pmatrix for 2nd branch length ");
                // }
                 Lj += get_transition_prob_bounded(pblj, nsk, sj, nstate) * L_sk_k[nj][sj];
            }
        }
    }
    if(debug) cout << "\tscoring: Lj\t" << Lj << endl;

    return Li * Lj;
}


// Find the most likely state for a node under each possible state, assuming current node has state nsk and parent node (connected by branch of length blen) has state np
double get_max_prob_children(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, evo_tree &rtree, double* pblen, int k, int nstate, int sp, int ni, int nj, int blen, int model){
    int debug = 0;

    vector<double> vec_li;
    double li = 0;
    // loop over possible si for a fixed state of parent (sp)
    for(int si=0; si<nstate; ++si){
        if (model == 0){
          li = get_transition_prob(rtree.mu, blen, sp, si);
        }
        else{
          li = pblen[sp + si * nstate];
        }

        if(debug){
            cout << "\tfor state " << si << endl;
            cout << "\t\tseparate likelihood "  << li << "\t" << L_sk_k[ni][si] << "\t" << L_sk_k[nj][si] << endl;
        }
        if(li > 0){
            li = log(li);
        }
        else li = SMALL_LNL;
        li += (L_sk_k[ni][si]);
        li += (L_sk_k[nj][si]);

        if(debug) cout << "\t\tscoring: Li\t" << li << endl;
        if(std::isnan(li) || li < SMALL_LNL) li = SMALL_LNL;
        if(debug) cout << "\t\trevised scoring: Li\t" << li << endl;
        vec_li.push_back(li);
    }
    int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
    double max_li = *max_element(vec_li.begin(), vec_li.end());
    assert(max_li == vec_li[max_i]);
    S_sk_k[k][sp] = max_i;

    if(debug){
        cout << "all likelihoods:";
        for(int i = 0; i < vec_li.size(); i++){
            cout << "\t" << vec_li[i];
        }
        cout << endl;
        cout << "for node: k " << k + 1 << ", parent state " << sp << ", max state is " << max_i << " with probability " << exp(max_li) << endl;
    }

    return max_li;
}



// Get possible combinations for a total copy number
void insert_tuple(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int i, int j, int k){
    int sum = pow(2, i + 1) + j + k;
    // cout << i << "\t" << j << "\t" << k << "\t" << sum << endl;
    if(sum >= 0 && sum <= cn_max){
        vector<int> c{i,j,k};
        decomp_table[sum].insert(c);
        comps.insert(c);
    }
}




// Get possible combinations for a total copy number with ordering of different types of events considered (at most 1 WGD)
void insert_tuple_order_withm(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int m1, int m2, int i, int j, int k, int j0, int k0){
    int sum = pow(2, i + 1) + m1 * j + k + 2 * m2 * j0 + 2 * k0;
    // cout << i << "\t" << j << "\t" << k << "\t" << sum << endl;
    if(sum >= 0 && sum <= cn_max){
        vector<int> c{i, j, k, j0, k0, m1, m2};
        decomp_table[sum].insert(c);
        comps.insert(c);
    }
}


// Get possible combinations for a total copy number with ordering of different types of events considered (at most 1 WGD)
void insert_tuple_order(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int i, int j, int k, int j0, int k0){
    for(int m = 0; m <= m_max; m++){
        for(int m = 0; m <= m_max; m++){
            int sum = pow(2, i + 1) + j + k + 2 * j0 + 2 * k0;
            // cout << i << "\t" << j << "\t" << k << "\t" << sum << endl;
            if(sum >= 0 && sum <= cn_max){
                vector<int> c{i, j, k, j0, k0};
                decomp_table[sum].insert(c);
                comps.insert(c);
            }
        }
    }
    // // Cases when there are multiple copies of a segment on a chromosome before chromosome gain/loss
    // if(j > 0){
    //     for(int m = 2; m <= m_max; m++){
    //         sum = pow(2, i + 1) + m * j + k + 2 * j0 + 2 * k0;
    //         if(sum >= 0 && sum <= cn_max){
    //             vector<int> c{i, j, k, j0, k0};
    //             decomp_table[sum].insert(c);
    //             comps.insert(c);
    //         }
    //         if(j0 > 0){
    //             for(int m2 = 2; m2 <= m_max; m2++){
    //                 sum = pow(2, i + 1) + m * j + k + 2 * m2 * j0 + 2 * k0;
    //                 if(sum >= 0 && sum <= cn_max){
    //                     vector<int> c{i, j, k, j0, k0};
    //                     decomp_table[sum].insert(c);
    //                     comps.insert(c);
    //                 }
    //             }
    //         }
    //     }
    // }
    // if(j0 > 0){
    //     for(int m2 = 2; m2 <= m_max; m2++){
    //         sum = pow(2, i + 1) + j + k + 2 * m2 * j0 + 2 * k0;
    //         if(sum >= 0 && sum <= cn_max){
    //             vector<int> c{i, j, k, j0, k0};
    //             decomp_table[sum].insert(c);
    //             comps.insert(c);
    //         }
    //     }
    // }
}


// Get possible allele-specific combinations for a total copy number
void insert_tuple_allele_specific(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int i, int j1, int j2, int k1, int k2){
    int sum = pow(2, i + 1) + j1 + j2 + k1 + k2;
    // cout << i << "\t" << j1 << "\t" << j2 << "\t" << k1 << "\t" << k2 << "\t" << sum << endl;
    if(sum >= 0 && sum <= cn_max){
        vector<int> c{i, j1, j2, k1, k2};
        decomp_table[sum].insert(c);
        comps.insert(c);
    }
}

// Adjust the value of m_max based on the maximum copy number for a sample
void adjust_m_max(){
    // For samples without WGD but with larger copy number
    int max_sum = pow(2, 1) + m_max * max_chr_change + max_site_change;
    int orig_m_max = m_max;
    int max_cn = max_sum;
    for(int i = 0; i < sample_max_cn.size(); i++){
        if(obs_num_wgd[i] > 0) continue;
        int cn = sample_max_cn[i];
        if(cn > max_cn){
            max_cn = cn;
        }
    }
    m_max = m_max + (max_cn - max_sum);
    int max_copy = 2 * (max_site_change + 1);    // at most 2 * (max_site_change + 1) copies of a segment before chromosome gain/loss
    m_max = m_max > max_copy ? max_copy : m_max;
    cout << "original m_max " << orig_m_max << endl;
    cout << "Max copy number for samples without WGD under original parameter setting is " << max_sum << endl;
    cout << "new m_max " << m_max << endl;
}


// All the input variables are global
map<int, set<vector<int>>> build_decomp_table()
{
    map<int, set<vector<int>>> decomp_table;
    for(int i=0; i<=cn_max; i++){
        set<vector<int>> comp;
        decomp_table[i] = comp;
    }
    if(is_total == 1){
        // for(int i=0; i<=max_wgd; i++){
        //     for(int j=0; j<=max_chr_change; j++){
        //         for(int k=0; k<=max_site_change; k++){
        //             insert_tuple(decomp_table, comps, cn_max, i, j, k);
        //             insert_tuple(decomp_table, comps, cn_max, i, -j, k);
        //             insert_tuple(decomp_table, comps, cn_max, i, j, -k);
        //             insert_tuple(decomp_table, comps, cn_max, i, -j, -k);
        //         }
        //     }
        // }
        // No need to consider addtional terms when no WGD occurs
        for(int j=0; j<=max_chr_change; j++){
            for(int k=0; k<=max_site_change; k++){
                insert_tuple_order(decomp_table, comps, cn_max, m_max, 0, j, k, 0, 0);
                insert_tuple_order(decomp_table, comps, cn_max, m_max, 0, -j, k, 0, 0);
                insert_tuple_order(decomp_table, comps, cn_max, m_max, 0, j, -k, 0, 0);
                insert_tuple_order(decomp_table, comps, cn_max, m_max, 0, -j, -k, 0, 0);
            }
        }

        // With WGD
        // For samples with WGD but with larger copy number
        // max_sum = pow(2, 2) + 2 * m_max * max_chr_change + 2 * max_site_change;
        // m_max = m_max;
        // for(int i = 0; i < sample_max_cn.size(); i++){
        //     if(obs_num_wgd[i] < 1) continue;
        //     int cn = sample_max_cn[i];
        //     if(cn > max_sum){
        //         m_max = m_max + (cn - max_sum);
        //     }
        // }

        for(int i=1; i<=max_wgd; i++){
            for(int j=0; j<=max_chr_change; j++){
                for(int k=0; k<=max_site_change; k++){
                    for(int j0=0; j0<=max_chr_change - j; j0++){
                        for(int k0=0; k0<=max_site_change - k; k0++){
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, k, j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, k, -j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, k, j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, k, -j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, k, -j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, k, j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, k, -j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, k, j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, -k, j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, -k, -j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, -k, j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, -k, -j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, -k, j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, -k, -j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, -k, j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, -k, -j0, -k0);
                        }
                    }
                }
            }
        }
    }else{
        for(int i=0; i<=max_wgd; i++){
            for(int j1=0; j1<=max_chr_change; j1++){
                for(int j2=0; j2<=max_chr_change; j2++){
                    for(int k1=0; k1<=max_site_change; k1++){
                        for(int k2=0; k2<=max_site_change; k2++){
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, -k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, -k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, -k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, -k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, -k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, -k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, -k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, -k1, -k2);
                        }
                    }
                }
            }
        }
    }
    return decomp_table;
}


// All the input variables are global, treating m_j as part of the vector
map<int, set<vector<int>>> build_decomp_table_withm()
{
    map<int, set<vector<int>>> decomp_table;
    for(int i=0; i<=cn_max; i++){
        set<vector<int>> comp;
        decomp_table[i] = comp;
    }
    if(is_total == 1){
        int m1_max, m2_max;  // No need to consider m_max when there is no chromosome change
        // No need to consider addtional terms when no WGD occurs.
        for(int j=0; j<=max_chr_change; j++){
            for(int k=0; k<=max_site_change; k++){
                if(j==0) m1_max = 0; else m1_max = m_max;
                for(int m1=0; m1<=m1_max; m1++){
                        insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, 0, 0, j, k, 0, 0);
                        insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, 0, 0, -j, k, 0, 0);
                        insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, 0, 0, j, -k, 0, 0);
                        insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, 0, 0, -j, -k, 0, 0);
                }
            }
        }

        for(int i=1; i<=max_wgd; i++){
            for(int j=0; j<=max_chr_change; j++){
                for(int j0=0; j0<=max_chr_change - j; j0++){
                    for(int k=0; k<=max_site_change; k++){
                        for(int k0=0; k0<=max_site_change - k; k0++){
                            if(j==0) m1_max = 0; else m1_max = m_max;
                            if(j0==0) m2_max = 0; else m2_max = m_max;
                            for(int m1=0; m1<=m1_max; m1++){
                                for(int m2=0; m2<=m2_max; m2++){
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, k, j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, k, -j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, k, j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, k, -j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, k, -j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, k, j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, k, -j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, k, j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, -k, j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, -k, -j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, -k, j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, -k, -j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, -k, j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, -k, -j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, -k, j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, -k, -j0, -k0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }else{
        for(int i=0; i<=max_wgd; i++){
            for(int j1=0; j1<=max_chr_change; j1++){
                for(int j2=0; j2<=max_chr_change; j2++){
                    for(int k1=0; k1<=max_site_change; k1++){
                        for(int k2=0; k2<=max_site_change; k2++){
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, -k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, -k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, -k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, -k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, -k1, k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, -k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, -k1, -k2);
                            insert_tuple_allele_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, -k1, -k2);
                        }
                    }
                }
            }
        }
    }
    return decomp_table;
}


// Assume the likelihood table is for each allele-specific copy number
double get_prob_children_decomp(vector<vector<double>>& L_sk_k, evo_tree &rtree, map<int, set<vector<int>>>& decomp_table, int sk, int cn_max, int nstate, double* pbli_wgd, double* pblj_wgd, double* pbli_chr, double* pblj_chr, double* pbli_seg, double* pblj_seg, int dim_wgd, int dim_chr, int dim_seg, int ni, int nj, int bli, int blj, int is_total){
    int debug = 0;
    int s_wgd, s_chr, s_seg;
    int e_wgd, e_chr, e_seg;
    // Each copy number is decomposed into a set of 3-tuples
    set<vector<int>> comp_start = decomp_table[sk];
    // get_decomposition(sk, s_wgd, s_chr, s_seg, decomp_table);

    double Li = 0;
    for(auto s : comp_start){
        s_wgd = s[0];
        s_chr = s[1] + s[3];
        s_seg = s[2] + s[4];
        // loop over possible si
        for(int si=0; si<nstate; ++si){
            // get the start and end state for each type
            int cn = state_to_total_cn(si, cn_max);
            set<vector<int>> comp_end = decomp_table[cn];
            // get_decomposition(si, e_wgd, e_chr, e_seg, decomp_table);
            for(auto e : comp_end){
                e_wgd = e[0];
                e_chr = e[1] + e[3];
                e_seg = e[2] + e[4];
                double prob_wgd = get_transition_prob_bounded(pbli_wgd, s_wgd, e_wgd, dim_wgd);
                double prob_chr = get_transition_prob_bounded(pbli_chr, s_chr, e_chr, dim_chr);
                double prob_seg = get_transition_prob_bounded(pbli_seg, s_seg, e_seg, dim_seg);
                double prob = prob_wgd * prob_chr * prob_seg * L_sk_k[ni][si];
                Li += prob;
                if(debug) cout << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob << "\n";
            }
          if(debug) cout << "\tscoring: Li\t" << Li << endl;
        }
    }

    double Lj = 0;
    for(auto s : comp_start){
        s_wgd = s[0];
        s_chr = s[1] + s[3];
        s_seg = s[2] + s[4];
        // loop over possible sj
        for(int sj=0; sj<nstate; ++sj){
            int cn = state_to_total_cn(sj, cn_max);
            set<vector<int>> comp_end = decomp_table[cn];
            // get the start and end state for each type
            // set<vector<int>> comp_end = decomp_table[sj];
            // get_decomposition(si, e_wgd, e_chr, e_seg, decomp_table);
            for(auto e : comp_end){
                e_wgd = e[0];
                e_chr = e[1] + e[3];
                e_seg = e[2] + e[4];
                double prob_wgd = get_transition_prob_bounded(pblj_wgd, s_wgd, e_wgd, dim_wgd);
                double prob_chr = get_transition_prob_bounded(pblj_chr, s_chr, e_chr, dim_chr);
                double prob_seg = get_transition_prob_bounded(pblj_seg, s_seg, e_seg, dim_seg);
                double prob = prob_wgd * prob_chr * prob_seg * L_sk_k[nj][sj];
                Lj += prob;
                if(debug) cout << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob << "\n";
            }
        }
        if(debug) cout << "\tscoring: Lj\t" << Lj << endl;
    }

    return Li * Lj;
}

// Assume the likelihood table is for each combination of states
double get_prob_children_decomp2(vector<vector<double>>& L_sk_k, evo_tree &rtree, set<vector<int>>& comps, int sk, int cn_max, int nstate, double* pbli_wgd, double* pblj_wgd, double* pbli_chr, double* pblj_chr, double* pbli_seg, double* pblj_seg, int dim_wgd, int dim_chr, int dim_seg, int ni, int nj, int bli, int blj, int is_total){
    int debug = 0;
    int s_wgd, s_chr, s_seg, s_chr2, s_seg2;
    int e_wgd, e_chr, e_seg, e_chr2, e_seg2;
    double prob_wgd, prob_chr, prob_seg, prob_chr2, prob_seg2, prob_chr_all, prob_seg_all;
    // Each copy number is decomposed into a set of 3-tuples
    set<vector<int>>::iterator iter = comps.begin();
    // It will move forward the passed iterator by passed value
    advance(iter, sk);
    vector<int> s = *iter;
    // get_decomposition(sk, s_wgd, s_chr, s_seg, decomp_table);
    s_wgd = s[0];
    s_chr = s[1];
    s_seg = s[2];
    s_chr2 = s[3];
    s_seg2 = s[4];
    // The indices for chromosome and segment matrix have to be ajusted
    int delta_chr = (dim_chr - 1)/2;
    int delta_seg = (dim_seg - 1)/2;

    if(debug){
        cout << "Starting state " << sk << "\t" << s_wgd << "\t" << s_chr << "\t" << s_seg << "\t" << s_chr2 << "\t" << s_seg2 << "\n";
        cout << "   Offset for chr and segment matrices " << delta_chr << "\t" << delta_seg << "\n";
    }


    double Li = 0;
    int si = 0;
    for(auto e : comps){
        double prob = 0;
        prob_wgd = 1;
        prob_chr = 1;
        prob_seg = 1;
        prob_chr2 = 1;
        prob_seg2 = 1;
        prob_chr_all = 1;
        prob_seg_all = 1;
        if(L_sk_k[ni][si] > 0){
            e_wgd = e[0];
            e_chr = e[1];
            e_seg = e[2];
            e_chr2 = e[3];
            e_seg2 = e[4];
            if(dim_wgd > 1) prob_wgd = pbli_wgd[s_wgd + e_wgd * dim_wgd];
            if(dim_chr > 1){
                prob_chr = pbli_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
                prob_chr2 = pbli_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
                prob_chr_all = prob_chr * prob_chr2;
            }
            if(dim_seg > 1){
                prob_seg = pbli_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
                prob_seg2 = pbli_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
                prob_seg_all = prob_seg * prob_seg2;
            }
            prob = prob_wgd * prob_chr_all * prob_seg_all * L_sk_k[ni][si];
            if(debug){
                cout << "End state " << si << "\t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2 << "\n";
                cout << "Prob for each event " << "\t" << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << prob << "\n";
            }
        }
        Li += prob;
        si++;
    }
    if(debug) cout << "\tscoring: Li\t" << Li << endl;

    double Lj = 0;
    int sj = 0;
    for(auto e : comps){
        double prob = 0;
        prob_wgd = 1;
        prob_chr = 1;
        prob_seg = 1;
        prob_chr2 = 1;
        prob_seg2 = 1;
        prob_chr_all = 1;
        prob_seg_all = 1;
        if(L_sk_k[nj][sj] > 0){
            e_wgd = e[0];
            e_chr = e[1];
            e_seg = e[2];
            e_chr2 = e[3];
            e_seg2 = e[4];
            if(dim_wgd > 1) prob_wgd = pblj_wgd[s_wgd + e_wgd * dim_wgd];
            if(dim_chr > 1){
                prob_chr = pblj_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
                prob_chr2 = pblj_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
                prob_chr_all = prob_chr * prob_chr2;
            }
            if(dim_seg > 1){
                prob_seg = pblj_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
                prob_seg2 = pblj_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
                prob_seg_all = prob_seg * prob_seg2;
            }
            prob = prob_wgd * prob_chr_all * prob_seg_all * L_sk_k[nj][sj];
            if(debug){
                cout << "End state " << sj << "\t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2 << "\n";
                cout << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << prob << "\n";
            }
        }
        Lj += prob;
        sj++;
    }
    if(debug) cout << "\tscoring: Lj\t" << Lj << endl;

    return Li * Lj;
}


// Assume the likelihood table is for each combination of states
double get_max_children_decomp2(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, evo_tree &rtree, set<vector<int>>& comps, int k, int nstate, double* pbli_wgd, double* pbli_chr, double* pbli_seg, int dim_wgd, int dim_chr, int dim_seg, int sp, int ni, int nj, int blen){
    int debug = 0;
    int s_wgd, s_chr, s_seg, s_chr2, s_seg2;
    int e_wgd, e_chr, e_seg, e_chr2, e_seg2;
    double prob_wgd, prob_chr, prob_seg, prob_chr2, prob_seg2, prob_chr_all, prob_seg_all;
    // Each copy number is decomposed into a set of 3-tuples
    set<vector<int>>::iterator iter = comps.begin();
    advance(iter, sp);
    vector<int> s = *iter;
    s_wgd = s[0];
    s_chr = s[1];
    s_seg = s[2];
    s_chr2 = s[3];
    s_seg2 = s[4];
    // The indices for chromosome and segment matrix have to be ajusted
    int delta_chr = (dim_chr - 1)/2;
    int delta_seg = (dim_seg - 1)/2;

    if(debug){
        cout << "Starting state " << sp << "\t" << s_wgd << "\t" << s_chr << "\t" << s_seg << "\t" << s_chr2 << "\t" << s_seg2 << "\n";
        cout << "   Offset for chr and segment matrices " << delta_chr << "\t" << delta_seg << "\n";
    }

    vector<double> vec_li;
    double li = 0;
    int si = 0;
    for(auto e : comps){
        prob_wgd = 1;
        prob_chr = 1;
        prob_seg = 1;
        prob_chr2 = 1;
        prob_seg2 = 1;
        prob_chr_all = 1;
        prob_seg_all = 1;

        e_wgd = e[0];
        e_chr = e[1];
        e_seg = e[2];
        e_chr2 = e[3];
        e_seg2 = e[4];
        if(dim_wgd > 1) prob_wgd = pbli_wgd[s_wgd + e_wgd * dim_wgd];
        if(dim_chr > 1){
            prob_chr = pbli_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
            prob_chr2 = pbli_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
            prob_chr_all = prob_chr * prob_chr2;
        }
        if(dim_seg > 1){
            prob_seg = pbli_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
            prob_seg2 = pbli_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
            prob_seg_all = prob_seg * prob_seg2;
        }
        li = prob_wgd * prob_chr_all * prob_seg_all;

        if(debug){
            cout << "End state " << si << "\t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2 << "\n";
            cout << "Prob for each event " << "\t" << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << li  << "\n";
        }

        if(li > 0){
            li = log(li);
        }
        else li = SMALL_LNL;
        li += (L_sk_k[ni][si]);
        li += (L_sk_k[nj][si]);

        if(debug) cout << "\t\tscoring: Li\t" << li << endl;
        if(std::isnan(li) || li < SMALL_LNL) li = SMALL_LNL;
        if(debug) cout << "\t\trevised scoring: Li\t" << li << endl;
        vec_li.push_back(li);

        si++;
     }
    int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
    double max_li = *max_element(vec_li.begin(), vec_li.end());
    assert(max_li == vec_li[max_i]);
    S_sk_k[k][sp] = max_i;

    if(debug){
        cout << "all likelihoods:";
        for(int i = 0; i < vec_li.size(); i++){
            cout << "\t" << vec_li[i];
        }
        cout << endl;
        cout << "for node: k " << k + 1 << ", parent state " << sp << ", max state is " << max_i << " with probability " << exp(max_li) << endl;
    }

    return max_li;
}


// Get the likelihood on one site of a chromosome (assuming higher level events on nodes)
// z: possible changes in copy number caused by chromosome gain/loss
void get_likelihood_site(vector<vector<double>>& L_sk_k, evo_tree &rtree, const vector<int>& knodes, map<double, double*>& pmats, const int has_wgd, const int z, const int model, int nstate){
  int debug = 0;
  if(debug){
      cout << "Computing likelihood for one site" << endl;
      cout << nstate << endl;
  }

  for(int kn=0; kn<knodes.size(); ++kn){
        int k = knodes[kn];
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
        double* pbli = pmats[bli];
        double* pblj = pmats[blj];

        if(debug) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni+1 << " , " << bli << "\t" <<  nj+1 << " , " << blj << endl;

        //loop over possible values of sk
        if(k == Ns + 1){    // root node is always normal
            if(debug) cout << "Getting likelihood for root node " << k << endl;
            int nsk = 2;
            if(model==2) nsk = 4;
            L_sk_k[k][nsk] = get_prob_children(L_sk_k, rtree, pbli, pblj, nsk, nstate, ni, nj, bli, blj, model);
        }
        else{
            for(int sk=0; sk<nstate; ++sk){
                int nsk = sk;  // state after changes by other large scale events
                if(has_wgd == 1) nsk = 2 * sk;
                nsk += z;
                if(debug) cout << "likelihood for state " << nsk << endl;
                if(nsk < 0 || nsk >= nstate) continue;
                // cout << " getting likelihood of children nodes " << endl;
                L_sk_k[k][nsk] = get_prob_children(L_sk_k, rtree, pbli, pblj, nsk, nstate, ni, nj, bli, blj, model);
            }
        }
  }
  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
  }
}


// Get the most likely state on one site of a chromosome (assuming higher level events on nodes)
void get_ancestral_states_site(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, evo_tree &rtree, const vector<int>& knodes, map<double, double*>& pmats, int nstate){
  int debug = 0;
  if(debug){
      cout << "Getting ancestral state for one site" << endl;
  }

  for(int kn=0; kn<knodes.size(); ++kn){
        int k = knodes[kn];
        int np = rtree.edges[rtree.nodes[k].e_in].start;
        double blen = rtree.edges[rtree.nodes[k].e_in].length;
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double* pblen = pmats[blen];

        if(debug) cout << "node:" << np + 1 << " -> " << rtree.nodes[k].id + 1 << " -> " << ni+1 << " , "  <<  nj+1 << " , " << blen << endl;

        //loop over possible values of sk
        if(k == 2 * Ns){    // root node is always normal, for edge (Ns+1, 2Ns)
            if(debug) cout << "Getting states for node MRCA " << k + 1 << endl;
            int sp = 2;
            if(model==2) sp = 4;
            if(debug) cout << "likelihood for state " << sp << endl;
            L_sk_k[k][sp] = get_max_prob_children(L_sk_k, S_sk_k, rtree, pblen, k, nstate, sp, ni, nj, blen, model);
        }
        else{
            for(int sp=0; sp<nstate; ++sp){  // looping over all possible states of its parent
                if(debug) cout << "likelihood for state " << sp << endl;
                L_sk_k[k][sp] = get_max_prob_children(L_sk_k, S_sk_k, rtree, pblen, k, nstate, sp, ni, nj, blen, model);
            }
        }
  }
  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
    print_tree_state(rtree, S_sk_k, nstate);
  }
}

// Get the likelihood on one site of a chromosome
// Assuming each observed copy number is composed of three type of events.
// Sum over all possible states for initial and final nodes
// Allow at most one WGD event along a branch
void get_likelihood_site_decomp(vector<vector<double>>& L_sk_k, evo_tree &rtree, map<int, set<vector<int>>>& decomp_table, const vector<int>& knodes, map<double, double*>& pmats_wgd, map<double, double*>& pmats_chr, map<double, double*>& pmats_seg, int dim_wgd, int dim_chr, int dim_seg, int cn_max, int nstate, int is_total){
  int debug = 0;
  if(debug){
      cout << "Computing likelihood for one site" << endl;
      cout << dim_wgd << "\t" << dim_chr << "\t"  << dim_seg << "\t"  << nstate << endl;
  }

  for(int kn=0; kn<knodes.size(); ++kn){
        int k = knodes[kn];
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

        double *pbli_wgd, *pblj_wgd;
        double *pbli_chr, *pblj_chr;
        double *pbli_seg, *pblj_seg;

        if(dim_wgd > 1){
            pbli_wgd = pmats_wgd[bli];
            pblj_wgd = pmats_wgd[blj];
        }
        if(dim_chr > 1){
            pbli_chr = pmats_chr[bli];
            pblj_chr = pmats_chr[blj];
        }
        if(dim_seg > 1){
            pbli_seg = pmats_seg[bli];
            pblj_seg = pmats_seg[blj];
        }

        if(debug) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni+1 << " , " << bli << "\t" <<  nj+1 << " , " << blj << endl;

        // loop over possible observed states of start nodes
        if(k == Ns + 1){    // root node is always normal
            // int sk = 4;
            // L_sk_k[k][sk] = get_prob_children_decomp(L_sk_k, rtree, decomp_table, sk, cn_max, nstate, pbli_wgd, pblj_wgd, pbli_chr, pblj_chr, pbli_seg, pblj_seg, dim_wgd, dim_chr, dim_seg, ni, nj, bli, blj, is_total);
            if(debug) cout << "Getting likelihood for root node " << k << endl;
            int sk = 0;
            for(auto v : comps){
                bool zeros = all_of(v.begin(), v.end(), [](int i) { return i==0; });
                if(zeros){
                    L_sk_k[k][sk] = get_prob_children_decomp2(L_sk_k, rtree, comps, sk, cn_max, nstate, pbli_wgd, pblj_wgd, pbli_chr, pblj_chr, pbli_seg, pblj_seg, dim_wgd, dim_chr, dim_seg, ni, nj, bli, blj, is_total);
                    break;
                }
                sk++;
            }
        }
        else{
            for(int sk=0; sk<nstate; ++sk){
                // cout << " getting likelihood of children nodes " << endl;
                // L_sk_k[k][sk] = get_prob_children_decomp(L_sk_k, rtree, decomp_table, sk, cn_max, nstate, pbli_wgd, pblj_wgd, pbli_chr, pblj_chr, pbli_seg, pblj_seg, dim_wgd, dim_chr, dim_seg, ni, nj, bli, blj, is_total);
                L_sk_k[k][sk] = get_prob_children_decomp2(L_sk_k, rtree, comps, sk, cn_max, nstate, pbli_wgd, pblj_wgd, pbli_chr, pblj_chr, pbli_seg, pblj_seg, dim_wgd, dim_chr, dim_seg, ni, nj, bli, blj, is_total);
            }
        }
  }
  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
  }
}

// Get the ancestral state on one site of a chromosome
// Assuming each observed copy number is composed of three type of events.
// Sum over all possible states for initial and final nodes
// Allow at most one WGD event along a branch
void get_ancestral_states_site_decomp(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, evo_tree &rtree, const vector<int>& knodes, map<double, double*>& pmats_wgd, map<double, double*>& pmats_chr, map<double, double*>& pmats_seg, int dim_wgd, int dim_chr, int dim_seg, int nstate){
  int debug = 0;
  if(debug){
      cout << "Getting ancestral state for one site under independent Markov chain model" << endl;
      cout << dim_wgd << "\t" << dim_chr << "\t"  << dim_seg << "\t"  << nstate << endl;
  }

  for(int kn=0; kn<knodes.size(); ++kn){
        int k = knodes[kn];
        int np = rtree.edges[rtree.nodes[k].e_in].start;
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blen = rtree.edges[rtree.nodes[k].e_in].length;

        double *pbli_wgd;
        double *pbli_chr;
        double *pbli_seg;

        if(dim_wgd > 1){
            pbli_wgd = pmats_wgd[blen];
        }
        if(dim_chr > 1){
            pbli_chr = pmats_chr[blen];
        }
        if(dim_seg > 1){
            pbli_seg = pmats_seg[blen];
        }

        if(debug) cout << "node:" << np + 1 << " -> " << rtree.nodes[k].id + 1 << " -> " << ni+1 << " , "  <<  nj+1 << " , " << blen << endl;

        // loop over possible observed states of start nodes
        if(k == Ns + 1){    // root node is always normal
            if(debug) cout << "Getting states for node MRCA " << k + 1 << endl;
            int sp = 0;
            for(auto v : comps){
                bool zeros = all_of(v.begin(), v.end(), [](int i) { return i==0; });
                if(zeros){
                    L_sk_k[k][sp] = get_max_children_decomp2(L_sk_k, S_sk_k, rtree, comps, k, nstate, pbli_wgd, pbli_chr, pbli_seg, dim_wgd, dim_chr, dim_seg, sp, ni, nj, blen);
                    break;
                }
                sp++;
            }
        }
        else{
            for(int sp=0; sp<nstate; ++sp){
                L_sk_k[k][sp] = get_max_children_decomp2(L_sk_k, S_sk_k, rtree, comps, k, nstate, pbli_wgd, pbli_chr, pbli_seg,  dim_wgd, dim_chr, dim_seg, sp, ni, nj, blen);
            }
        }
  }
  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
    print_tree_state(rtree, S_sk_k, nstate);
  }
}


// Infer the copy number of the MRCA given a tree at a site, assuming independent Markov chains
double reconstruct_marginal_ancestral_state_decomp(evo_tree &rtree, ofstream& fout, double min_asr){
    int debug = 0;
    if(debug) cout << "\treconstruct marginal ancestral state with independent chain model" << endl;

    if(!is_tree_valid(rtree, cons)){
        return SMALL_LNL;
    }

    vector<int> cn_mrca; // CNs for MRCA
    int nid = 2 * (Ns + 1) - 2;

    string header="node\tsite\tcn";
    for(int i = 0; i < comps.size(); i++){
        header += "\tprobablity " + to_string(i+1);
    }
    fout << header << endl;

    // For WGD model
    int dim_wgd = max_wgd + 1;
    int dim_chr = 2 * max_chr_change + 1;
    int dim_seg = 2 * max_site_change + 1;
    double *qmat_wgd, *qmat_chr, *qmat_seg;
    double *pmat_wgd, *pmat_chr, *pmat_seg;
    // Find the transition probability matrix for each branch
    map<double, double*> pmats_wgd;
    map<double, double*> pmats_chr;
    map<double, double*> pmats_seg;
    double *pmati_wgd, *pmatj_wgd;
    double *pmati_chr, *pmatj_chr;
    double *pmati_seg, *pmatj_seg;

    if(max_wgd > 0){
        qmat_wgd = new double[dim_wgd * dim_wgd];  // WGD
        memset(qmat_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
        get_rate_matrix_wgd(qmat_wgd, rtree.wgd_rate, max_wgd);
    }

    if(max_chr_change > 0){
        qmat_chr = new double[(dim_chr)*(dim_chr)];   // chromosome gain/loss
        memset(qmat_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
        get_rate_matrix_chr_change(qmat_chr, rtree.chr_gain_rate, rtree.chr_loss_rate, max_chr_change);
    }

    if(max_site_change > 0){
        qmat_seg = new double[(dim_seg)*(dim_seg)];  // segment duplication/deletion
        memset(qmat_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
        get_rate_matrix_site_change(qmat_seg, rtree.dup_rate, rtree.del_rate, max_site_change);
    }

    if(debug){
          cout << "Dimension of Qmat " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << "\n";
    }

    //create a list of nodes to loop over (only internal nodes), making sure the root is last
    vector<int> knodes;
    for(int k=Ns+2; k<rtree.ntotn; ++k) knodes.push_back( k );
    knodes.push_back(Ns+1);

    for(int kn=0; kn<knodes.size(); ++kn){
           int k = knodes[kn];
           double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
           double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

           // For WGD
           if(max_wgd > 0){
               pmati_wgd = new double[(dim_wgd)*(dim_wgd)];
               pmatj_wgd = new double[(dim_wgd)*(dim_wgd)];
               memset(pmati_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
               memset(pmatj_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
               if(pmats_wgd.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_wgd, pmati_wgd, bli, dim_wgd);
                   pmats_wgd[bli] = pmati_wgd;
               }
               if(pmats_wgd.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_wgd, pmatj_wgd, blj, dim_wgd);
                   pmats_wgd[blj] = pmatj_wgd;
               }
           }

           // For chr gain/loss
           if(max_chr_change > 0){
               pmati_chr = new double[(dim_chr)*(dim_chr)];
               pmatj_chr = new double[(dim_chr)*(dim_chr)];
               memset(pmati_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
               memset(pmatj_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
               if(pmats_chr.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_chr, pmati_chr, bli, dim_chr);
                   pmats_chr[bli] = pmati_chr;
               }
               if(pmats_chr.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_chr, pmatj_chr, blj, dim_chr);
                   pmats_chr[blj] = pmatj_chr;
               }
           }

           // For segment duplication/deletion
           if(max_site_change > 0){
               pmati_seg = new double[(dim_seg)*(dim_seg)];
               pmatj_seg = new double[(dim_seg)*(dim_seg)];
               memset(pmati_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
               memset(pmatj_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
               if(pmats_seg.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_seg, pmati_seg, bli, dim_seg);
                   pmats_seg[bli] = pmati_seg;
               }
               if(pmats_seg.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_seg, pmatj_seg, blj, dim_seg);
                   pmats_seg[blj] = pmatj_seg;
              }
          }
    }
    if(debug){
        for( auto it = pmats_wgd.begin(); it != pmats_wgd.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print( dim_wgd, dim_wgd, it->second, "  P-WGD matrix:" );
        }
        for( auto it = pmats_chr.begin(); it != pmats_chr.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print( dim_chr, dim_chr, it->second, "  P-CHR matrix:" );
        }
        for( auto it = pmats_seg.begin(); it != pmats_seg.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print( dim_seg, dim_seg, it->second, "  P-SEG matrix:" );
        }
    }

    int nstate = comps.size();
    double logL = 0;    // for all chromosmes

    for(int nchr=1; nchr<=vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "Computing likelihood on Chr " << nchr << endl;
      double site_logL = 0;   // log likelihood for all sites on a chromosome
      // Use a map to store computed log likelihood
      map<vector<int>, vector<vector<double>>> sites_lnl_map;
      // cout << " chromosome number change is " << 0 << endl;
      for(int nc=0; nc<vobs[nchr].size(); nc++){    // for each segment on the chromosome
          // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k;
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  L_sk_k = initialize_lnl_table_decomp(obs, nchr, rtree, comps, Ns, model, cn_max, nstate, is_total);
                  get_likelihood_site_decomp(L_sk_k, rtree, decomp_table, knodes, pmats_wgd, pmats_chr, pmats_seg, dim_wgd, dim_chr, dim_seg, cn_max, nstate, is_total);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  if(debug) cout << "\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              L_sk_k = initialize_lnl_table_decomp(obs, nchr, rtree, comps, Ns, model, cn_max, nstate, is_total);
              get_likelihood_site_decomp(L_sk_k, rtree, decomp_table, knodes, pmats_wgd, pmats_chr, pmats_seg, dim_wgd, dim_chr, dim_seg, cn_max, nstate, is_total);
          }
          // site_logL += extract_tree_lnl(L_sk_k, Ns, model);
          double lnl = extract_tree_lnl_decomp(L_sk_k, comps, Ns);
          site_logL += lnl;

          // Get the likelihood table of MRCA node (with largest ID) in the tree from likelihood table
          int state = distance(L_sk_k[nid].begin(), max_element(L_sk_k[nid].begin(), L_sk_k[nid].end()));

          set<vector<int>>::iterator iter = comps.begin();
          // It will move forward the passed iterator by passed value
          advance(iter, state);
          vector<int> c = *iter;
          // There may be multiple possible copy numbers given different coeffients for chromosome-level changes
          int cn = pow(2, c[0] + 1) +  c[1] + c[2] + 2 * c[3] + 2 * c[4];

          cn_mrca.push_back(cn);
          // Print the state of MRCA at this site
          string line = to_string(nid+1) + "\t" + to_string(nchr) + "_" + to_string(nc) + "\t" + to_string(cn);
          for(int i = 0; i < L_sk_k[nid].size(); i++){
              line += "\t" + to_string(L_sk_k[nid][i]);
          }
          fout << line << endl;

          if(debug){
              cout << "\t\tState with maximum likelihood " << state << " at chr " << nchr << " site " << nc << endl;
              cout << "\t\t" << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << c[3] << "\t" << c[4] << endl;
          }
      }
      logL += site_logL;
      if(debug){
          cout << "\nLikelihood for chromosome " << nchr << " is " << site_logL << endl;
      }
    } // for each chromosome
    if(debug){
        cout << "\nLikelihood without correcting acquisition bias: " << logL << endl;
        cout << "CNs at MRCA is: " << endl;
        for(int i = 0; i < cn_mrca.size(); i++){
            cout << cn_mrca[i] << endl;
        }
    }
    if(max_wgd > 0){
        free(qmat_wgd);
    }
    if(max_chr_change > 0){
        free(qmat_chr);
    }
    if(max_site_change > 0){
        free(qmat_seg);
    }
    for(auto m : pmats_wgd){
        free(m.second);
    }
    for(auto m : pmats_chr){
        free(m.second);
    }
    for(auto m : pmats_seg){
        free(m.second);
    }

    return logL;
}

// Infer the copy number of the MRCA given a tree at a site, assuming only segment duplication/deletion
double reconstruct_marginal_ancestral_state(evo_tree &rtree, ofstream& fout, double min_asr){
    int debug = 0;
    if(!is_tree_valid(rtree, cons)){
        return SMALL_LNL;
    }

    double logL = 0;    // for all chromosmes
    // For copy number instantaneous changes
    int nstate = cn_max + 1;
    if(model==2) nstate = (cn_max + 1) * (cn_max + 2) / 2;

    vector<int> cn_mrca; // CNs for MRCA
    int nid = 2 * (Ns + 1) - 2;
    if(debug) cout << "\treconstruct marginal ancestral state" << endl;

    string header="node\tsite\tcn\tmax probability";
    for(int i = 0; i < nstate; i++){
        header += "\tprobablity " + to_string(i+1);
    }
    fout << header << endl;

    vector<int> knodes;
    double *qmat = new double[(nstate)*(nstate)];
    memset(qmat, 0, (nstate)*(nstate)*sizeof(double));
    if(model > 0){
      if(debug){
          cout << "Getting rate matrix" << endl;
      }
      if(model==2){
          get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }else{
          get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }
    }

    //create a list of nodes to loop over (only internal nodes), making sure the root is last
    for(int k=Ns+2; k<rtree.ntotn; ++k) knodes.push_back( k );
    knodes.push_back(Ns+1);

    // Find the transition probability matrix for each branch
    map<double, double*> pmats;
    for(int kn=0; kn<knodes.size(); ++kn){
    int k = knodes[kn];
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
    if(model == 1 || model == 2){
        double *pmati = new double[(nstate)*(nstate)];
        double *pmatj = new double[(nstate)*(nstate)];
        memset(pmati, 0, (nstate)*(nstate)*sizeof(double));
        memset(pmatj, 0, (nstate)*(nstate)*sizeof(double));

        if(pmats.count(bli) == 0){
            get_transition_matrix_bounded(qmat, pmati, bli, nstate);
            pmats[bli] = pmati;
        }
        if(pmats.count(blj) == 0){
            get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
            pmats[blj] = pmatj;
        }
     }
    }

    if(debug){
      for( auto it = pmats.begin(); it != pmats.end(); ++it )
      {
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print( nstate, nstate, it->second, "  P matrix:" );
      }
    }

    for(int nchr=1; nchr<=vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "\tComputing likelihood on Chr " << nchr << endl;
      double site_logL = 0;   // log likelihood for all sites on a chromosome
      map<vector<int>, vector<vector<double>>> sites_lnl_map;

      for(int nc=0; nc<vobs[nchr].size(); nc++){    // for each segment on the chromosome
          if(debug) cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k;
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  L_sk_k = initialize_lnl_table(obs, rtree, Ns, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, pmats, 0, 0, model, nstate);
              }else{
                  // cout << "sites repeated" << end1;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              L_sk_k = initialize_lnl_table(obs, rtree, Ns, model, nstate, is_total);
              get_likelihood_site(L_sk_k, rtree, knodes, pmats, 0, 0, model, nstate);
          }
          double lnl = extract_tree_lnl(L_sk_k, Ns, model);
          site_logL += lnl;

          // Get the likelihood table of MRCA node (with largest ID) in the tree from likelihood table
          double max_ln = *max_element(L_sk_k[nid].begin(), L_sk_k[nid].end());
          int state = distance(L_sk_k[nid].begin(), max_element(L_sk_k[nid].begin(), L_sk_k[nid].end()));
          if(debug) cout << "\t\tState with maximum likelihood " << state << " at chr " << nchr << " site " << nc << endl;
          int cn = state;
          if(model == 2){
              cn = state_to_total_cn(state, cn_max);
          }
          cn_mrca.push_back(cn);

          string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc) + "\t" + to_string(cn) + "\t" + to_string(max_ln);
          for(int i = 0; i < L_sk_k[nid].size(); i++){
              line += "\t" + to_string(L_sk_k[nid][i]);
          }
          fout << setprecision(dbl::max_digits10) << line << endl;
      }
      logL += site_logL;
      if(debug){
          cout << "Site Likelihood for " << nchr << " is "  << site_logL << endl;
      }
    } // for each chromosome
    if(debug){
        cout << "\nLikelihood without correcting acquisition bias: " << logL << endl;
        cout << "CNs at MRCA is: " << endl;
        for(int i = 0; i < cn_mrca.size(); i++){
            cout << cn_mrca[i] << endl;
        }
    }

    free(qmat);
    for(auto m : pmats){
        free(m.second);
    }

    return logL;
}



// Infer the copy number of all internal nodes given a tree at a site, assuming independent Markov chains
void reconstruct_joint_ancestral_state_decomp(evo_tree &rtree, ofstream& fout, double min_asr){
    int debug = 0;
    if(debug) cout << "\treconstruct joint ancestral state with independent chain model" << endl;

    if(!is_tree_valid(rtree, cons)){
        return;
    }

    int max_id = 2 * Ns;

    // For WGD model
    int dim_wgd = max_wgd + 1;
    int dim_chr = 2 * max_chr_change + 1;
    int dim_seg = 2 * max_site_change + 1;
    double *qmat_wgd, *qmat_chr, *qmat_seg;
    double *pmat_wgd, *pmat_chr, *pmat_seg;
    // Find the transition probability matrix for each branch
    map<double, double*> pmats_wgd;
    map<double, double*> pmats_chr;
    map<double, double*> pmats_seg;
    double *pmati_wgd, *pmatj_wgd;
    double *pmati_chr, *pmatj_chr;
    double *pmati_seg, *pmatj_seg;

    if(max_wgd > 0){
        qmat_wgd = new double[dim_wgd * dim_wgd];  // WGD
        memset(qmat_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
        get_rate_matrix_wgd(qmat_wgd, rtree.wgd_rate, max_wgd);
    }

    if(max_chr_change > 0){
        qmat_chr = new double[(dim_chr)*(dim_chr)];   // chromosome gain/loss
        memset(qmat_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
        get_rate_matrix_chr_change(qmat_chr, rtree.chr_gain_rate, rtree.chr_loss_rate, max_chr_change);
    }

    if(max_site_change > 0){
        qmat_seg = new double[(dim_seg)*(dim_seg)];  // segment duplication/deletion
        memset(qmat_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
        get_rate_matrix_site_change(qmat_seg, rtree.dup_rate, rtree.del_rate, max_site_change);
    }

    if(debug){
          cout << "Dimension of Qmat " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << "\n";
    }

    //create a list of nodes to loop over (only internal nodes), making sure the root is last
    vector<int> knodes;
    for(int k=Ns+2; k<rtree.ntotn; ++k) knodes.push_back( k );
    knodes.push_back(Ns+1);

    for(int kn=0; kn<knodes.size(); ++kn){
           int k = knodes[kn];
           double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
           double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

           // For WGD
           if(max_wgd > 0){
               pmati_wgd = new double[(dim_wgd)*(dim_wgd)];
               pmatj_wgd = new double[(dim_wgd)*(dim_wgd)];
               memset(pmati_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
               memset(pmatj_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
               if(pmats_wgd.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_wgd, pmati_wgd, bli, dim_wgd);
                   pmats_wgd[bli] = pmati_wgd;
               }
               if(pmats_wgd.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_wgd, pmatj_wgd, blj, dim_wgd);
                   pmats_wgd[blj] = pmatj_wgd;
               }
           }

           // For chr gain/loss
           if(max_chr_change > 0){
               pmati_chr = new double[(dim_chr)*(dim_chr)];
               pmatj_chr = new double[(dim_chr)*(dim_chr)];
               memset(pmati_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
               memset(pmatj_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
               if(pmats_chr.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_chr, pmati_chr, bli, dim_chr);
                   pmats_chr[bli] = pmati_chr;
               }
               if(pmats_chr.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_chr, pmatj_chr, blj, dim_chr);
                   pmats_chr[blj] = pmatj_chr;
               }
           }

           // For segment duplication/deletion
           if(max_site_change > 0){
               pmati_seg = new double[(dim_seg)*(dim_seg)];
               pmatj_seg = new double[(dim_seg)*(dim_seg)];
               memset(pmati_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
               memset(pmatj_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
               if(pmats_seg.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_seg, pmati_seg, bli, dim_seg);
                   pmats_seg[bli] = pmati_seg;
               }
               if(pmats_seg.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_seg, pmatj_seg, blj, dim_seg);
                   pmats_seg[blj] = pmatj_seg;
              }
          }
    }
    if(debug){
        for( auto it = pmats_wgd.begin(); it != pmats_wgd.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print( dim_wgd, dim_wgd, it->second, "  P-WGD matrix:" );
        }
        for( auto it = pmats_chr.begin(); it != pmats_chr.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print( dim_chr, dim_chr, it->second, "  P-CHR matrix:" );
        }
        for( auto it = pmats_seg.begin(); it != pmats_seg.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print( dim_seg, dim_seg, it->second, "  P-SEG matrix:" );
        }
    }

    int nstate = comps.size();
    double logL = 0;    // for all chromosmes

    string header="node\tsite\tcn\tstate\tmax probability";
    for(int i = 0; i < nstate; i++){
        header += "\tprobablity " + to_string(i+1);
    }
    fout << header << endl;

    knodes.pop_back();  // no need to reconstruct root which is always normal

    for(int nchr=1; nchr<=vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with " << vobs[nchr].size() << " sites " << endl;
      map<vector<int>, vector<vector<double>>> sites_lnl_map;
      for(int nc=0; nc<vobs[nchr].size(); nc++){    // for each segment on the chromosome
          if(debug) cout << "\tfor site " << nc << " on chromosome " << nchr << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k(rtree.ntotn, vector<double>(nstate,0));
          vector<vector<int>> S_sk_k(rtree.ntotn, vector<int>(nstate,0));
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  initialize_asr_table_decomp(obs, rtree, pmats_wgd, pmats_chr, pmats_seg, L_sk_k, S_sk_k, Ns, model, nstate, is_total);
                  get_ancestral_states_site_decomp(L_sk_k, S_sk_k, rtree, knodes, pmats_wgd, pmats_chr, pmats_seg, dim_wgd, dim_chr, dim_seg, nstate);
              }else{
                  if(debug) cout << "\t\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              initialize_asr_table_decomp(obs, rtree, pmats_wgd, pmats_chr, pmats_seg, L_sk_k, S_sk_k, Ns, model, nstate, is_total);
              get_ancestral_states_site_decomp(L_sk_k, S_sk_k, rtree, knodes, pmats_wgd, pmats_chr, pmats_seg, dim_wgd, dim_chr, dim_seg, nstate);
          }

          map<int, int> asr_states;
          map<int, vector<int>> asr_cns = extract_tree_ancestral_state(rtree, L_sk_k, S_sk_k, Ns, model, asr_states);
          for(int nid = max_id; nid > Ns + 1; nid--){
              vector<int> cns = asr_cns[nid];
              int state = asr_states[nid];
              string line = to_string(nid+1) + "\t" + to_string(nchr) + "_" + to_string(nc);
              if(cns.size()>1){
                  line += "\t" + to_string(cns[0]) + "," + to_string(cns[1]);
              }else{
                line += "\t" + to_string(cns[0]);
              }
              line += "\t" + to_string(state) + "\t" + to_string(pow(10, L_sk_k[nid][state]));
              for(int i = 0; i < L_sk_k[nid].size(); i++){
                  line += "\t" + to_string(pow(10, L_sk_k[nid][i]));
              }
              fout << setprecision(dbl::max_digits10) << line << endl;
          }
      }
    } // for each chromosome

    if(max_wgd > 0){
        free(qmat_wgd);
    }
    if(max_chr_change > 0){
        free(qmat_chr);
    }
    if(max_site_change > 0){
        free(qmat_seg);
    }
    for(auto m : pmats_wgd){
        free(m.second);
    }
    for(auto m : pmats_chr){
        free(m.second);
    }
    for(auto m : pmats_seg){
        free(m.second);
    }

}


// Infer the copy number of all internal nodes given a tree at a site, assuming only segment duplication/deletion
void reconstruct_joint_ancestral_state(evo_tree &rtree, ofstream &fout, double min_asr){
    int debug = 0;
    if(debug) cout << "\treconstruct joint ancestral state" << endl;

    if(!is_tree_valid(rtree, cons)){
        return;
    }

    // For copy number instantaneous changes
    int nstate = cn_max + 1;
    if(model==2) nstate = (cn_max + 1) * (cn_max + 2) / 2;
    int max_id = 2 * Ns;

    string header="node\tsite\tcn\tstate\tmax probability";
    for(int i = 0; i < nstate; i++){
        header += "\tprobablity " + to_string(i+1);
    }
    fout << header << endl;

    vector<int> knodes;
    double *qmat = new double[(nstate)*(nstate)];
    memset(qmat, 0, (nstate)*(nstate)*sizeof(double));
    if(model > 0){
        if(debug){
            cout << "Getting rate matrix" << endl;
        }
        if(model==2){
            get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
        }else{
            get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
        }
    }

    //create a list of nodes to loop over (only internal nodes)
    for(int k=Ns+2; k<rtree.ntotn; ++k) knodes.push_back( k );
    knodes.push_back(Ns+1);     // needed for computing P-matrix

    // Find the transition probability matrix for each branch
    map<double, double*> pmats;
    for(int kn=0; kn<knodes.size(); ++kn){
      int k = knodes[kn];
      double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
      double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
      if(model == 1 || model == 2){
          double *pmati = new double[(nstate)*(nstate)];
          double *pmatj = new double[(nstate)*(nstate)];
          memset(pmati, 0, (nstate)*(nstate)*sizeof(double));
          memset(pmatj, 0, (nstate)*(nstate)*sizeof(double));

          if(pmats.count(bli) == 0){
              get_transition_matrix_bounded(qmat, pmati, bli, nstate);
              pmats[bli] = pmati;
          }
          if(pmats.count(blj) == 0){
              get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
              pmats[blj] = pmatj;
          }
       }
    }

    if(debug){
        for( auto it = pmats.begin(); it != pmats.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print( nstate, nstate, it->second, "  P matrix:" );
        }
    }

    knodes.pop_back();  // no need to reconstruct root which is always normal

    for(int nchr=1; nchr<=vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with " << vobs[nchr].size() << " sites " << endl;
      map<vector<int>, vector<vector<double>>> sites_lnl_map;
      for(int nc=0; nc<vobs[nchr].size(); nc++){    // for each segment on the chromosome
          if(debug) cout << "\tfor site " << nc << " on chromosome " << nchr << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k(rtree.ntotn, vector<double>(nstate,0));
          vector<vector<int>> S_sk_k(rtree.ntotn, vector<int>(nstate,0));
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  initialize_asr_table(obs, rtree, pmats, L_sk_k, S_sk_k, Ns, model, nstate, is_total);
                  get_ancestral_states_site(L_sk_k, S_sk_k, rtree, knodes, pmats, nstate);
              }else{
                  if(debug) cout << "\t\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              initialize_asr_table(obs, rtree, pmats, L_sk_k, S_sk_k, Ns, model, nstate, is_total);
              get_ancestral_states_site(L_sk_k, S_sk_k, rtree, knodes, pmats, nstate);
          }

          map<int, int> asr_states;     // The state ID used in rate matrix
          map<int, vector<int>> asr_cns = extract_tree_ancestral_state(rtree, L_sk_k, S_sk_k, Ns, model, asr_states);
          for(int nid = max_id; nid > Ns + 1; nid--){
              vector<int> cns = asr_cns[nid];
              int state = asr_states[nid];
              string line = to_string(nid+1) + "\t" + to_string(nchr) + "_" + to_string(nc);
              if(cns.size()>1){
                  line += "\t" + to_string(cns[0]) + "," + to_string(cns[1]);
              }else{
                line += "\t" + to_string(cns[0]);
              }
              line += "\t" + to_string(state) + "\t" + to_string(pow(10, L_sk_k[nid][state]));
              for(int i = 0; i < L_sk_k[nid].size(); i++){
                  line += "\t" + to_string(pow(10, L_sk_k[nid][i]));
              }
              fout << setprecision(dbl::max_digits10) << line << endl;
          }
      }
    } // for each chromosome

    free(qmat);
    for(auto m : pmats){
        free(m.second);
    }

}


// Get the likelihood on a set of chromosmes
double get_likelihood_chr(map<int, vector<vector<int>>>& vobs, evo_tree &rtree, const vector<int>& knodes, map<double, double*>& pmats, int has_wgd, int model, int nstate, int only_seg, int is_total=1){
    int debug = 0;
    double logL = 0;    // for all chromosmes
    double chr_gain = 0;
    double chr_loss = 0;

    for(int nchr=1; nchr<=vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "Computing likelihood on Chr " << nchr << endl;
      double chr_logL = 0;  // for one chromosome
      double chr_logL_normal = 0, chr_logL_gain = 0, chr_logL_loss = 0;
      double site_logL = 0;   // log likelihood for all sites on a chromosome
      int z = 0;    // no chr gain/loss
      // cout << " chromosome number change is " << 0 << endl;
      // Use a map to store computed log likelihood
      map<vector<int>, vector<vector<double>>> sites_lnl_map;

      for(int nc=0; nc<vobs[nchr].size(); nc++){    // for each segment on the chromosome
          // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k;
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  L_sk_k = initialize_lnl_table(obs, rtree, Ns, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, pmats, has_wgd, z, model, nstate);
              }else{
                  // cout << "sites repeated" << end1;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              L_sk_k = initialize_lnl_table(obs, rtree, Ns, model, nstate, is_total);
              get_likelihood_site(L_sk_k, rtree, knodes, pmats, has_wgd, z, model, nstate);
          }
          double lnl = extract_tree_lnl(L_sk_k, Ns, model);
          site_logL += lnl;

          if(debug){
              cout << "\nLikelihood for site " << nc << " is " << lnl << endl;
              print_tree_lnl(rtree, L_sk_k, nstate);
          }
      }

      double chr_normal = 1;
      if(!only_seg){
          chr_gain = rtree.chr_gain_rate;
          chr_loss = rtree.chr_loss_rate;
          if(debug){
              cout << "chromosome gain rate " << chr_gain << endl;
              cout << "chromosome loss rate " << chr_loss << endl;
              cout << "Number of chr so far " << vobs.size() << endl;
          }

          if(fabs(chr_loss-0) > SMALL_VAL){
             chr_normal -= chr_loss;
          }
          if(fabs(chr_gain-0) > SMALL_VAL){
             chr_normal -= chr_gain;
          }
      }

      chr_logL_normal = log(chr_normal) + site_logL;
      chr_logL += chr_logL_normal;
      if(debug){
          cout << "Likelihood without chr gain/loss for " << nchr << " is "  << chr_normal << endl;
          cout << "Site Likelihood for " << nchr << " is "  << site_logL << endl;
          cout << "Likelihood without chr gain/loss: " << chr_logL_normal << endl;
      }

      if(!only_seg){
          if(fabs(chr_loss-0) > SMALL_VAL){
              z = -1;
              double site_logL = 0;   // log likelihood for all sites on a chromosome
              // cout << " chromosome number change is " << z << endl;
              for(int nc=0; nc<vobs[nchr].size(); nc++){
                  // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
                  // for each site of the chromosome
                  vector<int> obs = vobs[nchr][nc];
                  vector<vector<double>> L_sk_k = initialize_lnl_table(obs, rtree, Ns, model, nstate, is_total);

                  get_likelihood_site(L_sk_k, rtree, knodes, pmats, has_wgd, z, model, nstate);
                  site_logL += extract_tree_lnl(L_sk_k, Ns, model);

                  if(debug){
                      print_tree_lnl(rtree, L_sk_k, nstate);
                  }
              } // for all sites on a chromosome

              chr_logL_loss = log(chr_loss) + site_logL;
              chr_logL +=  log(1 + exp(chr_logL_loss-chr_logL_normal));
              if(debug){
                  cout << "\nLikelihood before chr loss for " << nchr << " is " << site_logL << endl;
                  cout << "\nLikelihood after chr loss: " << chr_logL_loss << endl;
              }
          } // for all chromosome loss

          if(fabs(chr_gain-0) > SMALL_VAL){
              z = 1;
              double site_logL = 0;   // log likelihood for all sites on a chromosome
              // cout << " chromosome number change is " << z << endl;
              for(int nc=0; nc<vobs[nchr].size(); nc++){
                  // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
                  // for each site of the chromosome
                  vector<int> obs = vobs[nchr][nc];
                  vector<vector<double>> L_sk_k = initialize_lnl_table(obs, rtree, Ns, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, pmats, has_wgd, z, model, nstate);
                  site_logL += extract_tree_lnl(L_sk_k, Ns, model);

                  if(debug){
                      print_tree_lnl(rtree, L_sk_k, nstate);
                  }
              } // for all sites on a chromosome

              chr_logL_gain = log(chr_gain) + site_logL;
              if(chr_logL_loss > 0){
                  chr_logL += log(1 + 1 / (exp(chr_logL_normal-chr_logL_gain) + exp(chr_logL_loss-chr_logL_gain)));
              }
              else{
                  chr_logL += log(1 + exp(chr_logL_gain-chr_logL_normal));
              }

              if(debug){
                  cout << "\nLikelihood before chr gain for " << nchr << " is " << site_logL << endl;
                  cout << "\nLikelihood after chr gain: " << chr_logL_gain << endl;
              }
          } // for all chromosome loss
          // chr_logL = chr_logL_normal + log(1 + exp(chr_logL_loss-chr_logL_normal)) + log(1 + 1 / (exp(chr_logL_normal-chr_logL_gain) + exp(chr_logL_loss-chr_logL_gain)));
      }
      logL += chr_logL;
      if(debug){
          cout << "\nLikelihood after considering chr gain/loss for  " << nchr << " is " << logL << endl;
      }
    } // for each chromosome
    if(debug){
        cout << "\nLikelihood with chr gain/loss for all chromosmes: " << logL << endl;
    }
    return logL;
}

// Used when WGD is considered, dealing with mutations of different types at different levels
double get_likelihood_chr_decomp(map<int, vector<vector<int>>>& vobs, evo_tree &rtree, map<int, set<vector<int>>>& decomp_table, const vector<int>& knodes, map<double, double*>& pmats_wgd, map<double, double*>& pmats_chr, map<double, double*>& pmats_seg, int dim_wgd, int dim_chr, int dim_seg, int cn_max, int nstate, int is_total=1){
    int debug = 0;
    double logL = 0;    // for all chromosmes

    for(int nchr=1; nchr<=vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "Computing likelihood on Chr " << nchr << endl;
      double site_logL = 0;   // log likelihood for all sites on a chromosome
      // Use a map to store computed log likelihood
      map<vector<int>, vector<vector<double>>> sites_lnl_map;
      // cout << " chromosome number change is " << 0 << endl;
      for(int nc=0; nc<vobs[nchr].size(); nc++){    // for each segment on the chromosome
          // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k;
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  L_sk_k = initialize_lnl_table_decomp(obs, nchr, rtree, comps, Ns, model, cn_max, nstate, is_total);
                  get_likelihood_site_decomp(L_sk_k, rtree, decomp_table, knodes, pmats_wgd, pmats_chr, pmats_seg, dim_wgd, dim_chr, dim_seg, cn_max, nstate, is_total);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  if(debug) cout << "\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              L_sk_k = initialize_lnl_table_decomp(obs, nchr, rtree, comps, Ns, model, cn_max, nstate, is_total);
              get_likelihood_site_decomp(L_sk_k, rtree, decomp_table, knodes, pmats_wgd, pmats_chr, pmats_seg, dim_wgd, dim_chr, dim_seg, cn_max, nstate, is_total);
          }
          // site_logL += extract_tree_lnl(L_sk_k, Ns, model);
          double lnl = extract_tree_lnl_decomp(L_sk_k, comps, Ns);
          site_logL += lnl;

          if(debug){
              cout << "\nLikelihood for site " << nc << " is " << lnl << endl;
              cout << "CNs at this site: ";
              for(int i = 0; i < obs.size(); i++){
                  cout << "\t" << obs[i];
              }
              cout << endl;
              print_tree_lnl(rtree, L_sk_k, nstate);
          }
      }
      logL += site_logL;
      if(debug){
          cout << "\nLikelihood for chromosome " << nchr << " is " << site_logL << endl;
      }
    } // for each chromosome
    return logL;
}


// Compute the likelihood of dummy sites consisting entirely of 2s for the tree
double get_likelihood_invariant(evo_tree &rtree, const vector<int>& knodes, map<double, double*>& pmats, const int model, const int cn_max, int is_total=1){
    int debug = 0;
    int nstate = cn_max + 1;
    if(model == 2){
        nstate = (cn_max + 1) * (cn_max + 2) / 2;
    }
    double logL = 0;
    if(debug) cout << "Correcting for the skip of invariant sites" << endl;

    // Suppose the value is 2 for all samples
    int normal_cn  = 2;
    if(is_total == 0){
        normal_cn  = 4;
    }
    vector<int> obs(Ns, normal_cn);
    vector<vector<double>> L_sk_k = initialize_lnl_table(obs, rtree, Ns, model, nstate, is_total);

    for(int kn=0; kn<knodes.size(); ++kn){
        int k = knodes[kn];
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

        //loop over possible values of sk
        for(int sk=0; sk<nstate; ++sk){
          double Li = 0;
          // loop over possible si
          for(int si=0; si<nstate; ++si){
              if (model == 0){
                Li += get_transition_prob(rtree.mu, bli, sk, si) * L_sk_k[ni][si];
              }
              else{
                Li += get_transition_prob_bounded(pmats[bli], sk, si, nstate) * L_sk_k[ni][si];
              }
          }
          double Lj = 0;
          // loop over possible sj
          for(int sj=0; sj<nstate; ++sj){
              if (model == 0){
                   Lj += get_transition_prob(rtree.mu, blj, sk, sj) * L_sk_k[nj][sj];
              }
              else{
                   Lj += get_transition_prob_bounded(pmats[blj], sk, sj, nstate) * L_sk_k[nj][sj];
              }
          }

          L_sk_k[k][sk] = Li*Lj;
       }
    }

    logL = extract_tree_lnl(L_sk_k, Ns, model);
    // // assert(exp(logL) < 1);
    // long lkl = exp(logL);
    // double bias = 0;
    // if(lkl < 1)
    //     bias = log(1-lkl);
    if(debug){
        // cout << "The log likelihood of invariant sites is " << logL << endl;
        cout << "The likelihood of invariant sites is " << logL << endl;
        // cout << "The bias of invariant sites is " << bias << endl;
    }
    return logL;
}


// Compute the likelihood of dummy sites consisting entirely of 2s for the tree
double get_likelihood_invariant_decomp(evo_tree &rtree, map<int, set<vector<int>>>& decomp_table, const vector<int>& knodes, map<double, double*>& pmats_wgd, map<double, double*>& pmats_chr, map<double, double*>& pmats_seg, int dim_wgd, int dim_chr, int dim_seg, int cn_max, int nstate, int is_total=1){
    int debug = 0;
    double logL = 0;
    if(debug) cout << "Correcting for the skip of invariant sites" << endl;

    // Suppose the value is 2 for all samples
    int normal_cn  = 2;
    vector<int> obs(Ns, normal_cn);
    vector<vector<double>> L_sk_k = initialize_lnl_table_decomp(obs, 0, rtree, comps, Ns, model, cn_max, nstate, is_total);
    get_likelihood_site_decomp(L_sk_k, rtree, decomp_table, knodes, pmats_wgd, pmats_chr, pmats_seg, dim_wgd, dim_chr, dim_seg, cn_max, nstate, is_total);
    logL = extract_tree_lnl_decomp(L_sk_k, comps, Ns);

    if(debug){
        cout << "The likelihood of invariant sites is " << logL << endl;
    }
    return logL;
}


// Incorporate chromosome gain/loss and WGD
// Model 2: Treat total copy number as the observed data and the allele-specific information is missing
double get_likelihood_revised(evo_tree &rtree){
// double get_likelihood_revised(int Ns, int Nchar, int num_invar_bins, map<int, vector<vector<int>>>& vobs, evo_tree &rtree, int model, int cons, int cn_max, int only_seg, int correct_bias = 0, int is_total = 1){
  int debug = 0;
  if(debug) cout << "\tget_likelihood by matrix exponential" << endl;

  if(!is_tree_valid(rtree, cons)){
      return SMALL_LNL;
  }

  // For copy number instantaneous changes
  int nstate = cn_max + 1;
  if(model==2) nstate = (cn_max + 1) * (cn_max + 2) / 2;

  vector<int> knodes;

  double *qmat = new double[(nstate)*(nstate)];
  memset(qmat, 0, (nstate)*(nstate)*sizeof(double));
  if(model > 0){
      if(debug){
          cout << "Getting rate matrix" << endl;
      }
      if(model==2){
          get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }else{
          get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }
  }

  //create a list of nodes to loop over (only internal nodes), making sure the root is last
  for(int k=Ns+2; k<rtree.ntotn; ++k) knodes.push_back( k );
  knodes.push_back(Ns+1);

  // Find the transition probability matrix for each branch
  map<double, double*> pmats;
  for(int kn=0; kn<knodes.size(); ++kn){
    int k = knodes[kn];
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
    if(model == 1 || model == 2){
        double *pmati = new double[(nstate)*(nstate)];
        double *pmatj = new double[(nstate)*(nstate)];
        memset(pmati, 0, (nstate)*(nstate)*sizeof(double));
        memset(pmatj, 0, (nstate)*(nstate)*sizeof(double));

        if(pmats.count(bli) == 0){
            get_transition_matrix_bounded(qmat, pmati, bli, nstate);
            pmats[bli] = pmati;
        }
        if(pmats.count(blj) == 0){
            get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
            pmats[blj] = pmatj;
        }
     }
     // if(debug){
     //      cout << "Get Pmatrix for 1st branch length " << bli << endl;
     //      r8mat_print( nstate, nstate, pmats[bli], "  P matrix after change in map:" );
     //      r8mat_print( nstate, nstate, pmati, "  P matrix after change:" );
     //      cout << "Get Pmatrix for 2nd branch length " << blj << endl;
     //      r8mat_print( nstate, nstate, pmats[blj], "  P matrix after change in map:" );
     //      r8mat_print( nstate, nstate, pmatj, "  P matrix after change:" );
     // }
  }

  if(debug){
      for( auto it = pmats.begin(); it != pmats.end(); ++it )
      {
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print( nstate, nstate, it->second, "  P matrix:" );
      }
  }

  double logL = 0;

  // cout << "WGD rate is " << wgd << endl;
  if(only_seg){
      if(debug) cout << "Computing the likelihood without consideration of WGD" << endl;
      logL += get_likelihood_chr(vobs, rtree, knodes, pmats, 0, model, nstate, only_seg, is_total);
  }
  else{
      if(debug) cout << "Computing the likelihood with consideration of WGD" << endl;
      logL += (1 - rtree.wgd_rate) * get_likelihood_chr(vobs, rtree, knodes, pmats, 0, model, nstate, only_seg, is_total);
      logL += rtree.wgd_rate * get_likelihood_chr(vobs, rtree, knodes, pmats, 1, model, nstate, only_seg, is_total);
  }


  if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;
  if(correct_bias == 1){
      double lnl_invar = get_likelihood_invariant(rtree, knodes, pmats, model, cn_max, is_total);
      double bias = num_invar_bins * lnl_invar;
      logL = logL + bias;
      if(debug){
          cout << "Likelihood of an invariant bin " << lnl_invar << endl;
          cout << "Number of invariant bins " << num_invar_bins << endl;
          cout << "Bias to correct " << bias << endl;
          cout << "Final likelihood after correcting acquisition bias: " << logL << endl;
      }
  }

  if(std::isnan(logL) || logL < SMALL_LNL) logL = SMALL_LNL;
  if(debug){
      cout << "Final likelihood: " << logL << endl;
      cout << "Free memory" << endl;
  }

  free(qmat);
  for(auto m : pmats){
      free(m.second);
  }

  return logL;
}

// Computing likelihood when WGD and chr gain/loss are incorporated
// Assume likelihood is for allele-specific information
// double get_likelihood_decomp(int Ns, int Nchar, int num_invar_bins, map<int, vector<vector<int>>>& vobs, evo_tree &rtree, map<int, set<vector<int>>>& decomp_table, int cons, int cn_max, int max_wgd, int max_chr_change, int max_site_change, int correct_bias = 0, int is_total = 1){
double get_likelihood_decomp(evo_tree &rtree){
  int debug = 0;
  if(debug) cout << "\tget_likelihood from multiple chains" << endl;

  if(!is_tree_valid(rtree, cons)){
      return SMALL_LNL;
  }

  // For WGD model
  int dim_wgd = max_wgd + 1;
  int dim_chr = 2 * max_chr_change + 1;
  int dim_seg = 2 * max_site_change + 1;
  double *qmat_wgd, *qmat_chr, *qmat_seg;
  double *pmat_wgd, *pmat_chr, *pmat_seg;
  // Find the transition probability matrix for each branch
  map<double, double*> pmats_wgd;
  map<double, double*> pmats_chr;
  map<double, double*> pmats_seg;
  double *pmati_wgd, *pmatj_wgd;
  double *pmati_chr, *pmatj_chr;
  double *pmati_seg, *pmatj_seg;

  if(max_wgd > 0){
      qmat_wgd = new double[dim_wgd * dim_wgd];  // WGD
      memset(qmat_wgd, 0, (dim_wgd) * (dim_wgd) * sizeof(double));
      get_rate_matrix_wgd(qmat_wgd, rtree.wgd_rate, max_wgd);
  }

  if(max_chr_change > 0){
      qmat_chr = new double[(dim_chr) * (dim_chr)];   // chromosome gain/loss
      memset(qmat_chr, 0, (dim_chr) * (dim_chr) * sizeof(double));
      get_rate_matrix_chr_change(qmat_chr, rtree.chr_gain_rate, rtree.chr_loss_rate, max_chr_change);
  }

  if(max_site_change > 0){
      qmat_seg = new double[(dim_seg) * (dim_seg)];  // segment duplication/deletion
      memset(qmat_seg, 0, (dim_seg) * (dim_seg) * sizeof(double));
      get_rate_matrix_site_change(qmat_seg, rtree.dup_rate, rtree.del_rate, max_site_change);
  }

  if(debug){
        cout << "Dimension of Qmat " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << "\n";
  }

  //create a list of nodes to loop over (only internal nodes), making sure the root is last
  vector<int> knodes;
  for(int k = Ns+2; k < rtree.ntotn; ++k) knodes.push_back( k );
  knodes.push_back(Ns+1);

  for(int kn=0; kn<knodes.size(); ++kn){
         int k = knodes[kn];
         double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
         double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

         // For WGD
         if(max_wgd > 0){
             pmati_wgd = new double[(dim_wgd)*(dim_wgd)];
             pmatj_wgd = new double[(dim_wgd)*(dim_wgd)];
             memset(pmati_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
             memset(pmatj_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
             if(pmats_wgd.count(bli) == 0){
                 get_transition_matrix_bounded(qmat_wgd, pmati_wgd, bli, dim_wgd);
                 pmats_wgd[bli] = pmati_wgd;
             }
             if(pmats_wgd.count(blj) == 0){
                 get_transition_matrix_bounded(qmat_wgd, pmatj_wgd, blj, dim_wgd);
                 pmats_wgd[blj] = pmatj_wgd;
             }
         }

         // For chr gain/loss
         if(max_chr_change > 0){
             pmati_chr = new double[(dim_chr)*(dim_chr)];
             pmatj_chr = new double[(dim_chr)*(dim_chr)];
             memset(pmati_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
             memset(pmatj_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
             if(pmats_chr.count(bli) == 0){
                 get_transition_matrix_bounded(qmat_chr, pmati_chr, bli, dim_chr);
                 pmats_chr[bli] = pmati_chr;
             }
             if(pmats_chr.count(blj) == 0){
                 get_transition_matrix_bounded(qmat_chr, pmatj_chr, blj, dim_chr);
                 pmats_chr[blj] = pmatj_chr;
             }
         }

         // For segment duplication/deletion
         if(max_site_change > 0){
             pmati_seg = new double[(dim_seg)*(dim_seg)];
             pmatj_seg = new double[(dim_seg)*(dim_seg)];
             memset(pmati_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
             memset(pmatj_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
             if(pmats_seg.count(bli) == 0){
                 get_transition_matrix_bounded(qmat_seg, pmati_seg, bli, dim_seg);
                 pmats_seg[bli] = pmati_seg;
             }
             if(pmats_seg.count(blj) == 0){
                 get_transition_matrix_bounded(qmat_seg, pmatj_seg, blj, dim_seg);
                 pmats_seg[blj] = pmatj_seg;
            }
        }
  }
  if(debug){
      for( auto it = pmats_wgd.begin(); it != pmats_wgd.end(); ++it )
      {
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print( dim_wgd, dim_wgd, it->second, "  P-WGD matrix:" );
      }
      for( auto it = pmats_chr.begin(); it != pmats_chr.end(); ++it )
      {
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print( dim_chr, dim_chr, it->second, "  P-CHR matrix:" );
      }
      for( auto it = pmats_seg.begin(); it != pmats_seg.end(); ++it )
      {
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print( dim_seg, dim_seg, it->second, "  P-SEG matrix:" );
      }
  }

  double logL = 0;

  int nstate = comps.size();
  // cout << "Number of states is " << nstate << endl;
  logL = get_likelihood_chr_decomp(vobs, rtree, decomp_table, knodes, pmats_wgd, pmats_chr, pmats_seg, dim_wgd, dim_chr, dim_seg, cn_max, nstate, is_total);

  if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;
  if(correct_bias == 1){
      double lnl_invar = get_likelihood_invariant_decomp(rtree, decomp_table, knodes, pmats_wgd, pmats_chr, pmats_seg, dim_wgd, dim_chr, dim_seg, cn_max, nstate, is_total);
      double bias = num_invar_bins * lnl_invar;
      logL = logL + bias;
      if(debug){
          cout << "Likelihood of an invariant bin " << lnl_invar << endl;
          cout << "Number of invariant bins " << num_invar_bins << endl;
          cout << "Bias to correct " << bias << endl;
          cout << "Final likelihood after correcting acquisition bias: " << logL << endl;
      }
  }

  if(std::isnan(logL) || logL < SMALL_LNL) logL = SMALL_LNL;
  if(debug){
      cout << "Final likelihood: " << logL << endl;
      cout << "Free memory" << endl;
  }

  if(max_wgd > 0){
      free(qmat_wgd);
      for(auto m : pmats_wgd){
          free(m.second);
      }
  }
  if(max_chr_change > 0){
      free(qmat_chr);
      for(auto m : pmats_chr){
          free(m.second);
      }
  }
  if(max_site_change > 0){
      free(qmat_seg);
      for(auto m : pmats_seg){
          free(m.second);
      }
  }

  return logL;
}


// Compute the likelihood without grouping sites by chromosome, only considering segment duplication/deletion
/*
Ns: number of samples
Nchar: number of characters for each sample
vobs: the observed data matrix
rtree: the given tree
model: model of evolution
*/
double get_likelihood(int Ns, int Nchar, const vector<vector<int>>& vobs, evo_tree &rtree, int model, int cons, const int cn_max, int is_total=1){
  int debug = 0;
  if(debug) cout << "\tget_likelihood" << endl;
  int nstate = cn_max + 1;
  if(model==2) nstate = (cn_max + 1) * (cn_max + 2) / 2;

  // return 0 if the tree is not valid
  if(!is_tree_valid(rtree, cons)){
      return SMALL_LNL;
  }

  double *qmat = new double[(nstate)*(nstate)];
  memset(qmat, 0, (nstate)*(nstate)*sizeof(double));

  if(model > 0){
      if(debug){
          cout << "Getting rate matrix" << endl;
      }
      if(model==2){
          get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }else{
          get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }
  }

  double *pmati = new double[(nstate)*(nstate)];
  double *pmatj = new double[(nstate)*(nstate)];
  memset(pmati, 0, (nstate)*(nstate)*sizeof(double));
  memset(pmatj, 0, (nstate)*(nstate)*sizeof(double));

  double logL = 0;
  for(int nc=0; nc<Nchar; ++nc){
    vector<int> obs = vobs[nc];
    if(debug){
      cout << "char: " << nc << endl;
      for(int i=0; i<Ns; ++i) cout << "\t" << obs[i];
      cout << endl;
    }

    vector<vector<double>> L_sk_k = initialize_lnl_table(obs, rtree, Ns, model, nstate, is_total);
    if(debug){
        print_tree_lnl(rtree, L_sk_k, nstate);
    }

    //create a list of nodes to loop over, making sure the root is last
    vector<int> knodes;
    for(int k=Ns+2; k<rtree.ntotn; ++k) knodes.push_back( k );
    knodes.push_back(Ns+1);

    for(int kn=0; kn<knodes.size(); ++kn){
      int k = knodes[kn];
      int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
      double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
      int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
      double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
      if(model > 0){
           get_transition_matrix_bounded(qmat, pmati, bli, nstate);
           get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
           if(debug){
               cout << "Get Pmatrix for branch length " << bli << endl;
               r8mat_print( nstate, nstate, pmati, "  P matrix after change:" );
               cout << "Get Pmatrix for branch length " << blj << endl;
               r8mat_print( nstate, nstate, pmatj, "  P matrix after change:" );
           }
      }

      if(debug) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni+1 << " , " << bli << "\t" <<  nj+1 << " , " << blj << endl;

      //loop over possible values of sk
      for(int sk=0; sk<nstate; ++sk){
    	double Li = 0;
    	// loop over possible si
    	for(int si=0; si<nstate; ++si){
            if (model == 0){
                Li += get_transition_prob(rtree.mu, bli, sk, si) * L_sk_k[ni][si];
            }
            else{
                Li += get_transition_prob_bounded(pmati, sk, si, nstate) * L_sk_k[ni][si];
            }
    	  //cout << "\tscoring: Li\t" << li << "\t" << get_transition_prob(mu, bli, sk, si ) << "\t" << L_sk_k[ni][si] << endl;
        }
    	double Lj = 0;
    	// loop over possible sj
    	for(int sj=0; sj<nstate; ++sj){
            if (model == 0){
    	         Lj += get_transition_prob(rtree.mu, blj, sk, sj) * L_sk_k[nj][sj];
            }
            else{
                 Lj += get_transition_prob_bounded(pmatj, sk, sj, nstate) * L_sk_k[nj][sj];
            }
    	}
	    //cout << "scoring: sk" << sk << "\t" <<  Li << "\t" << Lj << endl;
	    L_sk_k[k][sk] = Li*Lj;
      }

      if(debug){
    	print_tree_lnl(rtree, L_sk_k, nstate);
      }
    }

    logL += extract_tree_lnl(L_sk_k, Ns, model);
  }

  if(debug) cout << "Final likelihood: " << logL << endl;

  free(qmat);
  free(pmati);
  free(pmatj);

  return logL;
}



//
// Maximum likelihood and optimisation functions
//

double my_f (const gsl_vector *v, void *params){
  evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current branch lengths
  vector<edge> enew;
  for(int i=0; i<(tree->nedge); ++i){
    enew.push_back( tree->edges[i] );
  }
  for(int i=0; i<(tree->nedge)-1; ++i){
    enew[i].length = exp( gsl_vector_get(v, i) );
  }
  evo_tree new_tree(Ns+1, enew);
  if(model==0){
      new_tree.mu = tree->mu;
  }
  else{
      new_tree.dup_rate = tree->dup_rate;
      new_tree.del_rate = tree->del_rate;
      new_tree.chr_gain_rate = tree->chr_gain_rate;
      new_tree.chr_loss_rate = tree->chr_loss_rate;
      new_tree.wgd_rate = tree->wgd_rate;
  }

  // return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, model, 0, cn_max);
  // return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, new_tree, model, 0, cn_max, only_seg, correct_bias, is_total);
  cons = 0;
  return -1.0*get_likelihood_revised(new_tree);
}

double my_f_mu (const gsl_vector *v, void *params){
  //cout << "my_f_mu" << endl;
  evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current branch lengths
  vector<edge> enew;
  for(int i=0; i<(tree->nedge); ++i){
    enew.push_back( tree->edges[i] );
  }
  for(int i=0; i<(tree->nedge)-1; ++i){
    enew[i].length = exp( gsl_vector_get(v, i) );
  }
  evo_tree new_tree(Ns+1, enew);
  if(model==0){
      new_tree.mu = exp( gsl_vector_get(v,tree->nedge-1) );  // 0 to nedge-2 are epars, nedge-1 is mu
  }
  else{
      new_tree.dup_rate = exp( gsl_vector_get(v,tree->nedge-1) );
      new_tree.del_rate = exp( gsl_vector_get(v,tree->nedge) );
      new_tree.chr_gain_rate = exp( gsl_vector_get(v,tree->nedge+1) );
      new_tree.chr_loss_rate = exp( gsl_vector_get(v,tree->nedge+2) );
      new_tree.wgd_rate = exp( gsl_vector_get(v,tree->nedge+3) );
  }

  // return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, model, 0, cn_max);
  // return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, new_tree, model, 0, cn_max, only_seg, correct_bias, is_total);
  cons = 0;
  return -1.0*get_likelihood_revised(new_tree);
}

double my_f_cons (const gsl_vector *v, void *params){
  evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current parameters observing timing constraints
  vector<edge> enew;
  for(int i=0; i<(tree->nedge); ++i){
    enew.push_back( tree->edges[i] );
  }

  // parameters coming in are internal branch lengths followed by total time
  int count = 0;
  for(int i=0; i<(tree->nedge)-1; ++i){
    if(enew[i].end > Ns){
      enew[i].length = exp( gsl_vector_get(v, count) );
      count++;
    }else{
      enew[i].length = 0;
    }
  }

  evo_tree new_tree(tree->nleaf, enew, exp( gsl_vector_get(v, count)), tree->tobs);
  if(model==0){
      new_tree.mu = tree->mu;
  }
  else{
      new_tree.dup_rate = tree->dup_rate;
      new_tree.del_rate = tree->del_rate;
      new_tree.chr_gain_rate = tree->chr_gain_rate;
      new_tree.chr_loss_rate = tree->chr_loss_rate;
      new_tree.wgd_rate = tree->wgd_rate;
  }

  // return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, model, 1, cn_max);
  // return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, new_tree, model, 1, cn_max, only_seg, correct_bias, is_total);
  cons = 1;
  return -1.0*get_likelihood_revised(new_tree);
}

double my_f_cons_mu (const gsl_vector *v, void *params){
  evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current parameters observing timing constraints
  vector<edge> enew;
  for(int i=0; i<(tree->nedge); ++i){
    enew.push_back( tree->edges[i] );
  }

  // parameters coming in are internal branch lengths followed by mutation rate
  int count = 0;
  for(int i=0; i<(tree->nedge)-1; ++i){
    if(enew[i].end > Ns){
      enew[i].length = exp( gsl_vector_get(v, count) );
      count++;
    }else{
      enew[i].length = 0;
    }
  }

  evo_tree new_tree(tree->nleaf, enew, exp( gsl_vector_get(v, count)), tree->tobs);
  if(model==0){
      new_tree.mu = exp( gsl_vector_get(v, count+1) );
  }
  else{
      new_tree.dup_rate = exp( gsl_vector_get(v,count+1) );
      new_tree.del_rate = exp( gsl_vector_get(v,count+2) );
      new_tree.chr_gain_rate = exp( gsl_vector_get(v,count+3) );
      new_tree.chr_loss_rate = exp( gsl_vector_get(v,count+4) );
      new_tree.wgd_rate = exp( gsl_vector_get(v,count+5) );
  }

  // return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, model, 1, cn_max);
  // return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, new_tree, model, 1, cn_max, only_seg, correct_bias, is_total);
  cons = 1;
  return -1.0*get_likelihood_revised(new_tree);
}

// Output error information without aborting
void my_err_handler(const char * reason, const char * file, int line, int gsl_errno){
    fprintf (stderr, "failed, %s, reason: %s, file %s line %d \n", gsl_strerror(gsl_errno), reason, file, line);
    return;
}


void get_variables(evo_tree &rtree, int model, int cons, int maxj, double *x){
    int debug = 0;
    // create a new tree from current value of parameters
    vector<edge> enew;
    for(int i=0; i<rtree.nedge; ++i){
      enew.push_back(rtree.edges[i]);
    }

    if(cons == 0){
      // The first element of x is not used for optimization
      // The index of x is added by 1 compared with index for simplex method
      for(int i=0; i<rtree.nedge-1; ++i){
          enew[i].length = x[i + 1];
      }
      evo_tree new_tree(Ns+1, enew);
      if(maxj==0){
          if(model == 0){
              new_tree.mu = rtree.mu;
          }
          else{
              new_tree.dup_rate = rtree.dup_rate;
              new_tree.del_rate = rtree.del_rate;
              new_tree.chr_gain_rate = rtree.chr_gain_rate;
              new_tree.chr_loss_rate = rtree.chr_loss_rate;
              new_tree.wgd_rate = rtree.wgd_rate;
              new_tree.mu = 0;
          }
      }else{
          if(model == 0){
              new_tree.mu = x[rtree.nedge];
              if(debug){
                  for(int i=0; i<rtree.nedge+1; i++){ cout << x[i] << '\n';}
                  cout << "mu value so far: " << new_tree.mu << endl;
              }
          }
          else{
              new_tree.dup_rate = x[rtree.nedge];
              new_tree.del_rate = x[rtree.nedge+1];
              new_tree.chr_gain_rate = x[rtree.nedge+2];
              new_tree.chr_loss_rate = x[rtree.nedge+3];
              new_tree.wgd_rate = x[rtree.nedge+4];
              new_tree.mu = 0;
              if(debug){
                  for(int i=0; i<rtree.nedge+2; i++){ cout << x[i] << '\n';}
                  cout << "dup_rate value so far: " << new_tree.dup_rate << endl;
                  cout << "del_rate value so far: " << new_tree.del_rate << endl;
                  cout << "chr_gain_rate value so far: " << new_tree.chr_gain_rate << endl;
                  cout << "chr_loss_rate value so far: " << new_tree.chr_loss_rate << endl;
                  cout << "wgd_rate value so far: " << new_tree.wgd_rate << endl;
              }
          }
      }
      rtree = evo_tree(new_tree);
    }else{
      int count = 0;
      for(int i=0; i<rtree.nedge-1; ++i){
        if(enew[i].end > Ns){
          	enew[i].length = x[count + 1];
          	count++;
        }else{
  	        enew[i].length = 0;
        }
      }
      if(debug){
          cout << "total height so far: " << x[count + 1] << endl;
      }
      evo_tree new_tree(rtree.nleaf, enew, x[count + 1], rtree.tobs);
      if(maxj==0){
        if(model == 0){
            new_tree.mu = rtree.mu;
        }
        else{
            new_tree.dup_rate = rtree.dup_rate;
            new_tree.del_rate = rtree.del_rate;
            new_tree.chr_gain_rate = rtree.chr_gain_rate;
            new_tree.chr_loss_rate = rtree.chr_loss_rate;
            new_tree.wgd_rate = rtree.wgd_rate;
            new_tree.mu = 0;
        }
      }else{
        if(model == 0){
            new_tree.mu = x[rtree.nintedge+2];
            if(debug){
                for(int i=0; i<=rtree.nintedge+2; i++){ cout << x[i] << '\n';}
                cout << "mu value so far: " << new_tree.mu << endl;
            }
        }
        else{
            new_tree.dup_rate = x[rtree.nintedge+2];
            new_tree.del_rate = x[rtree.nintedge+3];
            new_tree.chr_gain_rate = x[rtree.nintedge+4];
            new_tree.chr_loss_rate = x[rtree.nintedge+5];
            new_tree.wgd_rate = x[rtree.nintedge+6];
            new_tree.mu = 0;
            if(debug){
                for(int i=0; i<=rtree.nintedge+6; i++){ cout << x[i] << '\n';}
                cout << "dup_rate value so far: " << new_tree.dup_rate << endl;
                cout << "del_rate value so far: " << new_tree.del_rate << endl;
                cout << "chr_gain_rate value so far: " << new_tree.chr_gain_rate << endl;
                cout << "chr_loss_rate value so far: " << new_tree.chr_loss_rate << endl;
                cout << "wgd_rate value so far: " << new_tree.wgd_rate << endl;
            }
         }
      }
      rtree = evo_tree(new_tree);
    }
}


// Find the samples at the first sampling time point
// tobs stores sampling time difference relative to the first sample
int get_first_sample(vector<double> tobs){
    for (int i = 0; i < tobs.size(); i++){
        if(fabs(tobs[i] - 0) < SMALL_VAL){
            return i;
        }
    }
    return -1;
}

vector<double> get_blens_from_intervals(evo_tree &rtree, double *x){
    vector<double> intervals;
    // The estimated value may be nan
    for(int i=0; i < rtree.nleaf - 1; i++){
        double len = x[i+1];
        bool is_nan = std::isnan(len);
        if ( is_nan || (!is_nan && len < SMALL_BLEN)){
           len = SMALL_BLEN;
        }
        intervals.push_back(len);
    }
    // vector<double> intervals(x, x + sizeof x / sizeof x[0]);
    if(debug){
        rtree.print();
        cout << "Current length of intervals: " << endl;
        for(int i = 0; i<intervals.size(); i++){
            cout << i + 1 << "\t" << "\t" << intervals[i] << endl;
        }
        cout << "Corresponding " << rtree.nleaf << " nodes: ";
        for(int i = 0; i < rtree.top_tnodes.size(); i++){
            cout << rtree.top_tnodes[i] + 1 << "\t" << "Prevous node time " << rtree.node_times[rtree.top_tnodes[i]] << endl;
        }
    }

    vector<double> blens = rtree.get_edges_from_interval(intervals, rtree.top_tnodes);

    return blens;
}

// Update the tree after each iteration in the BFGS optimization
// Estimate time intervals rather than branch length in order to avoid negative terminal branch lengths
// Sort node times in increasing order and take the first Ns intervals
void get_variables_ntime(evo_tree &rtree, double *x){
    int debug = 0;
    int npar_ne;
    // create a new tree from current value of parameters
    if(debug){
        cout << "rtree before at address " << &rtree <<endl;
        rtree.print();
    }
    if(opt_one_branch==1){
        npar_ne = 1;
        if(debug) { cout << "update only one branch " << rtree.current_eid+1 << endl;    }
        if(cons == 0){
            rtree.edges[rtree.current_eid].length = x[1];
        }else{
            rtree.update_edge_from_ratio(x[1], rtree.current_eid);
        }
    }else{
        if(cons == 0){
            npar_ne = rtree.nedge - 1;
            // The first element of x is not used for optimization
            // The index of x is added by 1 compared with index for simplex method
            for(int i=0; i<rtree.nedge-1; i++){
              rtree.edges[i].length = x[i + 1];
            }
        }else{
            npar_ne = rtree.nintedge + 1;
            // vector<double> blens = get_blens_from_intervals(rtree, x);
            vector<double> ratios;
            // The estimated value may be nan
            // cout << "estimated x: ";
            for(int i=0; i < npar_ne; i++){
                double val = x[i + 1];
                // cout << "\t" << val;
                if(std::isnan(val)){
                    val = rtree.ratios[i];
                }
                ratios.push_back(val);
                rtree.ratios[i] = val;
            }
            // cout << endl;

            // vector<double> intervals(x, x + sizeof x / sizeof x[0]);
            if(debug){
                rtree.print();
                cout << "Current values of estimated variables: " << endl;
                for(int i = 0; i< ratios.size(); i++){
                    cout << i + 1 << "\t" << "\t" << ratios[i] << endl;
                }
            }
            // rtree.ratios = ratios;
            rtree.update_edges_from_ratios();
        }

        if(debug){
          cout << "New branch lengths: ";
          for(int i = 0; i < rtree.edges.size(); i++){
              cout << i + 1 << "\t" << "\t" << rtree.edges[i].length << endl;
          }
        }
        // update_tree(rtree, enew, x);
    }

    if(maxj == 1){
        if(model == 0){
            // rtree.nedge
            // rtree.nintedge+2 for constrained branches
            rtree.mu = x[npar_ne+1];
            if(debug){
                for(int i=0; i<=npar_ne+1; i++){ cout << x[i] << '\n';}
                cout << "mu value so far: " << rtree.mu << endl;
            }
        }
        else{
            rtree.dup_rate = x[npar_ne+1];
            rtree.del_rate = x[npar_ne+2];
            if(!only_seg){
                rtree.chr_gain_rate = x[npar_ne+3];
                rtree.chr_loss_rate = x[npar_ne+4];
                rtree.wgd_rate = x[npar_ne+5];
            }

            if(debug){
                for(int i=0; i<=npar_ne+1; i++){ cout << x[i] << '\n';}
                cout << "dup_rate value so far: " << rtree.dup_rate << endl;
                cout << "del_rate value so far: " << rtree.del_rate << endl;
                if(!only_seg){
                    cout << "chr_gain_rate value so far: " << rtree.chr_gain_rate << endl;
                    cout << "chr_loss_rate value so far: " << rtree.chr_loss_rate << endl;
                    cout << "wgd_rate value so far: " << rtree.wgd_rate << endl;
                }
            }
        }
    }

    if(debug){
        cout << "rtree after at address " << &rtree <<endl;
        rtree.print();
    }

}


/**
    the target function which needs to be optimized
    @param x the input vector x
    @return the function value at x
*/
double targetFunk(evo_tree &rtree, double x[]) {
// double targetFunk(evo_tree &rtree, int model, int cons, int maxj, int cn_max, int only_seg, int correct_bias, int is_total, double x[]) {
    // negative log likelihood
    // get_variables(rtree, model, cons, maxj, x);
    // get_variables_ntime(curr_tree, model, cons, maxj, only_seg, x);
    get_variables_ntime(rtree, x);
    if(model == 3){
        // cout << "Getting target function for optimization" << endl;
        // cout << max_wgd << "\t" << max_chr_change << "\t" << max_site_change << "\t" << decomp_table.size() << endl;
        // return -1.0*get_likelihood_decomp(Ns, Nchar, num_invar_bins, vobs, curr_tree, decomp_table, cons, cn_max, m_max, max_wgd, max_chr_change, max_site_change, correct_bias, is_total);
        return -1.0 * get_likelihood_decomp(rtree);
    }else{
        // return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, curr_tree, model, cons, cn_max, only_seg, correct_bias, is_total);
        return -1.0 * get_likelihood_revised(rtree);
    }
}


// Compute the likelihood of the tree with one parameter (one branch or one mutation rate)
double computeFunction(evo_tree& rtree, double value, int type = -1) {
    // Update the branch
    switch (type)
    {
        case -1:{
            rtree.current_it->length = value;
            rtree.current_it_back->length = value;
            // Need to update rtree.edges
            rtree.edges[rtree.current_eid].length = value;
            break;
        }
        case 0:{
            rtree.mu = value;
            break;
        }
        case 1:{
            rtree.dup_rate = value;
            break;
        }
        case 2:{
            rtree.del_rate = value;
            break;
        }
        case 3:{
            rtree.chr_gain_rate = value;
            break;
        }
        case 4:{
            rtree.chr_loss_rate = value;
            break;
        }
        case 5:{
            rtree.wgd_rate = value;
            break;
        }
        default:
            break;
    }

    double nlnl = 0;
    if(model == 3){
        nlnl = -get_likelihood_decomp(rtree);
    }else{
        nlnl = -get_likelihood_revised(rtree);
    }
    // cout << "\n result of computeFunction " << nlnl << endl;

    return nlnl;
}


// double optimFunc(evo_tree &rtree, int model, int cons, int maxj, int cn_max, int only_seg, int correct_bias, int is_total, int nvar, double *vars) {
//     return targetFunk(rtree, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, vars-1);
// }
//
// double optimGradient(evo_tree &rtree, int model, int cons, int maxj, int cn_max, int only_seg, int correct_bias, int is_total, int nvar, double *x, double *dfx) {
//     return derivativeFunk(rtree, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, nvar, x-1, dfx-1);
// }


// The number of parameters to estimate, different when the mutation rates are estimated
int get_ndim(int model, int maxj, int npar_ne, int only_seg){
    int ndim = 0;

    if(maxj==0){
      ndim = npar_ne;
    }else{
        if(model==0){
            ndim = npar_ne + 1;
        }
        else{
            if(only_seg){
                ndim = npar_ne + 2;
            }
            else{
                ndim = npar_ne + 5;
            }
        }
    }
    return ndim;
}



/*****************************************************
    One dimensional optimization with Brent method
*****************************************************/
/**
    This function calculate f(value) of the f() function, used by other general optimization method to minimize it.
    Please always override this function to adapt to likelihood or parsimony score.
    The default is for function f(x)=x.
    @param value x-value of the function
    @return f(value) of function f you want to minimize
*/
// double computeFunction(double value) { return value; }

/*****************************************************
	One dimensional optimization with Brent method
*****************************************************/

#define ITMAX 100
#define CGOLD 0.3819660
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Brents method in one dimension */
double brent_opt (evo_tree &rtree, int type, double ax, double bx, double cx, double tol,
                          double *foptx, double *f2optx, double fax, double fbx, double fcx) {
	int iter;
	double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double xw,wv,vx;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=bx;
	fx=fbx;
	if (fax < fcx) {
		w=ax;
		fw=fax;
		v=cx;
		fv=fcx;
	} else {
		w=cx;
		fw=fcx;
		v=ax;
		fv=fax;
	}

	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*foptx = fx;
			xw = x-w;
			wv = w-v;
			vx = v-x;
			*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
			          (v*v*xw + x*x*wv + w*w*vx);
			return x;
		}

		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0)
				p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}

		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=computeFunction(rtree, u, type);
		if (fu <= fx) {
			if (u >= x)
				a=x;
			else
				b=x;

			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x)
				a=u;
			else
				b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else
				if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
		}
	}

	*foptx = fx;
	xw = x-w;
	wv = w-v;
	vx = v-x;
	*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/(v*v*xw + x*x*wv + w*w*vx);

	return x;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN
#undef GOLD
#undef GLIMIT
#undef TINY

/**
    the brent method to find the value that minimizes the computeFunction().
    @return the x-value that minimize the function
    @param xmin lower bound
    @param xmax upper bound
    @param xguess first guess
    @param tolerance tolerance
    @param fx (OUT) function value at the minimum x found
    @param ferror (OUT) Dont know
*/
double minimizeOneDimen(evo_tree &rtree, int type, double xmin, double xguess, double xmax, double tolerance, double *fx, double *f2x) {
	double eps, optx, ax, bx, cx, fa, fb, fc;
	//int    converged;	/* not converged error flag */

	/* first attempt to bracketize minimum */
	if (xguess < xmin) xguess = xmin;
	if (xguess > xmax) xguess = xmax;
	eps = xguess*tolerance*50.0;
	ax = xguess - eps;
	if (ax < xmin) ax = xmin;
	bx = xguess;
	cx = xguess + eps;
	if (cx > xmax) cx = xmax;

	/* check if this works */
    // compute fb first to save some computation, if any
	fb = computeFunction(rtree, bx, type);
	fa = computeFunction(rtree, ax, type);
	fc = computeFunction(rtree, cx, type);

	/* if it works use these borders else be conservative */
	if ((fa < fb) || (fc < fb)) {
		if (ax != xmin) fa = computeFunction(rtree, xmin, type);
		if (cx != xmax) fc = computeFunction(rtree, xmax, type);
		ax = xmin;
		cx = xmax;
	}
	/*
	const int MAX_ROUND = 10;
	for (i = 0; ((fa < fb-tolerance) || (fc < fb-tolerance)) && (i<MAX_ROUND); i++) {
		// brent method require that fb is smaller than both fa and fc
		// find some random values until fb achieve this
			bx = (((double)rand()) / RAND_MAX)*(cx-ax) + ax;
			fb = computeFunction(bx);
	}*/

/*
	if ((fa < fb) || (fc < fb)) {
		if (fa < fc) { bx = ax; fb = fa; } else { bx = cx; fb = fc; }
		//cout << "WARNING: Initial value for Brent method is set at bound " << bx << endl;
	}*/
	//	optx = brent_opt(xmin, xguess, xmax, tolerance, fx, f2x, fa, fb, fc);
	//} else
	optx = brent_opt(rtree, type, ax, bx, cx, tolerance, fx, f2x, fa, fb, fc);
    if (*fx > fb) // if worse, return initial value
    {
        *fx = computeFunction(rtree, bx, type);
        return bx;
    }

	return optx; /* return optimal x */
}


/*****************************************************
    L-BFGS-B method
*****************************************************/

/**
	the approximated derivative function
	@param x the input vector x
	@param dfx the derivative at x
	@return the function value at x
*/
// double derivativeFunk(evo_tree &rtree, int model, int cons, int maxj, int cn_max, int only_seg, int correct_bias, int is_total, int ndim, double x[], double dfx[]) {
double derivativeFunk(evo_tree &rtree, int ndim, double x[], double dfx[]) {
    int debug = 0;
	double *h = new double[ndim+1];
    double temp;
    int dim;
	// double fx = targetFunk(rtree, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, x);
    double fx = targetFunk(rtree, x);
	for (dim = 1; dim <= ndim; dim++ ){
		temp = x[dim];
		h[dim] = ERROR_X * fabs(temp);
		if (h[dim] == 0.0) h[dim] = ERROR_X;
		x[dim] = temp + h[dim];
		h[dim] = x[dim] - temp;
		// dfx[dim] = (targetFunk(rtree, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, x));
        dfx[dim] = (targetFunk(rtree, x));
		x[dim] = temp;
	}
	for (dim = 1; dim <= ndim; dim++ ){
        dfx[dim] = (dfx[dim] - fx) / h[dim];
        if(debug){
            cout << "dfx[dim] " << dfx[dim] << endl;
        }
    }
    delete [] h;
	return fx;
}


double optimFunc(evo_tree &rtree, int nvar, double *vars) {
    return targetFunk(rtree, vars-1);
}

double optimGradient(evo_tree &rtree, int nvar, double *x, double *dfx) {
    return derivativeFunk(rtree, nvar, x-1, dfx-1);
}

// void lbfgsb(evo_tree &rtree, int model, int cons, int maxj, int cn_max, int only_seg, int correct_bias, int is_total, int n, int m, double *x, double *l, double *u, int *nbd,
void lbfgsb(evo_tree &rtree, int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, int *fail,
		double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT)
{
	char task[60];
	double f, *g, dsave[29], *wa;
	int tr = -1, iter = 0, *iwa, isave[44], lsave[4];

	/* shut up gcc -Wall in 4.6.x */

	for(int i = 0; i < 4; i++) lsave[i] = 0;

	if(n == 0) { /* not handled in setulb */
		*fncount = 1;
		*grcount = 0;
		// *Fmin = optimFunc(rtree, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, n, u);
        *Fmin = optimFunc(rtree, n, u);
		strcpy(msg, "NOTHING TO DO");
		*fail = 0;
		return;
	}
	if (nREPORT <= 0) {
		cerr << "REPORT must be > 0 (method = \"L-BFGS-B\")" << endl;
		exit(1);
	}
	switch(trace) {
	case 2: tr = 0; break;
	case 3: tr = nREPORT; break;
	case 4: tr = 99; break;
	case 5: tr = 100; break;
	case 6: tr = 101; break;
	default: tr = -1; break;
	}

	*fail = 0;
	g = (double*) malloc (n * sizeof(double));
	/* this needs to be zeroed for snd in mainlb to be zeroed */
	wa = (double *) malloc((2*m*n+4*n+11*m*m+8*m) * sizeof(double));
	iwa = (int *) malloc(3*n * sizeof(int));
	strcpy(task, "START");
	while(1) {
		/* Main workhorse setulb() from ../appl/lbfgsb.c : */
		setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task, tr, lsave, isave, dsave);
		/*    Rprintf("in lbfgsb - %s\n", task);*/
		if (strncmp(task, "FG", 2) == 0) {
			// f = optimGradient(rtree, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, n, x, g);
            f = optimGradient(rtree, n, x, g);
			if (!isfinite(f)) {
				cerr << "L-BFGS-B needs finite values of 'fn'" << endl;
				exit(1);
			}

		} else if (strncmp(task, "NEW_X", 5) == 0) {
			iter++;
			if(trace == 1 && (iter % nREPORT == 0)) {
				cout << "iter " << iter << " value " << f << endl;
			}
			if (iter > maxit) {
				*fail = 1;
				break;
			}
		} else if (strncmp(task, "WARN", 4) == 0) {
			*fail = 51;
			break;
		} else if (strncmp(task, "CONV", 4) == 0) {
			break;
		} else if (strncmp(task, "ERROR", 5) == 0) {
			*fail = 52;
			break;
		} else { /* some other condition that is not supposed to happen */
			*fail = 52;
			break;
		}
	}
	*Fmin = f;
	*fncount = *grcount = isave[33];
	if (trace) {
		cout << "final value " << *Fmin << endl;
		if (iter < maxit && *fail == 0)
			cout << "converged" << endl;
		else
			cout << "stopped after " << iter << " iterations\n";
	}
	strcpy(msg, task);
	free(g);
	free(wa);
	free(iwa);
}


/**
 Function to access the L-BFGS-B function, taken from IQ-TREE package which is further taken from HAL_HAS software package
 1. int nvar : The number of the variables
 2. double* vars : initial values of the variables
 3. double* lower : lower bounds of the variables
 4. double* upper : upper bounds of the variables
 5. double pgtol: gradient tolerance
 5. int maxit : max # of iterations
 @return minimized function value
 After the function is invoked, the values of x will be updated
*/
// double L_BFGS_B(evo_tree &rtree, int model, int cons, int maxj, int cn_max, int only_seg, int correct_bias, int is_total, int n, double* x, double* l, double* u, double pgtol, int maxit)
double L_BFGS_B(evo_tree &rtree, int n, double* x, double* l, double* u, double pgtol, int maxit) {
    int debug = 0;
    int i;
	double Fmin;
	int fail;
	int fncount;
	int grcount;
	char msg[100];

	int m = 10;          // number of BFGS updates retained in the "L-BFGS-B" method. It defaults to 5.

	int *nbd;           // 0: unbounded; 1: lower bounded; 2: both lower & upper; 3: upper bounded
	nbd = new int[n];
	for (i=0; i<n; i++)
		nbd[i] = 2;

	double factr = 1e+7; // control the convergence of the "L-BFGS-B" method.
	// Convergence occurs when the reduction in the object is within this factor
	// of the machine tolerance.
	// Default is 1e7, that is a tolerance of about 1e-8

//	double pgtol = 0;   // helps control the convergence of the "L-BFGS-B" method.
//    pgtol = 0.0;
	// It is a tolerance on the projected gradient in the current search direction.
	// Default is zero, when the check is suppressed

	int trace = 0;      // non-negative integer.
    if (debug)
        trace = 1;
	// If positive, tracing information on the progress of the optimization is produced.
	// Higher values may produce more tracing information.

	int nREPORT = 10;   // The frequency of reports for the "L-BFGS-B" methods if "trace" is positive.
	// Defaults to every 10 iterations.

/*#ifdef USE_OLD_PARAM
	lbfgsb(n, m, x, l, u, nbd, &Fmin, fn, gr1, &fail, ex,
			factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
#else*/

    // cout << "initial values in bfgs: ";
    // for(int mi = 0; mi < 3; mi++){
    //     cout << "\t" << x[mi];
    // }
    // cout << "\n";
    //
    // cout << "upper bound in bfgs: ";
    // for(int mi = 0; mi < 3; mi++){
    //     cout << "\t" << u[mi];
    // }
    // cout << "\n";

	// lbfgsb(rtree, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, n, m, x, l, u, nbd, &Fmin, &fail, factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
    // lbfgsb(n, m, x, l, u, nbd, &Fmin, &fail, factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
    lbfgsb(rtree, n, m, x, l, u, nbd, &Fmin, &fail, factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
//#endif

    if (fail == 51 || fail == 52) {
        cout << msg << endl;
    }

	delete[] nbd;

    return Fmin;
}


void init_branch_by_interval(double *variables, double *upper_bound, double *lower_bound, evo_tree &rtree, int npar_ne){
    int debug = 0;
    // estimate the top Ns intervals, which will be dynamic in different trees
    if(debug){
        cout << "Initializing the time intervals" << endl;
    }
    int sample1 = 0;
    rtree.get_ntime_interval(rtree.nleaf - 1, sample1);
    // If max_len is very small, this tree seems not feasible with current variables
    if( rtree.tobs[sample1] > 0 ){
        cout << "The first tip node in the proposed tree " << sample1 + 1 << " is not a first sample! " << endl;
    }

    for(int i=0; i< rtree.top_tinvls.size(); ++i){
      variables[i+1] = rtree.top_tinvls[i];
      lower_bound[i+1] = SMALL_BLEN;
      // cout << "Starting node of the interval " << rtree.top_tnodes[i] << "\t" << rtree.node_times[rtree.top_tnodes[i]] << endl;
      upper_bound[i+1] = age / 2;
      // When the interval is below the first sample, it should be smaller than the time difference between tip below it and 1st sample
      if(rtree.node_times[rtree.top_tnodes[i]] > rtree.node_times[sample1]){
          // Find the nodes below this node and gets its time
          // cout << "Tips below: " << endl;
          vector<int> tips = rtree.get_tips_below(rtree.top_tnodes[i]);
          vector<double> tips_obs;
          for(int i=0; i<tips.size(); i++){
              tips_obs.push_back(rtree.tobs[tips[i]]);
              // cout << "\t" << tips[i] << "\t" << tips_obs[i] << endl;
          }
          double max_len = *min_element(tips_obs.begin(), tips_obs.end());
          // If max_len is very small, this tree seems not feasible
          if( fabs(max_len - 0) < SMALL_VAL){
              cout << "The terminal branch length under " << rtree.top_tnodes[i] + 1 << " in the proposed tree will not be valid! " << endl;
              upper_bound[i+1] = SMALL_BLEN + SMALL_BLEN;
              // return rtree;
          }else{
              upper_bound[i+1] = max_len;
          }
      }
      if(debug){
          cout << "Top " << Ns << " time intervals: ";
          for(int i = 0; i < rtree.top_tinvls.size(); i++){
              cout << "\t" << rtree.top_tinvls[i];
          }
          cout<< endl;
          cout << "Corresponding " << Ns + 1 << " nodes: ";
          for(int i = 0; i < rtree.top_tnodes.size(); i++){
              cout << rtree.top_tnodes[i] + 1 << "\t" << "\t" << rtree.node_times[rtree.top_tnodes[i]] << endl;
          }
      }
}
}



// Using BFGS method to get the maximum likelihood with lower and upper bounds
// Note: the topology of rtree is fixed
void max_likelihood_BFGS(evo_tree &rtree, double &minL, double tolerance, int miter){
    int debug = 0;
    int npar_ne = 0;    // number of parameters to estimate
    int ndim = 0;
    double *variables, *upper_bound, *lower_bound;

    // Set variables
    if(opt_one_branch == 1){
        npar_ne = 1;
    }else{
        if(cons == 0){  // only estimate internal branches
            npar_ne = rtree.nedge-1;
        }else{
            npar_ne = rtree.nintedge + 1;
        }
    }
    if(debug) cout << "There are " << npar_ne << " parameters to estimate " << endl;

    ndim = get_ndim(model, maxj, npar_ne, only_seg);
    variables = new double[ndim+1];
    upper_bound = new double[ndim+1];
    lower_bound = new double[ndim+1];

    if(cons == 1){
        if(debug){
            cout << "Initializing variables related to branch length" << endl;
        }
        // init_branch_by_interval(variables, upper_bound, lower_bound, rtree, npar_ne);
        rtree.get_ratio_from_age(rtree.current_eid);
        if(debug){
            // cout << ratios.size() << "\t" << npar_ne << endl;
            cout << "time ratios obtained from node times: ";
            for(int i = 0; i < rtree.nnode; i++){
                cout << "\t" << rtree.ratios[i];
            }
            cout << endl;
        }
        // assert(rtree.nleaf == npar_ne);
        if(opt_one_branch == 1){
            int i = 0;
            variables[i+1] = rtree.ratios[i];
            lower_bound[i+1] = MIN_RATIO;
            upper_bound[i+1] = MAX_RATIO;
        }else{
            int i = 0;
            // root
            variables[i+1] = rtree.ratios[i];
            double max_obs = *max_element(rtree.tobs.begin(), rtree.tobs.end());
            lower_bound[i+1] = max_obs;
            upper_bound[i+1] = max_obs + age;

            for(int i=1; i < npar_ne; ++i){
              variables[i+1] = rtree.ratios[i];
              lower_bound[i+1] = MIN_RATIO;
              upper_bound[i+1] = MAX_RATIO;
            }
        }
    }else{
        if(opt_one_branch == 1){
            int i = 0;
            variables[i+1] = rtree.edges[rtree.current_eid].length;
            lower_bound[i+1] = SMALL_BLEN;
            upper_bound[i+1] = age;
        }else{
            for(int i=0; i<npar_ne; ++i){
              variables[i+1] = rtree.edges[i].length;
              lower_bound[i+1] = SMALL_BLEN;
              upper_bound[i+1] = age;
            }
        }
    }

    if(maxj==1){    // estimate mutation rates
        if(model == 0){
            int i = npar_ne;
            variables[i+1] = rtree.mu;
            lower_bound[i+1] = 0;
            upper_bound[i+1] = MAX_MRATE;
        }
        else{
            int i = npar_ne;
            variables[i+1] = rtree.dup_rate;
            lower_bound[i+1] = 0;
            upper_bound[i+1] = MAX_MRATE;

            i = npar_ne+1;
            variables[i+1] = rtree.del_rate;
            lower_bound[i+1] = 0;
            upper_bound[i+1] = MAX_MRATE;

            if(!only_seg){
                i = npar_ne+2;
                variables[i+1] = rtree.chr_gain_rate;
                lower_bound[i+1] = 0;
                upper_bound[i+1] = MAX_MRATE;

                i = npar_ne+3;
                variables[i+1] = rtree.chr_loss_rate;
                lower_bound[i+1] = 0;
                upper_bound[i+1] = MAX_MRATE;

                i = npar_ne+4;
                variables[i+1] = rtree.wgd_rate;
                lower_bound[i+1] = 0;
                upper_bound[i+1] = MAX_MRATE;
            }
        }
    }

    if(debug){
        cout << "Variables to estimate: " << endl;
        for(int i=0; i < ndim+1; i++){
            cout << "\t" << variables[i] << "\tLB:" << lower_bound[i] << "\tUB:" << upper_bound[i] << endl;
        }
    }

    // variables contains the parameters to estimate (branch length, mutation rate)
    minL = L_BFGS_B(rtree, ndim, variables+1, lower_bound+1, upper_bound+1, tolerance, miter);

    if (debug){
        cout << "lnL of current ML tree: " << minL << endl;
        rtree.print();
    }

    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables;

    // rtree has been updated in the optimization process
    // return rtree;
}


// given a tree, maximise the branch lengths (and optionally mu) assuming branch lengths are independent or constrained in time
evo_tree max_likelihood(evo_tree &rtree, int model, double& minL, const double ssize, const double tolerance, int miter, int cons=0, int maxj=0, int cn_max=4, int only_seg=0, int correct_bias=1, int is_total=1){
  int debug = 0;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  /* save original handler, install new handler */
  gsl_error_handler_t  *old_handler = gsl_set_error_handler (&my_err_handler);

  int npar_ne;
  int npar;
  gsl_vector *x;

  if(cons == 0){
    npar_ne = rtree.nedge-1;
    if(maxj==0){
      npar = npar_ne;
    }else{
      if(model==0){
          npar = npar_ne + 1;
      }
      if(model==1){
          npar = npar_ne + 5;
      }
    }

    // initialise the best guess branch length and mu if required
    x = gsl_vector_alloc (npar);
    for(int i=0; i<npar_ne; ++i){
      gsl_vector_set (x, i, log(rtree.edges[i].length));
    }
    if(maxj==1){
        if(model==0){
            gsl_vector_set (x, npar_ne, log(rtree.mu));
        }
        if(model==1){
            gsl_vector_set (x, npar_ne, log(rtree.dup_rate) );
            gsl_vector_set (x, npar_ne+1, log(rtree.del_rate) );
            gsl_vector_set (x, npar_ne+2, log(rtree.chr_gain_rate) );
            gsl_vector_set (x, npar_ne+3, log(rtree.chr_loss_rate) );
            gsl_vector_set (x, npar_ne+4, log(rtree.wgd_rate) );
        }
    }

    if(maxj==0){
      minex_func.f = my_f;
    }else{
      minex_func.f = my_f_mu;
    }

  }else{
    npar_ne = rtree.nintedge + 1;
    if(maxj==0){
      npar = npar_ne;
    }else{
        if(model==0){
            npar = npar_ne + 1;
        }
        if(model==1){
            npar = npar_ne + 5;
        }
    }

    x = gsl_vector_alloc (npar);
    // initialise with internal edges
    for(int i=0; i<rtree.nintedge; ++i){
      gsl_vector_set (x, i, log(rtree.intedges[i]->length));
    }

    // initialise with current total tree time
    gsl_vector_set (x, rtree.nintedge, log(rtree.get_total_time()) );
    if(maxj==1){
        if(model==0){
            gsl_vector_set (x, npar_ne, log(rtree.mu) );
        }
        if(model==1){
            gsl_vector_set (x, npar_ne, log(rtree.dup_rate) );
            gsl_vector_set (x, npar_ne+1, log(rtree.del_rate) );
            gsl_vector_set (x, npar_ne+2, log(rtree.chr_gain_rate) );
            gsl_vector_set (x, npar_ne+3, log(rtree.chr_loss_rate) );
            gsl_vector_set (x, npar_ne+4, log(rtree.wgd_rate) );
        }
    }

    if(maxj==0){
      minex_func.f = my_f_cons;
    }else{
      minex_func.f = my_f_cons_mu;
    }
  }
  evo_tree *p = &rtree;
  void * pv = p;

  // Set initial step sizes to 1
  //cout << "numbers: nedge / nintedge / npar / x->size: " << rtree.nedge << "\t" << rtree.nintedge << "\t" << npar << "\t" << x->size << endl;
  gsl_vector* ss = gsl_vector_alloc (npar);
  gsl_vector_set_all (ss, ssize);

  // Initialize method and iterate
  minex_func.n = npar;
  minex_func.params = pv;

  s = gsl_multimin_fminimizer_alloc (T, npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  size_t iter = 0;
  int status;
  double size;
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status) break;

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, tolerance);
    // gsl_set_error_handler_off();

    if(debug){
      printf ("%5lu, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n",
	      iter,
	      exp(gsl_vector_get (s->x, 0)),
	      exp(gsl_vector_get (s->x, s->x->size-1)),
	      s->fval, size);
    }

    if (status == GSL_SUCCESS){
      if(debug){
        	printf ("converged to minimum at\n");

        	printf ("%5lu, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n",
        		iter,
        		exp(gsl_vector_get (s->x, 0)),
        		exp(gsl_vector_get (s->x, s->x->size-1)),
        		s->fval, size);
        }
    }
  }
  while (status == GSL_CONTINUE && iter < miter);

  if( status == GSL_CONTINUE ){
    cout << "### WARNING: maximum likelihood did not converge" << endl;
  }

  // create a new tree
  vector<edge> enew;
  for(int i=0; i<rtree.nedge; ++i){
    enew.push_back( rtree.edges[i] );
  }

  if(cons == 0){
    for(int i=0; i<npar_ne; ++i){
      enew[i].length = exp(gsl_vector_get(s->x, i));
    }
    evo_tree new_tree(Ns+1, enew);
    if(maxj==0){
        if(model==0){
            new_tree.mu = rtree.mu;
        }
        if(model==1){
            new_tree.dup_rate = rtree.dup_rate;
            new_tree.del_rate = rtree.del_rate;
            new_tree.chr_gain_rate = rtree.chr_gain_rate;
            new_tree.chr_loss_rate = rtree.chr_loss_rate;
            new_tree.wgd_rate = rtree.wgd_rate;
            new_tree.mu = 0;
        }
    }else{
        if(model==0){
            new_tree.mu = exp(gsl_vector_get(s->x, npar_ne));
        }
        if(model==1){
            new_tree.dup_rate = exp(gsl_vector_get(s->x, npar_ne));
            new_tree.del_rate = exp(gsl_vector_get(s->x, npar_ne+1));
            new_tree.chr_gain_rate = exp(gsl_vector_get(s->x, npar_ne+2));
            new_tree.chr_loss_rate = exp(gsl_vector_get(s->x, npar_ne+3));
            new_tree.wgd_rate = exp(gsl_vector_get(s->x, npar_ne+4));
            new_tree.mu = 0;
        }
    }

    minL = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    /* restore original handler */
    gsl_set_error_handler (old_handler);

    return new_tree;

  }else{
    int count = 0;
    for(int i=0; i<rtree.nedge-1; ++i){
      if(enew[i].end > Ns){
        	enew[i].length = exp( gsl_vector_get(s->x, count) );
        	count++;
      }else{
	        enew[i].length = 0;
      }
    }
    // evo_tree new_tree(rtree.nleaf, enew, age, rtree.tobs);
    evo_tree new_tree(rtree.nleaf, enew, exp( gsl_vector_get(s->x, count)), rtree.tobs);
    if(maxj==0){
        if(model==0){
            new_tree.mu = rtree.mu;
        }
        if(model==1){
            new_tree.dup_rate = rtree.dup_rate;
            new_tree.del_rate = rtree.del_rate;
            new_tree.chr_gain_rate = rtree.chr_gain_rate;
            new_tree.chr_loss_rate = rtree.chr_loss_rate;
            new_tree.wgd_rate = rtree.wgd_rate;
            new_tree.mu = 0;
        }
    }else{
        if(model==0){
            new_tree.mu = exp(gsl_vector_get(s->x, npar_ne));
        }
        if(model==1){
            new_tree.dup_rate = exp(gsl_vector_get(s->x, npar_ne));
            new_tree.del_rate = exp(gsl_vector_get(s->x, npar_ne+1));
            new_tree.chr_gain_rate = exp(gsl_vector_get(s->x, npar_ne+2));
            new_tree.chr_loss_rate = exp(gsl_vector_get(s->x, npar_ne+3));
            new_tree.wgd_rate = exp(gsl_vector_get(s->x, npar_ne+4));
            new_tree.mu = 0;
        }
    }

    minL = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    /* restore original handler */
    gsl_set_error_handler (old_handler);

    // if (!is_tree_valid(new_tree)){
    //     // keep the original tree
    //     return rtree;
    // }

    return new_tree;
  }

}
