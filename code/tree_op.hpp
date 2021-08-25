#ifndef TREE_OP_HPP
#define TREE_OP_HPP

// This file contains functions related to evo_tree


#include "common.hpp"
#include "stats.hpp"
#include "evo_tree.hpp"


// using namespace std;


// Scaling tree height to 1/HEIGHT_SCALE if the current height is larger than the upper bound (patient age at last sample)
const int HEIGHT_SCALE = 3;
// // The difference from minmial height
// const int HEIGHT_OFFSET = 10;

// The time from beginning to the time of first sample
inline double get_total_time(const vector<double>& node_times, const double& max_tobs){
  return *max_element(node_times.begin(), node_times.end()) - max_tobs;
}


// The time from beginning to the time of last sample
inline double get_tree_height(const vector<double>& node_times){
  return *max_element(node_times.begin(), node_times.end());
}


// Check if there is negative branch length
inline bool is_blen_valid(const evo_tree& rtree){
  int nedge = 2 * rtree.nleaf - 2;

  for(int i = 0; i < nedge; ++i){
      if(rtree.edges[i].length < 0){
          return false;
      }
  }

  return true;
}


// update edge ending at i when its time is updated
inline void update_edge_len(evo_tree& rtree, int i){
  int nedge = 2 * rtree.nleaf - 2;

  for(int j = 0; j < nedge; j++){
      if (rtree.edges[j].end == i){
          rtree.edges[j].length = rtree.nodes[i].time - rtree.nodes[rtree.edges[j].start].time;
      }
      // assert(rtree.edges[j].length > 0);
  }
}


// Find the samples at the first sampling time point
// tobs stores sampling time difference relative to the first sample
inline int get_first_sample(const vector<double>& tobs){
    for (int i = 0; i < tobs.size(); i++){
        if(fabs(tobs[i] - 0) < SMALL_VAL){
            return i;
        }
    }
    return -1;
}


// check if the age of tip node is maintained
bool is_tip_age_valid(const vector<double>& node_ages, const vector<double>& tobs);

// check if node_ages and node_times are consistent
bool is_age_time_consistent(const vector<double>& node_times, const vector<double>& node_ages);

void check_node_age_ratio(evo_tree& tree, const vector<int>& knodes);

void check_nodes_transverse(evo_tree& tree);

// The string will be different for different labeled histories (topologies may be the same but timings are different)
string create_tree_string(const evo_tree& tree);


// The string will be unique for different topologies
// The key is to name the internal node based on tips
string create_tree_string_uniq(const evo_tree& tree);

// for different labeled histories
string order_tree_string(const string& tree);


// for different topologies
// Sort the parent nodes at first
string order_tree_string_uniq(const string& tree);


// Randomly swap two leaves
evo_tree perturb_tree(evo_tree& tree, long unsigned (*fp_myrng)(long unsigned));


// test if the tree is correctly built
void test_evo_tree(const evo_tree& tree);


// generate neutral coalescent trees
// here nsample is the number of cancer samples (not including germline node)
// here we directly calculate the edges in the tree
void generate_coal_tree(const int& nsample, gsl_rng* r, long unsigned (*fp_myrng)(long unsigned), vector<int>& edges, vector<double>& lengths, vector<double>& epoch_times, vector<double>& times, const int& Ne, const double& beta = 0, const double& gtime = 0.002739726);


// Scale the total time by given time
evo_tree generate_coal_tree(const int& nsample, gsl_rng* r, long unsigned (*fp_myrng)(long unsigned), int Ne = 1, double beta = 0, double gtime = 0.002739726);


void assign_tip_times(double delta_t, int Ns, gsl_rng* r, vector<double>& tobs, const vector<int>& edges, vector<double>& lengths);


// Simulate a random tree for mutataion generations
evo_tree generate_random_tree(int Ns, gsl_rng* r, long unsigned (*fp_myrng)(long unsigned), int Ne, int age, double beta, double gtime, double delta_t, int cons, int debug = 0);


// Create a new tree with the same topology as input tree but different branch lengths
evo_tree create_new_tree(gsl_vector* blens, evo_tree& rtree, const double& max_tobs, int cons);



// Adjust the time of sample1 so that it is not after all other tips. Then the subsequent adjustment will be addition, not introducing negative values
void adjust_sample1(evo_tree& rtree, const vector<double>& tobs, int sample1);


// Change the time and branch length for one tip
void adjust_one_tip(evo_tree& rtree, const vector<double>& tobs, int i, int sample1);


// Adjust all tips recursively so that tree height is smaller than age
void adjust_all_tips(evo_tree& rtree, const double& max_tobs, int age);


void adjust_blen(double& nx, double a, double b);


// Seem not work well
bool is_tip_valid(const evo_tree& rtree, const vector<double>& tobs, int sample1);

// Adjust the terminal branch lengths of a tree with fixed topolgy when the constraints are violated
// The initial or perturbed tree may not have tip nodes in accord with the sampling time information
// Ensure edge length keep positive
void adjust_tree_tips(evo_tree& rtree, const vector<double>& tobs, int age);


// The initial branch lengths are always positive.
// But after optimization, some branch lengths may become negative.
// Increase all branch lengths with the same value to keep tip timing difference
// The result tree may break the height constraint
void adjust_tree_blens_all(evo_tree& rtree);


// Only change negative branches. The tip differences will be adjusted later
void adjust_tree_blens(evo_tree& rtree);


// Scale the tree so that total tree height is in the lifetime of the patient
// This may cause tips nodes violating the given sampling time
void adjust_tree_height(evo_tree& rtree, gsl_rng* r, double min_height, double max_height, int scale = HEIGHT_SCALE);


// Check whether the set of branch lengths is valid under the time constraints
bool is_tree_valid(evo_tree& rtree, const double& max_tobs, int age, int cons);


// Save all branch lengths in a vector for optimization
// Assume lenvec is the size of all branches, using ID of end node to distinguish each edge
void save_branch_lengths(evo_tree& rtree, DoubleVector &lenvec, int startid = 0, Node* node = NULL, Node* dad = NULL);

// Restoring branch lengths from a vector
void restore_branch_lengths(evo_tree& rtree, DoubleVector &lenvec, int startid = 0, Node* node = NULL, Node* dad = NULL);


void save_mutation_rates(const evo_tree& rtree, DoubleVector &muvec);

void restore_mutation_rates(evo_tree& rtree, const DoubleVector &muvec);


// TODO: Build parsimony tree from copy number changes (breakpoints)
evo_tree build_parsimony_tree(int Ns, vector<vector<int>>& data);

// filename: a file with at least 3 columns: start, end, length
evo_tree read_tree_info(const string& filename, const int& Ns, int debug = 0);


// Read a newick tree
// evo_tree read_newick(const string& filename){
//     ifstream infile (filename.c_str());
//     if (infile.is_open()){
//       std::string line;
//       istringstream newickstream(incomingNewick);
//     }
//     else{
//         std::cerr << "Error: open of tree data unsuccessful: " <<  filename << std::endl;
//         exit(EXIT_FAILURE);
//     }
//
//     evo_tree new_tree(Ns+1, edges);
//     //new_tree.print();
//
//     return new_tree;
// }


// ajust time of tip nodes to be consistent with input time
// Given a random colescent tree, all the original tip times are equal
// Given a parsimony tree, the tip times may be different, which need to be made equal at first
void adjust_tip_time(evo_tree& rtree, const vector<double>& tobs, int Ns, int same_tip = 1, int debug = 0);


// Read parsimony trees built by other tools as starting trees, assign timings to tip nodes and initialize mutation rates
evo_tree read_parsimony_tree(const string& tree_file, const int& Ns, const vector<double>& rates, const vector<double>& tobs, gsl_rng* r, int age, int cons);


// vector<double> get_blens_from_intervals(evo_tree& rtree, double *x){
//     vector<double> intervals;
//     // The estimated value may be nan
//     for(int i = 0; i < rtree.nleaf - 1; i++){
//         double len = x[i + 1];
//         bool is_nan = std::isnan(len);
//         if ( is_nan || (!is_nan && len < BLEN_MIN)){
//            len = BLEN_MIN;
//         }
//         intervals.push_back(len);
//     }
//     // vector<double> intervals(x, x + sizeof x / sizeof x[0]);
//     if(debug){
//         rtree.print();
//         cout << "Current length of intervals: " << endl;
//         for(int i = 0; i < intervals.size(); i++){
//             cout << i + 1 << "\t" << "\t" << intervals[i] << endl;
//         }
//         cout << "Corresponding " << rtree.nleaf << " nodes: ";
//         for(int i = 0; i < rtree.top_tnodes.size(); i++){
//             cout << rtree.top_tnodes[i] + 1 << "\t" << "Prevous node time " << rtree.nodes[rtree.top_tnodes[i]].time << endl;
//         }
//     }
//
//     vector<double> blens = rtree.get_edges_from_interval(intervals, rtree.top_tnodes);
//
//     return blens;
// }


// Print branch lengths from time intervals
// void print_invl_blen(evo_tree& new_tree, string header){
//     cout << header << endl;
//     new_tree.print();
//
//     cout << "\tTop " << new_tree.nleaf - 1 << " time intervals: ";
//     for(int i = 0; i < new_tree.top_tinvls.size(); i++){
//         cout << "\t" << new_tree.top_tinvls[i];
//     }
//     cout<< endl;
//     cout << "\tCorresponding " << Ns + 1 << " nodes: ";
//     for(int i = 0; i < new_tree.top_tnodes.size(); i++){
//         cout << new_tree.top_tnodes[i] + 1 << "\t" << "\t" << new_tree.node_times[new_tree.top_tnodes[i]] << endl;
//     }
//
//     vector<double> blens = new_tree.get_edges_from_interval(new_tree.top_tinvls, new_tree.top_tnodes);
//     cout << "\tNew branch lengths: ";
//     for(int i = 0; i < blens.size(); i++){
//         cout << i + 1 << "\t" << "\t" << blens[i] << endl;
//     }
// }

#endif
