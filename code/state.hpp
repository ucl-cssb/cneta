#ifndef STATE_HPP
#define STATE_HPP


#include "likelihood.hpp"

// using namespace std;

void set_pmat(const evo_tree& rtree, int Ns, int nstate, int model, int cn_max, const vector<int>& knodes, vector<double>& blens, vector<double*>& pmat_per_blen, ofstream& fout);

void set_pmat_decomp(const evo_tree& rtree, MAX_DECOMP& max_decomp, int nstate, const vector<int>& knodes, DIM_DECOMP& dim_decomp, PMAT_DECOMP& pmat_decomp, ofstream& fout);

void print_tree_state(const evo_tree& rtree, const vector<vector<int>>& S_sk_k, int nstate);

// Get the states of the tree from likelihood table at one site, starting from MRCA
// Step 5 of algorithm of Pupko (2000)
void extract_tree_ancestral_state(const evo_tree& rtree, const set<vector<int>>& comps, const vector<vector<double>>& L_sk_k, const vector<vector<int>>& S_sk_k, int model, int cn_max, int is_total, int m_max, map<int, int> &asr_states);

// Create likelihood vectors and state vectors at the tip node for reconstructing joint ancestral state
// L_sk_k (S_sk_k) has one row for each tree node and one column for each possible state
// Step 1 of algorithm of Pupko (2000)
void initialize_asr_table(const vector<int>& obs, const evo_tree& rtree, map<double, double*>& pmats, vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, int model, int nstate, int is_total);

// Create likelihood vectors and state vectors at the tip node for reconstructing joint ancestral state, for independent chain model
// L_sk_k (S_sk_k) has one row for each tree node and one column for each possible state
void initialize_asr_table_decomp(const vector<int>& obs, const evo_tree& rtree, const set<vector<int>>& comps, MAX_DECOMP& max_decomp, PMAT_DECOMP& pmat_decomp, vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, int nstate, int is_total = 1);

// Find the most likely state for a node under each possible state, assuming current node has state nsk and parent node (connected by branch of length blen) has state np
// Step 2 of algorithm of Pupko (2000)
double get_max_prob_children(const vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, double* pblen, int k, int nstate, int sp, int ni, int nj, int blen, int model);

// Assume the likelihood table is for each combination of states
double get_max_children_decomp2(const vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, int k, int nstate, double* pbli_wgd, double* pbli_chr, double* pbli_seg, DIM_DECOMP& dim_decomp, int sp, int ni, int nj, int blen);

// Get the most likely state on one site of a chromosome (assuming higher level events on nodes)
// L_sk_k: the likelihood of the best reconstruction of the subtree at node sk provided that the father of sk has state k
// S_sk_k: the state of node sk provided that the father of sk has state k in the optimal conditional reconstruction
// Step 3 of algorithm of Pupko (2000)
void get_ancestral_states_site(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const vector<int>& knodes, const set<vector<int>>& comps, map<double, double*>& pmats, int nstate, int model);

// Get the probability of all possible copy numbers for one node at a site, for marginal reconstruction
string get_prob_line(const vector<vector<double>>& L_sk_k, int nid, int nchr, int nc, int is_total, int cn_max);

// Get the most likely copy number for each site at a specific node, for marginal reconstruction
void get_site_cnp(const vector<vector<double>>& L_sk_k, int nid, int nchr, int nc, int is_total, int cn_max, copy_number& cnp);

void print_node_cnp(ofstream& fout, const copy_number& cnp, int nid, int cn_max, int is_total);

// Get the ancestral state on one site of a chromosome
// Assuming each observed copy number is composed of three type of events.
// Sum over all possible states for initial and final nodes
// Allow at most one WGD event along a branch
void get_ancestral_states_site_decomp(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const vector<int>& knodes, const set<vector<int>>& comps, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int nstate);


// Infer the copy number of the MRCA given a tree at a site, assuming independent Markov chains
double reconstruct_marginal_ancestral_state_decomp(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, const vector<int>& knodes, const set<vector<int>>& comps, OBS_DECOMP& obs_decomp, int use_repeat, int infer_wgd, int infer_chr, int cn_max, string ofile, double min_asr, int is_total = 1);

// Infer the copy number of the MRCA given a tree at a site, assuming only segment duplication/deletion
double reconstruct_marginal_ancestral_state(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, const vector<int>& knodes, int model, int cn_max, int use_repeat, int is_total, string ofile, double min_asr);

// Infer the copy number of all internal nodes given a tree at a site, assuming independent Markov chains
void reconstruct_joint_ancestral_state_decomp(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, vector<int>& knodes, const set<vector<int>>& comps, MAX_DECOMP& max_decomp, int use_repeat, int cn_max, string ofile, double min_asr, int is_total = 1);

// Infer the copy number of all internal nodes given a tree at a site, assuming only segment duplication/deletion
void reconstruct_joint_ancestral_state(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, vector<int>& knodes, int model, int cn_max, int use_repeat, int is_total, int m_max, string ofile, double min_asr);


#endif
