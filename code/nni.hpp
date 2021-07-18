#ifndef NNI_HPP
#define NNI_HPP

#include "evo_tree.hpp"
#include "tree_op.hpp"
#include "model.hpp"
#include "likelihood.hpp"
#include "optimization.hpp"

// using namespace std;


inline bool is_a_branch(Node* node1, Node* node2){
    return (node1->findNeighbor(node2) != NULL && node2->findNeighbor(node1) != NULL);
}


inline bool is_inner_branch(Node* node1, Node* node2){
    return(node1->degree() >= 3 && node2->degree() >= 3 && is_a_branch(node1, node2));
}


// Two NNIs are considered conflicting if they operate on the same inner branch or adjacent branches.
void get_compatible_NNIs(vector<NNIMove> &nniMoves, vector<NNIMove> &compatibleNNIs);


// check age of the sibling of node
// find the top node u with one child v and another child c, only feasible when t(c) > t(v)
bool is_valid_NNI(const evo_tree& rtree, const Branch& curBranch);


// Get all branches where NNI is feasible, namely all internal branches
// NNI is only feasible when it does not violate parent-child age constraint
void get_NNI_branches(const evo_tree& rtree, Branches &nniBranches, Node* node, Node* dad);

// Find branches at most "depth" branches away from the tagged branch
void get_neighbor_inner_branches(const evo_tree& rtree, Node* node, Node* dad, int depth, Branches &surrBranches);


// Find branches at most two branches away from the tagged branch
void filter_NNI_branches(const evo_tree& rtree, vector<NNIMove> &appliedNNIs, Branches &nniBranches);


// Apply one NNI move (adjust branch length to satifisty time constaints when cons is true)
// Assume neighbors have been generated for the tree
void do_one_NNI(evo_tree& rtree, NNIMove& move, int cons);


// Update branch lengths related to the NNI move
// include four adjacent branches when nni5 is true
// used in do_all_NNIs
void change_NNI_Brans(evo_tree& rtree, NNIMove& nnimove, bool nni5);


// Simultaneously apply all NNIs, assigning new branch lengths to related branches
// Mutation rates are estimated at the same time if maxj = 1
void do_all_NNIs(evo_tree& rtree, vector<NNIMove> &compatibleNNIs, bool changeBran, bool nni5, int cons);


NNIMove get_random_NNI(Branch& branch, gsl_rng* r);


// do stochastic NNI
void do_random_NNIs(evo_tree& rtree, gsl_rng* r, int cons);


/**
   Search for the best NNI move corresponding to the chosen inner branch
   @return NNIMove the best NNI, this NNI could be worse than the current tree
   according to the evaluation scheme in use
   @param node1 1 of the 2 nodes on the branch
   @param node2 1 of the 2 nodes on the branch
   @param nniMoves (IN/OUT) detailed information of the 2 NNIs
   adapted from IQ-TREE package, phylotree.cpp
   need to compute likelihood
 */
NNIMove get_best_NNI_for_bran(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, NNIMove* nniMoves = NULL, bool nni5 = true);

// Find NNI increasing likelihood of current tree
void evaluate_NNIs(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, Branches &nniBranches, vector<NNIMove> &positiveNNIs, double curScore);


// Apply hill climbing perturbation to obtain a locally optimal tree (by NNI)
// score used in this function is log likelihood, the larger the better
// need to compute likelihood
void do_hill_climbing_NNI(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double loglh_epsilon, int speed_nni, bool nni5 = false);



/**************************** NNI operations used in MCMC approach ****************************/
// Use NNI for rooted (clock) trees
evo_tree do_random_NNI(evo_tree& rtree, gsl_rng* r, int debug);


/**************************** NNI operations used in test ****************************/
// Do NNI on an internal branch i, according to P293 (Yang, 2014)
vector<evo_tree> do_NNI_on_branch(evo_tree& rtree, int i, vector<NNIMove>& moves, int debug);

// Find all the NNI trees for a given tree
map<int, vector<evo_tree>> get_all_NNI_trees(evo_tree& rtree, vector<NNIMove>& moves, int debug);


#endif
