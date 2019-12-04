// run maximum likelihood inference

#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <sstream>
#include <ctime>
#include <map>
#include <unistd.h>


#ifdef _OPENMP
#include <omp.h>    // used for accelerating tree search
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
//#include <gsl/gsl_multimin.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "gzstream.h"

#include "evo_tree.hpp"
#include "genome.hpp"
#include "stats.hpp"
#include "utilities.hpp"


using namespace std;


// The number of tree shapes
static const int num_shapes[] = {1, 1, 1, 2, 3, 6, 11, 23, 46, 98, 207, 451, 983, 2179, 4850, 10905, 24631, 56011, 127912, 293547};
static const int num_trees[] = {1, 1, 3, 15, 105, 945, 10395, 135135, 2027025, 34459425, 654729075};
// The maximum number of samples in the tree which allows exhaustive search
const int LARGE_TREE = 11;
// The number of trees to search before terminating
const int MAX_TREE = 100;
// The maximum number of trees to perturb
const int MAX_TREE1 = 20;
// The maximum number of trees to refine
const int MAX_TREE2 = 5;
// The maximum number of times to refine the final set of trees
const int MAX_PERTURB = 100;
const double MAX_NLNL = 1e20;
const int MAX_OPT = 10; // max number of optimization for each tree
const int PRINT_PRECISION = 10;
int debug = 0;

// global value for tree search
map<string,int> searched_trees;
// map<string,int> searched_shapes;
bool exhausted_tree_search = false;


// Randomly pick a tree to perturb
evo_tree perturb_tree_set( const int& Ns, vector<evo_tree> trees ){
    int debug = 0;
    if(debug) cout << "\tperturb a set of trees" << endl;

    int count = 0;
    while(true){
        // randomly sample the fit population
        int ind = gsl_rng_uniform_int(r, trees.size());

        // generate a new tree
        evo_tree ttree = perturb_tree( Ns, trees[ind] );
        adjust_tree_blens(ttree);
        adjust_tree_height(ttree);
        adjust_tree_tips(ttree);

        string tstring = order_tree_string_uniq(create_tree_string_uniq( ttree ) );
        if ( searched_trees.find(tstring) == searched_trees.end() ) {
          // mark the tree
          searched_trees[ tstring ] = 0;
          return ttree;
        }
        else {
          //cout << "Tree already present" << endl;
          count++;
        }

        if(count > MAX_TREE){
          //cout << "\tperturb_tree cannot find new topologies" << endl;
          exhausted_tree_search = true;
          return ttree;
        }
    }
}

void print_invl_blen(evo_tree& new_tree, string header){
    cout << header << endl;
    new_tree.print();

    cout << "\tTop " << new_tree.nleaf - 1 << " time intervals: ";
    for(int i = 0; i < new_tree.top_tinvls.size(); i++){
        cout << "\t" << new_tree.top_tinvls[i];
    }
    cout<< endl;
    cout << "\tCorresponding " << Ns + 1 << " nodes: ";
    for(int i = 0; i < new_tree.top_tnodes.size(); i++){
        cout << new_tree.top_tnodes[i] + 1 << "\t" << "\t" << new_tree.node_times[new_tree.top_tnodes[i]] << endl;
    }

    vector<double> blens = new_tree.get_edges_from_interval(new_tree.top_tinvls, new_tree.top_tnodes);
    cout << "\tNew branch lengths: ";
    for(int i = 0; i<blens.size(); i++){
        cout << i + 1 << "\t" << "\t" << blens[i] << endl;
    }
}

// Apply hill climbing perturbation to obtain a locally optimal tree (by tip swapping)
evo_tree get_local_optimal_tree(evo_tree& tree, int Ngen, int max_perturb, int max_static, const double ssize, const double tolerance, const int miter, const int optim, const int model, const int cons, const int maxj, const int cn_max, const int only_seg, const int correct_bias, int is_total=1){
    int count_static = 0;
    int count = 0;
    int debug = 0;
    if(debug) cout << "Apply hill climbing perturbation to obtain a locally optimal tree" << endl;
    while(count <= Ngen && count_static < max_static ){
        // perturb T by NNI to get T'
        evo_tree ttree = perturb_tree(Ns, tree);
        string tstring = order_tree_string_uniq(create_tree_string_uniq(ttree));
        // perturb the tree until it is new
        // int count_tree = 1;
        // while (searched_trees.find(tstring) != searched_trees.end() && count_tree < max_perturb){
        //     ttree = perturb_tree(Ns, ttree);
        //     tstring = order_tree_string_uniq(create_tree_string_uniq(ttree));
        //     count_tree += 1;
        // }
        if ( searched_trees.find(tstring) == searched_trees.end() ) {
            // mark the tree
            searched_trees[ tstring ] = 0;
        }
        // cout << "Perturb " << count_tree << " times for tree " << index[i] << endl;
        if(cons){
            adjust_tree_blens(ttree);
            adjust_tree_height(ttree);
            adjust_tree_tips(ttree);
        }

        // optimize T' to get T"
        evo_tree otree;
        double Lf = 0;
        if(optim == 0){
            while(!(Lf>0)){
                Lf = 0;
                otree = max_likelihood(ttree, model, Lf, ssize, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
            }
        }
        if(optim == 1){
            while(!(Lf>0)){
                Lf = 0;
                otree = max_likelihood_BFGS(ttree, Lf, tolerance, miter);
            }
        }
        otree.score = Lf;
        assert(tree.score > 0);
        if(Lf < tree.score){
            // cout << "tree score before " << tree.score << endl;
            tree = otree;
            // cout << "tree score after " << tree.score << endl;
        }
        else{
            count_static += 1;
        }
        // A tree has been maximized at least five times
        count += 1;
        // stop if maximum value does not improve for several times
        if( max_static > 0 && count_static == max_static ){
          cout << "\t### static likelihood value. Finishing optimizing tree on ngen = " << count << endl;
        }
    }
    return tree;
}

// Build parsimony tree from copy number changes (breakpoints)
evo_tree build_parsimony_tree(int Ns, vector<vector<int>> data){
    vector<int> nodes;
    vector<edge> edges;

    evo_tree ptree = evo_tree(Ns + 1, edges);

    return ptree;
}


// Read parsimony trees built by other tools
evo_tree read_parsimony_tree(const string& tree_file, int Ns, const vector<double>& rates, vector<double>& tobs){
    int debug = 0;
    if(debug)   cout << tree_file << endl;
    evo_tree rtree = read_tree_info(tree_file, Ns);
    // for(int i=0; i<tobs.size();i++){
    //     cout << tobs[i] << endl;
    // }
    rtree.tobs = tobs;

    // The branch lengths in parsimony tree may be very large
    adjust_tree_blens(rtree);
    adjust_tree_height(rtree);
    adjust_tree_tips(rtree);
    if(debug) rtree.print();

    if(rates.size()>1){
      rtree.dup_rate = rates[0];
      rtree.del_rate = rates[1];
      if(rates.size() > 2){
          rtree.chr_gain_rate = rates[2];
          rtree.chr_loss_rate = rates[3];
          rtree.wgd_rate = rates[4];
      }
    }
    else{
      rtree.mu = rates[0];
    }

    return rtree;
}


// Generate initial set of unique trees, at most Npop trees
vector<evo_tree> get_initial_trees(int init_tree, string dir_itrees, int Npop, const vector<double>& rates, int max_tree_num, int cons, int Ne = 1, double beta = 0, double gtime=1){
    int debug = 0;
    vector<evo_tree> trees;

    if(init_tree == 1){
        assert(dir_itrees != "");
        string fname;
        boost::filesystem::path p(dir_itrees);
        for (auto&& x : boost::filesystem::directory_iterator(p)){
            fname = x.path().string();
            evo_tree rtree = read_parsimony_tree(fname, Ns, rates, tobs);
            string tstring = order_tree_string_uniq(create_tree_string_uniq(rtree));

            if ( searched_trees.find(tstring) == searched_trees.end() ) {
                // mark the tree
                searched_trees[ tstring ] = 0;
            }
            rtree.score = MAX_NLNL;
            trees.push_back(rtree);
        }
    }
    else{
        int num_tree = 0;
        int n =  (max_tree_num < Npop) ? max_tree_num: Npop;
        while(num_tree < n){
            evo_tree rtree;
            rtree = generate_coal_tree(Ns, Ne, beta, gtime);
            // string tstring = rtree.make_newick(0);
            // tstring.erase(remove_if(tstring.begin(), tstring.end(), [](char c) { return !(c == '(' || c == ')'); }), tstring.end());
            string tstring = order_tree_string_uniq(create_tree_string_uniq(rtree));

            if ( searched_trees.find(tstring) == searched_trees.end() ) {
                // mark the tree
                searched_trees[ tstring ] = 0;
                num_tree += 1;
            }
            else{
                continue;
            }

            if(cons){
                rtree.tobs = tobs;
                adjust_tree_blens(rtree);
                adjust_tree_height(rtree);
                adjust_tree_tips(rtree);
            }

            assert(is_tree_valid(rtree, cons));

            if(rates.size()>1){
              rtree.dup_rate = rates[0];
              rtree.del_rate = rates[1];
              if(!only_seg){
                  rtree.chr_gain_rate = rates[2];
                  rtree.chr_loss_rate = rates[3];
                  rtree.wgd_rate = rates[4];
              }
            }
            else{
              rtree.mu = rates[0];
            }
            rtree.score = MAX_NLNL;
            trees.push_back( rtree );
        }
    }

    if(debug){
        ofstream out_tree("./inital_trees.txt");
        for(int i=0; i < trees.size(); ++i){
            int precision = 5;
            string newick = trees[i].make_newick(precision);
            out_tree << newick << endl;
            out_tree << order_tree_string_uniq(create_tree_string_uniq(trees[i])) << endl;
        }
        out_tree.close();
    }

    return trees;
}


// There may be multiple trees with the same likelihood, getting the one with minimal score (negative log likelihood)
vector<evo_tree> find_best_trees(const vector<evo_tree>& trees, const vector<double>& lnLs, vector<int>& index, int n){
    vector<evo_tree> btrees;
    int x=0;

    iota( index.begin(), index.end(), x++);
    sort( index.begin(), index.end(), [&](int i,int j){ return lnLs[i] < lnLs[j]; } );

    for(int i=0; i < n; ++i){
        btrees.push_back(trees[index[i]]);
    }

    return btrees;
}


void check_matrix_row_sum(double *mat, int nstate){
    for(int i = 0; i < nstate; i++){
        double sum = 0;
        for(int j = 0; j < nstate; j++){    // jth column
            sum += mat[i + j * nstate];
        }
        cout << "Total probability of changing from " << i << " is " << sum << endl;
    }
}

// Check the relationship among number of segments, mutation rates and branch lengths
void check_chain_by_branch_m2(int num_seg, const vector<double>& rates, double blen, int cn_max, ofstream& fout){
    double dup_rate = rates[0];
    double del_rate = rates[1];
    int nstate = (cn_max + 1) * (cn_max + 2) / 2;
    double *qmat = new double[(nstate)*(nstate)];
    memset(qmat, 0, (nstate)*(nstate)*sizeof(double));
    get_rate_matrix_allele_specific(qmat, dup_rate, del_rate, cn_max);
    // r8mat_print( nstate, nstate, qmat, "  Q matrix:" );
    cout << "Checking row sum for Q matrix" << endl;
    check_matrix_row_sum(qmat, nstate);

    double *pmat = new double[(nstate)*(nstate)];
    memset(pmat, 0, (nstate)*(nstate)*sizeof(double));
    get_transition_matrix_bounded(qmat, pmat, blen, nstate);
    // r8mat_print( nstate, nstate, pmat, "  P matrix:" );
    cout << "Checking row sum for P matrix" << endl;
    check_matrix_row_sum(pmat, nstate);

    // Print the probability of changing to each state
    map<int, double> state_prob;
    for(int i = 0; i < cn_max; i++){
        state_prob[i] = 0;
    }
    // for(int i = 0; i < nstate; i++){    // ith column
    //     int cn = state_to_total_cn(i, cn_max);
    //     double total_prob = 0;
    //     for(int j = 0; j < nstate; j++){    // jth row
    //         int idx = j + i * nstate;
    //         double prob = pmat[idx];
    //         total_prob += prob;    // assume each state is equally likely
    //         // cout << "Probability of changing from state " << j << " to state " << i << " is " << prob << endl;
    //     }
    //     state_prob[cn] += total_prob;
    // }
    //


    // Suppose the starting state is normal
    int s = 4;
    for(int j = 0; j < nstate; j++){    // jth column
        int cn = state_to_total_cn(j, cn_max);
        int idx = s + j * nstate;
        double prob = pmat[idx];
        state_prob[cn] += prob;
    }

    fout << cn_max << "\t"  << blen;
    for(int i = 0; i < rates.size(); i++){
        fout << "\t" << rates[i];
    }
    for(auto sp : state_prob){
        int cn = sp.first;
        double prob = sp.second;
        fout <<  "\t" << prob;
        cout << "Probability of changing from normal to state " << cn << " is " << prob << endl;
    }
    fout << endl;

    free(qmat);
    free(pmat);
}


// Check the relationship among number of segments, mutation rates and branch lengths
void check_chain_by_branch_m3(int is_total, int num_seg, const vector<double>& rates, double blen, int cn_max, int m_max, int max_wgd, int max_chr_change, int max_site_change, ofstream& fout){
    double dup_rate = rates[0];
    double del_rate = rates[1];
    double chr_gain_rate = rates[2];
    double chr_loss_rate = rates[3];
    double wgd_rate = rates[4];

    int nstate = (cn_max + 1) * (cn_max + 2) / 2;
    int dim_wgd = max_wgd + 1;
    int dim_chr = 2 * max_chr_change + 1;
    int dim_seg = 2 * max_site_change + 1;
    double *qmat_wgd, *qmat_chr, *qmat_seg;
    double *pmat_wgd, *pmat_chr, *pmat_seg;
    int delta_chr, delta_seg;  // The indices for chromosome and segment matrix have to be ajusted

    cout << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << "\n";
    cout << "Checking row sum for Q matrix" << endl;

    // For WGD model
    if(max_wgd > 0){
        qmat_wgd = new double[dim_wgd * dim_wgd];  // WGD
        memset(qmat_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
        get_rate_matrix_wgd(qmat_wgd, wgd_rate, max_wgd);
        check_matrix_row_sum(qmat_wgd, dim_wgd);

        pmat_wgd = new double[(dim_wgd)*(dim_wgd)];
        memset(pmat_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
        get_transition_matrix_bounded(qmat_wgd, pmat_wgd, blen, dim_wgd);
        check_matrix_row_sum(pmat_wgd, dim_wgd);
    }

    if(max_chr_change>0){
        qmat_chr = new double[(dim_chr)*(dim_chr)];   // chromosome gain/loss
        memset(qmat_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
        get_rate_matrix_chr_change(qmat_chr, chr_gain_rate, chr_loss_rate, max_chr_change);
        check_matrix_row_sum(qmat_chr, dim_chr);
        pmat_chr = new double[(dim_chr)*(dim_chr)];
        memset(pmat_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
        get_transition_matrix_bounded(qmat_chr, pmat_chr, blen, dim_chr);
        check_matrix_row_sum(pmat_chr, dim_chr);
        delta_chr = (dim_chr - 1)/2;
    }

    if(max_site_change>0){
        qmat_seg = new double[(dim_seg)*(dim_seg)];  // segment duplication/deletion
        memset(qmat_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
        get_rate_matrix_site_change(qmat_seg, dup_rate, del_rate, max_site_change);
        check_matrix_row_sum(qmat_seg, dim_seg);
        pmat_seg = new double[(dim_seg)*(dim_seg)];
        memset(pmat_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
        get_transition_matrix_bounded(qmat_seg, pmat_seg, blen, dim_seg);
        check_matrix_row_sum(pmat_seg, dim_seg);
        delta_seg = (dim_seg - 1)/2;
    }

    // Print the probability of changing to each state
    map<int, double> state_prob;
    for(int i = 0; i < cn_max; i++){
        state_prob[i] = 0;
    }

    decomp_table = build_decomp_table();
    int s_wgd = 0;
    int s_chr = 0;
    int s_seg = 0;
    int s_chr2 = 0;
    int s_seg2 = 0;

    cout << "Start state: " << s_wgd << "\t" << s_chr << "\t" << s_seg << "\n";

    for(int cn = 0; cn <= cn_max; cn++){    // jth column
        set<vector<int>> comp = decomp_table[cn];
        cout << "Copy number " << cn << endl;
        for(auto e : comp){
            double prob = 1;
            double prob_wgd = 0;
            double prob_chr = 0, prob_chr2 = 0;
            double prob_seg = 0, prob_seg2 = 0;

            int e_wgd = e[0];
            int e_chr = e[1];
            int e_seg = e[2];
            int e_chr2 = 0;
            int e_seg2 = 0;

            if(max_wgd > 0){
                prob_wgd = pmat_wgd[s_wgd + e_wgd * dim_wgd];
                prob = prob_wgd;
            }
            if(max_chr_change > 0){
                prob_chr = pmat_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
            }
            if(max_site_change > 0){
                prob_seg = pmat_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
            }

            if(is_total == 1){
                if(max_chr_change > 0){
                    prob = prob * prob_chr;
                }
                if(max_site_change > 0){
                    prob = prob * prob_seg;
                }
                cout << "   End state: " << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << prob << endl;
            }else{
                e_seg = e[3];
                e_chr2 = e[2];
                e_seg2 = e[4];

                if(max_chr_change > 0){
                    prob_chr2 = pmat_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
                    prob = prob * prob_chr * prob_chr2;
                }
                if(max_site_change > 0){
                    prob_seg = pmat_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
                    prob_seg2 = pmat_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
                    prob = prob * prob_seg * prob_seg2;
                }
                cout << "   End state: " << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2 << "\t" << prob << endl;
            }
            state_prob[cn] += prob;
        }
        if(comp.size()==0) state_prob[cn] = 0;
    }

    fout << cn_max << "\t"  << blen;
    for(int i = 0; i < rates.size(); i++){
        fout << "\t" << rates[i];
    }
    for(auto sp : state_prob){
        int cn = sp.first;
        double prob = sp.second;
        fout <<  "\t" << prob;
        cout << "Probability of changing from normal to state " << cn << " is " << prob << endl;
    }
    fout << endl;

    if(max_wgd>0){
        free(qmat_wgd);
        free(pmat_wgd);
    }
    if(max_chr_change>0){
        free(qmat_chr);
        free(pmat_chr);
    }
    if(max_site_change>0){
        free(qmat_seg);
        free(pmat_seg);
    }
}


void get_compatible_NNIs(vector<NNIMove> &nniMoves, vector<NNIMove> &compatibleNNIs) {
    int debug = 0;
    compatibleNNIs.clear();
	for (vector<NNIMove>::iterator it1 = nniMoves.begin(); it1 != nniMoves.end(); it1++) {
		bool select = true;
		for (vector<NNIMove>::iterator it2 = compatibleNNIs.begin(); it2 != compatibleNNIs.end(); it2++) {
			if ((*it1).node1 == (*(it2)).node1
					|| (*it1).node2 == (*(it2)).node1
					|| (*it1).node1 == (*(it2)).node2
					|| (*it1).node2 == (*(it2)).node2) {
				select = false;
                break;
            }
        }
		if (select) {
            compatibleNNIs.push_back(*it1);
        }
    }
    if(debug){
        cout << "There are " << compatibleNNIs.size() << " compatible NNIs" << endl;
        for(int i = 0; i < compatibleNNIs.size(); i++){
            cout << compatibleNNIs[i].node1->id + 1 << "\t" << compatibleNNIs[i].node2->id + 1 << endl;
        }
    }
}


bool is_a_branch(Node *node1, Node *node2) {
    return (node1->findNeighbor(node2) != NULL && node2->findNeighbor(node1) != NULL);
}

bool is_inner_branch(Node *node1, Node *node2){
    return(node1->degree() >= 3 && node2->degree() >= 3 && is_a_branch(node1, node2));
}

void get_NNI_branches(Branches &nniBranches, Node *node, Node *dad) {
    assert(node != NULL);
    FOR_NEIGHBOR_IT(node, dad, it) {
            if (is_inner_branch((*it)->node, node)) {
                Branch curBranch;
                curBranch.first = (*it)->node;
                curBranch.second = node;
                int branchID = pairInteger(curBranch.first->id, curBranch.second->id);
                nniBranches.insert(pair<int, Branch>(branchID, curBranch));
            }
            get_NNI_branches(nniBranches, (*it)->node, node);
        }
}


// Find branches at most "depth" branches away from the tagged branch
void get_neighbor_inner_branches(Node *node, Node *dad, int depth, Branches &surrBranches) {
    cout << "Find branches at most certain branches away from the tagged branch" << endl;
    if (depth == 0)
      return;
      FOR_NEIGHBOR_IT(node, dad, it) {
          if (!(*it)->node->is_leaf()) {
              Branch curBranch;
              curBranch.first = node;
              curBranch.second = (*it)->node;
              int branchID = pairInteger(node->id, (*it)->node->id);
              if (surrBranches.find(branchID) == surrBranches.end())
                  surrBranches.insert(pair<int,Branch>(branchID, curBranch));
              get_neighbor_inner_branches((*it)->node, node, depth-1, surrBranches);
          }
      }
}

// Find branches at most two branches away from the tagged branch
void filter_NNI_branches(vector<NNIMove> &appliedNNIs, Branches &nniBranches) {
    cout << "Find branches at most two branches away from the tagged branch" << endl;
    for (vector<NNIMove>::iterator it = appliedNNIs.begin(); it != appliedNNIs.end(); it++) {
        Branch curBranch;
        curBranch.first = it->node1;
        curBranch.second = it->node2;
        int branchID = pairInteger(it->node1->id, it->node2->id);
        if (nniBranches.find(branchID) == nniBranches.end())
            nniBranches.insert(pair<int,Branch>(branchID, curBranch));
        get_neighbor_inner_branches(it->node1, it->node2, 2, nniBranches);
        get_neighbor_inner_branches(it->node2, it->node1, 2, nniBranches);
    }
}

// Apply one NNI move
void do_one_NNI(NNIMove &move) {
    int debug = 0;

    if(debug) cout << " Apply one NNI move" << endl;
    Node *node1 = move.node1;
    Node *node2 = move.node2;
    assert(node1->degree() == 3 && node2->degree() == 3);

    NeighborVec::iterator node1Nei_it = move.node1Nei_it;
    NeighborVec::iterator node2Nei_it = move.node2Nei_it;
    Neighbor *node1Nei = *(node1Nei_it);
    Neighbor *node2Nei = *(node2Nei_it);

    NeighborVec::iterator it;
    if(debug){
        cout << " old neighbors of node " << node1->id + 1 << endl;
        FOR_NEIGHBOR(node1, node2, it)
            cout << "\t" << (*it)->node->id + 1;
        cout << endl;
        cout << " old neighbors of node " << node2->id + 1 << endl;
        FOR_NEIGHBOR(node2, node1, it)
            cout << "\t" << (*it)->node->id + 1;
        cout << endl;
    }

    // do the NNI swap
    node1->updateNeighbor(node1Nei_it, node2Nei);
    node2Nei->node->updateNeighbor(node2, node1);

    node2->updateNeighbor(node2Nei_it, node1Nei);
    node1Nei->node->updateNeighbor(node1, node2);

    if(debug){
        cout << " new neighbors of node " << node1->id + 1 << endl;
        FOR_NEIGHBOR(node1, node2, it)
            cout << "\t" << (*it)->node->id + 1;
        cout << endl;
        cout << " new neighbors of node " << node2->id + 1 << endl;
        FOR_NEIGHBOR(node2, node1, it)
            cout << "\t" << (*it)->node->id + 1;
        cout << endl;
    }
}


// Update branch lengths related to the NNI move
void change_NNI_Brans(NNIMove nnimove, bool nni5) {
    int debug = 0;
    if(debug) cout << "  Update branch lengths related to the NNI move" << endl;
	Node *node1 = nnimove.node1;
	Node *node2 = nnimove.node2;

    NeighborVec::iterator it;
    if(debug)
    {
        cout << " old neighbors of node " << node1->id + 1 << endl;
        FOR_NEIGHBOR(node1, node2, it)
            cout << "\t" << (*it)-> node->id + 1 << "," << (*it)->id + 1 << "," << (*it)->length;
        cout << endl;
        cout << " old neighbors of node " << node2->id + 1 << endl;
        FOR_NEIGHBOR(node2, node1, it)
            cout << "\t" << (*it)-> node->id + 1 << "," << (*it)->id + 1 << "," << (*it)->length;
        cout << endl;
    }


	Neighbor *node1_node2_nei = (Neighbor*) node1->findNeighbor(node2);
	Neighbor *node2_node1_nei = (Neighbor*) node2->findNeighbor(node1);
	node1_node2_nei->setLength(nnimove.newLen[0]);
	node2_node1_nei->setLength(nnimove.newLen[0]);

    if(debug)
    {
        cout << " new neighbors of node " << node1->id + 1 << endl;
        FOR_NEIGHBOR(node1, node2, it)
            cout << "\t" << (*it)->node->id + 1 << "\t" << (*it)->id + 1 << "\t" << (*it)->length;
        cout << endl;
        cout << " new neighbors of node " << node2->id + 1 << endl;
        FOR_NEIGHBOR(node2, node1, it)
            cout << "\t" << (*it)->node->id + 1 << "\t" << (*it)->id + 1 << "\t" << (*it)->length;
        cout << endl;
    }

	if (nni5) {
		int i = 1;
		Neighbor* nei;
		Neighbor* nei_back;
		NeighborVec::iterator it;
		FOR_NEIGHBOR(node1, node2, it)
		{
			nei = (*it)->node->findNeighbor(node1);
			nei_back = (node1)->findNeighbor((*it)->node);
			nei->setLength(nnimove.newLen[i]);
			nei_back->setLength(nnimove.newLen[i]);
			i++;
		}
		FOR_NEIGHBOR(node2, node1, it)
		{
			nei = (*it)->node->findNeighbor(node2);
			nei_back = (node2)->findNeighbor((*it)->node);
			nei->setLength(nnimove.newLen[i]);
			nei_back->setLength(nnimove.newLen[i]);
			i++;
		}
	}

}

// Simultaneously apply all NNIs, assigning new branch lengths to related branches
void do_all_NNIs(vector<NNIMove> &compatibleNNIs, bool changeBran, bool nni5) {
    // cout << "Simultaneously apply all NNIs, assigning new branch lengths to related branches" << endl;
    for (vector<NNIMove>::iterator it = compatibleNNIs.begin(); it != compatibleNNIs.end(); it++) {
		do_one_NNI(*it);
        if (changeBran) {
            // apply new branch lengths
			change_NNI_Brans(*it, nni5);
        }
    }
}


bool isABranch(Node* node1, Node* node2) {
    return (node1->findNeighbor(node2) != NULL && node2->findNeighbor(node1) != NULL);
}

bool isInnerBranch(Node* node1, Node* node2) {
    return(node1->degree() >= 3 && node2->degree() >= 3 && isABranch(node1, node2));
}

NNIMove get_random_NNI(Branch &branch) {
    assert(isInnerBranch(branch.first, branch.second));
    // for rooted tree
    if (((Neighbor*)branch.first->findNeighbor(branch.second))->direction == TOWARD_ROOT) {
        // swap node1 and node2 if the direction is not right, only for nonreversible models
        Node *tmp = branch.first;
        branch.first = branch.second;
        branch.second = tmp;
    }
    NNIMove nni;
    nni.node1 = (Node*) branch.first;
    nni.node2 = (Node*) branch.second;

    FOR_NEIGHBOR_IT(branch.first, branch.second, node1NeiIt)
        if (((Neighbor*)*node1NeiIt)->direction != TOWARD_ROOT) {
            nni.node1Nei_it = node1NeiIt;
            break;
        }
    // int randInt = random_int(branch.second->neighbors.size()-1);
    int randInt = myrng(branch.second->neighbors.size()-1);
    int cnt = 0;
    FOR_NEIGHBOR_IT(branch.second, branch.first, node2NeiIt) {
        // if this loop, is it sure that direction is away from root because node1->node2 is away from root
        if (cnt == randInt) {
            nni.node2Nei_it = node2NeiIt;
            break;
        } else {
            cnt++;
        }
    }
	assert(*nni.node1Nei_it != NULL && *nni.node2Nei_it != NULL);
    assert(((Neighbor*)*nni.node1Nei_it)->direction != TOWARD_ROOT && ((Neighbor*)*nni.node2Nei_it)->direction != TOWARD_ROOT);
    nni.newloglh = 0.0;
    return nni;
}

// do stochastic NNI
evo_tree do_random_NNIs(evo_tree &rtree) {
    int debug = 0;
    if(debug) cout << "doing stochastic NNIs" << endl;

    int cntNNI = 0;
    int numRandomNNI;
    Branches nniBranches;
    numRandomNNI = floor((rtree.nleaf - 3) * 0.5);
    if (rtree.nleaf >= 4 && numRandomNNI == 0)
        numRandomNNI = 1;

    evo_tree ntree(rtree);

    while (cntNNI < numRandomNNI) {
        nniBranches.clear();
        get_NNI_branches(nniBranches, &ntree.nodes[rtree.nleaf], NULL);
        if (nniBranches.size() == 0) break;
        // Convert the map data structure Branches to vector of Branch
        vector<Branch> vectorNNIBranches;
        for (Branches::iterator it = nniBranches.begin(); it != nniBranches.end(); ++it) {
            vectorNNIBranches.push_back(it->second);
        }
        // int randInt = random_int((int) vectorNNIBranches.size());
        int randInt = myrng(vectorNNIBranches.size());
        NNIMove randNNI = get_random_NNI(vectorNNIBranches[randInt]);
        do_one_NNI(randNNI);
        cntNNI++;
    }
    if (debug)
	    cout << "Tree perturbation: number of random NNI performed = " << cntNNI << endl;

    return ntree;
}


// Optimizing the branch (node1, node2)
void optimize_one_branch(evo_tree& rtree, Node *node1, Node *node2) {
    int debug = 0;
    if(debug) cout << "\tOptimizing the branch " << node1->id + 1 << ", " << node2->id + 1 << endl;
    if ((node1->id == Ns + 1 || node2->id == Ns + 1))
        return; // does not optimize virtual branch from root

    // current_it and current_it_back point to the same branch
    rtree.current_it = (Neighbor*) node1->findNeighbor(node2);
    assert(rtree.current_it);
    rtree.current_it_back = (Neighbor*) node2->findNeighbor(node1);
    assert(rtree.current_it_back);
    assert(rtree.current_it->length == rtree.current_it_back->length);

    // cout << "node pointed by rtree.current_it " << rtree.current_it->node->id + 1 << ", edge id " << rtree.current_it->id + 1 << ", direction " << rtree.current_it->direction<< endl;
    // cout << "node pointed by rtree.current_it_back " << rtree.current_it_back->node->id + 1 << ", edge id " << rtree.current_it_back->id + 1 << ", direction " << rtree.current_it->direction<< endl;

    double current_len = rtree.current_it->length;
    assert(current_len >= 0.0);
    /** minimum branch length for optimization, default 0.000001 */
    double min_branch_length = 0.000001;
    /** maximum branch length for optimization, default 100 */
    double max_branch_length = 100;
    double negative_lh = 0;
    double ferror, optx;

    // Brent method
    optx = minimizeOneDimen(rtree, min_branch_length, current_len, max_branch_length, min_branch_length, &negative_lh, &ferror);
    if (debug) {
        cout << "\tUsing Brent method to optimize the likelihood" << endl;
        cout << "\tmax Brent logl: " << -negative_lh << " optimized branch length " << optx << endl;
    }
    // cout << "\n rtree.current_it->length before assigning optx " << current_len << endl;
    rtree.current_it->length = optx;
    rtree.current_it_back->length = optx;
    rtree.edges[rtree.current_it->id].length = optx;
    // cout << "rtree.current_it->length after assigning optx " << rtree.current_it->length << endl;
}


// Find the approximate likelihood of an NNI-tree
NNIMove get_best_NNI_for_bran(evo_tree& rtree, Node *node1, Node *node2, NNIMove* nniMoves = NULL, bool nni5 = true) {
    int debug = 0;

    if(debug) cout << "Computing the approximate likelihood of an NNI-tree" << endl;

	assert(!node1->is_leaf() && !node2->is_leaf());
    assert(node1->degree() == 3 && node2->degree() == 3);

    if (((Neighbor*)node1->findNeighbor(node2))->direction == TOWARD_ROOT) {
        // swap node1 and node2 if the direction is not right, only for nonreversible models
        // cout << "swapping node " << node1->id + 1 << ", " << node2->id + 1 << endl;
        Node *tmp = node1;
        node1 = node2;
        node2 = tmp;
    }

	int IT_NUM = (nni5) ? 6 : 2;
	NeighborVec::iterator it;
	NeighborVec::iterator saved_it[6];
	int id = 0;
	saved_it[id++] = node1->findNeighborIt(node2);
	saved_it[id++] = node2->findNeighborIt(node1);

	if (nni5) {
		FOR_NEIGHBOR(node1, node2, it)
			saved_it[id++] = (*it)->node->findNeighborIt(node1);

		FOR_NEIGHBOR(node2, node1, it)
			saved_it[id++] = (*it)->node->findNeighborIt(node2);
	}
	assert(id == IT_NUM);

	Neighbor *saved_nei[6];
    int mem_id = 0;
	// save Neighbor and allocate new Neighbor pointer
	for (id = 0; id < IT_NUM; id++) {
		saved_nei[id] = (*saved_it[id]);
        *saved_it[id] = new Neighbor(saved_nei[id]->node, saved_nei[id]->length, saved_nei[id]->id);
        ((Neighbor*)*saved_it[id])->direction = ((Neighbor*)saved_nei[id])->direction;
	}

    // cout << "saved neighbors: " << endl;
    // for(int i = 0; i < IT_NUM; i++){
    //     Neighbor * nei = saved_nei[i];
    //     cout << nei->node->id + 1 << "\t" << nei->length << endl;
    // }

	// get the Neighbor again since it is replaced for saving purpose
	Neighbor* node12_it = (Neighbor*) node1->findNeighbor(node2);
	Neighbor* node21_it = (Neighbor*) node2->findNeighbor(node1);

    int cnt;
    bool newNNIMoves = false;
    if (!nniMoves) {
		// Initialize the 2 NNI moves
    	newNNIMoves = true;
    	nniMoves = new NNIMove[2];
    	nniMoves[0].node1 = NULL;
    }

    if (nniMoves[0].node1) {
    	// assuming that node1Nei_it and node2Nei_it is defined in nniMoves structure
    	for (cnt = 0; cnt < 2; cnt++) {
    		// sanity check
    		if (!node1->findNeighbor((*nniMoves[cnt].node1Nei_it)->node)) cout << "Error in neighbors" << endl;
    		if (!node2->findNeighbor((*nniMoves[cnt].node2Nei_it)->node)) cout << "Error in neighbors" << endl;
    	}
    } else {
        // cout << "creating NNI move" << endl;
        cnt = 0;
        FOR_NEIGHBOR_IT(node1, node2, node1_it)
        if (((Neighbor*)*node1_it)->direction != TOWARD_ROOT)
        {
			cnt = 0;
			FOR_NEIGHBOR_IT(node2, node1, node2_it) {
				//   Initialize the 2 NNI moves
				nniMoves[cnt].node1Nei_it = node1_it;
				nniMoves[cnt].node2Nei_it = node2_it;
				cnt++;
			}
			break;
        }
        assert(cnt == 2);
    }

    // Initialize node1 and node2 in nniMoves
	nniMoves[0].node1 = nniMoves[1].node1 = node1;
	nniMoves[0].node2 = nniMoves[1].node2 = node2;
    nniMoves[0].newloglh = nniMoves[1].newloglh = -DBL_MAX;

    if(debug){
        cout << " NNI move 1 " << nniMoves[0].node1->id + 1  << "\t" << nniMoves[0].node2->id + 1 << endl;
        cout << " NNI move 2 " << nniMoves[1].node1->id + 1  << "\t" << nniMoves[1].node2->id + 1 << endl;
    }


    for (cnt = 0; cnt < 2; cnt++)
    {
        // do the NNI swap
    	NeighborVec::iterator node1_it = nniMoves[cnt].node1Nei_it;
    	NeighborVec::iterator node2_it = nniMoves[cnt].node2Nei_it;
        Neighbor *node1_nei = *node1_it;
        Neighbor *node2_nei = *node2_it;

        if(debug){
            cout << " \ndoing NNI move" << cnt + 1 << endl;
            cout << " old neighbors of node " << node1->id + 1 << endl;
            FOR_NEIGHBOR(node1, node2, it)
                cout << "\t" << (*it)->id + 1;
            cout << endl;
            cout << " old neighbors of node " << node2->id + 1 << endl;
            FOR_NEIGHBOR(node2, node1, it)
                cout << "\t" << (*it)->id + 1;
            cout << endl;
        }

        node1->updateNeighbor(node1_it, node2_nei);
        node2_nei->node->updateNeighbor(node2, node1);
        node2->updateNeighbor(node2_it, node1_nei);
        node1_nei->node->updateNeighbor(node1, node2);

        if(debug){
            cout << " new neighbors of node " << node1->id + 1<< endl;
            for (it = node1->neighbors.begin(); it != node1->neighbors.end(); it++)
                cout << "\t" << (*it)->node->id + 1;
            cout << endl;
            cout << " new neighbors of node " << node2->id + 1 << endl;
            for (it = node2->neighbors.begin(); it != node2->neighbors.end(); it++)
                cout << "\t" << (*it)->node->id + 1;
            cout << endl;
        }

        int nni5_num_eval = 1;
        for (int step = 0; step < nni5_num_eval; step++) {
    		int i=1;
            if (nni5) {
                // cout << "\noptimizing neighbors of node 1 (excluding node2)" << endl;
    			FOR_NEIGHBOR(node1, node2, it)
    			{
    				optimize_one_branch(rtree, node1, (Node *) (*it)->node);
    				node1->findNeighbor((*it)->node)->getLength(nniMoves[cnt].newLen[i]);
    				i++;
    			}
            }

            // cout << "\noptimizing (node 1, node2)" << endl;
    		optimize_one_branch(rtree, node1, node2);
            node1->findNeighbor(node2)->getLength(nniMoves[cnt].newLen[0]);

            if (nni5) {
                // cout << "\noptimizing neighbors of node 2 (excluding node1)" << endl;
    			FOR_NEIGHBOR(node2, node1, it)
    			{
    				optimize_one_branch(rtree, node2, (Node *) (*it)->node);
    				node2->findNeighbor((*it)->node)->getLength(nniMoves[cnt].newLen[i]);
    				i++;
    			}
    		}
        }
		double score;
        if(model == 3){
            score = get_likelihood_decomp(rtree);
        }else{
            score = get_likelihood_revised(rtree);
        }

        if(debug){
            cout << "\nNNI " << node1->id << " - " << node2->id << ": " << score << endl;
        }

		nniMoves[cnt].newloglh = score;

        // else, swap back, also recover the branch lengths
		node1->updateNeighbor(node1_it, node1_nei);
		node1_nei->node->updateNeighbor(node2, node1);
		node2->updateNeighbor(node2_it, node2_nei);
		node2_nei->node->updateNeighbor(node1, node2);
		// ONLY FOR CHECKING WITH OLGA's PLEN MODEL
		//node1_nei->length = node2_nei->length = saved_len;
    }

	 // restore the Neighbor*
     if(debug)  cout << "\nrestore the neighbors" << endl;
	 for (id = IT_NUM-1; id >= 0; id--) {
		 if (*saved_it[id] == rtree.current_it) rtree.current_it = (Neighbor*) saved_nei[id];
		 if (*saved_it[id] == rtree.current_it_back) rtree.current_it_back = (Neighbor*) saved_nei[id];

		 delete (*saved_it[id]);
		 (*saved_it[id]) = saved_nei[id];
	 }
	 // restore the length of 4 branches around node1, node2
	 FOR_NEIGHBOR(node1, node2, it)
		 (*it)->setLength((*it)->node->findNeighbor(node1));
	 FOR_NEIGHBOR(node2, node1, it)
		 (*it)->setLength((*it)->node->findNeighbor(node2));

	 NNIMove res;
	 if (nniMoves[0].newloglh > nniMoves[1].newloglh) {
		 res = nniMoves[0];
	 } else {
		 res = nniMoves[1];
	 }
	if (newNNIMoves) {
		delete [] nniMoves;
	}
	return res;
}


// Find NNI increasing likelihood of current tree
void evaluate_NNIs(evo_tree &rtree, Branches &nniBranches, vector<NNIMove>  &positiveNNIs, double curScore) {
    int debug = 0;
    if(debug)
        cout << "\nSearching NNIs increasing likelihood of current tree with score " << curScore << endl;
    for (Branches::iterator it = nniBranches.begin(); it != nniBranches.end(); it++) {
        NNIMove nni = get_best_NNI_for_bran(rtree, (Node *) it->second.first, (Node *) it->second.second);
        if(debug) cout << "\n   one positive NNI found: " << nni.node1->id + 1 << ", " << nni.node2->id + 1 <<" with score " << nni.newloglh << endl;
        if (nni.newloglh > curScore) {
            positiveNNIs.push_back(nni);
        }
    }
    if(debug){
        cout << "\nThere are " << positiveNNIs.size() << " positive NNIs" << endl;
        for(int i = 0; i < positiveNNIs.size(); i++){
            cout << "\t" << positiveNNIs[i].node1->id << ", " << positiveNNIs[i].node2->id << endl;
        }
        cout << endl;
    }
}

// Assume lenvec is the size of all branches, using ID of end node to distinguish each edge
void save_branch_lengths(evo_tree &rtree, DoubleVector &lenvec, int startid = 0, Node *node = NULL, Node *dad = NULL) {
    int debug = 0;
    if (!node) {
        node = &(rtree.nodes[rtree.root_node_id]);   // root
        int branchNum = rtree.edges.size();
        if(debug){
            rtree.print();
            cout << "Root is " << node->id + 1 << " with " << branchNum << " edges "  << node->neighbors.size() << " neighbors " << endl;
        }
        if (lenvec.empty()) lenvec.resize(branchNum);
    }

    if(debug){
        cout << "saving neighboring branch lengths of node " << node->id + 1 << " to a vector" << endl;
        cout << "initial size of all branches:";
        for(int i = 0; i < lenvec.size(); i++){
            cout << "\t" << lenvec[i];
        }
        cout << endl;
        cout << " traversing neighbors of node " << node->id + 1 << endl;
    }

    FOR_NEIGHBOR_IT(node, dad, it){
        // cout << " neighbor " <<  (*it)->id + 1 << " with length " <<  (*it)->length << endl;
    	(*it)->getLength(lenvec, (*it)->id);
        // cout  << " length of neighbor " << (*it)->id + 1 << " is " << lenvec[(*it)->id] << endl;
    	save_branch_lengths(rtree, lenvec, startid, (Node *) (*it)->node, node);
    }

    if(debug){
        cout << "size of all branches:";
        for(int i = 0; i < lenvec.size(); i++){
            cout << "\t" << lenvec[i];
        }
        cout << endl;
    }
}


void restore_branch_lengths(evo_tree &rtree, DoubleVector &lenvec, int startid = 0, Node *node = NULL, Node *dad = NULL) {
    int debug = 0;
    if(debug) cout << "\nrestoring branch lengths from a vector" << endl;
    if (!node) {
        node = &(rtree.nodes[rtree.root_node_id]);   // root
        cout << "Root is " << node->id + 1 << endl;
        assert(lenvec.size() == rtree.edges.size());
    }
    FOR_NEIGHBOR_IT(node, dad, it){
    	(*it)->setLength(lenvec, (*it)->id + startid);
        (*it)->node->findNeighbor(node)->setLength(lenvec, (*it)->id + startid);
    	restore_branch_lengths(rtree, lenvec, startid, (Node *) (*it)->node, node);
    }

    if(debug){
        cout << "size of all branches after restoring:";
        for(int i = 0; i < lenvec.size(); i++){
            cout << "\t" << lenvec[i];
        }
        cout << endl;
    }
}


void compute_best_traversal(evo_tree& rtree, NodeVector &nodes, NodeVector &nodes2) {
    Node *farleaf = rtree.find_farthest_leaf();
//    Node *farleaf = root;

    // double call to farthest leaf to find the longest path on the tree
    rtree.find_farthest_leaf(farleaf);
    // cout << "Tree diameter: " << farleaf->height << endl;

    rtree.get_preorder_branches(nodes, nodes2, farleaf);
}

double optimize_all_branches(evo_tree& rtree, int my_iterations, double tolerance) {
    int debug = 0;
    if (debug)
        cout << "\nOptimizing branch lengths (max " << my_iterations << " loops)..." << endl;

    NodeVector nodes, nodes2;
    compute_best_traversal(rtree, nodes, nodes2);

    // cout << "nodes in best traversal: ";
    // for(auto n : nodes){
    //     cout << "\t" << n->id + 1;
    // }
    // cout << endl;
    // cout << "nodes2 in best traversal: ";
    // for(auto n : nodes2){
    //     cout << "\t" << n->id + 1;
    // }
    // cout << endl;

    double tree_lh = 0;
    if(model == 3){
        tree_lh = get_likelihood_decomp(rtree);
    }else{
        tree_lh = get_likelihood_revised(rtree);
    }
    if (debug) {
        cout << "Initial tree log-likelihood: " << tree_lh << endl;
    }

    DoubleVector lenvec;
    for (int i = 0; i < my_iterations; i++) {
    	save_branch_lengths(rtree, lenvec, 0);
        for (int j = 0; j < nodes.size(); j++) {
            optimize_one_branch(rtree, (Node *)nodes[j], (Node *)nodes2[j]);
            if (debug) {
                cout << " optimizing branch " << nodes[j]->id << " " << nodes2[j]->id << endl;
            }
        }
        double new_tree_lh;
        if(model == 3){
            new_tree_lh = get_likelihood_decomp(rtree);
        }else{
            new_tree_lh = get_likelihood_revised(rtree);
        }
        //cout<<"After opt  log-lh = "<<new_tree_lh<<endl;
        if (debug) {
            cout << "tree log likelihood after iteration " << i + 1 << " : ";
            cout << new_tree_lh << endl;
        }

        if (new_tree_lh < tree_lh - tolerance*0.1) {
        	// IN RARE CASE: tree log-likelihood decreases, revert the branch length and stop
        	if (debug)
                cout << "tree log-likelihood decreases" << endl;
        		cout << "NOTE: Restoring branch lengths as tree log-likelihood decreases after branch length optimization: "
        			<< tree_lh << " -> " << new_tree_lh << endl;

        	restore_branch_lengths(rtree, lenvec);

            double max_delta_lh = 1.0;

            if(model == 3){
                new_tree_lh = get_likelihood_decomp(rtree);
            }else{
                new_tree_lh = get_likelihood_revised(rtree);
            }

            if (fabs(new_tree_lh-tree_lh) > max_delta_lh) {
                cout << endl;
                cout << "new_tree_lh: " << new_tree_lh << "   tree_lh: " << tree_lh << endl;
            }
        	assert(fabs(new_tree_lh-tree_lh) < max_delta_lh);
        	return new_tree_lh;
        }

        // only return if the new_tree_lh >= tree_lh! (in rare case that likelihood decreases, continue the loop)
        if (tree_lh <= new_tree_lh && new_tree_lh <= tree_lh + tolerance) {
            if(debug) cout << "tree log-likelihood increases" << endl;
        	return new_tree_lh;
        }

        tree_lh = new_tree_lh;
    }

    if(debug) cout << "current score " << tree_lh << endl;

    return tree_lh;
}




// Apply hill climbing perturbation to obtain a locally optimal tree (by NNI)
// score used in this function is log likelihood, the larger the better
evo_tree do_hill_climbing_NNI(evo_tree& rtree, int speed_nni = 0, bool nni5 = true) {
    int debug = 0;

    double loglh_epsilon = 0.001;
    unsigned int totalNNIApplied = 0;
    unsigned int numSteps = 0;
    const int MAXSTEPS = Ns;

    Branches nniBranches;
    vector<NNIMove> positiveNNIs;
    vector<NNIMove> appliedNNIs;
    bool changeBran = true;

    evo_tree ntree(rtree);  // new tree after NNI
    double curScore;
    if(model == 3){
        curScore = get_likelihood_decomp(rtree);
    }else{
        curScore = get_likelihood_revised(rtree);
    }
    double curBestScore = curScore;
    double originalScore = curScore;

    // ntree.generate_neighbors();     // Get neighbor information for NNI applications
    for (numSteps = 1; numSteps <= MAXSTEPS; numSteps++) {
        if(debug){
            cout << "Step " << numSteps << endl;
            ntree.print();
        }

        double oldScore = curScore;
        // save all current branch lengths
        DoubleVector lenvec;
        save_branch_lengths(ntree, lenvec, 0);
        nniBranches.clear();

        // Do reduced NNI search.
        if (speed_nni) {
            if(debug) cout << "Doing reduced NNI search." << endl;
            // speedNNI option: only evaluate NNIs that are 2 branches away from the previously applied NNI
            Branches filteredNNIBranches;
            filter_NNI_branches(appliedNNIs, filteredNNIBranches);
            for (Branches::iterator it = filteredNNIBranches.begin(); it != filteredNNIBranches.end(); it++) {
                Branch curBranch = it->second;
                Neighbor* nei = (Neighbor*) curBranch.first->findNeighbor(curBranch.second);
                int branchID =  pairInteger(curBranch.first->id, curBranch.second->id);
                nniBranches.insert(pair<int, Branch>(branchID, curBranch));
            }
        } else {
            if(debug) cout << "Doing complete NNI search." << endl;
            get_NNI_branches(nniBranches, &ntree.nodes[rtree.nleaf], NULL);
        }

        if(debug){
            cout << " There are " <<  nniBranches.size() << " NNI branches" << endl;
            for(auto p: nniBranches){
                int branchID = p.first;
                Branch curBranch = p.second;
                cout << branchID << "\t" << curBranch.first->id + 1 << "," << curBranch.second->id + 1 << endl;
            }
        }

        // Only consider NNIs that increase the likelihood of current tree
        positiveNNIs.clear();
        evaluate_NNIs(rtree, nniBranches, positiveNNIs, curScore);
        if (positiveNNIs.size() == 0) {
            break;
        }

        /* sort all positive NNI moves (ASCENDING) */
        sort(positiveNNIs.begin(), positiveNNIs.end());
        /* remove conflicting NNIs */
        appliedNNIs.clear();
        if(debug) cout << "Getting compatible NNIs" << endl;
        get_compatible_NNIs(positiveNNIs, appliedNNIs);

        // do non-conflicting positive NNIs
        do_all_NNIs(appliedNNIs, changeBran, nni5);
        curScore = optimize_all_branches(ntree, 1, loglh_epsilon);
        double bfgs_score = 0;
        // max_likelihood_BFGS(ntree, bfgs_score, loglh_epsilon, 2000);
        // cout << "BRENT score " << curScore << ", BFGS score " << -bfgs_score << endl;

        if (curScore < appliedNNIs.at(0).newloglh - loglh_epsilon) {
            if(debug) cout << "Tree getting worse: curScore = " << curScore << " / best score = " <<  appliedNNIs.at(0).newloglh << endl;
            // tree cannot be worse if only 1 NNI is applied
            if (appliedNNIs.size() > 1) {
                // revert all applied NNIs
                do_all_NNIs(appliedNNIs, changeBran, nni5);
                restore_branch_lengths(ntree, lenvec);
                // only do the best NNI
                if(debug) cout << "only do the best NNI" << endl;
                appliedNNIs.resize(1);
                do_all_NNIs(appliedNNIs, changeBran, nni5);
                curScore = optimize_all_branches(ntree, 1, loglh_epsilon);
                // ntree = max_likelihood_BFGS(ntree, curScore, tolerance, miter);
                // cout << "current score " << curScore << ", best NNI score " << appliedNNIs.at(0).newloglh << endl;
                assert(curScore > appliedNNIs.at(0).newloglh - 0.1);
            } else{
                if(debug) cout << "Applied all NNIs successfully " << endl;
                //  cout << "current score " << curScore << ", best NNI score " << appliedNNIs.at(0).newloglh << endl;
                assert(curScore > appliedNNIs.at(0).newloglh - 0.1 && "Using one NNI reduces LogL");
            }
            totalNNIApplied++;
        } else {
            totalNNIApplied += appliedNNIs.size();
        }

        if(curScore < oldScore - loglh_epsilon){    // error
            cout << "$$$$$$$$: " << curScore << "\t" << oldScore << "\t" << curScore - oldScore << endl;

        }
        if (curScore - oldScore <  loglh_epsilon)   // no improvement
            break;
        if (curScore > curBestScore + 0.1) {
            curBestScore = curScore;
        }
    }

    if (totalNNIApplied == 0 && debug) {
        cout << "NOTE: Input tree is already NNI-optimal" << endl;
    }
    if (numSteps == MAXSTEPS) {
        cout << "WARNING: NNI search needs unusual large number of steps (" << numSteps << ") to converge!" << endl;
    }
    if(curScore < originalScore - loglh_epsilon){ // error
        cout << "AAAAAAAAAAAAAAAAAAA: " << curScore << "\t" << originalScore << "\t" << curScore - originalScore << endl;
    }

    // Get the new tree (already changed when applying NNIs)
    // ntree.score = curScore;
    if(debug){
        cout << "new tree after hill climbing NNI" << endl;
        ntree.print();
    }

    return ntree;
}


// Only feasible for trees with fewer than 12 samples
// Do maximization multiple times (determined by Ngen), since numerical optimizations are local hill-climbing algorithms and may converge to a local peak
evo_tree do_exhaustive_search(string real_tstring, int Ngen, const int init_tree, const string& dir_itrees, const int& max_static, const vector<double>& rates, const double ssize, const double tolerance, const int miter, const int optim, int Ne = 1, double beta = 0, double gtime=1){
    // initialize candidate tree set
    if(Ns > LARGE_TREE){
        cout << "\nFor data with larger than " << LARGE_TREE << " samples, it is very slow!" << endl;
        // return evo_tree();
    }
    // int max_tree_num = fact(Ns) * fact(Ns - 1) / exp2(Ns-1); // For trees
    int max_tree_num = num_trees[Ns - 1];
    cout << "\nMaximum number of possible trees to explore " << max_tree_num << endl;
    vector<evo_tree> init_trees = get_initial_trees(init_tree, dir_itrees, max_tree_num, rates, max_tree_num, cons, Ne, beta, gtime);
    vector<evo_tree> best_trees(init_trees);
    assert(max_tree_num == init_trees.size());
    vector<double> lnLs(max_tree_num,0);
    vector<int> index(max_tree_num);
    cout << "Initial number of trees " << max_tree_num << endl;
    if(debug) cout << "String for real tree is " << real_tstring << endl;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i < max_tree_num; ++i){
        string tstring = order_tree_string_uniq(create_tree_string_uniq(init_trees[i]));
        if(debug) cout << "String for tree " << i << " is " << tstring << endl;
        string tid = to_string(i);
        if(tstring == real_tstring){
            tid = tid + "(real)";
        }
        // if(debug) cout << "Maximization for tree " << i << endl;
        evo_tree best_tree;
        double Lf = MAX_NLNL;
        int count_static = 0;
        double min_lnL = MAX_NLNL;
        int count = 0;

        while(count < Ngen){
            // if(debug) cout << "Maximization in iteration " << count << endl;
            if(optim == 0){
                best_tree = max_likelihood(init_trees[i], model, Lf, ssize, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
            }else{
                // init_trees[i] has already been the same as best_tree
                best_tree = max_likelihood_BFGS(init_trees[i], Lf, tolerance, miter);
            }
            count++;
            // cout << Lf << " \t " << min_lnL << endl;
            if(Lf < min_lnL){   // Find a better tree
                min_lnL = Lf;
                lnLs[i] = Lf;
                best_trees[i] = best_tree;
                best_trees[i].score = Lf;
            }else{
                count_static++;
            }
            if(count_static >= max_static){
                if(debug) cout << "\tstatic likelihood " << lnLs[i] << " for tree " << i << endl;
                break;
            }
            // Change branch lengths randomly
            // gsl_vector* rblens = gsl_vector_alloc(init_trees[i].nedge);
            // for(int i=0; i<init_trees[i].nedge; ++i){
            //   gsl_vector_set(rblens, i, runiform(r, 1, age/2));
            // }
            // init_trees[i] = create_new_tree(rblens, init_trees[i], 0);
        }
        cout.precision(dbl::max_digits10);
        cout << "Score for tree " << tid << " is: " << lnLs[i] << endl;
    }

    // Output best tree in C
    cout << "The number of trees searched is " << searched_trees.size() << endl;
    // cout << "Output best tree so far" << endl;
    // evo_tree min_lnL_tree = trees[0];
    evo_tree min_lnL_tree = find_best_trees(best_trees, lnLs, index, 1)[0];

    // Weird that the score of a tree may be zero
    // min_lnL_tree.print();
    // cout << index[0] << endl;
    // for(int i=0; i<lnLs.size(); i++){
    //     cout << lnLs[i] << endl;
    //     cout << trees[i].score << endl;
    // }
    cout << "FINISHED. MIN -ve logL = " << lnLs[index[0]] << endl;
    cout << "The best tree reported is tree " << index[0] << endl;
    // ofstream out_tree("./example/searched_trees.txt");
    // for(int i=0; i < max_tree_num; ++i){
    // // for(auto it : searched_trees){
    //     // out_tree << it.first << endl;
    //     out_tree << order_tree_string_uniq(create_tree_string_uniq(trees[i])) << endl;
    // }
    // out_tree.close();

    return min_lnL_tree;
}


// Npop determines the maximum number of unique trees to try
evo_tree do_hill_climbing(const int Npop, const int Ngen, const int init_tree, const string& dir_itrees, const int& max_static, const vector<double>& rates, double ssize, double tolerance, int miter, int optim, int Ne = 1, double beta = 0, double gtime=1){
    int debug = 0;
    // initialize candidate tree set
    int max_tree_num = INT_MAX;
    // cout << "Maximum number of possible trees to explore " << max_tree_num << endl;
    vector<evo_tree> trees = get_initial_trees(init_tree, dir_itrees, Npop, rates, max_tree_num, cons, Ne, beta, gtime);
    int num2init = trees.size();
    vector<double> lnLs(num2init,0);
    vector<int> index(num2init);
    cout << "\tInitial number of trees " << num2init << endl;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i < num2init; ++i){
        evo_tree rtree;
        double Lf = 0;
        if(optim == 0){
            while(!(Lf>0)){
                Lf = 0;
                rtree = max_likelihood(trees[i], model, Lf, ssize, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
            }
        }
        if(optim == 1){
            while(!(Lf>0)){
                Lf = 0;
                // rtree = max_likelihood_BFGS(trees[i], model, Lf, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
                rtree = max_likelihood_BFGS(trees[i], Lf, tolerance, miter);
            }
        }
        rtree.score = Lf;
        trees[i] = rtree;
        lnLs[i] = Lf;
        // cout << "Score for tree " << i << " is " << Lf << endl;
    }

    // Select top MAX_TREE trees for hill climbing NNI to obtain locally optimal ML trees
    int num2perturb = (trees.size() < MAX_TREE1) ? trees.size() : MAX_TREE1;
    cout << "\tNumber of trees to perturb " << num2perturb << endl;
    cout << "\tNumber of maximum times to do perturbation " << Ngen << endl;
    cout << "\tNumber of times with no improvement before stopping perturbation " << max_static << endl;
    // int max_perturb = fact(Ns)/(2*fact(Ns-2));
    // cout << "Maximum number of reachable neighbors for one perturbation " << max_perturb << endl;
    vector<evo_tree> trees2 = find_best_trees(trees, lnLs, index, num2perturb);
    vector<double> lnLs2(num2perturb,0);
    vector<int> index2(num2perturb);
    // Perturb trees randomly
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i < num2perturb; ++i){
        // trees2[i] = get_local_optimal_tree(trees2[i], Ngen, max_perturb, max_static, ssize, tolerance, miter, optim, model, cons, maxj, cn_max, only_seg, correct_bias, is_total);
        trees2[i] = do_hill_climbing_NNI(trees2[i]);
        lnLs2[i] = trees2[i].score;
        // cout << "Score " << lnLs2[i] << endl;
    }

    // Keep top 5 trees for further optimization to escape from local optima
    int num2refine = (trees2.size() < MAX_TREE2) ? trees2.size(): MAX_TREE2;
    cout << "\tNumber of trees to refine " << num2refine << endl;
    vector<evo_tree> trees3 = find_best_trees(trees2, lnLs2, index2, num2refine);
    vector<double> lnLs3(num2refine, 0);
    vector<int> index3(num2refine);
    for(int i=0; i < num2refine; ++i){
        lnLs3[i] = trees3[i].score;
        // cout << lnLs3[i] << endl;
    }
    // int x=0;
    // iota( index3.begin(), index3.end(), x++);
    // sort( index3.begin(), index3.end(), [&](int i,int j){ return lnLs3[i] < lnLs3[j]; } );

    int count = 0;
    // Perturb trees randomly
    while(count < MAX_PERTURB){
        int i = gsl_rng_uniform_int(r, trees3.size());
        count += 1;
        if(debug) cout << "\t\tPerturb tree " << i << endl;
        // Randomly perturb a tree to allow escape from local optima, doing 0.5(n-2) random NNIs
        evo_tree ttree(trees3[i]);
        ttree = do_random_NNIs(ttree);
        // for(int i = 0; i < trees3[i].intedges.size()/2; i++){
        //     ttree = do_random_NNI(ttree, debug);
        // }
        // evo_tree ttree = perturb_tree(Ns, trees3[i]);
        // int c = 0;
        // int rnum_perturb = gsl_rng_uniform_int(r, max_perturb);
        // string tstring = order_tree_string_uniq(create_tree_string_uniq(ttree));
        // while( c < rnum_perturb && searched_trees.find(tstring) != searched_trees.end()){
        //     ttree = perturb_tree(Ns, ttree);
        //     tstring = order_tree_string_uniq(create_tree_string_uniq(ttree));
        //     c += 1;
        // }
        // if ( searched_trees.find(tstring) == searched_trees.end() ) {
        //     // mark the tree
        //     searched_trees[ tstring ] = 0;
        // }
        // if(cons){
        //     adjust_tree_blens(ttree);
        //     adjust_tree_height(ttree);
        //     adjust_tree_tips(ttree);
        // }
        // cout << "Optimized score before " << trees3[i].score << endl;
        // ttree.score = trees3[i].score;
        // evo_tree otree = get_local_optimal_tree(ttree, Ngen, max_perturb, max_static, ssize, tolerance, miter, optim, model, cons, maxj, cn_max, only_seg, correct_bias, is_total);
        evo_tree otree = do_hill_climbing_NNI(ttree);
        // cout << "Optimized score after " << otree.score << endl;
        evo_tree btree = find_best_trees(trees3, lnLs3, index3, 1)[0];
        assert(btree.score == lnLs3[index3[0]]);
        double max_lnl = lnLs3[index3[0]];
        double min_lnl = lnLs3[index3[index3.size()-1]];
        // for(int i=0; i < num2refine; ++i){
        //     cout << lnLs3[index3[i]] << endl;
        // }
        // cout << "Best score " << lnLs3[index3[0]] << endl;
        // cout << "Worst score " << lnLs3[index3.size()-1] << endl;
        if(otree.score < max_lnl){  // better than best tree in C
            // cout << "Replace worst tree" << endl;
            trees3[index3[index3.size()-1]] = otree;
            lnLs3[index3[index3.size()-1]] = otree.score;
            count = 0;
        }
        else{
            // Replace the worst tree
            if(otree.score < min_lnl){ // better than worst tree in C
                // cout << "Replace worst tree" << endl;
                trees3[index3[index3.size()-1]] = otree;
                lnLs3[index3[index3.size()-1]] = otree.score;
            }
        }
    }

    // Output best tree in C
    cout << "\tThe number of trees searched is " << searched_trees.size() << endl;
    // cout << "Output best tree so far" << endl;
    evo_tree min_lnL_tree = find_best_trees(trees3, lnLs3, index3, 1)[0];
    double min_lnL = min_lnL_tree.score;
    // min_lnL_tree.print();
    cout << "\nFINISHED. MIN -ve logL = " << min_lnL << endl;

    // print out searched trees
    if(debug){
        ofstream out_tree("./searched_trees.txt");
        for(auto it : searched_trees){
            out_tree << it.first << endl;
        }
        out_tree.close();
    }

    return min_lnL_tree;
}



evo_tree do_evolutionary_algorithm(const int& Npop, const int& Ngen, const int init_tree, const string& dir_itrees, const int& max_static, const vector<double>& rates, const double ssize, const double tolerance, const int miter, const int optim, const int model, const int cons, const int maxj, const int cn_max, const int only_seg, const int correct_bias, int is_total=1, int Ne = 1, double beta = 0, double gtime=1){
  //cout << "Running evolutionary algorithm" << endl;
  // create initial population of trees. Sample from coalescent trees
  vector<evo_tree> trees = get_initial_trees(init_tree, dir_itrees, Npop, rates, Npop, cons, Ne, beta, gtime);
  vector<double> lnLs(2 * Npop,0);
  double min_lnL = MAX_NLNL;
  evo_tree min_lnL_tree;
  double Lf = 0;
  int count_static = 0;

  for(int g=0; g<Ngen; ++g){
    // Growth stage: create Npop copies + Npop copies with mutation (topolgy changes)
    // The top scoring trees have already been scored and optimised
    vector<evo_tree> new_trees;
    vector<evo_tree> opt_trees;

    if( g == 0 ){
      for(int i=0; i<Npop; ++i){
	         new_trees.push_back( trees[i] );
      }
      for(int i=0; i<Npop; ++i){
          evo_tree ntree = perturb_tree_set(Ns, trees);
          ntree.score = - get_likelihood_revised(ntree);
    	  //new_trees.push_back( perturb_tree(Ns, Nchar, trees[i]) );
    	  new_trees.push_back(ntree);
      }

      // Selection: score all the trees in new_trees (all got maximized)
      for(int i=0; i<2*Npop; ++i){
        evo_tree otree;
        if(optim == 0){
          otree = max_likelihood(new_trees[i], model, Lf, ssize, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
        }
        if(optim == 1){
          // new_trees[i].print();
          // otree = max_likelihood_BFGS(new_trees[i], model, Lf, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
          otree = max_likelihood_BFGS(new_trees[i], Lf, tolerance, miter);
        }
        otree.score = Lf;
        // cout << "otree tobs " << otree.tobs[0] << endl;
        lnLs[i] = Lf;

        string tstring = order_tree_string_uniq(create_tree_string_uniq( otree ) );
        searched_trees[ tstring ] += 1;

        opt_trees.push_back( otree );
        //cout << "g/i/lnL:\t" << g << "\t" << i << "\t" << lnLs[i] << endl;
      }
    }else{
          // Leave this subpopulation unchanged
          for(int i=0; i<Npop; ++i){
        	new_trees.push_back( trees[i] );
        	opt_trees.push_back( trees[i] );
        	lnLs[i] = new_trees[i].score;
          }

          // cout << "Size of new_trees " << new_trees.size() << endl;
          // cout << "Size of opt_trees " << new_trees.size() << endl;

          // Perturb this subpopulation
          for(int i=0; i<Npop; ++i){
            evo_tree ntree = perturb_tree_set(Ns, trees);
            ntree.score = - get_likelihood_revised(ntree);
      	    //new_trees.push_back( perturb_tree(Ns, Nchar, trees[i]) );
      	    new_trees.push_back(ntree);
            // new_trees of size 2 Npop
            evo_tree otree;
            if(optim == 0){
        	    otree = max_likelihood(new_trees[Npop + i], model, Lf, ssize, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
            }
            if(optim == 1){
                // new_trees[Npop + i].print();
                // otree = max_likelihood_BFGS(new_trees[Npop + i], model, Lf, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
                otree = max_likelihood_BFGS(new_trees[Npop + i], Lf, tolerance, miter);
            }
        	otree.score = Lf;
            // cout << "otree tobs " << otree.tobs[0] << endl;
        	lnLs[Npop + i] = Lf;

            string tstring = order_tree_string_uniq(create_tree_string_uniq( otree ) );
            searched_trees[ tstring ] += 1;

        	opt_trees.push_back( otree );
        	//cout << "g/i/lnL:\t" << g << "\t" << Npop+i << "\t" << lnLs[i] << endl;
        }
      // Randomly sample again?
    }

    vector<int> index(2*Npop);
    int x=0;
    iota(index.begin(), index.end(), x++);
    sort(index.begin(), index.end(), [&](int i,int j){ return lnLs[i]<lnLs[j];});

    // Selection: calculate mean fitness of top half of population
    double meand = 0;
    for(int k=0; k<Npop; ++k){
      //cout << "\t\t" << lnLs[ index[k] ] << endl;
      meand += lnLs[ index[k] ];
    }
    meand = meand/Npop;
    if( g%1 == 0) cout << "g / av dist / top dist / trees searched \t" << g << "\t" << meand << "\t" << lnLs[ index[0] ] << "\t" << searched_trees.size() << endl;

    // Selection: select top half
    for(int i=0; i<Npop; ++i){
      trees[i] = opt_trees[ index[i] ];
    }

    // Selection: record the best (lowest) scoring tree
    if( lnLs[ index[0] ] < min_lnL ){
      min_lnL = lnLs[ index[0] ];
      min_lnL_tree = opt_trees[ index[0] ];
      count_static = 0;
      min_lnL_tree.print();
    }else{
      cout << "min static" << endl;
      count_static += 1;
    }

    if( max_static > 0 && count_static == max_static ){
      cout << "\t### static likelihood function. Finishing on ngen = " << g << endl;
      break;
    }

    fill(lnLs.begin(), lnLs.end(), 0);

    // A tree has been maximized at least five times
    int sum_max_num = 0;
    for(auto it : searched_trees){
        sum_max_num += it.second;
    }
    // cout << "Total times of maximization " << sum_max_num << endl;

    if( exhausted_tree_search == true && sum_max_num > MAX_OPT * searched_trees.size()){
      cout << "\tperturb_tree struggling to find new topologies. Either exhausted possible trees or local minimum" << endl;
      break;
    }
  }

  cout << "FINISHED. MIN -ve logL = " << min_lnL << endl;

  // print out searched trees
  // ofstream out_tree("./test/searched_trees.txt");
  // for(auto it : searched_trees){
  //     out_tree << it.first << endl;
  // }
  // out_tree.close();

  return min_lnL_tree;
}

/*
vector<edge> create_edges_from_nodes( const vector<node>& nodes, const vector<double>& node_times ){
  vector<edge> enew;

  // Add leaf and internal nodes
  int root_id = 0;
  int id = 0;
  for(int i=0; i<nodes.size(); ++i){
    if(nodes[i].parent == -1) root_id = nodes[i].id;
  }

  for(int i=0; i<nodes.size(); ++i){
    if(nodes[i].id != root_id && nodes[i].id != root_id-1 && nodes[i].parent != root_id){
      enew.push_back( edge(id, nodes[i].parent, nodes[i].id, node_times[nodes[i].id] - node_times[nodes[i].parent]) );
      id++;
    }
  }
  // Add root nodes
  for(int i=0; i<nodes.size(); ++i){
  if( nodes[i].parent == root_id && nodes[i].id != root_id-1){
    enew.push_back( edge(id, nodes[i].parent, nodes[i].id, node_times[nodes[i].id] - node_times[nodes[i].parent]) );
    id++;
  }
  }
  enew.push_back( edge(id, root_id, Ns, 0) );

  return enew;
}
*/


evo_tree init_tree_from_file(const string& tree_file, int Ns, int Nchar, int model, int only_seg, vector<double> tobs, vector<double> rates){
    evo_tree test_tree = read_tree_info(tree_file, Ns);
    // test_tree.print();
    if(model==0){
        test_tree.mu = 1.0/Nchar;
    }
    else{
        test_tree.dup_rate = rates[0];
        test_tree.del_rate = rates[1];
        if(!only_seg){
            test_tree.chr_gain_rate = rates[2];
            test_tree.chr_loss_rate = rates[3];
            test_tree.wgd_rate = rates[4];
        }
    }
    test_tree.tobs = tobs;

    return test_tree;
}


// Get the mutation rates in year along each branch
// num_total_bins: number of total segments in the (haploid) genome
// unit of mutation rate: duplication/deletion rate -- bin/allele/year, chromosome gain/loss rate -- chromosome/year, WGD date -- year
vector<int> compute_mutation_rates(evo_tree tree, int only_seg, int num_total_bins){
    double mu_est;  // mutation rate per year
    vector<double> mu_all;
    vector<int> nmuts;

    if(!only_seg){
        mu_est = NORM_PLOIDY * NUM_CHR * (tree.chr_gain_rate + tree.chr_loss_rate) + NORM_PLOIDY * num_total_bins * (tree.dup_rate + tree.del_rate) + tree.wgd_rate;
    }
    else{
        mu_est = NORM_PLOIDY * num_total_bins * (tree.dup_rate + tree.del_rate);
    }

    mu_all.clear();
    for(int i = 0; i<tree.lengths.size(); i++){
        mu_all.push_back(mu_est);
    }

    nmuts = tree.get_nmuts(mu_all);

    return nmuts;
}

void check_single_blen(int model){
    // Check how different mutation rates affect the rates of reaching different states
    // For model 2, only first two rates matter
    // vector<vector<double>> rate_range{
    //     {0.0001, 0.0001, 0.0001, 0.0001, 0.0001},
    //     // {0.0001, 0.0001, 0.001, 0.001, 0.01},
    //     {0.001, 0.001, 0.001, 0.001, 0.001},
    //     {0.01, 0.01, 0.01, 0.01, 0.01},
    //     {0.1, 0.1, 0.1, 0.1, 0.1},
    //     {1, 1, 1, 1, 1}};

   // For model 3
    vector<vector<double>> rate_range{
        {0.001, 0.001, 0, 0, 0.001},
        {0.001, 0.001, 0, 0, 0.01},
        {0.001, 0.001, 0, 0, 0.1},
        {0.01, 0.01, 0, 0, 0.001},
        {0.01, 0.01, 0, 0, 0.01},
        {0.01, 0.01, 0, 0, 0.1}};

    vector<vector<int>> change_range{
        {1, 0, 2},
        {1, 0, 4},
        {1, 0, 8},
        {1, 0, 2},
        {1, 0, 4},
        {1, 0, 8}};
    // failed when adding branch length 50 for model 3
    vector<double> blen_range{1, 10, 20, 30, 40};
    // vector<double> blen_range{1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    // vector<int> cn_range{4, 6, 8, 10, 12};
    vector<int> cn_range{4, 6, 8};
    int num_seg = 100;
    int m_max = 1;
    int max_wgd = 2;
    int max_chr_change = 2;
    int max_site_change = 4;
    for (int c = 0; c < cn_range.size(); c++){
        int cn_max = cn_range[c];
        string fname = "chain_state_cn" + to_string(cn_max) + ".txt";
        ofstream fout(fname);
        for(auto rates : rate_range){
            // cout << "Mutation rates: " << rate_range[i] << endl;
            for(int j = 0; j < blen_range.size(); j++){
                double blen = blen_range[j];
                cout << "\tBranch lengths: " << blen << endl;
                if(model == 2){
                    check_chain_by_branch_m2(num_seg, rates, blen, cn_max, fout);
                }else{
                    int is_total = 1;
                    vector<int> num_changes = change_range[j];
                    int max_wgd = num_changes[0];
                    int max_chr_change = num_changes[1];
                    int max_site_change = num_changes[2];
                    check_chain_by_branch_m3(is_total, num_seg, rates, blen, cn_max, m_max, max_wgd, max_chr_change, max_site_change, fout);
                }
            }
        }
    }
}



void check_NNI(const string& tree_file){
    // original tree
    tobs = vector<double>(Ns, 0);
    vector<double> rates(5, 0);
    int model = 2;
    int only_seg = 1;
    int Ns = 5;
    int debug = 1;

    evo_tree real_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    real_tree.print();

    // tree after NNI
    // do_random_NNI(real_tree, debug);

    // do_NNI_on_branch(real_tree, 0, debug);
    // do_NNI_on_branch(real_tree, 1, debug);
    // do_NNI_on_branch(real_tree, 2, debug);

    // vector<evo_tree> nni_trees_all = get_all_NNI_trees(real_tree, debug);
    // cout << "There are " << nni_trees_all.size() << " NNI trees" << endl;
    vector<NNIMove> moves;
    map<int, vector<evo_tree>> nni_trees_all = get_all_NNI_trees(real_tree, moves, debug);
    cout << "There are " << moves.size() << " NNI moves" << endl;

    real_tree.print();
}


// Run the program on a given tree with different modes of estimation (branch length constrained or not, mutation rate estimated or not)
void run_test(const string& tree_file, int Ns, int num_total_bins, int Nchar, int model, int cn_max, int only_seg, int correct_bias, int is_total, const vector<double>& tobs, const vector<vector<int>>& vobs0, int Nchar0, const vector<double>& rates, double ssize, double tolerance, double miter){
    // MLE testing
    //static const int arr1[] = {8,5, 8,1, 9,2, 9,3, 10,9, 10,8, 11,4, 11,10, 7,11, 7,6 };
    //vector<int> e (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
    //for(int i=0; i<e.size();++i) e[i] = e[i] - 1;
    //static const double arr2[] = {18.49, 38.49, 51.71, 31.71, 0.51, 3.73, 22.2, 0.013, 0.99, 0};
    //vector<double> l (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

    // read in true tree
    evo_tree test_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    double Ls = 0.0;
    double Lf = 0;  // Used in maximization
    string fname = "";
    ofstream out_tree;
    evo_tree min_tree;
    double mu_est;
    vector<double> mu_all;
    vector<int> nmuts;

    cout << "Computing likelihood of the given tree" << endl;
    Ls = get_likelihood(Ns, Nchar0, vobs0, test_tree, model, 0, cn_max);
    cout << "\nOriginal tree -ve likelihood with original method: " << -Ls << endl;

    Ls = get_likelihood_revised(test_tree);
    cout << "\nOriginal tree -ve likelihood with revised method: " << -Ls << endl;

    cout << "\n\n### Running optimisation: branches free, mu fixed" << endl;
    // min_tree = max_likelihood(test_tree, model, Lf, ssize, tolerance, miter, 0, 0, cn_max, only_seg, correct_bias);
    // cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree.mu <<  endl;
    // min_tree.print();
    // fname = tree_file + ".opt00_gsl.txt";
    // out_tree.open(fname);
    // min_tree.write(out_tree);
    // out_tree.close();
    // The original tree gets changed during maximization
    test_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    Lf = 0;
    cons = 0;
    maxj = 0;
    // min_tree = max_likelihood_BFGS(test_tree, model, Lf, tolerance, miter, 0, 0, cn_max, only_seg, correct_bias);
    min_tree = max_likelihood_BFGS(test_tree, Lf, tolerance, miter);
    cout << "\nMinimised tree likelihood / mu by BFGS : " << Lf << "\t" << min_tree.dup_rate << "\t" << min_tree.del_rate;
    if(!only_seg){
        cout << "\t" << min_tree.chr_gain_rate << "\t" << min_tree.chr_loss_rate << "\t" << min_tree.wgd_rate;
    }
    cout << endl;
    min_tree.print();
    fname = tree_file + ".opt00_bfgs.txt";
    out_tree.open(fname);
    min_tree.write_with_mut(out_tree, compute_mutation_rates(min_tree, only_seg, num_total_bins));
    out_tree.close();


    cout << "\n\n### Running optimisation: branches free, mu free" << endl;
    // test_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    // Lf = 0;
    // min_tree = max_likelihood(test_tree, model, Lf, ssize, tolerance, miter, 0, 1, cn_max, only_seg, correct_bias);
    // cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree.mu << endl;
    // min_tree.print();
    // fname = tree_file + ".opt01_gsl.txt";
    // out_tree.open(fname);
    // min_tree.write(out_tree);
    // out_tree.close();
    test_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    Lf = 0;
    cons = 0;
    maxj = 1;
    // min_tree = max_likelihood_BFGS(test_tree, model, Lf, tolerance, miter, 0, 1, cn_max, only_seg, correct_bias);
    min_tree = max_likelihood_BFGS(test_tree, Lf, tolerance, miter);
    cout << "\nMinimised tree likelihood / mu by BFGS : " << Lf << "\t" << min_tree.dup_rate << "\t" << min_tree.del_rate;
    if(!only_seg){
        cout << "\t" << min_tree.chr_gain_rate << "\t" << min_tree.chr_loss_rate << "\t" << min_tree.wgd_rate;
    }
    cout << endl;
    min_tree.print();
    fname = tree_file + ".opt01_bfgs.txt";
    out_tree.open(fname);
    // min_tree.write(out_tree
    min_tree.write_with_mut(out_tree, compute_mutation_rates(min_tree, only_seg, num_total_bins));
    out_tree.close();


    cout << "\n\n### Running optimisation: branches constrained, mu fixed" << endl;
    // test_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    // Lf = 0;
    // min_tree = max_likelihood(test_tree, model, Lf, ssize, tolerance, miter, 1, 0, cn_max, only_seg, correct_bias);
    // cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree.mu <<endl;
    // min_tree.print();
    // fname = tree_file + ".opt10_gsl.txt";
    // out_tree.open(fname);
    // min_tree.write(out_tree);
    // out_tree.close();
    test_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    Lf = 0;
    cons = 1;
    maxj = 0;
    min_tree = max_likelihood_BFGS(test_tree, Lf, tolerance, miter);
    // min_tree = max_likelihood_BFGS(test_tree, model, Lf, tolerance, miter, 1, 0, cn_max, only_seg, correct_bias);
    cout << "\nMinimised tree likelihood / mu by BFGS: " << Lf << "\t" << min_tree.dup_rate << "\t" << min_tree.del_rate;
    if(!only_seg){
        cout << "\t" << min_tree.chr_gain_rate << "\t" << min_tree.chr_loss_rate << "\t" << min_tree.wgd_rate;
    }
    cout << endl;
    min_tree.print();
    fname = tree_file + ".opt10_bfgs.txt";
    out_tree.open(fname);
    // min_tree.write(out_tree);
    min_tree.write_with_mut(out_tree, compute_mutation_rates(min_tree, only_seg, num_total_bins));
    out_tree.close();


    cout << "\n\n### Running optimisation: branches constrained, mu free" << endl;
    // test_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    // Lf = 0;
    // min_tree = max_likelihood(test_tree, model, Lf, ssize, tolerance, miter, 1, 1, cn_max, only_seg, correct_bias);
    // cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree.mu*Nchar <<endl;
    // min_tree.print();
    // fname = tree_file + ".opt11_gsl.txt";
    // out_tree.open(fname);
    // min_tree.write(out_tree);
    // out_tree.close();
    test_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    Lf = 0;
    cons = 1;
    maxj = 1;
    min_tree = max_likelihood_BFGS(test_tree, Lf, tolerance, miter);
    // min_tree = max_likelihood_BFGS(test_tree, model, Lf, tolerance, miter, 1, 1, cn_max, only_seg, correct_bias);
    cout << "\nMinimised tree likelihood / mu by BFGS : " << Lf << "\t" << min_tree.dup_rate << "\t" << min_tree.del_rate;
    if(!only_seg){
        cout << "\t" << min_tree.chr_gain_rate << "\t" << min_tree.chr_loss_rate << "\t" << min_tree.wgd_rate;
    }
    cout << endl;
    min_tree.print();
    fname = tree_file + ".opt11_bfgs.txt";
    out_tree.open(fname);
    // min_tree.write(out_tree);
    min_tree.write_with_mut(out_tree, compute_mutation_rates(min_tree, only_seg, num_total_bins));
    out_tree.close();

}


// Compute the likelihood of a tree given the observed copy number profile
double compute_tree_likelihood(const string& tree_file, map<int, set<vector<int>>>& decomp_table, int Ns, int Nchar, int num_invar_bins, map<int, vector<vector<int>>>& vobs, vector<double>& tobs, vector<double>& rates, int model, int cons, int cn_max, int max_wgd, int max_chr_change, int max_site_change, int only_seg, int correct_bias, int is_total=1){
    evo_tree tree;
    // string format = tree_file.substr(tree_file.find_last_of(".") + 1);
    // cout << "Format of the input tree file is " << format << endl;

    // if (format == "txt"){
    tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
    // }
    // if (format == "nex"){
    //     tree = read_nexus(tree_file);
    // }

    cout << "Node times: ";
    for(int i = 0; i<tree.node_times.size(); i++){
        cout << "\t" << tree.node_times[i];
    }
    cout<< endl;

    int sample1 = 0;
    tree.get_ntime_interval(Ns, sample1);
    cout << "Top " << Ns << " time intervals: ";
    for(int i = 0; i < tree.top_tinvls.size(); i++){
        cout << "\t" << tree.top_tinvls[i];
    }
    cout<< endl;
    cout << "Corresponding " << Ns + 1 << " nodes: ";
    for(int i = 0; i < tree.top_tnodes.size(); i++){
        cout << tree.top_tnodes[i] + 1 << "\t" << "\t" << tree.node_times[tree.top_tnodes[i]] << endl;
    }

    vector<double> blens = tree.get_edges_from_interval(tree.top_tinvls, tree.top_tnodes);
    cout << "Branch lengths computed from intervals: ";
    for(int i = 0; i<blens.size(); i++){
        cout << i + 1 << "\t" << "\t" << blens[i] << endl;
    }

    double Ls = 0;
    if(only_seg){
        // Ls = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, tree, model, cons, cn_max, only_seg, correct_bias, is_total);
        Ls = get_likelihood_revised(tree);
    }else{
        cout << "\nComputing the likelihood based on independent Markov chain model " << endl;
        // Ls = get_likelihood_decomp(Ns, Nchar, num_invar_bins, vobs, tree, decomp_table, cons, cn_max, m_max, max_wgd, max_chr_change, max_site_change, correct_bias, is_total);
        Ls = get_likelihood_decomp(tree);
    }


    return Ls;
}

void maximize_tree_likelihood(const string& tree_file, const string& ofile, int Ns, int Nchar, vector<double>& tobs, vector<double>& rates, double ssize, double tolerance, int miter, int model, int optim, int cons, int maxj, int cn_max, int only_seg, int correct_bias, int is_total){
    evo_tree tree;
    // string format = tree_file.substr(tree_file.find_last_of(".") + 1);
    // cout << "Format of the input tree file is " << format << endl;

    // if (format == "txt"){
        tree = read_tree_info(tree_file, Ns);
        tree.print();
        tree.tobs = tobs;
    // }
    // if (format == "nex"){
    //     tree = read_nexus(tree_file);
    // }
    if(cons){
        adjust_tree_blens(tree);
        adjust_tree_height(tree);
        adjust_tree_tips(tree);
    }

    if(model==0){
        tree.mu = 1.0/Nchar;
    }
    else{
        tree.dup_rate = rates[0];
        tree.del_rate = rates[1];
        if(!only_seg){
            assert(rates.size() > 2);
            tree.chr_gain_rate = rates[2];
            tree.chr_loss_rate = rates[3];
            tree.wgd_rate = rates[4];
        }
    }

    double Lf = 0;
    evo_tree min_tree;
    if(optim == 1){
        // min_tree = max_likelihood_BFGS(tree, model, Lf, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
        min_tree = max_likelihood_BFGS(tree, Lf, tolerance, miter);
    }
    if(optim == 0){
        min_tree = max_likelihood(tree, model, Lf, ssize, tolerance, miter, cons, maxj, cn_max, only_seg, correct_bias, is_total);
    }
    cout << "\nMinimised tree likelihood: " << Lf << endl;
    if(maxj == 1){
        cout << "Estimated mutation rates: " << endl;
        cout << "\t" << min_tree.dup_rate << "\t" << min_tree.del_rate;
        if(!only_seg){
            cout << "\t" << min_tree.chr_gain_rate << "\t" << min_tree.chr_loss_rate << "\t" << min_tree.wgd_rate;
        }
        cout << endl;
    }
    min_tree.print();

    ofstream out_tree(ofile);
    min_tree.write(out_tree);
    out_tree.close();
}


void find_ML_tree(string real_tstring, int total_chr, int num_total_bins, string ofile, int tree_search, int Npop, int Ngen, int init_tree, string dir_itrees, int max_static, double ssize, double tolerance, int miter, int optim, const vector<double>& rates, int model, int cons, int maxj, int cn_max, int only_seg, int correct_bias, int is_total=1, int Ne = 1, double beta = 0, double gtime=1){
    evo_tree min_lnL_tree;
    if(tree_search == 0){
        cout << "\nSearching tree space with evolutionary algorithm" << endl;
        min_lnL_tree = do_evolutionary_algorithm(Npop, Ngen, init_tree, dir_itrees, max_static, rates, ssize, tolerance, miter, optim, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, Ne, beta, gtime);
    }
    else if(tree_search == 1){
        cout << "\nSearching tree space with hill climbing algorithm" << endl;
        min_lnL_tree = do_hill_climbing(Npop, Ngen, init_tree, dir_itrees, max_static, rates, ssize, tolerance, miter, optim, Ne, beta, gtime);
    }
    else{
        cout << "\nSearching tree space exhaustively" << endl;
        // cout << "Parameters: " << Ngen << "\t" << Ns << "\t" << Nchar << "\t" << num_invar_bins << "\t" << model << "\t" << cons << "\t" << cn_max << "\t" << only_seg << "\t" << correct_bias << "\t" << is_total << endl;
        min_lnL_tree = do_exhaustive_search(real_tstring, Ngen, init_tree, dir_itrees, max_static, rates, ssize, tolerance, miter, optim, Ne, beta, gtime);
    }

    if(maxj==1){
        if(model==0){
            cout << "Estimated mutation rate:  " << min_lnL_tree.mu << endl;
        }
        else{
            cout << "Estimated duplication rate:  " << min_lnL_tree.dup_rate << endl;
            cout << "Estimated deletion rate:  " << min_lnL_tree.del_rate << endl;
            if(!only_seg){
                cout << "Estimated chromosome gain rate (year):  " << min_lnL_tree.chr_gain_rate << endl;
                cout << "Estimated chromosome loss rate (year):  " << min_lnL_tree.chr_loss_rate << endl;
                cout << "Estimated whole genome doubling rate (year):  " << min_lnL_tree.wgd_rate << endl;
            }
        }
    }
    // Write out the top tree
    //cout << "Best fitting tree, -ve lnL = " << global_min << endl;
    cout.precision(PRINT_PRECISION);
    min_lnL_tree.print();

    // Recompute the likelihood to check it is not affected by tree adjustment
    double nlnL = 0;
    if(model == 3){
        // nlnL = -get_likelihood_decomp(Ns, Nchar, num_invar_bins, vobs, min_lnL_tree, decomp_table, cons, cn_max, m_max, max_wgd, max_chr_change, max_site_change, correct_bias, is_total);
        nlnL = -get_likelihood_decomp(min_lnL_tree);
    }else{
        nlnL = -get_likelihood_revised(min_lnL_tree);
        // nlnL = -get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, min_lnL_tree, model, cons, cn_max, only_seg, correct_bias, is_total);
    }
    cout.precision(dbl::max_digits10);
    cout << "Recomputed negative log likelihood " << nlnL << endl;

    // Check the validity of the tree
    if(cons==1){
        if(!is_tree_valid(min_lnL_tree, cons)){
            cout << "The final tree is not valid!" << endl;
        }
    }

    double mu_est;
    if(!only_seg){
        mu_est = NORM_PLOIDY * total_chr * (min_lnL_tree.chr_gain_rate + min_lnL_tree.chr_loss_rate) + NORM_PLOIDY * num_total_bins * (min_lnL_tree.dup_rate + min_lnL_tree.del_rate) + min_lnL_tree.wgd_rate;
    }
    else{
        mu_est = NORM_PLOIDY * num_total_bins * (min_lnL_tree.dup_rate + min_lnL_tree.del_rate);
    }

    cout << "Estimated total mutation rate per year " << mu_est << endl;
    vector<double> mu_all;
    for(int i = 0; i<min_lnL_tree.lengths.size(); i++){
        mu_all.push_back(mu_est);
    }
    vector<int> nmuts = min_lnL_tree.get_nmuts(mu_all);

    ofstream out_tree(ofile);
    // min_lnL_tree.write(out_tree);
    min_lnL_tree.write_with_mut(out_tree, nmuts);
    out_tree.close();

    string ofile_nex = ofile + ".nex";
    ofstream nex_tree(ofile_nex);
    int precision = 5;
    string newick = min_lnL_tree.make_newick(precision);
    min_lnL_tree.write_nexus(newick, nex_tree);
    nex_tree.close();

    ofile_nex =  ofile + ".nmut.nex";
    ofstream nex_tree2(ofile_nex);
    precision = 0;
    newick = min_lnL_tree.make_newick_nmut(precision, nmuts);
    min_lnL_tree.write_nexus(newick, nex_tree2);
    nex_tree2.close();
}



int main (int argc, char ** const argv) {
    int Npop, Ngen, max_static, miter, bootstrap, seed, optim, mode, tree_search, init_tree;
    double tolerance, ssize, mu, dup_rate, del_rate, chr_gain_rate, chr_loss_rate, wgd_rate, max_rate, beta, gtime;
    string datafile, timefile, ofile, tree_file, dir_itrees;
    int is_bin, incl_all, Ne;
    int infer_marginal_state, infer_joint_state, speed_nni;
    double min_asr;

    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ;
    po::options_description required("Required parameters");
    required.add_options()
     ("cfile,c", po::value<string>(&datafile)->required(), "input copy number profile file")
     ;
    po::options_description optional("Optional parameters");
    optional.add_options()
     ("tfile,t", po::value<string>(&timefile)->required(), "input time information file")

    ("nsample,s", po::value<int>(&Ns)->default_value(5), "number of samples or regions")
    ("ofile,o", po::value<string>(&ofile)->default_value("maxL-tree.txt"), "output tree file with maximum likelihood")
    // ("ofile_nex,o", po::value<string>(&ofile_nex)->default_value("maxL-tree.nex"), "output tree file with maximum likelihood in NEXUS format")

    ("tree_file", po::value<string>(&tree_file)->default_value(""), "input tree file ")

    ("cn_max", po::value<int>(&cn_max)->default_value(4), "maximum copy number of a segment")
    ("is_total", po::value<int>(&is_total)->default_value(1), "whether or not the input is total copy number")
    ("is_bin", po::value<int>(&is_bin)->default_value(1), "whether or not the input copy number is for each bin. If not, the input copy number is read as it is. Or else, consecutive bins will be merged")
    ("incl_all", po::value<int>(&incl_all)->default_value(1), "whether or not to include all the input copy numbers without further propressing")
    ("infer_wgd", po::value<int>(&infer_wgd)->default_value(0), "whether or not to infer WGD status of a sample from its copy numbers")
    ("infer_chr", po::value<int>(&infer_chr)->default_value(0), "whether or not to infer chromosome gain/loss status of a sample from its copy numbers")
    ("use_repeat", po::value<int>(&use_repeat)->default_value(1), "whether or not to use repeated site patterns when computing the likelihood")

    ("infer_marginal_state", po::value<int>(&infer_marginal_state)->default_value(1), "whether or not to infer marginal ancestral state of MRCA")
    ("infer_joint_state", po::value<int>(&infer_joint_state)->default_value(1), "whether or not to infer joint ancestral state of all internal nodes")
    ("min_asr", po::value<double>(&min_asr)->default_value(0.5), "minimum posterior probability to determine the best ancestral state")

    ("m_max", po::value<int>(&m_max)->default_value(1), "maximum number of copies of a segment in a chromosome")
    ("max_wgd", po::value<int>(&max_wgd)->default_value(1), "maximum number of WGD")
    ("max_chr_change", po::value<int>(&max_chr_change)->default_value(1), "maximum number of chromosome changes")
    ("max_site_change", po::value<int>(&max_site_change)->default_value(2), "maximum number of segment changes")

    ("epop", po::value<int>(&Ne)->default_value(2), "effective population size of cell populations")
    ("gtime", po::value<double>(&gtime)->default_value(1), "generation time in year")
    ("beta", po::value<double>(&beta)->default_value(0), "population growth rate")

    ("npop,p", po::value<int>(&Npop)->default_value(100), "number of population in genetic algorithm or maximum number of initial trees")
    ("ngen,g", po::value<int>(&Ngen)->default_value(50), "number of generation in genetic algorithm or maximum number of times to perturb/optimize a tree")
    ("nstop,e", po::value<int>(&max_static)->default_value(10), "Stop after this number of times that the tree does not get improved")
    ("tolerance,r", po::value<double>(&tolerance)->default_value(1e-2), "tolerance value in maximization methods")
    ("miter,m", po::value<int>(&miter)->default_value(2000), "maximum number of iterations in maximization")
    ("ssize,z", po::value<double>(&ssize)->default_value(0.01), "initial step size")

    ("only_seg", po::value<int>(&only_seg)->default_value(0), "Whether or not to only consider segment-level mutations (0: include chromosome gain/loss and whole genome doubling, 1: only consider segment-level mutations)")
    ("mu,x", po::value<double>(&mu)->default_value(0.02), "mutation rate")
    ("dup_rate", po::value<double>(&dup_rate)->default_value(0.01), "duplication rate")
    ("del_rate", po::value<double>(&del_rate)->default_value(0.01), "deletion rate")
    ("chr_gain_rate", po::value<double>(&chr_gain_rate)->default_value(0), "chromosome gain rate")
    ("chr_loss_rate", po::value<double>(&chr_loss_rate)->default_value(0), "chromosome loss rate")
    ("wgd_rate", po::value<double>(&wgd_rate)->default_value(0), "WGD (whole genome doubling) rate")
    // ("max_rate", po::value<double>(&max_rate)->default_value(0.05), "The maximum rate of a mutation event (upper bound of BFGS optimization)")

    ("model,d", po::value<int>(&model)->default_value(2), "model of evolution (0: Mk, 1: one-step bounded (total), 2: one-step bounded (allele-specific, 3: independent Markov chains)")
    ("constrained", po::value<int>(&cons)->default_value(1), "constraints on branch length (0: none, 1: fixed total time)")
    ("fixm", po::value<int>(&maxj)->default_value(0), "estimation of mutation rate (0: mutation rate fixed to be the given value, 1: estimating mutation rate)")
    ("optim", po::value<int>(&optim)->default_value(1), "method of optimization (0: Simplex, 1: L-BFGS-B)")
    ("tree_search", po::value<int>(&tree_search)->default_value(1), "method of searching tree space (0: Genetic algorithm, 1: Random-restart hill climbing, 2: Exhaustive search)")
    ("speed_nni", po::value<int>(&speed_nni)->default_value(1), "whether or not to do reduced NNI while doing hill climbing NNIs")
    ("init_tree", po::value<int>(&init_tree)->default_value(0), "method of building inital tree (0: Random coalescence tree, 1: Maximum parsimony tree)")
    ("dir_itrees", po::value<string>(&dir_itrees)->default_value(""), "directory containing inital trees")
    ("correct_bias", po::value<int>(&correct_bias)->default_value(1), "correct ascertainment bias")

    ("bootstrap,b", po::value<int>(&bootstrap)->default_value(0), "doing bootstrap or not")
    ("mode", po::value<int>(&mode)->default_value(1), "running mode of the program (0: Compute maximum likelihood tree from copy number profile; 1: Test on example data; 2: Compute the likelihood of a given tree with branch length; 3: Compute the maximum likelihood of a given tree)")
    ("verbose", po::value<int>(&debug)->default_value(0), "verbose level (0: default, 1: debug)")
    ("seed", po::value<int>(&seed)->default_value(0), "seed used for generating random numbers")
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
          cout << "svtreeml [version 0.1], a program to build a maximum likelihood phylogenetic tree from copy number profile" << endl;
          cout << "This program can also be used to compute the likelihood of a tree with fixed branch lengths or maximum likelihood of a tree with fixed topology (and/or mutation rates), given the copy number profile." << endl;
          return 1;
      }
      po::notify(vm);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    setup_rng(seed);

    if (cons==0){
      cout << "Assuming the tree is unconstrained when doing optimization" << endl;
    }
    else{
      cout << "Assuming the tree is constrained by age at sampling time" << endl;
    }

    if (maxj==0){
      cout << "\nAssuming mutation rate is fixed " << endl;
    }
    else{
      cout << "\nEstimating mutation rate" << endl;
      if (only_seg==0){
        cout << "\tEstimating mutation rates for segment duplication/deletion, chromosome gain/loss, and whole genome doubling " << endl;
      }
      else{
        cout << "\tEstimating mutation rates for segment duplication/deletion " << endl;
      }
    }

    if (correct_bias==0){
      cout << "\nNot correcting acquisition bias in likelihood computation " << endl;
    }
    else{
      cout << "\nCorrecting acquisition bias in likelihood computation " << endl;
    }

    if (optim==0){
      cout << "\nUsing Simplex method for optimization" << endl;
    }
    else{
      cout << "\nUsing L-BFGS-B method for optimization" << endl;
    }

    // tobs already defined globally
    if(timefile != ""){
        tobs = read_time_info(timefile, Ns, age);
    }
    else{
        tobs = vector<double>(Ns, 0);
    }

    if(cons){
        cout << "\nThe age of patient at the first sampling time: " << age << endl;
    }
    else{
        age = MAX_AGE;
    }

    vector<double> rates;
    if(model==0){
        rates.push_back(mu);
        cout << "\nAssuming Mk model " << endl;
    }
    else{
        if(model == 1) cout << "\nAssuming One-step bounded model " << endl;
        if(model == 2) cout << "\nAssuming One-step allele-specific model " << endl;
        if(model == 3){
            cout << "\nAssuming independent Markov chains" << endl;
            if(is_total == 0){
                cout << "This model only supports total copy number for now!" << endl;
                exit(1);
            }
            only_seg = 0;
        }
        rates.push_back(dup_rate);
        rates.push_back(del_rate);
        if(!only_seg){
            rates.push_back(chr_gain_rate);
            rates.push_back(chr_loss_rate);
            rates.push_back(wgd_rate);
        }
    }

    int num_total_bins = 0;
    Nchar = 0;
    num_invar_bins = 0;
    map<int, vector<vector<int>>> data;
    cout << "\nReading input copy numbers" << endl;
    if(is_bin==1){
        data = read_data_var_regions_by_chr(datafile, Ns, cn_max, num_invar_bins, num_total_bins, Nchar, obs_num_wgd, obs_change_chr, is_total);
        cout << "   Merging consecutive bins in the input" << endl;
    }else{
        if(incl_all==1){
            cout << "   Using all input segments " << endl;
            correct_bias = 0;
        }else{
            cout << "   Using variable input segments " << endl;
        }
        data = read_data_regions_by_chr(datafile, Ns, cn_max, num_invar_bins, num_total_bins, Nchar, obs_num_wgd, obs_change_chr, incl_all, is_total);
    }
    // cout << "Number of invariant bins " << num_invar_bins << endl;
    if(num_total_bins == num_invar_bins){
        cout << "There are no variant segments in the input data!" << endl;
        exit(1);
    }
    if(use_repeat == 1){
        cout << "   Using site repeats to spped up likelihood computation " << endl;
    }

    vobs = get_obs_vector_by_chr(data);
    if(infer_wgd == 1){
        for(int i=0; i<Ns; ++i){
            cout << "Sample " << i+1 << " probably has " << obs_num_wgd[i] << " WGD events" << endl;
        }
    }

    // Build the table after reading input file
    if(model == 3){
        adjust_m_max();
        cout << "maximum number of WGD events is " << max_wgd << endl;
        cout << "maximum number of chromosme gain/loss events on one chromosme is " << max_chr_change << endl;
        cout << "maximum number of segment duplication/deletion events is " << max_site_change << endl;
        decomp_table = build_decomp_table();
        // decomp_table = build_decomp_table_withm();
        cout << "\tNumber of states is " << comps.size() << endl;
        print_decomp_table(decomp_table);
        print_comps(comps);
    }

    string real_tstring = "";
    if(tree_file != ""){
      evo_tree real_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);
      real_tstring = order_tree_string_uniq(create_tree_string_uniq(real_tree));
      double Ls = 0;
      if(model == 3){
          // Ls = -get_likelihood_decomp(Ns, Nchar, num_invar_bins, vobs, real_tree, decomp_table, cons, cn_max, m_max, max_wgd, max_chr_change, max_site_change, correct_bias, is_total);
          Ls = -get_likelihood_decomp(real_tree);
      }else{
          // Ls = -get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, real_tree, model, cons, cn_max, only_seg, correct_bias, is_total);
          if(debug){
              cout << "Parameters for likelihood computation: " << Ns << "\t" << Nchar << "\t" << num_invar_bins << "\t" << model << "\t" << cons << "\t" << cn_max << "\t" << only_seg << "\t" << correct_bias << "\t" << is_total << endl;
          }
          Ls = -get_likelihood_revised(real_tree);
      }
      cout.precision(dbl::max_digits10);
      cout << "Score for tree real is: " << Ls << endl;
    }

    if (mode == 0){
      cout << "\nBuilding maximum likelihood tree from copy number profile " << endl;
      if (init_tree==0){
        cout << "   Using random coalescence trees as initial trees " << endl;
      }
      else{
        cout << "   Using maximum parsimony trees as initial trees " << endl;
      }

      if(bootstrap){
          cout << "\nDoing bootstapping " << endl;
          get_bootstrap_vector_by_chr(data, vobs);
          // cout << endl;
          if(debug){
              cout << " Copy number matrix after bootstapping" << endl;
              for(auto it : vobs){
                  cout << it.first << "\t" << it.second.size() << endl;
              }
          }
      }
      cout << "\nNumber of invariant bins after reading input is: " << num_invar_bins << endl;
      int total_chr = data.rbegin()->first;
      find_ML_tree(real_tstring, total_chr, num_total_bins, ofile, tree_search, Npop, Ngen, init_tree, dir_itrees, max_static, ssize, tolerance, miter, optim, rates, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, Ne, beta, gtime);
  }
  else if(mode == 1){
      cout << "Running test on tree " << tree_file << endl;
      // int num_invar_bins0 = 0;
      // int num_total_bins0 = 0;
      // int Nchar0 = 0;
      // cout << "Running data without grouping by chromosome" << endl;
      // vector<vector<int>> data0 = read_data_var_regions(datafile, Ns, cn_max, num_invar_bins0, num_total_bins0, Nchar0, is_total);
      // // Construct the CN matrix
      // // cout << "The number of sites used in vobs0: " << data0.size() << endl;
      // vector<vector<int>> vobs0;
      // for(int nc=0; nc<data0.size(); ++nc){
      //   vector<int> obs;
      //   for(int i=0; i<Ns; ++i){
      //     obs.push_back( data0[nc][i+3] );
      //   }
      //   vobs0.push_back( obs );
      // }
      //
      // run_test(tree_file, Ns, num_total_bins, Nchar, model, cn_max, only_seg, correct_bias, is_total, tobs, vobs0, Nchar0, rates, ssize, tolerance, miter);

      check_NNI(tree_file);
  }
  else if(mode == 2){
      cout << "Computing the likelihood of a given tree from copy number profile " << endl;
      double lnl = compute_tree_likelihood(tree_file, decomp_table, Ns, Nchar, num_invar_bins, vobs, tobs, rates, model, cons, cn_max, max_wgd, max_chr_change, max_site_change, only_seg, correct_bias, is_total);
      cout << "The negative log likelihood of the input tree is " << -lnl << endl;
  }
  else if(mode == 3){
      cout << "Computing maximum likelihood of a given tree from copy number profile " << endl;
      maximize_tree_likelihood(tree_file, ofile, Ns, Nchar, tobs, rates, ssize, tolerance, miter, model, optim, cons, maxj, cn_max, only_seg, correct_bias, is_total);
  }
  else{
      cout << "Inferring ancestral states of a given tree from copy number profile " << endl;
      // maximize_tree_likelihood(tree_file, ofile, Ns, Nchar, tobs, rates, ssize, tolerance, miter, model, optim, cons, maxj, cn_max, only_seg, correct_bias, is_total);
      evo_tree real_tree = init_tree_from_file(tree_file, Ns, Nchar, model, only_seg, tobs, rates);

      if(infer_marginal_state==1){
          cout << "\tInferring marginal ancestral states" << endl;
          string ofile_mrca = ofile + ".mrca.state";
          ofstream fout(ofile_mrca);
          double lnl = 0;
          if(model==3){
             lnl = reconstruct_marginal_ancestral_state_decomp(real_tree, fout, min_asr);
         }else{
             lnl = reconstruct_marginal_ancestral_state(real_tree, fout, min_asr);
         }
         cout << "The negative log likelihood of the input tree is " << -lnl << endl;
      }

      if(infer_joint_state == 1){
          cout << "\tInferring joint ancestral states" << endl;
          string ofile_state = ofile + ".joint.state";
          ofstream fout_state(ofile_state);
          if(model==3){

             reconstruct_joint_ancestral_state_decomp(real_tree, fout_state, min_asr);
         }else{
             reconstruct_joint_ancestral_state(real_tree, fout_state, min_asr);
         }
      }
  }
}
