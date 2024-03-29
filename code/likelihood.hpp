#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

//
// Likelihood computation
//

#include "parse_cn.hpp"
#include "evo_tree.hpp"
#include "tree_op.hpp"
#include "model.hpp"

// using namespace std;


// from book Effective STL
struct DeleteObject{
  template<typename T>
  void operator()(const T* ptr) const{
    delete [] ptr;
  }
};


// information derived from input data for DECOMP model
struct OBS_DECOMP{
  int m_max;   // maximum copy of a segment before chr-level events, used in likelihood table initialization

  // used to get dimension of Pmatrix at different levels
  int max_wgd;
  int max_chr_change;
  int max_site_change;

  vector<int> obs_num_wgd;  // possible number of WGD events, used in likelihood table initialization
  vector<vector<int>> obs_change_chr; // possible number of chr-level events, used in likelihood table initialization
  // vector<int> sample_max_cn;  // not used in likelihood computation
};


// information derived from input data
struct MAX_DECOMP{
  int m_max;   // maximum copy of a segment before chr-level events, used in likelihood table initialization

  // used to get dimension of Pmatrix at different levels
  int max_wgd;
  int max_chr_change;
  int max_site_change;
};



// information derived from input data
struct DIM_DECOMP{
  int dim_wgd;
  int dim_chr;
  int dim_seg;
};


struct PMAT_DECOMP{
  map<double, double*> pmats_wgd;
  map<double, double*> pmats_chr;
  map<double, double*> pmats_seg;
};


// used in computing likelihood of children nodes
struct PROB_DECOMP{
  double* pbli_wgd;
  double* pblj_wgd;

  double* pbli_chr;
  double* pblj_chr;

  double* pbli_seg;
  double* pblj_seg;
};


// variables used in likelihood computation
struct LNL_TYPE{
  int model;
  int cn_max;
  int is_total; // whether or not the input is total copy number

  // used to check tree validity
  int cons;
  double max_tobs;  // maximum sampling time
  int patient_age;   // used in BFGS constraints

  int use_repeat;   // whether or not to use repeated site patterns, used in get_likelihood_chr*
  int correct_bias; // Whether or not to correct acquisition bias, used in get_likelihood_*
  int num_invar_bins; // used in get_likelihood_*

  int only_seg; // Whether or not to only consider segment-level mutations, used in get_likelihood_revised

  int infer_wgd; // whether or not to infer WGD status of a sample, called in initialize_lnl_table_decomp
  int infer_chr; // whether or not to infer chromosome gain/loss status of a sample, called in initialize_lnl_table_decomp

  vector<int> knodes;
};

const double LARGE_LNL = -1e9;
const double SMALL_LNL = -1e20;
const double MAX_NLNL = 1e20;


/****************** common functions *******************/
void print_lnl_at_tips(const evo_tree& rtree, const vector<int>& obs, const vector<vector<double>>& L_sk_k, int nstate);


void print_tree_lnl(const evo_tree& rtree, vector<vector<double>>& L_sk_k, int nstate);

// From https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
    const std::vector<T>& vec,
    const Compare& compare)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
    return p;
}

template <typename T>
std::vector<T> apply_permutation(
    const std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](std::size_t i){ return vec[i]; });
    return sorted_vec;
}



/****************** functions for non DECOMP model *******************/
// Create likelihood vectors at the tip node, one table for each site
// obs: a vector of CNs for one site
// L_sk_k has one row for each tree node and one column for each possible state
// only used in get_likelihood_chr
// extracted as a function to avoid duplication in selection statement
void initialize_lnl_table(vector<vector<double>>& L_sk_k, const vector<int>& obs, const evo_tree& rtree, int model, int nstate, int is_total);


// nstate = cn_max + 1
// only used in get_likelihood_site
// extracted as a function to avoid duplication in selection statement
inline double get_prob_children(vector<vector<double>>& L_sk_k, const evo_tree& rtree, double* pbli, double* pblj, int nsk, int ni, int nj, int bli, int blj, int model, int nstate){
    double Li = 0.0;
    for(int si = 0; si < nstate; ++si){
      if(L_sk_k[ni][si] > 0){
        if(model == MK){
          Li += get_transition_prob(rtree.mu, bli, nsk, si) * L_sk_k[ni][si];
        }else{
          Li += pbli[nsk + si * nstate] * L_sk_k[ni][si];
        }
      }
    }

    double Lj = 0.0;
    for(int sj = 0; sj < nstate; ++sj){
      if(L_sk_k[nj][sj] > 0){
        if(model == MK){
          Lj += get_transition_prob(rtree.mu, blj, nsk, sj) * L_sk_k[nj][sj];
        }else{
          Lj += pblj[nsk + sj * nstate] * L_sk_k[nj][sj];
        }
      }
    }

    return Li * Lj;
}



// Get the likelihood on one site of a chromosome (assuming higher level events on nodes)
// z: possible changes in copy number caused by chromosome gain/loss
// only used in get_likelihood_chr
// extracted as a function to avoid duplication in selection statement
void get_likelihood_site(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, const int& has_wgd, const int& z, const int& model, const int& nstate);


// Get the likelihood on a set of chromosmes
// only used in get_likelihood_revised
// extracted as a function to avoid duplication in selection statement
double get_likelihood_chr(map<int, vector<vector<int>>>& vobs, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, const int& has_wgd, const int& only_seg, const int& use_repeat, const int& model, const int& nstate, const int& is_total);



// Incorporate chromosome gain/loss and WGD
// Model 2: Treat total copy number as the observed data and the haplotype-specific information is missing
// additional parameters moved to lnl_type to facilitate calling in optimization
double get_likelihood_revised(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, LNL_TYPE& lnl_type);


/* Compute the likelihood without grouping sites by chromosome, only considering site duplication/deletion (not used)
Precondition: the tree is valid
Ns: number of samples
Nchar: number of characters for each sample
vobs: the observed data matrix
rtree: the given tree
model: model of evolution
*/
double get_likelihood(const vector<vector<int>>& vobs, evo_tree& rtree, const vector<int>& knodes, int model, int cn_max, int is_total);


// Get the likelihood of the tree from likelihood table
double extract_tree_lnl(vector<vector<double>>& L_sk_k, int Ns, int model);


/************** functions for model DECOMP **************/

// L_sk_k has one row for each tree node and one column for each possible state; chr starting from 1
// This function is critical in obtaining correct likelihood. If one tip is not initialized, the final likelihood will be 0.
// nstate = comps.size();
vector<vector<double>> initialize_lnl_table_decomp(vector<int>& obs, OBS_DECOMP& obs_decomp, int chr, const evo_tree& rtree, const set<vector<int>>& comps, int infer_wgd, int infer_chr, int cn_max, int is_total = 1);


// Assume the likelihood table is for each total copy number (no WGD order considered, deprecated)
double get_prob_children_decomp(vector<vector<double>>& L_sk_k, const evo_tree& rtree, map<int, set<vector<int>>>& decomp_table, int sk, int cn_max, int nstate, PROB_DECOMP& prob_decomp, DIM_DECOMP& dim_decomp, int ni, int nj, int bli, int blj, int is_total);


// Assume the likelihood table is for each combination of states
double get_prob_children_decomp2(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, int sk, PROB_DECOMP& prob_decomp, DIM_DECOMP& dim_decomp, int ni, int nj, int bli, int blj, int cn_max, int nstate, int is_total);


// Get the likelihood on one site of a chromosome
// Assuming each observed copy number is composed of three type of events.
// Sum over all possible states for initial and final nodes
// Allow at most one WGD event along a branch
void get_likelihood_site_decomp(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, const vector<int>& knodes, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int cn_max, int is_total);

// Used when WGD is considered, dealing with mutations of different types at different levels (no WGD order considered, deprecated)
double get_likelihood_chr_decomp(map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const evo_tree& rtree, const set<vector<int>>& comps, const vector<int>& knodes, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int infer_wgd, int infer_chr, int use_repeat, int cn_max, int is_total);

// Computing likelihood when WGD and chr gain/loss are incorporated
// Assume likelihood is for haplotype-specific information
double get_likelihood_decomp(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type);

// Get the likelihood of the tree from likelihood table of state combinations
double extract_tree_lnl_decomp(vector<vector<double>>& L_sk_k, const set<vector<int>>& comps, int Ns);

#endif
