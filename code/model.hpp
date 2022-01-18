#ifndef MODEL_HPP
#define MODEL_HPP

// This file implements the model of copy number evolution

#include "common.hpp"   // for model enum

// using namespace std;  // must before r8lib.hpp


#include "matexp/matrix_exponential.hpp"
#include "matexp/r8lib.hpp"



typedef vector<double>::iterator DBIter;  // convenience typedefs
typedef pair<DBIter, DBIter> DBIterPair;


// to validate the rate matrix (row sum should be 0)
bool check_matrix_row_sum(double *mat, int nstate);

void check_pmats_blen2(int nstate, const vector<double>& blens, const vector<double*> pmat_per_blen);

/*************** MK model *****************/
// Mk model based on JC model
inline double get_transition_prob(const double& mu, const double& blength, const int& sk, const int& sj){
  //if(debug) cout << "\tget_transition_prob" << endl;

  // assume a basic Markov model here cf JC
  // Qij = u/5 when i!=j
  // we have five states 0, 1, 2, 3, 4
  // Pij(t) = exp(-ut) k_ij + (1 - exp(-ut))*pi_j
  // pi_0 = pi_1 = pi_2 = pi_3 = pi_4 = 1/5

  double prob = 0.0;

  if(sk == sj){
    prob = exp(-mu * blength) + (1 - exp(-mu * blength)) * 0.2;
  }else{
    prob = (1 - exp(-mu * blength)) * 0.2;
  }

  return prob;
}

/*************** BOUND model *****************/

// model when copy number is bounded by 0 and cn_max
void get_rate_matrix_bounded(double* m, const double& dup_rate, const double& del_rate, const int& cn_max);

// A matrix with dimension (cn_max + 1) * (cn_max + 2) / 2
// Suppose copy number configuration is in specific order such as:  0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2	 2/1	3/0	0/4	 1/3	  2/2	 3/1	4/0
void get_rate_matrix_allele_specific(double* m, const double& dup_rate, const double& del_rate, const int& cn_max);

// n = cn_max + 1 for model 1 (total copy number)
void get_transition_matrix_bounded(double* q, double* p, const double& t, const int& n);

// not used in practice to save effeorts in function call
double get_transition_prob_bounded(double* p, const int& sk, const int& sj, const int& n);


/*************** DECOMP model *****************/

// model which decompose copy numbers into different levels

// rate matrix for segment-level CNAs
void get_rate_matrix_site_change(double* m, const double& dup_rate, const double& del_rate, const int& site_change_max);

// rate matrix for chr-level CNAs
void get_rate_matrix_chr_change(double* m, const double& chr_gain_rate, const double& chr_loss_rate, const int& chr_change_max);

// Maximum allowed number of WGD events over a time interval: 2
void get_rate_matrix_wgd(double* m, const double& wgd_rate, const int& wgd_max = 2);

// Get possible combinations for a total copy number (n = (alpha, beta, gamma))
// deprecated due to excluding of WGD
void insert_tuple(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int i, int j, int k);

// Get possible combinations for a total copy number with ordering of different types of events considered (at most 1 WGD), mainly before and after WGD
void insert_tuple_order(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int i, int j, int k, int j0, int k0);

// Get possible combinations for a total copy number with ordering of different types of events considered (at most 1 WGD), mainly before and after WGD (seem problematic)
void insert_tuple_order_withm(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int m1, int m2, int i, int j, int k, int j0, int k0);

// Get possible haplotype-specific combinations for a total copy number (not used in practice for now)
void insert_tuple_allele_specific(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int i, int j1, int j2, int k1, int k2);

// Adjust the value of m_max based on the maximum copy number for a sample, not used for now
void adjust_m_max(const vector<int>& obs_num_wgd, const vector<int>& sample_max_cn, int m_max, int max_chr_change, int max_site_change);

// list all the possible decomposition of a copy number, treating m_j as part of the vector
void build_decomp_table(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int max_wgd, int max_chr_change, int max_site_change, int is_total = 1);

// list all the possible decomposition of a copy number, treating m_j before and after WGD differently
void build_decomp_table_withm(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int max_wgd, int max_chr_change, int max_site_change, int is_total = 1);

void print_comps(const set<vector<int>>& comps);

void print_decomp_table(const map<int, set<vector<int>>>& decomp_table);


#endif
