#ifndef PARSE_CN_HPP
#define PARSE_CN_HPP


#include "common.hpp"   // for enum model
#include "stats.hpp"


// collection of functions for reading/writing

// const int num_total_bins = 4401;



// template <typename T>
// vector<size_t> sort_indexes(const vector<T> &v){
//   // initialize original index locations
//   vector<size_t> idx(v.size());
//   iota(idx.begin(), idx.end(), 0);
//
//   // sort indexes based on comparing values in v
//   sort(idx.begin(), idx.end(),
//        [&v](size_t i1, size_t i2){return v[i1] < v[i2];});
//
//   return idx;
// }


inline bool is_equal_vector(const vector<int>& bin1, const vector<int>& bin2){
    assert(bin1.size() == bin2.size());

    for(int i = 0; i < bin1.size(); i++){
        if(bin1[i] != bin2[i]){
            return false;
        }
    }
    return true;
}


// Convert the allele-specific state (ordered by 0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2 2/1	3/0	 0/4 1/3 2/2 3/1	4/0) to total copy number
int state_to_total_cn(int state, int cn_max);

// Convert the allele-specific state (ordered by 0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2 2/1	3/0	 0/4 1/3 2/2 3/1	4/0) to total copy number
int state_to_allele_cn(int state, int cn_max, int& cnA, int& cnB);

// Change the allele specific copy number to the state used in substitution rate matrix
int allele_cn_to_state(int cnA, int cnB);


// Find the potential number of WGDs for each sample
void get_num_wgd(const vector<vector<vector<int>>>& s_info, int cn_max, vector<int>& obs_num_wgd, int is_total, int debug = 0);


// Find the potential number of chromosome changes for each sample
void get_change_chr(const vector<vector<vector<int>>>& s_info, vector<vector<int>>& obs_change_chr, int cn_max, int is_total, int debug = 0);


// Find the largest copy number in a sample
void get_sample_mcn(const vector<vector<vector<int>>>& s_info, vector<int>& sample_max_cn, int cn_max, int is_total, int debug = 0);


// Distinguish invariable and variable sites;
// Combine adjacent invariable sites
// TODO: Distinguish different types of invariant sites: 2 vs 4, which have different likelihoods
vector<vector<int>> get_invar_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int is_total = 1, int debug = 0);

vector<vector<int>> get_all_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int incl_all, int is_total = 1, int debug = 0);

map<int, vector<vector<int>>>  group_segs_by_chr(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug = 0);

vector<vector<int>> group_segs(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug = 0);

// Get the input matrix of copy numbers by chromosome
map<int, vector<vector<int>>> get_obs_vector_by_chr(map<int, vector<vector<int>>>& data, const int& Ns);

void get_bootstrap_vector_by_chr(map<int, vector<vector<int>>>& data, map<int, vector<vector<int>>>& vobs, gsl_rng* r);


/******************* read input file ***********************/
// Read the samping time and patient age of each sample
// Ns: number of samples
// age: age at earliest sample
vector<double> read_time_info(const string& filename, const int& Ns, int& age, int debug = 0);

// Format of input total copy number: sample, chr, seg, copy_number
// Format of input allele specific copy number: sample, chr, seg, copy_number A, copy_number B -> converted into a specific number by position
vector<vector<vector<int>>> read_cn(const string& filename, int Ns, int &num_total_bins, int cn_max, int is_total = 1, int debug = 0);


// Read the input copy numbers
vector<vector<int>> read_data_var_regions(const string& filename, const int& Ns, const int& cn_max, int& num_invar_bins, int& num_total_bins, int& seg_size, vector<int>&  obs_num_wgd, vector<vector<int>>& obs_change_chr, vector<int>& sample_max_cn, int model, int is_total = 1, int debug = 0);


// Read the input copy numbers and group them by chromosome
map<int, vector<vector<int>>> read_data_var_regions_by_chr(const string& filename, const int& Ns, const int& cn_max, int& num_invar_bins, int& num_total_bins, int &seg_size, vector<int>&  obs_num_wgd, vector<vector<int>>& obs_change_chr, vector<int>& sample_max_cn, int model, int is_total = 1, int is_bin = 1, int incl_all = 1, int debug = 0);



#endif
