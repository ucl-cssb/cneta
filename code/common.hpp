#ifndef COMMON_HPP
#define COMMON_HPP

#include <cstdio>
#include <iostream> // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <fstream>
#include <sstream>

#include <string>
#include <cstring>

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <algorithm>
#include <random>

#include <map>
#include <vector>
#include <set>
#include <stack>

#include <boost/format.hpp>   // used in newick format
#include <boost/algorithm/string.hpp>

#include "gzstream.h"


using namespace std;

// MK model was used for test, not suitable for CN evolution
// INFINITE was used for simulating Poisson process
enum MODEL {MK, BOUNDT, BOUNDA, DECOMP, INFINITE = 3};

const int MAX_AGE = 100;

const int NUM_CHR = 22;
const int NORM_PLOIDY = 2;
const int NORM_ALLElE_STATE = 4;

const int PRINT_PRECISION = 10;

const double SMALL_VAL = 1.0e-10;   // used to compare floats

const int WGD_CUTOFF = 3;    // genome ploidy to determine WGD

// key: chr, seg, copy_number
typedef map<int, map<int, int>> copy_number;

const string VERSION = "1.0";


// read-only
struct INPUT_PROPERTY{
  int Ns;
  int cn_max;
  int model;

  int is_total;
  int is_rcn;
  int is_bin;
  int incl_all;
};

// obtained from input file
struct INPUT_DATA{
  int num_invar_bins;
  int num_total_bins;
  int seg_size;  // Nchar

  vector<int> obs_num_wgd;
  vector<vector<int>> obs_change_chr;
  vector<int> sample_max_cn;
};


// Read strings separated by space into a vector.
// num: the number of element in the string. If not given, assume it is the first number in the string
template <typename T>
void get_vals_from_str(vector<T>& vals, string str_vals, int num = 0){
    assert(str_vals != "");
    stringstream ss(str_vals);
    if(num==0)   ss >> num;
    for(int i = 0; i < num; i++){
        T d1;
        ss >> d1;
        vals.push_back(d1);
    }
    // cout << "vector from " << str_vals << ": ";
    // for(auto v : vals){
    //     cout << "\t" << v;
    // }
    // cout << endl;
}


#endif
