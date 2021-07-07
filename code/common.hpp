#ifndef COMMON_HPP
#define COMMON_HPP

#include <cstdio>
#include <iostream> // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <fstream>
#include <sstream>

#include <string>
#include <cstring>

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


enum MODEL {MK, BOUNDT, BOUNDA, DECOMP};

const int MAX_AGE = 100;

const int NUM_CHR = 22;
const int NORM_PLOIDY = 2;

const int PRINT_PRECISION = 10;

const double SMALL_VAL = 1.0e-10;   // used to compare floats

#endif
