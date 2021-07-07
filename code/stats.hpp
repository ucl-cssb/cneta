#ifndef STATS_HPP
#define STATS_HPP


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include <cstdio>
#include <iostream>

#include <numeric>  // for accumulate, partial_sum
#include <unistd.h>  // for getpid

#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <cmath>

#include <vector>


using namespace std;


// Pairing function: see https://en.wikipedia.org/wiki/Pairing_function
inline int pairInteger(int int1, int int2){
    if(int1 <= int2){
        return ((int1 + int2) * (int1 + int2 + 1)/2 + int2);
    }else{
        return ((int1 + int2) * (int1 + int2 + 1)/2 + int1);
    }
}

void setup_rng(gsl_rng* r, unsigned set_seed);


// factorial for choose(n,k)
inline int fact(int n){
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}


// wrapper for uniform
double runiform(gsl_rng* r, double a, double b);



// sample an element according to a probability vector
// gsl_ran_multinomial?
int rchoose(gsl_rng* r, const vector<double>& rates);


#endif
