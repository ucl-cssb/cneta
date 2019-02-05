#include <stdio.h>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

using namespace std;


/*
This file uses integer linear programing (ILP) to find the minimal number of copy number events that transform one copy number profile to another copy number profile.
It uses IBM ILOG CPLEX Optimization Studio 12.8, which is free for academic use.
To install CPLEX Studio, you can download the software from https://www.ibm.com/products/ilog-cplex-optimization-studio.
Please change the variables according to your local installation of CPLEX Studio.
*/

// Compute the probability for vectors of different sizes according to distributions
const int MEAN_DUP_SIZE = 10;
const int MEAN_DEL_SIZE = 10;
const int MAX_DUP_SIZE = 30;
const int MAX_DEL_SIZE = 30;
// The maximum size of a segment
const int BIN_SIZE = 385;
// const int MAX_DUP_SIZE = 3;
// const int MAX_DEL_SIZE = 3;
// const int BIN_SIZE = 10;


int get_num_vec();
void create_bvectors(int num_vec, int** bvecs, double* bprobs, int* bsizes);
void create_model(IloModel model, IloNumVarArray cols, IloRangeArray range, const int* cnv_diff, int num_vec, int** bvecs, const double* bprobs, const int* bsizes);
int solve_ilp(int* profile1, int* profile2);


int main(int argc, char **argv)
{
    int profile1[BIN_SIZE] = {0};
    int profile2[BIN_SIZE] = {0};
    // create two example CNV profiles
    profile1[BIN_SIZE-1] = 1;
    for (int i = BIN_SIZE-1; i>=BIN_SIZE-5; i--){
        profile2[i] = 2;
    }

    solve_ilp(profile1, profile2);
}

// Find the number of total vectors
int get_num_vec(){
  // The number of vectors with minimal size
  int num_min_vec = BIN_SIZE;
  int num_max_vec = BIN_SIZE + 1 - MAX_DUP_SIZE;
  int num_vec_pos = 0;
  for (int i = num_max_vec; i <= num_min_vec; i++){
      num_vec_pos += i;
  }
  int num_vec_neg = num_vec_pos;
  int num_vec = num_vec_pos + num_vec_neg + 1;
  cout << "The number of base vectors is: " << num_vec << endl;
  return num_vec;
}


// create a 2D array of base vectors
void create_bvectors(int num_vec, int** bvecs, double* bprobs, int* bsizes)
{
    int i,j,k;
    int last_pos;
    int n=0;    // The count of base vectors
    // create positive vectors of increasing size
    for(i=1; i<=MAX_DUP_SIZE; i++){
        // cout<<"positive chunks of size "<< i << "\n";
        // the last starting position for a specific duplication size
        last_pos = BIN_SIZE - i + 1;
        for (j=1; j<=last_pos; j++){
            for(k=0; k<i; k++){
                bvecs[n][j+k-1]=1;
            }
            // get the weight of each vector
            double prob = gsl_ran_exponential_pdf(i, MEAN_DUP_SIZE);
            bprobs[n] = -log(prob);
            bsizes[n] = i;
            n = n + 1;
        }
    }
    // create negative vectors of increasing size
    for(i=1; i<=MAX_DEL_SIZE; i++){
        // cout<<"negative chunks of size "<< i << "\n";
        // the last starting position for a specific duplication size
        last_pos = BIN_SIZE - i + 1;
        for (j=1; j<=last_pos; j++){
            for(k=0; k<i; k++){
                bvecs[n][j+k-1]=-1;
            }
            // get the weight of each vector
            double prob = gsl_ran_exponential_pdf(i, MEAN_DEL_SIZE);
            bprobs[n] = -log(prob);
            bsizes[n] = i;
            n = n + 1;
        }
    }
    // create the vector for whole genome doubling
    for(i=0; i<BIN_SIZE; i++){
        bvecs[n][i] = 1;
    }
    double prob = gsl_ran_exponential_pdf(BIN_SIZE, MEAN_DUP_SIZE);
    bprobs[n] = -log(prob);
    bsizes[n] = BIN_SIZE;

    // cout << "The probability of each vector: " << endl;
    // for(i=0; i<num_vec; i++){
    //     cout << bprobs[i] << " ";
    // }
    // cout << endl;
    //
    // cout << "The size of each vector: " << endl;
    // for(i=0; i<num_vec; i++){
    //     cout << bsizes[i] << " ";
    // }
    // cout << endl;
    //
    // cout << "The content of each vector: " << endl;
    // for(i=0; i<num_vec; i++){
    //     for (j=0; j<BIN_SIZE;j++){
    //         cout << bvecs[i][j] << " ";
    //     }
    //     cout << endl;
    // }
}

// Create ILP model by column
// To populate by column, we first create the range constraints and the
// objective, and then create the variables and add them to the ranges and
// objective using column expressions.
void create_model(IloModel model, IloNumVarArray cols, IloRangeArray range, const int* cnv_diff, int num_vec, int** bvecs, const double* bprobs, const int* bsizes)
{
    int i, j;
    IloEnv env = model.getEnv();
    IloObjective obj = IloMinimize(env);

    // Add the difference of two CNV profiles b1...bm
    // cout << "Adding the difference of two CNV profiles" << endl;
    for (i=0; i<BIN_SIZE; i++)
    {
        range.add(IloRange(env, cnv_diff[i], cnv_diff[i]));
    }

    // cout << "Adding the base vectors" << endl;
    // Add the base vectors v1...vn
    for (j=0; j<num_vec; j++)
    {
        // cout << "Adding vector weight" << endl;
        IloNumColumn col = obj(bprobs[j]);
        // IloNumColumn col = obj(bsizes[j]);
        // IloNumColumn col = obj(1);
        // cout << "Adding vector value" << endl;
        for (i = 0; i < BIN_SIZE; i++) {
            col += range[i](bvecs[j][i]);
        }
        // cout << "Apending vector" << endl;
        cols.add(IloNumVar(col, 0, IloInfinity, ILOINT));
        col.end();
    }

    model.add(obj);
    model.add(range);
}


/*
Input: Two CNV profiles (int arrays)
Output: The optimal transformations from one profile to another
*/
int solve_ilp(int* profile1, int* profile2){
    int i,j;
    int num_vec = get_num_vec();

    // The probability of each base vector
    double* bprobs = new double[num_vec];
    // The size of each base vector
    int* bsizes = new int[num_vec];
    // The list of base vectors
    int** bvecs = new int*[num_vec];
    for(i=0; i<num_vec; i++){
        bvecs[i] = new int[BIN_SIZE];
    }
    for(i=0; i<num_vec; i++){
        for(j=0; j<BIN_SIZE; j++)
        {
            bvecs[i][j] = 0;
        }
    }

    create_bvectors(num_vec, bvecs, bprobs, bsizes);

    // Right hand side for the constraints
    int* cnv_diff = new int[BIN_SIZE];
    // create a sample CNV profile difference
    for(i=0; i<BIN_SIZE; i++){
        cnv_diff[i] = profile2[i] - profile1[i];
    }
    cout << "The difference between two CNV profiles: " << endl;
    for(i=0; i<BIN_SIZE; i++){
        cout << cnv_diff[i] << " ";
    }
    cout << endl;

    IloEnv env;
    try {
       IloModel model(env);
       IloNumVarArray cols(env);
       IloRangeArray range(env);

       cout << "Filling the constraint matrix ..." << endl;
       create_model(model, cols, range, cnv_diff, num_vec, bvecs, bprobs, bsizes);

       cout << "Solving the ILP ..." << endl;
       IloCplex cplex(model);
       cplex.exportModel("mcomp_ilp.lp");

       // Optimize the problem and obtain solution.
       if ( !cplex.solve() ) {
          env.error() << "Failed to optimize LP" << endl;
          throw(-1);
       }

       IloNumArray vals(env);
       env.out() << "Solution status = " << cplex.getStatus() << endl;
       env.out() << "Solution value  = " << cplex.getObjValue() << endl;
       cplex.getValues(vals, cols);
       // env.out() << "Values        = " << vals << endl;
       // Finding the composite base vectors
       cout<< "The index of composite base vectors (starting from 0): " << endl;
       for (i=0; i<num_vec; i++){
           if(vals[i] > 0){
               cout<<i<<" ";
           }
       }
       cout<< endl;
       cout<< "The content of composite base vectors: " << endl;
       for (i=0; i<num_vec; i++){
           if(vals[i] > 0){
               cout << "The " << i + 1 << "th base vector: " << endl;
               for(j=0; j<BIN_SIZE; j++){
                   if(bvecs[i][j]!=0)
                   {
                       cout << j << ": " << bvecs[i][j] << "; ";
                   }
               }
               cout<<endl;
           }
       }
       cout << endl;
       // cplex.getSlacks(vals, range);
       // env.out() << "Slacks        = " << vals << endl;
       // cplex.getDuals(vals, range);
       // env.out() << "Duals         = " << vals << endl;
       // cplex.getReducedCosts(vals, cols);
       // env.out() << "Reduced Costs = " << vals << endl;
   }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }

   env.end();
   return 0;
}
