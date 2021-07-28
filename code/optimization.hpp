#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

//
// Maximum likelihood and optimisation functions
//

#include <gsl/gsl_multimin.h>

#include "lbfgsb/lbfgsb_new.h"

#include "common.hpp"
#include "likelihood.hpp"
// #include "nni.hpp"

// using namespace std;

const double ERROR_X = 1.0e-4;
// The minimum mutation rates allowed
const double MIN_MRATE = 1.0e-20;
// The maximum mutation rates allowed
const double MAX_MRATE = 1;
// The minimum age ratio allowed
const double MIN_RATIO = 1e-2;
const double MAX_RATIO = 0.99;


// bundle of variables used in optimization (values from input)
struct OPT_TYPE{
  int maxj;   // used to determine optimization function

  double tolerance;
  int miter;

  int opt_one_branch;

  vector<double> tobs;  // sample times, used to get constaints in optimization

  double scale_tobs;   // scale factor to get lower limit of root age when doing constrained optimization (BFGS) based on maximimum sample time difference
};

struct GSL_PARAM{
  evo_tree rtree;
  map<int, vector<vector<int>>> vobs;
  LNL_TYPE lnl_type;
};

/*****************************************************
    optimization with GSL simplex methods (deprecated)
*****************************************************/

double my_f(const gsl_vector *v, void *params);
double my_f_mu(const gsl_vector *v, void *params);
double my_f_cons(const gsl_vector *v, void *params);
double my_f_cons_mu(const gsl_vector *v, void *params);


// Output error information without aborting
inline void my_err_handler(const char * reason, const char * file, int line, int gsl_errno){
    fprintf (stderr, "failed, %s, reason: %s, file %s line %d \n", gsl_strerror(gsl_errno), reason, file, line);
    return;
}

// given a tree, maximise the branch lengths (and optionally mu) assuming branch lengths are independent or constrained in time
// use GSL simplex optimization
void max_likelihood(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, const vector<double>& tobs, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double& min_nlnl, const double& ssize);


/*****************************************************
    One dimensional optimization with Brent method
*****************************************************/

// Compute the likelihood of the tree with one parameter (one branch or one mutation rate)
// Used in Brent optimization
// return negative likelihood function for minimalization
double computeFunction(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double value, int type = -1);



/* Brents method in one dimension */
/**
    This function calculate f(value) of the f() function, used by other general optimization method to minimize it.
    Please always override this function to adapt to likelihood or parsimony score.
    The default is for function f(x)=x.
    @param value x-value of the function
    @return f(value) of function f you want to minimize
*/
double brent_opt(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int type, double ax, double bx, double cx, double tol,
                          double *foptx, double *f2optx, double fax, double fbx, double fcx);


/**
    the brent method to find the value that minimizes the computeFunction().
    @return the x-value that minimize the function
    @param xmin lower bound
    @param xmax upper bound
    @param xguess first guess
    @param tolerance tolerance
    @param fx (OUT) function value at the minimum x found
    @param ferror (OUT) Dont know
*/
double minimizeOneDimen(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int type, double xmin, double xguess, double xmax, double tolerance, double *fx, double *f2x);


// Using Brent method to optimize the likelihood of mutation rate
double optimize_mutation_rates(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double tolerance);


// Optimizing the branch length of (node1, node2) with Brent method
// Optimize mutation rates if necessary
double optimize_one_branch(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double tolerance, Node* node1, Node* node2);


// used in optimize_all_branches
void compute_best_traversal(evo_tree& rtree, NodeVector &nodes, NodeVector &nodes2);


// optimize each branch
double optimize_all_branches(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int my_iterations, double tolerance);



/*****************************************************
    L-BFGS-B method
*****************************************************/

// The number of parameters to estimate, different when the mutation rates are estimated
inline int get_ndim(int maxj, int npar_ne, int model, int only_seg){
    int ndim = 0;

    if(!maxj){
      ndim = npar_ne;
    }else{
        if(model == MK){
            ndim = npar_ne + 1;
        }else{
            if(only_seg){
                ndim = npar_ne + 2;
            }else{
                ndim = npar_ne + 5;
            }
        }
    }
    return ndim;
}


// Update the tree after each iteration in the BFGS optimization (deprecated)
void update_variables(evo_tree& rtree, int model, int cons, int maxj, double *x);


// Update the tree after each iteration in the BFGS optimization
// Estimate time intervals rather than branch length in order to avoid negative terminal branch lengths
// Sort node times in increasing order and take the first Ns intervals
void update_variables_transformed(evo_tree& rtree, double *x, LNL_TYPE& lnl_type, OPT_TYPE& opt_type);


/**
    the target function which needs to be optimized
    @param x the input vector x
    @return the function value at x
*/
double targetFunk(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double x[]);


/**
	the approximated derivative function
	@param x the input vector x
	@param dfx the derivative at x
	@return the function value at x
*/
double derivativeFunk(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int ndim, double x[], double dfx[]);

double optimFunc(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int nvar, double *vars);

double optimGradient(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int nvar, double *x, double *dfx);


void lbfgsb(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, int *fail,
		double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT);


/**
 Function to access the L-BFGS-B function, taken from IQ-TREE package which is further taken from HAL_HAS software package
 1. int nvar : The number of the variables
 2. double* vars : initial values of the variables
 3. double* lower : lower bounds of the variables
 4. double* upper : upper bounds of the variables
 5. double pgtol: gradient tolerance
 5. int maxit : max # of iterations
 @return minimized function value
 After the function is invoked, the values of x will be updated
*/
double L_BFGS_B(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int n, double* x, double* l, double* u);

// Using BFGS method to get the maximum likelihood with lower and upper bounds (minimalize negative likelihood function)
// Note: the topology of rtree is fixed, yet its branch lengths may be updated in the optimization process
void max_likelihood_BFGS(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double &min_nlnl,int debug = 0);


// Optimizing the branch length of (node1, node2) with BFGS to incorporate constraints imposed by patient age and tip timings.
// Optimize mutation rates if necessary
double optimize_one_branch_BFGS(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, Node* node1, Node* node2);



#endif
