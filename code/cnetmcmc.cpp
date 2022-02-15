// run Bayesian MCMC inference

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "gzstream.h"

#include "evo_tree.hpp"
#include "stats.hpp"
#include "parse_cn.hpp"
#include "nni.hpp"

// using namespace std;

const double LOG_MINUS_INFINITY = numeric_limits<double>::lowest();

int debug = 0;

gsl_rng* r;
unsigned seed;


/********* input parameters ***********/
int Ns;   // number of samples
int age;

// Set max_* as global variables to avoid adding more parameters in maximization
int m_max;
int max_wgd;
int max_chr_change;
int max_site_change;

int cn_max;
int is_total; // whether or not the input is total copy number

int model;
int cons;
int maxj;

int use_repeat;   // whether or not to use repeated site patterns, used in get_likelihood_chr*
int correct_bias; // Whether or not to correct acquisition bias, used in get_likelihood_*
int num_invar_bins;   // number of invariant sites

int only_seg; // Whether or not to only consider segment-level mutations, used in get_likelihood_revised

int infer_wgd; // whether or not to infer WGD status of a sample, called in initialize_lnl_table_decomp
int infer_chr; // whether or not to infer chromosome gain/loss status of a sample, called in initialize_lnl_table_decomp

int nstate;

/********* derived from input ***********/
map<int, vector<vector<int>>> vobs;   // CNP for each site, grouped by chr
vector<double> tobs; // input times for each sample, should be the same for all trees, defined when reading input, used in likelihood computation
double max_tobs;

int Nchar;  // number of sites


vector<int> obs_num_wgd;  // possible number of WGD events
vector<vector<int>> obs_change_chr;
vector<int> sample_max_cn;

map<int, set<vector<int>>> decomp_table;  // possible state combinations for observed copy numbers
set<vector<int>> comps;

LNL_TYPE lnl_type;
OBS_DECOMP obs_decomp;

//create a list of nodes to loop over, making sure the root is last
vector<int> knodes;

// unary function and pointer to unary function
// allows use of gsl rng for standard template algorithms
inline long unsigned myrng(long unsigned n){
  return gsl_rng_uniform_int(r, n);
}

long unsigned (*fp_myrng)(long unsigned);


// Compute the total lengths of branches to estimate
double get_total_blens(gsl_vector* blens, int num_branch){
    double total_blens = 0;
    for(int i = 0; i < num_branch; i++){
        total_blens += gsl_vector_get(blens, i);
    }
    // cout << "Total branch length is: " << total_blens << endl;
    return total_blens;
}


/*
 Tree height: unifrom(lbe-fbe, lbe-dob)
 Prior on both tree topology and node times
*/
double get_prior_tree_height(vector<double> prior_parameters)
{
    double min_height = prior_parameters[0];
    double max_height = prior_parameters[1];
    double log_prior_on_TH = log(1/(max_height-min_height));

    return log_prior_on_TH;
}

/*
 coalescent tree prior on both tree topology and node times
*/
double get_prior_pop_size(int pop_size)
{
    double val = 1/pop_size;
    double log_prior = log(1/val);

    return log_prior;
}

/*
 exponential prior for one branch length (branch connecting MRCA and LUCA)
*/
double get_prior_indivual_blen(double blen, double lambda)
{
    double log_prior = log(lambda) - lambda * blen;

    return log_prior;
}


/*
 The log prior probability of branch length
 Use compound Dirichlet distribution
 Namely, gama distribution on total tree length, Dirichlet prior on branch length proportions
*/
double get_prior_blen(gsl_vector* blens, int num_branch, vector<double> prior_parameters, vector<double> alphas)
{
    assert(prior_parameters.size() > 1);
    double a = prior_parameters[0];    // shape of Gamma prior on total_blens
    double b = prior_parameters[1];    // scale of Gamma prior on total_blens

    // Calculate Gamma prior on tree length (total_blens)
    double total_blens = get_total_blens(blens, num_branch);
    double log_gamma_prior_on_TL = (a - 1.0) * log(total_blens) - total_blens/b - a*log(b) - lgamma(a);
    double log_edge_length_proportions_prior = 0;
    // Calculate Dirichlet prior on edge length proportions
    //
    // Note that, for n edges, the Dirichlet prior density is
    //
    // p1^{a1-1} p2^{a2-1} ... pn^{an-1}
    // ------------------------------
    //    n*Gamma(c) / Gamma(n*c)
    //
    // where n = num_branch, pk = edge length k / total_blens and Gamma is the Gamma function.
    // If c == 1, then both numerator and denominator equal 1, so it is pointless
    // do loop over edge lengths.
    if(prior_parameters.size() > 2){
        double c = prior_parameters[2];    // parameter of Dirichlet prior on edge length proportions
        log_edge_length_proportions_prior = lgamma(num_branch*c);
        if(c != 1.0)
        {
            for(int i = 0; i < num_branch; i++)
            {
                double edge_length_proportion = gsl_vector_get(blens, i)/total_blens;
                log_edge_length_proportions_prior += (c - 1.0) * log(edge_length_proportion);
            }
            log_edge_length_proportions_prior -= lgamma(c) * num_branch;
        }
    }
    else{
        assert(alphas.size()==num_branch);
        double sum_alpha = 0;
        // cout << "   each alpha " << endl;
        for(int i = 0; i < num_branch; i++)
        {
            double alpha = alphas[i];
            sum_alpha += alpha;
            log_edge_length_proportions_prior -= lgamma(alpha);
            double edge_length_proportion = gsl_vector_get(blens, i)/total_blens;
            log_edge_length_proportions_prior += (alpha - 1.0) * log(edge_length_proportion);
            // cout << "\t" << i  << "\t" << alpha << "\t" << lgamma(alpha) <<  "\t" << (alpha - 1.0) * log(edge_length_proportion) << endl;
        }
        // if(debug){
        //     cout << "   sum of alpha " << sum_alpha << "\t" << lgamma(sum_alpha) << endl;
        //     cout << "   log_edge_length_proportions_prior " << log_edge_length_proportions_prior << endl;
        // }
        log_edge_length_proportions_prior += lgamma(sum_alpha);
    }

    double log_prior = log_gamma_prior_on_TL + log_edge_length_proportions_prior;
    return log_prior;
}

// Assume a uniform prior: all labeled rooted topologies are equally likely
double get_prior_topology(int Ns){
    int  n = Ns + 1;    // number of tips
    // number of rooted labeled bifurcating trees
    // double log_num_topologies = lgamma(2.0*n - 3.0 + 1.0) - (n - 2.0)*log(2.0) - lgamma(n - 2.0 + 1.0);
    // number of possible labeled histories
    double log_num_topologies = lgamma(n + 1.0) + lgamma(n - 1.0 + 1.0) - (n - 1.0)*log(2.0);
    return -log_num_topologies;
}

double get_prior_mutation_gamma(double mu, vector<double> prior_parameters){
    double a = prior_parameters[0];    // shape of Gamma prior on mutation rate
    double b = prior_parameters[1];    // scale of Gamma prior on mutation rate
    // Calculate Gamma prior on mutation rate
    double log_mut = (a - 1.0) * log(mu) - mu/b - a*log(b) - lgamma(a);
    return log_mut;
}


/*---------------------------------------------------------------------------------
|
|   CdfNormal
|
|   Calculates the cumulative density distribution (CDF) for the normal using:
|
|   Hill, I. D.  1973.  The normal integral.  Applied Statistics, 22:424-427.
|      (AS66)
|
---------------------------------------------------------------------------------*/
// Extracted from MrBayes code in utils.c
double CdfNormal( double x)
{
    int             invers = 0;
     double          p, limit = 10.0, t = 1.28, y = x*x/2.0;

    if(x < 0.0)
        {
        invers = 1;
        x  *= -1.0;
        }
    if(x > limit)
        return (invers ? 0 : 1);
    if(x < t)
        p = 0.5 - x * (0.398942280444 - 0.399903438504 * y /
            (y + 5.75885480458 - 29.8213557808 /
            (y + 2.62433121679 + 48.6959930692 /
            (y + 5.92885724438))));
    else
        p = 0.398942280385 * exp(-y) /
            (x - 3.8052e-8 + 1.00000615302 /
            (x + 3.98064794e-4 + 1.98615381364 /
            (x - 0.151679116635 + 5.29330324926 /
            (x + 4.8385912808 - 15.1508972451 /
            (x + 0.742380924027 + 30.789933034 /
            (x + 3.99019417011))))));

    return (invers ? p : 1-p);
}

// normal for log10 mu
double get_prior_mutation_lnormal(double lmu, vector<double> prior_parameters){
    double mu = prior_parameters[0];
    double sigma = prior_parameters[1];
    // Calculate Normal prior on log of mutation rate
    // double log_mut =  -(0.5* log(2*M_PI) +  log(sigma)) - pow(lmu - mu,2) / (2*sigma*sigma);

    // Calculate truncated Normal prior on log of mutation rate
    double z = (lmu - mu) / sigma;
    double z_0 = (RATE_MIN_LOG - mu) /sigma;
    double z_1 = (RATE_MAX_LOG - mu) /sigma;
    double normConst_0 = CdfNormal(z_0);
    double normConst_1 = CdfNormal(z_1);

    double log_mut = -(0.5* log(2*M_PI) +  log(sigma)) - z * z / 2.0 - log(normConst_1 - normConst_0);

    return log_mut;
}


gsl_vector* initialize_branch_length(evo_tree& rtree, int cons=0)
{
    int npar;
    gsl_vector *x;
    double blen;
    if(!cons){
      npar = rtree.edges.size() - 1;
      // initialise the best guess branch length and mu if required
      x = gsl_vector_alloc (npar);
      for(int i = 0; i < npar; ++i){
          blen = rtree.edges[i].length;
          gsl_vector_set (x, i, blen);
      }
    }
    else{   // only estimate internal edges for all tumor samples
      npar = rtree.nleaf - 2;;
      x = gsl_vector_alloc (npar);
      // initialise with internal edges
      if(debug){
          cout << "Initialise with internal edges" << endl;
      }
      vector<edge*> intedges = rtree.get_internal_edges();
      for(int i = 0; i < npar; ++i){
          blen = intedges[i]->length;
          gsl_vector_set (x, i, blen);
      }
    }
    // double num_branch = x -> size;
    // cout << "There are " << num_branch <<  " parameters for estimating branch lengths " << endl;
    return x;
}


evo_tree create_new_coal_tree(int nsample, int epop){
  //cout << "GENERATING COAL TREE" << endl;
  vector<int> edges;
  vector<double> lengths;
  vector<double> epoch_times;
  vector<double> times;

  int nlin = nsample;
  vector<int> nodes;

  // For compatibility with ape in R, root node must be labelled
  // 0 to nsample-1 are leaf, nsample is germline, nsample+1 is root
  // all other nodes ids start from here
  int node_count = nsample+2;

  // create leaf nodes
  for(int i = 0; i < nsample; ++i) nodes.push_back(i);

  // create vector of event times
  // total nodes = 2*(nsample+1) -1
  for(int i = 0; i < (2*nsample+1); ++i) times.push_back(0.0);

  double t_tot = 0;
  while(nlin > 1){
    // sample a time from Exp( combinations(k,2) )
    double lambda = fact(nlin)/( 2*fact(nlin-2) );
    double t = gsl_ran_exponential(r, 1/lambda) * epop;
    t_tot += t ;

    // choose two random nodes from available list
    random_shuffle(nodes.begin(), nodes.end(), fp_myrng);

    // edge node_count -> node
    edges.push_back(node_count);
    edges.push_back(nodes[nodes.size()-1]);
    lengths.push_back(t_tot - times[ nodes[nodes.size()-1] ] );

    edges.push_back(node_count);
    edges.push_back(nodes[nodes.size()-2]);
    lengths.push_back(t_tot - times[ nodes[nodes.size()-2] ] );

    // update time for this node
    //cout << "\t" << node_count+1 << "\t" << t_tot << endl;
    epoch_times.push_back(t_tot);
    times[ node_count ] = t_tot;

    nodes.pop_back();
    nodes.pop_back();

    nodes.push_back(node_count);
    node_count++;
    nlin--;
  }

  // create the root and germline nodes and edges
  double lambda = 1;
  double t = gsl_ran_exponential(r, 1/lambda) * epop;
  t_tot += t;
  epoch_times.push_back(t_tot);

  // add in the time for the root and germline nodes
  times[nsample+1] = t_tot;
  times[nsample] = t_tot;

  edges.push_back(nsample+1);
  edges.push_back(node_count-1);
  lengths.push_back(t);

  edges.push_back(nsample+1);
  edges.push_back(nsample);
  lengths.push_back(0);

  // invert the times
  for(int i = 0; i < times.size(); ++i){
    times[i] = t_tot - times[i];
  }

  // cout << "total time of tree: " << t_tot << " : ";
  // for(int i = 0; i < epoch_times.size(); ++i) cout << "\t" << epoch_times[i];
  // cout << endl;
  //
  // // invert the times
  // cout << "times of nodes:" << endl;
  // for(int i = 0; i < times.size(); ++i){
  //   cout << i+1 << "\t" << times[i] << endl;
  // }

  evo_tree ret(nsample+1, edges, lengths);
  return ret;
}


// Use NNI for clock trees
int move_topology(double& log_hastings_ratio, evo_tree& rtree){
    // int debug = 0;
    if(debug){
        printf ("Before:\n");
        rtree.print();
        // Output to a file to facilitate ploting
        ofstream out_tree("./sim-data-mcmc-tree-before.txt");
        rtree.write(out_tree);
        out_tree.close();
    }

    /* pick an interior branch (excluding the top edge connecting to root Ns+1), around which it is possible to make an NNI */
    double tv;
    double tc;
    int i, u, v, c, a, b;
    edge *e_uv, *e_uc, *e_va;
    vector<int> children;
    double rnum;
    int count = 0;
    double blen_uc;
    double blen_va;
    // const int COUNT_SAMPLE = 20;
    tc = -2;
    tv = -1;
    set<double> chosen; // Stop if all the internal edges have been tested

    vector<edge*> intedges = rtree.get_internal_edges();
    int nintedge = intedges.size();
    while (tc <= tv && chosen.size() < nintedge - 1){
        count++;
        // Exclude the edges connecting to unaltered genome
        assert(intedges[nintedge-1]->start == Ns + 1);
        i = gsl_rng_uniform_int(r, nintedge-1);
        chosen.insert(i);
        if(debug){
            cout << "Picking " << i+1 << "th internal edge: " << intedges[i]->start+1 << ", " << intedges[i]->end+1 << endl;
        }
        e_uv = intedges[i];
        // Find the children of the start node
        u = e_uv->start;
        children = rtree.nodes[u].daughters;
        // Select the child with subtree as v
        v = *max_element(children.begin(), children.end());
        c = *min_element(children.begin(), children.end());

        assert(v > Ns+1);

        tv = rtree.nodes[v].time;
        tc = rtree.nodes[c].time;

        // Find the children of v
        children = rtree.nodes[v].daughters;
        assert(children.size()==2);
        rnum = runiform(r, 0, 1);
        if(rnum > 0.5){
            a = children[0];
        }
        else{
            a = children[1];
        }
    }

    if(tc <= tv){
        if(debug) cout << "Cannot update topology!" << endl;
        return 0;
    }
    if(debug){
        cout << "Picked an edge after " << count << " iteration(s)" << endl;
        cout << "u, v, a, c " << u+1 << "\t" << v+1 << "\t" << a+1 << "\t" << c+1 << endl;
    }

    /* record branch lengths */
    // oldALength = e_va->length;
    // oldCLength = e_uc->length;

    /* make topology change */
    if(debug){
        cout << "Making topology change" << endl;
    }
    // (u,c) --> (u,a); (v,a) --> (v,c)
    // move the subtree below
    // Find the length of edge uc and va
    for(int j=0; j < rtree.edges.size(); j++){
        edge *e = &rtree.edges[j];
        if(e->start == u && e->end == c){
            blen_uc = e->length;
            e->end = a;
            e_uc = e;
        }
        if(e->start == v && e->end == a){
            blen_va = e->length;
            e->end = c;
            e_va = e;
        }
    }

    // update branch length to keep node ages/times unchanged
    e_va->length = blen_uc - e_uv->length;
    e_uc->length = e_uv->length + blen_va;
    adjust_blen(e_va->length, BLEN_MIN, BLEN_MAX);
    adjust_blen(e_uc->length, BLEN_MIN, BLEN_MAX);

    rtree.generate_nodes();

    /* adjust branch lengths */
    // e_uc->length = rtree.node_times[a] - rtree.node_times[u];
    // e_va->length = rtree.node_times[c] - rtree.node_times[v];

    if(debug){
        printf("After:\n");
        rtree.print();
        // Output to a file to facilitate ploting
        ofstream out_tree("./test/sim-data-mcmc-tree-after.txt");
        rtree.write(out_tree);
        out_tree.close();
    }
    return 1;
}


// use normal proposal
double move_mrates_normal(double& log_hastings_ratio, double prev_mu, double sigma){
    // double a = RATE_MIN;
    // double b = RATE_MAX;
    double a = RATE_MIN_LOG;
    double b = RATE_MAX_LOG;
    double nx =  prev_mu + gsl_ran_gaussian(r, sigma);
    double n, e;
    if(nx <= a){
        e = a - nx;
        n = floor(e/(b-a));
        nx = a + (e - n*(b-a));
    }
    if(nx >= b){
        e = nx - b;
        n = floor(e/(b-a));
        nx = b - (e - n*(b-a));
    }
    assert(nx > a && nx < b);
    return nx;
}

// use normal proposal
double move_normal(double& log_hastings_ratio, double prev, double a, double b, double sigma){
    double nx =  prev + gsl_ran_gaussian(r, sigma);
    double n, e;
    if(nx <= a){
        e = a - nx;
        n = floor(e/(b-a));
        nx = a + (e - n*(b-a));
    }
    if(nx >= b){
        e = nx - b;
        n = floor(e/(b-a));
        nx = b - (e - n*(b-a));
    }
    assert(nx > a && nx < b);
    return nx;
}


// Propose a new branch length with multiplier proposal
gsl_vector* move_blens_multiplier(double& log_hastings_ratio, double a, double b, double lambda, gsl_vector* prev_blens, int num_branch, int i){
    if(debug) cout << "Update branch length " << i << endl;
    gsl_vector *blens = gsl_vector_alloc(num_branch);
    for (int i = 0; i < num_branch; i++){
        gsl_vector_set(blens, i, gsl_vector_get(prev_blens, i));
    }
    double x, nx, y, ny, c;

    x = gsl_vector_get(prev_blens, i);
    y = (b - x)/(x - a);
    c = exp(lambda * (runiform(r, 0, 1) - 0.5));
    ny = y * c;
    nx = (b + a * ny)/(ny + 1);
    // adjust_blen(nx, a, b);
    gsl_vector_set(blens, i, nx);

    if(debug){
        cout << "old x y:" << x << "," << y << "\tnew x y:" << nx << "," << ny << endl;
    }

    log_hastings_ratio += log(c) + 2 * log(abs((nx-a)/(x-a)));

    return blens;
}


// Propose a new branch length with multiplier proposal for all branches
gsl_vector* move_blens_multiplier_all(double& log_hastings_ratio, double a, double b, double lambda, gsl_vector* prev_blens, int num_branch){
    // cout << "Update branch length" << endl;
    gsl_vector *blens = gsl_vector_alloc(num_branch);
    double x, nx, y, ny, c;

    for (int i = 0; i < num_branch; i++){
        x = gsl_vector_get(prev_blens, i);
        y = (b - x)/(x - a);
        c = exp(lambda * (runiform(r, 0, 1) - 0.5));
        ny = y * c;
        nx = (b + a * ny)/(ny + 1);
        // adjust_blen(nx, a, b);
        gsl_vector_set(blens, i, nx);

        if(debug){
            cout << "branch " << i << "\told x y:" << x << "," << y << "\tnew x y:" << nx << "," << ny << endl;
        }
        log_hastings_ratio += log(c) + 2 * log(abs((nx-a)/(x-a)));
    }
    return blens;
}


// Propose a new branch length with normal proposal
gsl_vector* move_blens_normal(double& log_hastings_ratio, double sigma, gsl_vector* prev_blens, int num_branch){
    // cout << "Update branch length" << endl;
    gsl_vector *blens = gsl_vector_alloc(num_branch);
    for (int i = 0; i < num_branch; i++){
        double x = gsl_vector_get(prev_blens, i);
        double nx =  x + gsl_ran_gaussian(r, sigma);
        adjust_blen(nx, BLEN_MIN, BLEN_MAX);
        gsl_vector_set(blens, i, nx);
    }
    return blens;
}

// scale the whole tree by applying a multiplier to all branch lengths
gsl_vector* move_blens_scale(double& log_hastings_ratio, double lambda, gsl_vector* prev_blens, int num_branch){
    // cout << "Update branch length" << endl;
    gsl_vector *blens = gsl_vector_alloc(num_branch);
    double m = exp(lambda * (runiform(r, 0, 1) - 0.5));
    // cout << "Multiplier is " << m << endl;
    log_hastings_ratio = log(m) * num_branch;
    // Update all branches one by one
    for(int j=0; j<num_branch; j++){
        gsl_vector_set(blens, j, gsl_vector_get(prev_blens, j) * m);
    }
    return blens;
}

// scale a single branch by applying a multiplier
double move_blens_single(double& log_hastings_ratio, double lambda, double prev_blen){
    double m = exp(lambda * (runiform(r, 0, 1) - 0.5));
    // cout << "Multiplier is " << m << endl;
    log_hastings_ratio = log(m);
    double  new_blen =  prev_blen * m;
    return new_blen;
}

// Update branch lengths with Bactrian proposal
gsl_vector* move_blens_bactrian(double& log_hastings_ratio, double m, double sigma, gsl_vector* prev_blens, int num_branch){
    // cout << "Update branch length" << endl;
    gsl_vector *blens = gsl_vector_alloc(num_branch);
    double x, nx, y, z, p;
    // Update all branches at once
    for(int j=0; j<num_branch; j++){
        x = gsl_vector_get(prev_blens, j);
        z = gsl_ran_gaussian(r, 1);
        p = runiform(r, 0, 1);
        if(p < 0.5){
            y = z * sqrt(1 - m * m) + m;
        }
        else{
            y = z * sqrt(1 - m * m) - m;
        }
        nx = x + y * sigma;
        adjust_blen(nx, BLEN_MIN, BLEN_MAX);
        gsl_vector_set(blens, j, nx);
    }
    return blens;
}


void update_topology(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, double lambda_topl, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1){
    double prev_log_prior = 1, log_prior = 1, prev_log_likelihood;
    double log_hastings_ratio = 0;
    bool accept = false;
    evo_tree ntree;
    // int debug = 0;

    double rn = runiform(r, 0, 1);
    // Change topology with probability lambda_topl
    if(rn > lambda_topl){
        cout << "Not updating tree topology this time" << endl;
        return;
    }

    if(debug){
        cout << "Updating tree topology " << endl;
    }

    // prev_log_prior = get_prior_topology(Ns);
    if(sample_prior){
        prev_log_likelihood = 1;
    }else{
      // if(!is_tree_valid(rtree, tobs, age, cons)){
      //    prev_log_likelihood = SMALL_LNL;
      // }else{
        // rtree.print();
      if(model == DECOMP){
          prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
      }else{
          prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
      }
      // }
    }

    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    // rtree will be changed after move
    ntree = evo_tree(rtree);
    int changed = move_topology(log_hastings_ratio, ntree);
    if(!changed){
        if(debug) cout << "Canot propose a new tree topology this time" << endl;
        if(n_draw > n_burnin)  nrejects++;
        return;
    }
    if(sample_prior){
        log_likelihood = 1;
    }else{
        // ntree.print();
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    // log_prior = get_prior_topology(Ns);
    if(debug){
        cout << "   log hastings ratio of topology proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    // cout << "Accept or reject the proposal" << endl;
    double lnl_ratio = log_likelihood - prev_log_likelihood;
    if(lnl_ratio == 0) cout << "Same likelihood after proposal!" << endl;

    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;

    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);

    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;

    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    // cout << "Print the accepted state" << endl;
    if(accept)
    {
        if(debug) cout << "accept tree topolgy" << endl;
        rtree = evo_tree(ntree);
        if(n_draw > n_burnin){
            naccepts++;
        }
    }
    else{
        if(debug) cout << "reject tree topolgy" << endl;
        if(n_draw > n_burnin)  nrejects++;
    }
}


void update_blen(evo_tree& rtree, int branch_i, int model, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_blen, vector<double> alphas, double lambda, double lambda_all, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    gsl_vector *prev_blens, *blens;
    int num_branch;
    double prev_log_prior, prev_log_likelihood, log_prior;
    double log_hastings_ratio = 0;
    bool accept = false;
    evo_tree ntree;
    // int debug = 0;

    if(debug){
        cout << "Updating branch lengths " << endl;
    }

    prev_blens = initialize_branch_length(rtree, cons);
    num_branch = prev_blens->size;
    if(debug){
        cout << "Previous branch length: " << endl;
        for(int k = 0; k < num_branch; k++){
            cout << "\t" << gsl_vector_get(prev_blens, k);
        }
        cout << endl;
        // cout << "There are " << num_branch <<  " branches to estimate" << endl;
    }
    prev_log_prior = get_prior_blen(prev_blens, num_branch, prior_parameters_blen, alphas);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double max_blen = BLEN_MAX;
    if(cons){
        max_blen = get_total_time(rtree.get_node_times(), max_tobs);
    }
    // cout << "Max branch length to propose " << max_blen << endl;
    if(branch_i == -1){
        // cout << "update all branches " << endl;
        blens = move_blens_multiplier_all(log_hastings_ratio, BLEN_MIN, max_blen, lambda_all, prev_blens, num_branch);
    }else{
        blens = move_blens_multiplier(log_hastings_ratio, BLEN_MIN, max_blen, lambda, prev_blens, num_branch, branch_i);
    }
    // blens = move_blens_bactrian(log_hastings_ratio, lambda_all, sigma, prev_blens, num_branch);
    // gsl_vector* blens = move_blens_normal(sigma, prev_blens, num_branch, log_hastings_ratio);

    // Create a new tree with the proposed tree length
    if(debug){
        cout << "Proposed branch length: " << endl;
        for(int k = 0; k < num_branch; k++){
            cout << "\t" << gsl_vector_get(blens, k);
        }
        cout << endl;
        cout << "   Create a new tree with the proposed tree length" << endl;
        // cout << "number of estimated branches " << blens->size << endl;
    }

    ntree = create_new_tree(blens, rtree, max_tobs, cons);
    if(cons){
        adjust_tree_blens(ntree);
        adjust_tree_tips(ntree, tobs, age);
    }

    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    log_prior = get_prior_blen(blens, num_branch, prior_parameters_blen, alphas);
    if(debug){
        cout << "   log hastings ratio of branch length proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    // cout << "Accept or reject the proposal" << endl;
    accept = false;
    double lnl_ratio = log_likelihood - prev_log_likelihood;
    if(lnl_ratio == 0) cout << "Same likelihood after proposal!" << endl;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    // cout << "Print the accepted state" << endl;
    if(accept)
    {
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        if(n_draw > n_burnin)  nrejects++;
        gsl_vector_memcpy(blens, prev_blens);
    }
}


void update_tree_height(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters, double sigma, int sample_prior, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    // uniform prior
    double prev_log_prior = 1, prev_log_likelihood, log_prior = 1;
    double log_hastings_ratio = 0;
    bool accept = false;
    double old_val, new_val;
    evo_tree ntree;

    double min_height = prior_parameters[0];
    double max_height = prior_parameters[1];

    if(debug){
        cout << "Updating tree height " << endl;
    }

    old_val = get_tree_height(rtree.get_node_times());
    // prev_log_prior = get_prior_tree_height(prior_parameters);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    new_val = move_normal(log_hastings_ratio, old_val, min_height, max_height, sigma);
    ntree = evo_tree(rtree);
    // Set the new tree height
    // ntree.tree_height = new_val;
    // ntree.total_time = new_val - *max_element(rtree.tobs.begin(), rtree.tobs.end());
    // Rescale the tree
    double ratio =  new_val/old_val;
    ntree.scale_time(ratio);
    // // Update terminal branches
    // for(int i = 0; i < ntree.nleaf-1; ++i){
    //   vector<int> es = ntree.get_ancestral_edges( ntree.nodes[i].id );
    //   reverse(es.begin(),es.end());
    //   ntree.edges[ es.back() ].length = ntree.total_time + ntree.tobs[ ntree.nodes[i].id ];
    //   for(int j=0; j<es.size()-1; ++j){
    //     ntree.edges[ es.back() ].length -= ntree.edges[es[j]].length;
    //     if(ntree.edges[ es.back() ].length < 0){
    //         cout << "Cannot update tree height this time" << endl;
    //         return;
    //     }
    //   }
    // }
    // ntree.calculate_node_times();
    // ntree.get_internal_edges();
    // ntree.lengths.clear();

    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    // log_prior = get_prior_tree_height(prior_parameters);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    double lnl_ratio = log_likelihood - prev_log_likelihood;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    if(accept)
    {
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        if(n_draw > n_burnin)  nrejects++;
    }
}


// Update the effective population size for coalescent model
void update_pop_size(evo_tree& rtree, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    // uniform prior
    double prev_log_prior = 1, prev_log_likelihood, log_prior = 1;
    double log_hastings_ratio = 0;
    bool accept = false;
    double old_val, new_val;
    evo_tree ntree;

    if(debug){
        cout << "Updating model parameters " << endl;
    }

    old_val = get_tree_height(rtree.get_node_times());
    // prev_log_prior = get_prior_tree_height(prior_parameters);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    new_val = move_normal(log_hastings_ratio, old_val, 1, 100, sigma);
    ntree = evo_tree(rtree);
    // Set the new tree height
    double ratio = new_val/old_val;
    ntree.scale_time(ratio);

    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    // log_prior = get_prior_tree_height(prior_parameters);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    double lnl_ratio = log_likelihood - prev_log_likelihood;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    if(accept)
    {
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        if(n_draw > n_burnin)  nrejects++;
    }
}

void update_mutation_rates(evo_tree& rtree, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_rate, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    double prev_log_prior, prev_log_likelihood, log_prior;
    double log_hastings_ratio = 0;
    bool accept = false;
    double new_mu;
    evo_tree ntree;

    if(debug){
        cout << "Updating model parameters " << endl;
    }

    assert(rtree.mu > 0);
    prev_log_prior = get_prior_mutation_gamma(rtree.mu, prior_parameters_rate);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    new_mu = move_mrates_normal(log_hastings_ratio, rtree.mu, sigma);
    ntree = evo_tree(rtree);
    ntree.mu = new_mu;
    // if(debug){
    //     cout << "Old mu " << rtree.mu << endl;
    //     cout << "epopw mu " << ntree.mu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    log_prior = get_prior_mutation_gamma(ntree.mu, prior_parameters_rate);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    double lnl_ratio = log_likelihood - prev_log_likelihood;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    if(accept)
    {
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        if(n_draw > n_burnin)  nrejects++;
    }
}

//
void update_mutation_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    double prev_log_prior, prev_log_likelihood, log_prior;
    double log_hastings_ratio = 0.0;
    bool accept = false;
    evo_tree ntree;

    if(debug){
        cout << "Updating model parameters " << endl;
    }

    assert(rtree.mu > 0);
    double lmu = log10(rtree.mu);
    prev_log_prior = get_prior_mutation_lnormal(lmu, prior_parameters_mut);
    if(sample_prior){
        prev_log_likelihood = 1;
    }else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.mu = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old mu " << rtree.mu << "\t" << lmu << endl;
    //     cout << "   epopw mu " << ntree.mu << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    double lnl_ratio = log_likelihood - prev_log_likelihood;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    if(accept){
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }else{
        if(n_draw > n_burnin)  nrejects++;
    }
}



void update_deletion_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    double prev_log_prior, prev_log_likelihood, log_prior;
    double log_hastings_ratio = 0;
    bool accept = false;
    evo_tree ntree;

    if(debug){
        cout << "Updating model parameters " << endl;
    }

    assert(rtree.del_rate > 0);
    double lmu = log10(rtree.del_rate);
    prev_log_prior = get_prior_mutation_lnormal(lmu, prior_parameters_mut);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.del_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old deletion rate " << rtree.del_rate << "\t" << lmu << endl;
    //     cout << "   epopw deletion rate " << ntree.del_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    double lnl_ratio = log_likelihood - prev_log_likelihood;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    if(accept)
    {
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        if(n_draw > n_burnin)  nrejects++;
    }
}


void update_duplication_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    double prev_log_prior, prev_log_likelihood, log_prior;
    double log_hastings_ratio = 0;
    bool accept = false;
    evo_tree ntree;

    if(debug){
        cout << "Updating model parameters " << endl;
    }

    assert(rtree.dup_rate > 0);
    double lmu = log10(rtree.dup_rate);
    prev_log_prior = get_prior_mutation_lnormal(lmu, prior_parameters_mut);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.dup_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old duplication rate " << rtree.dup_rate << "\t" << lmu << endl;
    //     cout << "   epopw duplication rate " << ntree.dup_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    double lnl_ratio = log_likelihood - prev_log_likelihood;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    if(accept)
    {
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        if(n_draw > n_burnin)  nrejects++;
    }
}

void update_cgain_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    double prev_log_prior, prev_log_likelihood, log_prior;
    double log_hastings_ratio = 0;
    bool accept = false;
    evo_tree ntree;

    if(debug){
        cout << "Updating model parameters " << endl;
    }

    assert(rtree.chr_gain_rate > 0);
    double lmu = log10(rtree.chr_gain_rate);
    prev_log_prior = get_prior_mutation_lnormal(lmu, prior_parameters_mut);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.chr_gain_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old duplication rate " << rtree.dup_rate << "\t" << lmu << endl;
    //     cout << "   epopw duplication rate " << ntree.dup_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    double lnl_ratio = log_likelihood - prev_log_likelihood;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    if(accept)
    {
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        if(n_draw > n_burnin)  nrejects++;
    }
}


void update_closs_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    double prev_log_prior, prev_log_likelihood, log_prior;
    double log_hastings_ratio = 0;
    bool accept = false;
    evo_tree ntree;

    if(debug){
        cout << "Updating model parameters " << endl;
    }

    assert(rtree.chr_loss_rate > 0);
    double lmu = log10(rtree.chr_loss_rate);
    prev_log_prior = get_prior_mutation_lnormal(lmu, prior_parameters_mut);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.chr_loss_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old duplication rate " << rtree.dup_rate << "\t" << lmu << endl;
    //     cout << "   epopw duplication rate " << ntree.dup_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    double lnl_ratio = log_likelihood - prev_log_likelihood;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    if(accept)
    {
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        if(n_draw > n_burnin)  nrejects++;
    }
}


void update_wgd_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, int& nrejects, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias, int is_total=1) {
    double prev_log_prior, prev_log_likelihood, log_prior;
    double log_hastings_ratio = 0;
    bool accept = false;
    evo_tree ntree;

    if(debug){
        cout << "Updating model parameters " << endl;
    }

    assert(rtree.wgd_rate > 0);
    double lmu = log10(rtree.wgd_rate);
    prev_log_prior = get_prior_mutation_lnormal(lmu, prior_parameters_mut);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            prev_log_likelihood = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            prev_log_likelihood = get_likelihood_revised(rtree, vobs, lnl_type);
        }
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.wgd_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old duplication rate " << rtree.dup_rate << "\t" << lmu << endl;
    //     cout << "   epopw duplication rate " << ntree.dup_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        if(model == DECOMP){
            log_likelihood = get_likelihood_decomp(ntree, vobs, obs_decomp, comps, lnl_type);
        }else{
            log_likelihood = get_likelihood_revised(ntree, vobs, lnl_type);
        }
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    double lnl_ratio = log_likelihood - prev_log_likelihood;
    double prior_ratio = log_prior - prev_log_prior;
    double log_diff = log_hastings_ratio + lnl_ratio + prior_ratio;
    double accept_prob = 0;
    if(log_diff < -100.0)  accept_prob = 0.0;
    else if(log_diff > 0.0)    accept_prob = 1.0;
    else    accept_prob = exp(log_diff);
    if(debug) cout << "ratios are " << log_hastings_ratio << "\t" << lnl_ratio << "\t" << prior_ratio << "\t"  << log_diff << endl;
    if(runiform(r, 0, 1) < accept_prob){
        accept = true;
    }
    if(accept)
    {
        if(n_draw > n_burnin){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        if(n_draw > n_burnin)  nrejects++;
    }
}


// Given a tree, find the MAP estimation of the branch lengths (and optionally mu) assuming branch lengths are independent or constrained in time
void run_mcmc(evo_tree& rtree, int model, const int n_draws, const int n_burnin, const int n_gap, vector<double> proposal_parameters, vector<double> prior_parameters_blen, vector<double> prior_parameters_height, vector<double> alphas, vector<double> prior_parameters_mut, double lambda_topl, string ofile, string tfile, const ITREE_PARAM& itree_param, int sample_prior=0, int fix_topology=0, int cons=0, int maxj=0, int cn_max = 4, int only_seg = 1, int correct_bias=0, int is_total=1){
        ofstream fout_trace(ofile);
        ofstream fout_tree(tfile);

        int precision = 6;
        int naccepts_topology = 0, nrejects_topology = 0, nsel_topology = 0;
        int naccepts_blen = 0, nrejects_blen = 0, nsel_blen = 0;
        int naccepts_height = 0, nrejects_height = 0, nsel_height = 0;
        int naccepts_mrate = 0, nrejects_mrate = 0, nsel_mrate = 0;
        int naccepts_dup = 0, nrejects_dup = 0, nsel_dup = 0;
        int naccepts_del = 0, nrejects_del = 0, nsel_del = 0;
        int naccepts_gain = 0, nrejects_gain = 0, nsel_gain = 0;
        int naccepts_loss = 0, nrejects_loss = 0, nsel_loss = 0;
        int naccepts_wgd = 0, nrejects_wgd = 0, nsel_wgd = 0;
        int nedge = 2 * rtree.nleaf - 2;
        int nintedge = rtree.nleaf - 2;
        vector<int> naccepts_bli(nedge-1, 0), nrejects_bli(nedge-1, 0), nsel_bli(nedge-1, 0);
        vector<int> naccepts_bli_cons(nintedge, 0), nrejects_bli_cons(nintedge, 0), nsel_bli_cons(nintedge, 0);

        double log_likelihood;
        double mu_lmut, sigma_lmut, sigma_mut;
        double mu_ldup, sigma_ldup, sigma_dup;
        double mu_ldel, sigma_ldel, sigma_del;
        double mu_lgain, sigma_lgain, sigma_gain;
        double mu_lloss, sigma_lloss, sigma_loss;
        double mu_lwgd, sigma_lwgd, sigma_wgd;

        double lambda = proposal_parameters[0];     // multiplier proposal
        double lambda_all = proposal_parameters[1];     // multiplier proposal
        // double lambda_all = proposal_parameters[1];     // Bactrian proposal
        double sigma_blen = proposal_parameters[2];     // normal proposal for branch length
        double sigma_height;

        if(model == MK){
            mu_lmut = prior_parameters_mut[0];
            sigma_lmut = prior_parameters_mut[1];

            sigma_mut = proposal_parameters[3];
        }
        else{
            assert(prior_parameters_mut.size()>9);
            assert(proposal_parameters.size()>7);

            mu_ldup = prior_parameters_mut[0];
            sigma_ldup = prior_parameters_mut[1];

            mu_ldel = prior_parameters_mut[2];
            sigma_ldel = prior_parameters_mut[3];

            mu_lgain = prior_parameters_mut[4];
            sigma_lgain = prior_parameters_mut[5];

            mu_lloss = prior_parameters_mut[6];
            sigma_lloss = prior_parameters_mut[7];

            mu_lwgd = prior_parameters_mut[8];
            sigma_lwgd = prior_parameters_mut[9];

            sigma_dup = proposal_parameters[3];
            sigma_del = proposal_parameters[4];
            sigma_gain = proposal_parameters[5];
            sigma_loss = proposal_parameters[6];
            sigma_wgd = proposal_parameters[7];
        }
        if(cons){
            sigma_height = proposal_parameters[proposal_parameters.size()-1];
        }

        string header;

        if(model == 0){
            header = "state\tlnl";
            if(maxj == 1){
                header += "\tmu";
            }
        }
        else{
            header = "state\tlnl";
            if(maxj == 1){
                header += "\tdup_rate\tdel_rate\tgain_rate\tloss_rate\twgd_rate";
            }
        }

        if(cons){
            header += "\theight";
            for(int j=1; j<=nintedge; j++){
                header += "\tl" + to_string(j);
            }
        }else{
            for(int j=1; j<=nedge-1; j++){
                header += "\tl" + to_string(j);
            }
        }

        fout_trace << "# Parameters" << endl;   // Add one line on top for compatibility with MrBayes format
        fout_trace << header << endl;

        fout_tree << "#nexus" << endl;
        fout_tree << "begin trees;" << endl;

        // select each operator stochasticly.
        // possible operators when the tree is constrained and mutation rates are estimated
        // 0: topology, 1: all branch lengths, 2: some branch lengths, 3: one branch length; 4: mutation rate; 5: all rates in model 3; 6: duplication rate; 7: deletion rate; 8: chromosome gain rate; 9: chromosome loss rate; 10: WGD rate; 11: tree height
        int sel_type = 0;
        gsl_ran_discrete_t*  dis;
        double prob_move_topol = 0;
        double prob_move_blen = 0;
        double prob_move_rate = 0;
        for(int i = 1; i <= n_draws; ++i){
            // Randomly choose one type of operator: topology, branch length, mutation rate
            if(maxj){
                if(fix_topology){
                    prob_move_topol = 0;
                    prob_move_blen = 0.5;
                    prob_move_rate = 0.5;
                }
                else{
                    prob_move_topol = 1/3;
                    prob_move_blen = 1/3;
                    prob_move_rate = 1/3;
                }
            }else{
                if(fix_topology){
                    prob_move_topol = 0;
                    prob_move_blen = 1;
                    prob_move_rate = 0;
                }else{
                    prob_move_topol = 0.5;
                    prob_move_blen = 0.5;
                    prob_move_rate = 0;
                }
            }
            double probs[3] = {prob_move_topol, prob_move_blen, prob_move_rate};
            dis = gsl_ran_discrete_preproc(3, probs);
            sel_type = gsl_ran_discrete(r, dis);

            if(sel_type == 0){
                nsel_topology++;
                update_topology(rtree, model, log_likelihood, naccepts_topology, nrejects_topology, i, n_burnin, n_gap, lambda_topl, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
            }
            else if(sel_type == 1){
                if(cons){
                    double prob_move_height = 0.5;
                    double prob_move_blen_all = 0.5;
                    double probs[2] = {prob_move_height, prob_move_blen_all};
                    dis = gsl_ran_discrete_preproc(2, probs);
                    int blen_update_tpye = gsl_ran_discrete(r, dis);

                    if(blen_update_tpye==0){
                        nsel_height++;
                        update_tree_height(rtree, model, log_likelihood, naccepts_height, nrejects_height, i, n_burnin, n_gap, prior_parameters_height, sigma_height, sample_prior, cn_max, only_seg, correct_bias, is_total);
                    }else{
                        // either update one branch or all branches
                        int nblen = nintedge + 1;
                        double p = (double) 1/(nblen);
                        double probs[2] = {p, 1-p};
                        dis = gsl_ran_discrete_preproc(2, probs);
                        int blen_update = gsl_ran_discrete(r, dis);
                        if(blen_update == 0){
                            nsel_blen++;
                            update_blen(rtree, -1, model, log_likelihood, naccepts_blen, nrejects_blen, i, n_burnin, n_gap, prior_parameters_blen, alphas, lambda, lambda_all, sigma_blen, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
                        }else{
                            // randomly update one branch
                            int bli = gsl_rng_uniform_int(r, nintedge);
                            nsel_bli_cons[bli]++;
                            update_blen(rtree, bli, model, log_likelihood, naccepts_bli_cons[bli], nrejects_bli_cons[bli], i, n_burnin, n_gap, prior_parameters_blen, alphas, lambda, lambda_all, sigma_blen, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
                        }

                    }
                }else{
                    double p = (double) 1/(nedge);
                    double probs[2] = {p, 1-p};
                    dis = gsl_ran_discrete_preproc(2, probs);
                    int blen_update = gsl_ran_discrete(r, dis);
                    if(blen_update == 0){
                        nsel_blen++;
                        update_blen(rtree, -1, model, log_likelihood, naccepts_blen, nrejects_blen, i, n_burnin, n_gap, prior_parameters_blen, alphas, lambda, lambda_all, sigma_blen, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
                    }else{
                        // randomly update one branch
                        int bli = gsl_rng_uniform_int(r, nedge-1);
                        nsel_bli[bli]++;
                        update_blen(rtree, bli, model, log_likelihood, naccepts_bli[bli], nrejects_bli[bli], i, n_burnin, n_gap, prior_parameters_blen, alphas, lambda, lambda_all, sigma_blen, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
                    }

                }
            }else{
                if(model == MK){
                    nsel_mrate++;
                    vector<double> prior_parameters_mu({mu_lmut, sigma_lmut});
                    update_mutation_rates_lnormal(rtree, model, log_likelihood, naccepts_mrate, nrejects_mrate, i, n_burnin, n_gap, prior_parameters_mu, sigma_mut, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);

                }else{
                    double prob_move_dup = 0.2;
                    double prob_move_del = 0.2;
                    double prob_move_gain = 0.2;
                    double prob_move_loss = 0.2;
                    double prob_move_wgd = 0.2;
                    if(only_seg==1){
                        prob_move_gain = 0;
                        prob_move_loss = 0;
                        prob_move_wgd = 0;
                    }
                    // if(model == DECOMP){
                    //     if(max_wgd == 0) prob_move_wgd = 0;
                    //     if(max_chr_change == 0){
                    //         prob_move_gain = 0;
                    //         prob_move_loss = 0;
                    //     }
                    // }
                    double probs[5] = {prob_move_dup, prob_move_del, prob_move_gain, prob_move_loss, prob_move_wgd};
                    dis = gsl_ran_discrete_preproc(5, probs);
                    int rate_update_tpye = gsl_ran_discrete(r, dis);
                    switch(rate_update_tpye){
                        case 0:{
                            nsel_dup++;
                            vector<double> prior_parameters_dup({mu_ldup, sigma_ldup});
                            update_duplication_rates_lnormal(rtree, model, log_likelihood, naccepts_dup, nrejects_dup, i, n_burnin, n_gap, prior_parameters_dup, sigma_dup, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
                            break;
                        }

                        case 1:{
                            nsel_del++;
                            vector<double> prior_parameters_del({mu_ldel, sigma_ldel});
                            update_deletion_rates_lnormal(rtree, model, log_likelihood, naccepts_del, nrejects_del, i, n_burnin, n_gap, prior_parameters_del, sigma_del, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
                            break;
                        }

                        case 2:{
                            nsel_gain++;
                            vector<double> prior_parameters_gain({mu_lgain, sigma_lgain});
                            update_cgain_rates_lnormal(rtree, model, log_likelihood, naccepts_gain, nrejects_gain, i, n_burnin, n_gap, prior_parameters_gain, sigma_gain, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
                            break;
                        }

                        case 3:{
                            nsel_loss++;
                            vector<double> prior_parameters_loss({mu_lloss, sigma_lloss});
                            update_closs_rates_lnormal(rtree, model, log_likelihood, naccepts_loss, nrejects_loss, i, n_burnin, n_gap, prior_parameters_loss, sigma_loss, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
                            break;
                        }

                        case 4:{
                            nsel_wgd++;
                            vector<double> prior_parameters_wgd({mu_lwgd, sigma_lwgd});
                            update_wgd_rates_lnormal(rtree, model, log_likelihood, naccepts_wgd, nrejects_wgd, i, n_burnin, n_gap, prior_parameters_wgd, sigma_wgd, sample_prior, cons, cn_max, only_seg, correct_bias, is_total);
                            break;
                        }

                        default:{
                        cout << "Wrong mutation rate type!" << endl;
                        break;
                        }
                    }
                }
            }

            // Discard burn in samples
            if(i <= n_burnin){
                continue;
            }
            if(i % n_gap == 0){
                // print out the accepted proposal
                // string str_tree = order_tree_string(create_tree_string(rtree));
                // fout << i - n_burnin << "\t" << str_tree << "\t"  << log_likelihood << "\t" << rtree.mu ;
                // fout_trace << (i - n_burnin)/n_gap << "\t"  << log_likelihood;
                fout_trace << i << "\t"  << log_likelihood;

                if(maxj){
                    if(model == MK){
                        fout_trace << "\t" << rtree.mu ;
                    }else{
                        if(!only_seg){
                            fout_trace << "\t" << rtree.dup_rate << "\t" << rtree.del_rate << "\t" << rtree.chr_gain_rate << "\t" << rtree.chr_loss_rate << "\t" << rtree.wgd_rate;
                        }else{
                            fout_trace << "\t" << rtree.dup_rate << "\t" << rtree.del_rate;
                        }
                    }
                }

                // print out the branch lengths
                if(cons){
                    vector<edge*> intedges = rtree.get_internal_edges();
                    fout_trace << "\t" << get_tree_height(rtree.get_node_times());
                    for(int k = 0; k < intedges.size(); ++k){
                        fout_trace << "\t" << intedges[k]->length;
                    }
                }else{
                    for(int k = 0; k < nedge-1; k++){
                        fout_trace << "\t" << rtree.edges[k].length;
                    }
                }
                fout_trace << endl;

                string newick = rtree.make_newick(precision);
                // fout_tree << "tree " << (i - n_burnin)/n_gap << " = " << newick << ";" << endl;
                fout_tree << "tree " << i << " = " << newick << ";" << endl;
            }
        }

        fout_tree << "end;" << endl;

        fout_trace.close();
        fout_tree.close();

        // double n_keep = (n_draws - n_burnin)/n_gap;

        cout << "move\tTuning\t#accept\t#reject\t#nsel\tPr(m)\tPr(acc|m)\n";
        if(!fix_topology){
            double accept_rate_topology = (double) naccepts_topology / (nrejects_topology + naccepts_topology);
            double sel_rate_topology = (double) nsel_topology / n_draws;
            cout << "topology proposal (NNI - narrow)\t - \t" << naccepts_topology << "\t" << nrejects_topology << "\t" << nsel_topology << "\t" << sel_rate_topology << "\t" <<  accept_rate_topology << endl;
        }

        double accept_rate_blen = (double) naccepts_blen / (nrejects_blen + naccepts_blen);
        double sel_rate_blen = (double) nsel_blen / n_draws;
        cout << "branch length all " << lambda_all << "\t"  << naccepts_blen << "\t" << nrejects_blen << "\t"<< nsel_blen << "\t" << sel_rate_blen << "\t" <<  accept_rate_blen << endl;

        if(!cons){
            for(int j = 0; j < naccepts_bli.size(); j++){
                double accept_rate_blen = (double) naccepts_bli[j] / (nrejects_bli[j] + naccepts_bli[j]);
                double sel_rate_blen = (double) nsel_bli[j] / n_draws;
                cout << "branch length " << j+1 << "\t" << lambda << "\t"  << naccepts_bli[j] << "\t" << nrejects_bli[j] << "\t"<< nsel_blen << "\t" << sel_rate_blen << "\t" <<  accept_rate_blen << endl;
            }
        }
        else{
            double accept_rate_height = (double) naccepts_height / (nrejects_height + naccepts_height);
            double sel_rate_height = (double) nsel_height / n_draws;
            cout << "tree height\t" << sigma_height << "\t" << naccepts_height << "\t" << nrejects_height << "\t" << nsel_height << "\t"<< sel_rate_height << "\t" <<  accept_rate_height << endl;

            for(int j = 0; j < naccepts_bli_cons.size(); j++){
                double accept_rate_blen = (double) naccepts_bli_cons[j] / (nrejects_bli_cons[j] + naccepts_bli_cons[j]);
                double sel_rate_blen = (double) nsel_bli[j] / n_draws;
                cout << "branch length " << j+1 << "\t" << lambda << "\t"  << naccepts_bli_cons[j] << "\t" << nrejects_bli_cons[j] << "\t"<< nsel_blen << "\t" << sel_rate_blen << "\t" <<  accept_rate_blen << endl;
            }
        }

        if(maxj){
            if(model == MK){
                double accept_rate_mrate = (double) naccepts_mrate / (nrejects_mrate + naccepts_mrate);
                double sel_rate_mrate = (double) nsel_mrate / n_draws;
                cout << "mutation rate\t" << sigma_mut << "\t" << naccepts_mrate << "\t" << nrejects_mrate << "\t" << nsel_mrate << "\t" << sel_rate_mrate << "\t" <<  accept_rate_mrate << endl;
            }
            else{
                double accept_rate_dup = (double) naccepts_dup / (nrejects_dup + naccepts_dup);
                double sel_rate_dup = (double) nsel_dup / n_draws;
                cout << "duplication rate\t" << sigma_dup << "\t"  << naccepts_dup<< "\t" << nrejects_dup << "\t" << nsel_dup << "\t" << sel_rate_dup << "\t" <<  accept_rate_dup << endl;

                double accept_rate_del =(double) naccepts_del / (nrejects_del + naccepts_del);
                double sel_rate_del = (double) nsel_del / n_draws;
                cout << "deletion rate\t" << sigma_del << "\t"  << naccepts_del << "\t" << nrejects_del << "\t" << nsel_del << "\t" << sel_rate_del << "\t" <<  accept_rate_del << endl;

                if(!only_seg){
                    double accept_rate_gain = (double) naccepts_gain / (nrejects_gain + naccepts_gain);
                    double sel_rate_gain = (double) nsel_gain / n_draws;
                    cout << "chromosome gain rate\t" << sigma_gain << "\t"  << naccepts_gain << "\t" << nrejects_gain << "\t" << nsel_gain << "\t" << sel_rate_gain << "\t" <<  accept_rate_gain << endl;

                    double accept_rate_loss = (double) naccepts_loss / (nrejects_loss + naccepts_loss);
                    double sel_rate_loss = (double) nsel_loss / n_draws;
                    cout << "chromosome loss rate\t" << sigma_loss << "\t"  << naccepts_loss << "\t" << nrejects_loss << "\t" << nsel_loss << "\t" << sel_rate_loss << "\t" <<  accept_rate_loss << endl;

                    double accept_rate_wgd = (double) naccepts_wgd / (nrejects_wgd + naccepts_wgd);
                    double sel_rate_wgd = (double) nsel_wgd / n_draws;
                    cout << "whole genome doubling rate\t" << sigma_wgd << "\t" << naccepts_wgd << "\t" << nrejects_wgd << "\t" << nsel_wgd << "\t" << sel_rate_wgd << "\t" <<  accept_rate_wgd << endl;
                }
            }
        }

        cout << "\nTuning: The value of the operator’s tuning parameter, or ’-’ if the operator can’t be optimized." << endl;
        cout << "#accept: The total number of times a proposal by this operator has been accepted (excluding burnin samples)." << endl;
        cout << "#reject: The total number of times a proposal by this operator has been rejected (excluding burnin samples)." << endl;
        cout << "#nsel: The total number of times a proposal is selected." << endl;
        cout << "Pr(m): The probability this operator is chosen in a step of the MCMC (i.e. the normalized weight)." << endl;
        cout << "Pr(acc|m): The acceptance probability (#accept as a fraction of the total proposals for this operator, excluding burnin samples)." << endl;
}


// https://stackoverflow.com/questions/25444449/how-do-i-convert-a-stdstring-containing-doubles-to-a-vector-of-doubles
vector<double> get_vertex_indices(string const& pointLine)
{
  istringstream iss(pointLine);

  return vector<double>{
    istream_iterator<double>(iss),
    istream_iterator<double>()
  };
}


// Read parsimony trees built by other tools
evo_tree read_reference_tree(const string& tree_file, int Ns, const vector<double>& rates, const vector<double>& tobs){
    // int debug = 0;
    if(debug)   cout << tree_file << endl;
    evo_tree rtree = read_tree_info(tree_file, Ns);
    // for(int i = 0; i < tobs.size();i++){
    //     cout << tobs[i] << endl;
    // }
    // rtree.tobs = tobs;
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

    cout << "Tree height " << get_tree_height(rtree.get_node_times()) << endl;
    cout << "Total time " << get_total_time(rtree.get_node_times(), max_tobs) << endl;

    return rtree;
}


// Assign mutation rates to the tree
void revise_init_tree(evo_tree& rtree, const vector<double> rates, const vector<double>& tobs, int cons){
    if(rates.size() > 1){
      rtree.dup_rate = rates[0];
      rtree.del_rate = rates[1];
      if(rates.size() > 2){
          rtree.chr_gain_rate = rates[2];
          rtree.chr_loss_rate = rates[3];
          rtree.wgd_rate = rates[4];
      }
    }else{
      rtree.mu = rates[0];
    }

    if(cons){
        double min_height = *max_element(tobs.begin(), tobs.end());
        double max_height = age + min_height;
        adjust_tree_height(rtree, r, min_height, max_height);
        adjust_tree_tips(rtree, tobs, age);
        adjust_tree_blens(rtree);
    }
    // cout << "Total time of random tree " << rtree.total_time << endl;
}


void run_with_reference_tree(string rtreefile, int Ns, int Nchar, int num_invar_bins, int fix_topology, int model, int cons, int maxj, int cn_max, int only_seg, int correct_bias, int is_total, const vector<double>& ref_rates, const vector<double>& tobs, const vector<double>& rates, int n_draws, int n_burnin, int n_gap, const vector<double>& proposal_parameters, const vector<double>& prior_parameters_blen, const vector<double>& prior_parameters_height, const vector<double>& alphas, const vector<double>& prior_parameters_mut, double lambda_topl, string trace_param_file, string trace_tree_file, int sample_prior, const ITREE_PARAM& itree_param){
    // MLE testing
    // read in true tree
    evo_tree test_tree = read_reference_tree(rtreefile, Ns, ref_rates, tobs);
    double Ls = 0.0;
    Ls = get_likelihood_revised(test_tree, vobs, lnl_type);
    // cout << "\nOriginal tree likelihood: " << Ls << endl;
    cout << 0 << "\t" << Ls  << "\t" << test_tree.mu << "\t" << test_tree.dup_rate << "\t" << test_tree.del_rate;
    // Output the original tree length
    double avg_blen = 0;
    for(int i = 0; i < test_tree.edges.size(); ++i){
      avg_blen += test_tree.edges[i].length;
      cout << "\t" << test_tree.edges[i].length;
    }
    avg_blen = avg_blen / test_tree.edges.size();
    cout << "\t" << avg_blen;
    cout << endl;
    cout << "Internal EDGES:" << endl;
    cout << "\tid\tstart\tend\tlength" << endl;
    vector<edge*> intedges = test_tree.get_internal_edges();
    for(int i = 0; i < intedges.size(); ++i){
      cout << "\t" << intedges[i]->id+1 << "\t" << intedges[i]->start+1 << "\t" << intedges[i]->end+1 << "\t" << intedges[i]->length << endl;
    }

    double Lf = 0;
    stringstream sstm;
    ofstream out_tree;

    cout << "Generate the start tree" << endl;
    evo_tree rtree;
    if(fix_topology){
        cout << "   starting with a random tree with the same topology as the reference tree" << endl;
        int nedge = test_tree.edges.size();
        gsl_vector* rblens = gsl_vector_alloc(nedge);
        for(int i = 0; i < nedge; ++i){
          gsl_vector_set(rblens, i, runiform(r, 1, age));
        }
        rtree = create_new_tree(rblens, test_tree, max_tobs, 0);
    }
    else{
        cout << "   starting with a random coalescence tree" << endl;
        rtree = generate_coal_tree(Ns, r, fp_myrng, itree_param);
    }
    revise_init_tree(rtree, rates, tobs, cons);

    // sstm << "./test1/sim-data-" << cons << maxj << "-mcmc-tree-start.txt";
    // out_tree.open(sstm.str());
    // rtree.write(out_tree);
    // out_tree.close();
    // sstm.str("");

    Ls = get_likelihood_revised(rtree, vobs, lnl_type);
    // Ls = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, 0, cn_max, only_seg, correct_bias, is_total);
    cout << "\nRandom tree likelihood: " << Ls << endl;

    // Estimate branch length with MCMC
    cout << "\n\n### Running MCMC" << endl;
    run_mcmc(rtree, model, n_draws, n_burnin, n_gap, proposal_parameters, prior_parameters_blen, prior_parameters_height, alphas, prior_parameters_mut, lambda_topl, trace_param_file, trace_tree_file, itree_param, sample_prior, fix_topology, cons, maxj, cn_max, only_seg, correct_bias, is_total);
    // cout << "\nMinimised tree likelihood by MCMC / mu : " << Lf << "\t" << min_tree.mu*Nchar <<  endl;
}



int main (int argc, char ** const argv) {
    int n_draws, n_burnin, n_gap, clock;
    double lambda, sigma_blen, lambda_all, lambda_topl, sigma_height;
    double dirichlet_param, tlen_shape, tlen_scale;
    double rate_shape, rate_scale;
    string datafile, timefile, trace_param_file, trace_tree_file, rtreefile, config_file, seg_file;
    string dirichlet_alpha;
    int sample_prior, fix_topology;
    double rmu, rdup_rate, rdel_rate, rgain_rate, rloss_rate, rwgd_rate;
    double mu, dup_rate, del_rate, chr_gain_rate, chr_loss_rate, wgd_rate;
    double sigma_mut, sigma_dup, sigma_del, sigma_gain, sigma_loss, sigma_wgd;
    double sigma_lmut, sigma_ldup, sigma_ldel, sigma_lgain, sigma_lloss, sigma_lwgd;
    int is_bin, incl_all;
    int epop;
    double gtime, beta;
    // int is_stochastic;
    int init_tree;
    string file_itree;

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
            ("tfile,t", po::value<string>(&timefile)->default_value(""), "input time information file")
            ("config_file", po::value<string>(&config_file)->default_value(""), "configuration file of input parameters")
            ("Ns,s", po::value<int>(&Ns)->default_value(5), "number of samples or regions")
            ("cn_max", po::value<int>(&cn_max)->default_value(4), "maximum copy number of a segment")
            ("is_bin", po::value<int>(&is_bin)->default_value(1), "whether or not the input copy number is for each bin. If not, the input copy number is read as it is")
            ("incl_all", po::value<int>(&incl_all)->default_value(1), "whether or not to include all the input copy numbers without further propressing")
            ("is_total", po::value<int>(&is_total)->default_value(1), "whether or not the input is total copy number")

            // Parameters used for generating initial coalescence tree
            ("epop", po::value<int>(&epop)->default_value(2), "effective population size of cell populations")
            ("gtime", po::value<double>(&gtime)->default_value(1), "generation time in year")
            ("beta", po::value<double>(&beta)->default_value(0), "population growth rate")

            ("trace_param_file", po::value<string>(&trace_param_file)->default_value("trace-mcmc-params.txt"), "output trace file of parameter values")
            ("trace_tree_file", po::value<string>(&trace_tree_file)->default_value("trace-mcmc-trees.txt"), "output trace file of trees")
            ("seg_file", po::value<string>(&seg_file)->default_value(""), "output file with the postprocessed copy number matrix for tree building ")
            ("n_burnin,r", po::value<int>(&n_burnin)->default_value(9000), "number of burnin samples")
            ("n_draws,n", po::value<int>(&n_draws)->default_value(10000), "number of posterior draws to keep")
            ("n_gap", po::value<int>(&n_gap)->default_value(1), "sampling every kth samples ")
            ("sample_prior", po::value<int>(&sample_prior)->default_value(0), "whether or not to sample from the prior only")

            ("model,d", po::value<int>(&model)->default_value(2), "model of evolution (0: Mk, 1: one-step bounded (total), 2: one-step bounded (haplotype-specific), 3: independent Markov chains")
            ("clock", po::value<int>(&clock)->default_value(0), "model of molecular clock (0: strict global, 1: random local clock)")
            ("fix_topology", po::value<int>(&fix_topology)->default_value(0), "whether or not to fix the topology of the tree")
            ("cons", po::value<int>(&cons)->default_value(1), "constraints on branch length (0: none, 1: fixed total time)")
            ("maxj", po::value<int>(&maxj)->default_value(1), "estimation of mutation rate (0: mutation rate fixed to be the given value, 1: estimating mutation rate)")
            ("correct_bias", po::value<int>(&correct_bias)->default_value(1), "correct ascertainment bias")

            ("init_tree", po::value<int>(&init_tree)->default_value(0), "method of building inital tree for MCMC sampling (0: Random coalescence tree, 1: Provided tree, 2: Tree with the same topology as the true tree)")
            ("file_itree", po::value<string>(&file_itree)->default_value(""), "path of inital tree for MCMC sampling")

            ("m_max", po::value<int>(&m_max)->default_value(1), "maximum number of copies of a segment in a chromosome")
            ("max_wgd", po::value<int>(&max_wgd)->default_value(1), "maximum number of WGD")
            ("max_chr_change", po::value<int>(&max_chr_change)->default_value(1), "maximum number of chromosome changes")
            ("max_site_change", po::value<int>(&max_site_change)->default_value(2), "maximum number of segment changes")

            ("rtreefile", po::value<string>(&rtreefile)->default_value(""), "reference tree file")
            ("rmu", po::value<double>(&rmu)->default_value(0.02), "mutation rate of the reference tree")
            ("rdup_rate", po::value<double>(&rdup_rate)->default_value(0.01), "site duplication rate of the reference tree")
            ("rdel_rate", po::value<double>(&rdel_rate)->default_value(0.01), "site deletion rate of the reference tree")
            ("rgain_rate", po::value<double>(&rgain_rate)->default_value(0.01), "chromosome gain rate of the reference tree")
            ("rloss_rate", po::value<double>(&rloss_rate)->default_value(0.01), "chromosome loss rate of the reference tree")
            ("rwgd_rate", po::value<double>(&rwgd_rate)->default_value(0.01), "whole genome doubling rate of the reference tree")

            ("sigma_blen,g", po::value<double>(&sigma_blen)->default_value(2.5), "sigma for proposal of branch length")
            ("sigma_height", po::value<double>(&sigma_height)->default_value(2.5), "sigma for proposal of tree height")

            ("only_seg", po::value<int>(&only_seg)->default_value(0), "Whether or not to only consider segment-level mutations (0: include chromosome gain/loss and whole genome doubling, 1: only consider segment-level mutations)")
            ("mu,x", po::value<double>(&mu)->default_value(0.025), "mean of mutation rate")
            ("dup_rate", po::value<double>(&dup_rate)->default_value(0.01), "mean of site duplication rate")
            ("del_rate", po::value<double>(&del_rate)->default_value(0.01), "mean of site deletion rate")
            ("chr_gain_rate", po::value<double>(&chr_gain_rate)->default_value(0.001), "mean of chromosome gain rate")
            ("chr_loss_rate", po::value<double>(&chr_loss_rate)->default_value(0.001), "mean of chromosome loss rate")
            ("wgd_rate", po::value<double>(&wgd_rate)->default_value(0.0001), "mean of WGD (whole genome doubling) rate")

            ("sigma_mut", po::value<double>(&sigma_mut)->default_value(0.05), "sigma for proposal of mutation rate")
            ("sigma_dup", po::value<double>(&sigma_dup)->default_value(0.05), "sigma for proposal of site duplication rate")
            ("sigma_del", po::value<double>(&sigma_del)->default_value(0.05), "sigma for proposal of site deletion rate")
            ("sigma_gain", po::value<double>(&sigma_gain)->default_value(0.05), "sigma for proposal of chromosome gain rate")
            ("sigma_loss", po::value<double>(&sigma_loss)->default_value(0.05), "sigma for proposal of chromosome loss rate")
            ("sigma_wgd", po::value<double>(&sigma_wgd)->default_value(0.05), "sigma for proposal of whole genome doubling rate")

            ("lambda_topl", po::value<double>(&lambda_topl)->default_value(1), "lambda for topology proposal")
            ("lambda,l", po::value<double>(&lambda)->default_value(0.5), "lambda for multiplier proposal of individual branch length")
            ("lambda_all,m", po::value<double>(&lambda_all)->default_value(0.5), "lambda for multiplier proposal of all branch lengths")

            ("sigma_lmut", po::value<double>(&sigma_lmut)->default_value(0.05), "sigma for prior of log10 of mutation rate")
            ("sigma_ldup", po::value<double>(&sigma_ldup)->default_value(0.05), "sigma for prior of log10 of site duplication rate")
            ("sigma_ldel", po::value<double>(&sigma_ldel)->default_value(0.05), "sigma for prior of log10 of site deletion rate")
            ("sigma_lgain", po::value<double>(&sigma_lgain)->default_value(0.05), "sigma for prior of log10 of chromosome gain rate")
            ("sigma_lloss", po::value<double>(&sigma_lloss)->default_value(0.05), "sigma for prior of log10 of chromosome loss rate")
            ("sigma_lwgd", po::value<double>(&sigma_lwgd)->default_value(0.05), "sigma for prior of log10 of whole genome doubling rate")
            ("rate_shape", po::value<double>(&rate_shape)->default_value(1), "shape of gamma prior for mutation rate")
            ("rate_scale", po::value<double>(&rate_scale)->default_value(0.1), "scale of gamma prior for mutation rate")

            ("dirichlet_alpha", po::value<string>(&dirichlet_alpha)->default_value(""), "Dirichlet prior for branch length proportion")
            ("dirichlet_param", po::value<double>(&dirichlet_param)->default_value(0), "Symmetric Dirichlet prior for branch length proportion")
            ("tlen_shape", po::value<double>(&tlen_shape)->default_value(10), "shape of gamma prior for total branch length")
            ("tlen_scale", po::value<double>(&tlen_scale)->default_value(1), "scale of gamma prior for total branch length")

            // ("output_dir", po::value<string>(&output_dir)->default_value("./"), "output directory")
            // ("is_stochastic", po::value<int>(&is_stochastic)->default_value(1), "type of move in each sampling step (0: sequential, 1: stochastic)")
            ("verbose", po::value<int>(&debug)->default_value(0), "verbose level (0: default, 1: debug)")
            ("seed", po::value<unsigned>(&seed)->default_value(0), "seed used for generating random numbers")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(required).add(optional);
    po::options_description config_file_options;
    config_file_options.add(generic).add(required).add(optional);
    po::variables_map vm;

    try {
            po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
            if(vm.count("help")) {
                    cout << cmdline_options << endl;
                    return 1;
            }

            if(vm.count("version")) {
                    cout << "cnetmcmc [version 0.1], a program to build a phylogenetic tree  from total or haplotype-specific copy number profiles of multiple (spatio-temporal) samples for a single patient with MCMC approach" << endl;
                    cout << "This program can also be used to estimate parameters given a tree of fixed topology." << endl;
                    return 1;
            }
            po::notify(vm);
            cout << "configuration file: " << config_file << endl;
            if(config_file!=""){
                ifstream ifs(config_file.c_str());
                po::store(po::parse_config_file(ifs, config_file_options), vm);
                po::notify(vm);
            }
    } catch (const exception& e) {
            cerr << e.what() << endl;
            return 1;
    }

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    setup_rng(r, seed);

    fp_myrng = &myrng;

    ITREE_PARAM itree_param{epop, beta, gtime};

    if(model == MK){
        cout << "Assuming Mk model " << endl;
    }
    if(model == BOUNDT) cout << "Assuming One-step Bounded model " << endl;
    if(model == BOUNDA) cout << "Assuming One-step haplotype-specific model " << endl;

    if(!cons){
        cout << "Assuming the tree is unconstrained " << endl;
    }else{
        cout << "Assuming the tree is constrained by age at sampling time" << endl;
    }

    if(!maxj){
        cout << "Assuming mutation rate is fixed " << endl;
    }else{
        cout << "Estimating mutation rate" << endl;
    }

    if(!correct_bias){
        cout << "Not correcting acquisition bias in likelihood computation " << endl;
    }else{
        cout << "Correcting acquisition bias in likelihood computation " << endl;
    }

    if(maxj){
        if(!only_seg){
            cout << "Estimating mutation rates for site duplication/deletion, chromosome gain/loss, and whole genome doubling " << endl;
        }else{
            cout << "Estimating mutation rates for site duplication/deletion " << endl;
        }
    }

    int num_total_bins = 0;
    Nchar = 0;
    num_invar_bins = 0;
    map<int, vector<vector<int>>> data;
    cout << "\nReading input copy numbers" << endl;
    if(is_bin){
        cout << "   Merging consecutive bins in the input" << endl;
    }else{
        if(incl_all){
            cout << "   Using all input segments " << endl;
            correct_bias = 0;
        }else{
            cout << "   Using variable input segments " << endl;
        }
    }

    INPUT_PROPERTY input_prop{Ns, cn_max, model, is_total, 0, is_bin, incl_all};
    INPUT_DATA input_data{num_invar_bins, num_total_bins, Nchar, obs_num_wgd, obs_change_chr, sample_max_cn};
    data = read_data_var_regions_by_chr(datafile, input_prop, input_data, seg_file, debug);

    // assign variables back for those changed during input parsing
    num_invar_bins = input_data.num_invar_bins;
    num_total_bins = input_data.num_total_bins;
    Nchar = input_data.seg_size;
    obs_num_wgd = input_data.obs_num_wgd;
    obs_change_chr = input_data.obs_change_chr;
    sample_max_cn = input_data.sample_max_cn;

    vobs = get_obs_vector_by_chr(data, Ns);

    int nleaf = Ns + 1;
    for(int k= (nleaf + 1); k < (2 * nleaf - 1); ++k) knodes.push_back(k);
    knodes.push_back(nleaf);

    // tobs already defined globally
    tobs = read_time_info(timefile, Ns, age);
    if(cons){
        cout << "The age of patient at the first sampling time: " << age << endl;
    }
    else{
        age = MAX_AGE;
        tobs = vector<double>(Ns, 0);
    }


    // parameters for proposal distribution
    vector<double> proposal_parameters({lambda, lambda_all, sigma_blen});
    if(model == MK){
        proposal_parameters.push_back(sigma_mut);
    }
    else{
        proposal_parameters.push_back(sigma_dup);
        proposal_parameters.push_back(sigma_del);
        proposal_parameters.push_back(sigma_gain);
        proposal_parameters.push_back(sigma_loss);
        proposal_parameters.push_back(sigma_wgd);
    }
    if(cons){
        proposal_parameters.push_back(sigma_height);
    }

    vector<double> prior_parameters_blen;
    prior_parameters_blen.push_back(tlen_shape);
    prior_parameters_blen.push_back(tlen_scale);
    if(dirichlet_param > 0){
        prior_parameters_blen.push_back(dirichlet_param);
    }

    vector<double> alphas = get_vertex_indices(dirichlet_alpha);
    if(debug){
        cout << "Alphas:";
        for(int i = 0; i < alphas.size(); i++){
            cout << "\t" << alphas[i];
        }
        cout << endl;
    }
    vector<double> prior_parameters_rate({rate_shape, rate_scale});
    vector<double> prior_parameters_mut;
    if(model == MK){ // Mk model
        prior_parameters_mut.push_back(log10(mu));
        prior_parameters_mut.push_back(sigma_lmut);
    }
    else{ // bounded model
        prior_parameters_mut.push_back(log10(dup_rate));
        prior_parameters_mut.push_back(sigma_ldup);
        prior_parameters_mut.push_back(log10(del_rate));
        prior_parameters_mut.push_back(sigma_ldel);
        // if(!only_seg){
            prior_parameters_mut.push_back(log10(chr_gain_rate));
            prior_parameters_mut.push_back(sigma_lgain);
            prior_parameters_mut.push_back(log10(chr_loss_rate));
            prior_parameters_mut.push_back(sigma_lloss);
            prior_parameters_mut.push_back(log10(wgd_rate));
            prior_parameters_mut.push_back(sigma_lwgd);
        // }
    }

    double min_height = *max_element(tobs.begin(), tobs.end());
    double max_height = age + min_height;
    vector<double> prior_parameters_height;
    prior_parameters_height.push_back(min_height);
    prior_parameters_height.push_back(max_height);

    vector<double> rates;
    if(model == MK){
        rates.push_back(mu);
    }else{  //bouned model
        rates.push_back(dup_rate);
        rates.push_back(del_rate);
        // if(!only_seg){
            rates.push_back(chr_gain_rate);
            rates.push_back(chr_loss_rate);
            rates.push_back(wgd_rate);
        // }
    }

    // Build the table after reading input file
    if(model == DECOMP){
        // adjust_m_max();
        cout << "maximum number of WGD events is " << max_wgd << endl;
        cout << "maximum number of chromosome gain/loss events on one chromosome is " << max_chr_change << endl;
        cout << "maximum number of site duplication/deletion events is " << max_site_change << endl;
        build_decomp_table(decomp_table, comps, cn_max, m_max, max_wgd, max_chr_change, max_site_change, is_total);
        // build_decomp_table_withm(decomp_table, comps, cn_max, m_max, max_wgd, max_chr_change, max_site_change, is_total);
        cout << "\tNumber of states is " << comps.size() << endl;
        print_decomp_table(decomp_table);
    }

    max_tobs = *max_element(tobs.begin(), tobs.end());
    lnl_type = {model, cn_max, is_total, cons, max_tobs, age, use_repeat, correct_bias, num_invar_bins, only_seg, infer_wgd, infer_chr, knodes};

    obs_decomp = {m_max, max_wgd, max_chr_change, max_site_change, obs_num_wgd, obs_change_chr};

    vector<double> ref_rates;
    if(model == MK){
        ref_rates.push_back(rmu);
    }else{
        ref_rates.push_back(rdup_rate);
        ref_rates.push_back(rdel_rate);
        ref_rates.push_back(rgain_rate);
        ref_rates.push_back(rloss_rate);
        ref_rates.push_back(rwgd_rate);
    }

    if(rtreefile != "") {
        run_with_reference_tree(rtreefile, Ns, Nchar, num_invar_bins, fix_topology, model, cons, maxj, cn_max, only_seg, correct_bias, is_total, ref_rates, tobs, rates, n_draws,  n_burnin,  n_gap, proposal_parameters, prior_parameters_blen, prior_parameters_height, alphas, prior_parameters_mut, lambda_topl, trace_param_file, trace_tree_file, sample_prior, itree_param);
    }
    else{
        cout << "\nGenerate the start tree" << endl;
        evo_tree rtree;
        if(init_tree == 0){
          cout << "\tUsing random coalescence tree " << endl;
          // start with a random coalescence tree
          rtree = generate_coal_tree(Ns, r, fp_myrng, itree_param);
        }else if(init_tree == 1){
          cout << "\tUsing provided tree " << endl;
          rtree = read_tree_info(file_itree, Ns);
        }else{
          cout << "\tUsing random tree with the same topology as the real tree " << endl;
          evo_tree test_tree = read_reference_tree(rtreefile, Ns, ref_rates, tobs);
          int nedge = test_tree.edges.size();
          gsl_vector* rblens = gsl_vector_alloc(nedge);
          for(int i = 0; i < nedge; ++i){
            gsl_vector_set(rblens, i, runiform(r, 1, age));
          }
          rtree = create_new_tree(rblens, test_tree, max_tobs, 0);
        }

        cout << "\tAssigning mutation rates to the initial tree" << endl;
        revise_init_tree(rtree, rates, tobs, cons);

        // double Ls = get_likelihood_revised(rtree, vobs, lnl_type);
        cout << "\tGetting start tree likelihood" << endl;
        double Ls = 0;
        if(model == DECOMP){
            Ls = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
        }else{
            Ls = get_likelihood_revised(rtree, vobs, lnl_type);
        }
        cout << "   Random start tree likelihood: " << Ls << endl;

        // Estimate branch length with MCMC
        cout << "\n\n### Running MCMC" << endl;
        run_mcmc(rtree, model, n_draws, n_burnin, n_gap, proposal_parameters, prior_parameters_blen, prior_parameters_height, alphas, prior_parameters_mut, lambda_topl, trace_param_file, trace_tree_file, itree_param, sample_prior, fix_topology, cons, maxj, cn_max, only_seg, correct_bias, is_total);
    }
}
