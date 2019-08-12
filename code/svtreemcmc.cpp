// run Bayesian MCMC inference

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
#include <limits>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "gzstream.h"

#include "evo_tree.hpp"
#include "genome.hpp"
#include "stats.hpp"
#include "utilities.hpp"

using namespace std;


const double BLEN_MIN = 1e-10;
const double BLEN_MAX = 100;
const double RATE_MIN = 1e-10;
const double RATE_MAX = 1;
const double RATE_MIN_LOG = -10;
const double RATE_MAX_LOG = 0;
const double LOG_MINUS_INFINITY = numeric_limits<double>::lowest();

int debug = 0;

// Compute the total lengths of branches to estimate
double get_total_blens(gsl_vector* blens, int num_branch){
    double total_blens = 0;
    for(int i=0; i < num_branch; i++){
        total_blens += gsl_vector_get(blens, i);
    }
    // cout << "Total branch length is: " << total_blens << endl;
    return total_blens;
}

void adjust_blen(double& nx, double a, double b){
    // double a = BLEN_MIN;
    // double b = BLEN_MAX;
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
        if (c != 1.0)
        {
            for(int i=0; i<num_branch; i++)
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
        for(int i=0; i<num_branch; i++)
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

    if (x < 0.0)
        {
        invers = 1;
        x  *= -1.0;
        }
    if (x > limit)
        return (invers ? 0 : 1);
    if (x < t)
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
    if(cons == 0){
      npar = rtree.nedge-1;
      // initialise the best guess branch length and mu if required
      x = gsl_vector_alloc (npar);
      for(int i=0; i<npar; ++i){
          blen = rtree.edges[i].length;
          gsl_vector_set (x, i, blen);
      }
    }
    else{   // only estimate internal edges for all tumor samples
      npar = rtree.nintedge;
      x = gsl_vector_alloc (npar);
      // initialise with internal edges
      if(debug){
          cout << "Initialise with internal edges" << endl;
      }
      for(int i=0; i<rtree.nintedge; ++i){
          blen = rtree.intedges[i]->length;
          gsl_vector_set (x, i, blen);
      }
    }
    // double num_branch = x -> size;
    // cout << "There are " << num_branch <<  " parameters for estimating branch lengths " << endl;
    return x;
}


// Create a new tree with the same topology as input tree but different branch lengths
evo_tree create_new_tree(gsl_vector* blens, evo_tree& rtree, int cons){
    // create a new tree
    vector<edge> enew;
    // cout << "copy original edges" << endl;
    for(int i=0; i<rtree.nedge; ++i){
      enew.push_back( rtree.edges[i] );
    }
    if(cons){
        // cout << "branches constrained" << endl;
        int count = 0;
        for(int i=0; i<rtree.nedge-1; ++i){
            if(enew[i].end > Ns){
                // cout << "count " << count << endl;
                enew[i].length = gsl_vector_get(blens, count);
                count++;
            }else{
                enew[i].length = 0;
            }
        }

        // Time to first sample
        // cout << "time to 1st sample: " << total_time << endl;
        evo_tree new_tree(rtree.nleaf, enew, rtree.get_total_time(), rtree.tobs);
        new_tree.mu = rtree.mu;
        new_tree.dup_rate = rtree.dup_rate;
        new_tree.del_rate = rtree.del_rate;

        return new_tree;
    }
    else{
        // cout << "branches unconstrained" << endl;
        for(int i=0; i<rtree.nedge-1; ++i){
            enew[i].length = gsl_vector_get(blens,i);
        }

        evo_tree new_tree(Ns+1, enew);
        new_tree.mu = rtree.mu;
        new_tree.dup_rate = rtree.dup_rate;
        new_tree.del_rate = rtree.del_rate;

        return new_tree;
    }
}

evo_tree create_new_coal_tree(int nsample, int Ne){
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
  for(int i=0; i<nsample; ++i) nodes.push_back(i);

  // create vector of event times
  // total nodes = 2*(nsample+1) -1
  for(int i=0; i<(2*nsample+1); ++i) times.push_back(0.0);

  double t_tot = 0;
  while(nlin > 1){
    // sample a time from Exp( combinations(k,2) )
    double lambda = fact(nlin)/( 2*fact(nlin-2) );
    double t = gsl_ran_exponential(r, 1/lambda) * Ne;
    t_tot += t ;

    // choose two random nodes from available list
    random_shuffle(nodes.begin(), nodes.end(),fp);

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
  double t = gsl_ran_exponential(r, 1/lambda) * Ne;
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
  for(int i=0; i<times.size(); ++i){
    times[i] = t_tot - times[i];
  }

  // cout << "total time of tree: " << t_tot << " : ";
  // for(int i=0; i<epoch_times.size(); ++i) cout << "\t" << epoch_times[i];
  // cout << endl;
  //
  // // invert the times
  // cout << "times of nodes:" << endl;
  // for(int i=0; i<times.size(); ++i){
  //   cout << i+1 << "\t" << times[i] << endl;
  // }

  evo_tree ret(nsample+1, edges, lengths);
  return ret;
}


// Use NNI for clock trees
int move_topology(double& log_hastings_ratio, evo_tree& rtree){
    if(debug){
        printf ("Before:\n");
        rtree.print();
        // Output to a file to facilitate ploting
        ofstream out_tree("./test/sim-data-mcmc-tree-before.txt");
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

    while (tc <= tv && chosen.size() < rtree.nintedge - 1){
        count++;
        // Exclude the edges connecting to unaltered genome
        assert(rtree.intedges[rtree.nintedge-1]->start == Ns+1);
        i = myrng(rtree.nintedge-1);
        chosen.insert(i);
        if(debug){
            cout << "Picking " << i+1 << "th internal edge: " << rtree.intedges[i]->start+1 << ", " << rtree.intedges[i]->end+1 << endl;
        }
        e_uv = rtree.intedges[i];
        // Find the children of the start node
        u = e_uv->start;
        children = rtree.nodes[u].daughters;
        // Select the child with subtree as v
        v = *max_element(children.begin(), children.end());
        c = *min_element(children.begin(), children.end());

        assert(v > Ns+1);

        tv = rtree.node_times[v];
        tc = rtree.node_times[c];

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
        cout << "Cannot update topology!" << endl;
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

    // update branch length to keep node ages unchanged
    e_va->length = blen_uc - e_uv->length;
    e_uc->length = e_uv->length + blen_va;
    adjust_blen(e_va->length, BLEN_MIN, BLEN_MAX);
    adjust_blen(e_uc->length, BLEN_MIN, BLEN_MAX);

    rtree.nodes.clear();
    rtree.generate_nodes();

    rtree.intedges.clear();
    rtree.generate_int_edges();

    /* adjust branch lengths */
    // e_uc->length = rtree.node_times[a] - rtree.node_times[u];
    // e_va->length = rtree.node_times[c] - rtree.node_times[v];

    if(debug){
        printf ("After:\n");
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
gsl_vector* move_blens_multiplier(double& log_hastings_ratio, double a, double b, double lambda, gsl_vector* prev_blens, int num_branch){
    // cout << "Update branch length" << endl;
    gsl_vector *blens = gsl_vector_alloc(num_branch);
    // double a = BLEN_MIN;
    // double b = BLEN_MAX;
    double x, nx, y, ny, c;
    for (int i=0; i < num_branch; i++){
        x = gsl_vector_get(prev_blens, i);
        y = (b - x)/(x - a);
        c = exp(lambda * (runiform(r, 0, 1) - 0.5));
        ny = y * c;
        nx = (b + a * ny)/(ny + 1);
        adjust_blen(nx, a, b);
        gsl_vector_set(blens, i, nx);
        log_hastings_ratio += log(c) + 2 * log(abs((nx-a)/(x-a)));
    }
    return blens;
}


// Propose a new branch length with normal proposal
gsl_vector* move_blens_normal(double& log_hastings_ratio, double sigma, gsl_vector* prev_blens, int num_branch){
    // cout << "Update branch length" << endl;
    gsl_vector *blens = gsl_vector_alloc(num_branch);
    for (int i=0; i < num_branch; i++){
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


void update_topology(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, double lambda_topl, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias){
    double prev_log_prior = 1, prev_log_likelihood, log_prior = 1, log_diff;
    double log_hastings_ratio = 0;
    bool accept = false;
    evo_tree ntree;

    double rn = runiform(r, 0, 1);
    // Change topology with probability lambda_topl
    if (rn > lambda_topl){
        cout << "Not updating tree topology this time" << endl;
        return;
    }

    if(debug){
        cout << "Updating tree topology " << endl;
    }

    prev_log_prior = get_prior_topology(Ns);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }

    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    // rtree will be changed after move
    ntree = evo_tree(rtree);
    int changed = move_topology(log_hastings_ratio, ntree);
    if(changed==0){
        cout << "Canot propose a new tree topology this time" << endl;
        return;
    }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    // log_prior = get_prior_topology(Ns);
    if(debug){
        cout << "   log hastings ratio of topology proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    // cout << "Accept or reject the proposal" << endl;
    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    // cout << "Print the accepted state" << endl;
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}


void update_blen(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_blen, vector<double> alphas, double lambda, double mbactrian, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias) {
    gsl_vector *prev_blens, *blens;
    int num_branch;
    double prev_log_prior, prev_log_likelihood, log_prior, log_diff;
    double log_hastings_ratio = 0;
    bool accept;
    evo_tree ntree;

    if(debug){
        cout << "Updating branch lengths " << endl;
    }

    prev_blens = initialize_branch_length(rtree, cons);
    num_branch = prev_blens->size;
    if(debug){
        cout << "Previous branch length: " << endl;
        for(int k=0; k<num_branch; k++){
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
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double min_blen = BLEN_MIN;
    double max_blen = BLEN_MAX;
    if(cons){
        max_blen = rtree.get_total_time();
    }
    // cout << "Max branch length to propose " << max_blen << endl;
    blens = move_blens_multiplier(log_hastings_ratio, min_blen, max_blen, lambda, prev_blens, num_branch);
    // blens = move_blens_bactrian(log_hastings_ratio, mbactrian, sigma, prev_blens, num_branch);
    // gsl_vector* blens = move_blens_normal(sigma, prev_blens, num_branch, log_hastings_ratio);
    // Create a new tree with the proposed tree length
    if(debug){
        cout << "   Create a new tree with the proposed tree length" << endl;
        // cout << "number of estimated branches " << blens->size << endl;
    }

    ntree = create_new_tree(blens, rtree, cons);
    adjust_tree_blens(ntree);
    adjust_tree_tips(ntree);

    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    log_prior = get_prior_blen(blens, num_branch, prior_parameters_blen, alphas);
    if(debug){
        cout << "   log hastings ratio of branch length proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    // cout << "Accept or reject the proposal" << endl;
    accept = false;
    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        // cout << "log difference " << log_diff << endl;
        if (log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    // cout << "Print the accepted state" << endl;
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
    else{
        gsl_vector_memcpy(blens, prev_blens);
    }
}


void update_tree_height(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters, double sigma, int sample_prior, int cn_max, int only_seg, int correct_bias) {
    // uniform prior
    double prev_log_prior = 1, prev_log_likelihood, log_prior = 1, log_diff;
    double log_hastings_ratio = 0;
    bool accept = false;
    double old_val, new_val;
    evo_tree ntree;

    double min_height = prior_parameters[0];
    double max_height = prior_parameters[1];

    if(debug){
        cout << "Updating tree height " << endl;
    }

    old_val = rtree.get_tree_height();
    // prev_log_prior = get_prior_tree_height(prior_parameters);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, 1, cn_max, only_seg, correct_bias);
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    // Min height must be larger than the largest time of internal nodes
    // for(int i=rtree.nleaf; i< rtree.ntotn; i++){
    //     if(rtree.node_times[i] > min_height){
    //         min_height = rtree.node_times[i];
    //     }
    // }
    new_val = move_normal(log_hastings_ratio, old_val, min_height, max_height, sigma);
    ntree = evo_tree(rtree);
    // Set the new tree height
    // ntree.tree_height = new_val;
    // ntree.total_time = new_val - *max_element(rtree.tobs.begin(), rtree.tobs.end());
    // Rescale the tree
    double ratio =  new_val/old_val;
    ntree.scale_time(ratio);
    // // Update terminal branches
    // for(int i=0; i<ntree.nleaf-1; ++i){
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
    // ntree.generate_int_edges();
    // ntree.lengths.clear();
    // for(int i=0; i<ntree.nedge; ++i){
    //     ntree.lengths.push_back(ntree.edges[i].length);
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, ntree, model, 1, cn_max, only_seg, correct_bias);
    }
    // log_prior = get_prior_tree_height(prior_parameters);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}


// Update the effective population size for coalescent model
void update_pop_size(evo_tree& rtree, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias) {
    // uniform prior
    double prev_log_prior = 1, prev_log_likelihood, log_prior = 1, log_diff;
    double log_hastings_ratio = 0;
    bool accept = false;
    double old_val, new_val;
    evo_tree ntree;

    if(debug){
        cout << "Updating model parameters " << endl;
    }

    old_val = rtree.get_tree_height();
    // prev_log_prior = get_prior_tree_height(prior_parameters);
    if(sample_prior){
        prev_log_likelihood = 1;
    }
    else{
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
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
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    // log_prior = get_prior_tree_height(prior_parameters);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}

void update_mutation_rates(evo_tree& rtree, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_rate, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias) {
    double prev_log_prior, prev_log_likelihood, log_prior, log_diff;
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
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    new_mu = move_mrates_normal(log_hastings_ratio, rtree.mu, sigma);
    ntree = evo_tree(rtree);
    ntree.mu = new_mu;
    // if(debug){
    //     cout << "Old mu " << rtree.mu << endl;
    //     cout << "New mu " << ntree.mu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    log_prior = get_prior_mutation_gamma(ntree.mu, prior_parameters_rate);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}

//
void update_mutation_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias) {
    double prev_log_prior, prev_log_likelihood, log_prior, log_diff;
    double log_hastings_ratio = 0;
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
    }
    else{
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.mu = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old mu " << rtree.mu << "\t" << lmu << endl;
    //     cout << "   New mu " << ntree.mu << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}



void update_deletion_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias) {
    double prev_log_prior, prev_log_likelihood, log_prior, log_diff;
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
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.del_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old deletion rate " << rtree.del_rate << "\t" << lmu << endl;
    //     cout << "   New deletion rate " << ntree.del_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}


void update_duplication_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias) {
    double prev_log_prior, prev_log_likelihood, log_prior, log_diff;
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
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.dup_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old duplication rate " << rtree.dup_rate << "\t" << lmu << endl;
    //     cout << "   New duplication rate " << ntree.dup_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}

void update_cgain_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias) {
    double prev_log_prior, prev_log_likelihood, log_prior, log_diff;
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
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.chr_gain_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old duplication rate " << rtree.dup_rate << "\t" << lmu << endl;
    //     cout << "   New duplication rate " << ntree.dup_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}


void update_closs_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias) {
    double prev_log_prior, prev_log_likelihood, log_prior, log_diff;
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
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.chr_loss_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old duplication rate " << rtree.dup_rate << "\t" << lmu << endl;
    //     cout << "   New duplication rate " << ntree.dup_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}


void update_wgd_rates_lnormal(evo_tree& rtree, int model, double& log_likelihood, int& naccepts, const int n_draw, const int n_burnin, const int n_gap, vector<double> prior_parameters_mut, double sigma, int sample_prior, int cons, int cn_max, int only_seg, int correct_bias) {
    double prev_log_prior, prev_log_likelihood, log_prior, log_diff;
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
        prev_log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    if(debug){
        cout << "   Previous prior and likelihood " << prev_log_prior << "\t" << prev_log_likelihood << endl;
    }

    double new_lmu = move_mrates_normal(log_hastings_ratio, lmu, sigma);
    ntree = evo_tree(rtree);
    ntree.wgd_rate = pow(10, new_lmu);
    // if(debug){
    //     cout << "   Old duplication rate " << rtree.dup_rate << "\t" << lmu << endl;
    //     cout << "   New duplication rate " << ntree.dup_rate << "\t" << new_lmu << endl;
    // }
    if(sample_prior){
        log_likelihood = 1;
    }
    else{
        log_likelihood = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
    }
    log_prior = get_prior_mutation_lnormal(new_lmu, prior_parameters_mut);
    if(debug){
        cout << "   log hastings ratio of rate proposal " << log_hastings_ratio << endl;
        cout << "   Current prior and likelihood " << log_prior << "\t" << log_likelihood << endl;
    }

    if (log_prior > LOG_MINUS_INFINITY)
    {
        log_diff = log_hastings_ratio;
        log_diff += (log_likelihood + log_prior) - (prev_log_likelihood + prev_log_prior);
        if ( log(runiform(r, 0, 1)) < log_diff){
            accept = true;
        }
    }
    if (accept)
    {
        if(n_draw > n_burnin && n_draw % n_gap == 0){
            naccepts++;
        }
        rtree = evo_tree(ntree);
    }
}



// Given a tree, find the MAP estimation of the branch lengths (and optionally mu) assuming branch lengths are independent or constrained in time
void run_mcmc(evo_tree& rtree, int model, const int n_draws, const int n_burnin, const int n_gap, vector<double> proposal_parameters, vector<double> prior_parameters_blen, vector<double> prior_parameters_height, vector<double> alphas, vector<double> prior_parameters_mut, double lambda_topl, string ofile, string tfile, int sample_prior=0, int fix_topology=0, int cons=0, int maxj=0, int cn_max = 4, int only_seg = 1, int correct_bias=0){
        ofstream fout_trace(ofile);
        ofstream fout_tree(tfile);

        int precision = 6;
        int naccepts_toplogy = 0;
        int naccepts_blen = 0;
        int naccepts_height = 0;
        int naccepts_mrate = 0;
        int naccepts_dup = 0;
        int naccepts_del = 0;
        int naccepts_gain = 0;
        int naccepts_loss = 0;
        int naccepts_wgd = 0;

        double log_likelihood;
        double mu_lmut, sigma_lmut, sigma_mut;
        double mu_ldup, sigma_ldup, sigma_dup;
        double mu_ldel, sigma_ldel, sigma_del;
        double mu_lgain, sigma_lgain, sigma_gain;
        double mu_lloss, sigma_lloss, sigma_loss;
        double mu_lwgd, sigma_lwgd, sigma_wgd;

        double lambda = proposal_parameters[0];     // multiplier proposal
        double mbactrian = proposal_parameters[1];     // Bactrian proposal
        double sigma_blen = proposal_parameters[2];     // normal proposal for branch length
        double sigma_height;

        if(model==0){
            mu_lmut = prior_parameters_mut[0];
            sigma_lmut = prior_parameters_mut[1];

            sigma_mut = proposal_parameters[3];
        }
        if(model==1){
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
        if(model==0){
            header = "state\tlnl\tmu";
        }
        if(model==1){
            header = "state\tlnl\tdup_rate\tdel_rate\tgain_rate\tloss_rate\twgd_rate";
        }
        if(cons){
            header += "\theight";
            for(int j=1; j<=rtree.nintedge; j++){
                header += "\tl" + to_string(j);
            }
        }
        else{
            for(int j=1; j<=rtree.nedge-1; j++){
                header += "\tl" + to_string(j);
            }
        }

        fout_trace << header << endl;

        fout_tree << "#nexus" << endl;
        fout_tree << "begin trees;" << endl;

        for (int i = 1; i <= n_draws; ++i)
        {
            if(fix_topology==0){
                update_topology(rtree, model, log_likelihood, naccepts_toplogy, i, n_burnin, n_gap, lambda_topl, sample_prior, cons, cn_max, only_seg, correct_bias);
            }

            update_blen(rtree, model, log_likelihood, naccepts_blen, i, n_burnin, n_gap, prior_parameters_blen, alphas, lambda, mbactrian, sigma_blen, sample_prior, cons, cn_max, only_seg, correct_bias);

            if(cons){
                update_tree_height(rtree, model, log_likelihood, naccepts_height, i, n_burnin, n_gap, prior_parameters_height, sigma_height, sample_prior, cn_max, only_seg, correct_bias);
            }

            if(maxj){
                // update_mutation_rates(rtree, log_likelihood, naccepts_mrate, i, n_burnin, n_gap, prior_parameters_mut, sigma_mut, cons, cn_max, only_seg, correct_bias);
                if(model==0){
                    vector<double> prior_parameters_mu({mu_lmut, sigma_lmut});
                    update_mutation_rates_lnormal(rtree, model, log_likelihood, naccepts_mrate, i, n_burnin, n_gap, prior_parameters_mu, sigma_mut, sample_prior, cons, cn_max, only_seg, correct_bias);
                }
                if(model==1){
                    vector<double> prior_parameters_dup({mu_ldup, sigma_ldup});
                    update_duplication_rates_lnormal(rtree, model, log_likelihood, naccepts_dup, i, n_burnin, n_gap, prior_parameters_dup, sigma_dup, sample_prior, cons, cn_max, only_seg, correct_bias);

                    vector<double> prior_parameters_del({mu_ldel, sigma_ldel});
                    update_deletion_rates_lnormal(rtree, model, log_likelihood, naccepts_del, i, n_burnin, n_gap, prior_parameters_del, sigma_del, sample_prior, cons, cn_max, only_seg, correct_bias);

                    if(only_seg==0){
                        vector<double> prior_parameters_gain({mu_lgain, sigma_lgain});
                        update_cgain_rates_lnormal(rtree, model, log_likelihood, naccepts_gain, i, n_burnin, n_gap, prior_parameters_gain, sigma_gain, sample_prior, cons, cn_max, only_seg, correct_bias);

                        vector<double> prior_parameters_loss({mu_lloss, sigma_lloss});
                        update_closs_rates_lnormal(rtree, model, log_likelihood, naccepts_loss, i, n_burnin, n_gap, prior_parameters_loss, sigma_loss, sample_prior, cons, cn_max, only_seg, correct_bias);

                        vector<double> prior_parameters_wgd({mu_lwgd, sigma_lwgd});
                        update_wgd_rates_lnormal(rtree, model, log_likelihood, naccepts_wgd, i, n_burnin, n_gap, prior_parameters_wgd, sigma_wgd, sample_prior, cons, cn_max, only_seg, correct_bias);
                    }
                }
            }
            // Discard burn in samples
            if(i <= n_burnin){
                continue;
            }
            if (i % n_gap == 0){
                // print out the accepted proposal
                // string str_tree = order_tree_string(create_tree_string(rtree));
                // fout << i - n_burnin << "\t" << str_tree << "\t"  << log_likelihood << "\t" << rtree.mu ;
                fout_trace << (i - n_burnin)/n_gap << "\t"  << log_likelihood;
                if(model==0){
                    fout_trace << "\t" << rtree.mu ;
                }
                if(model==1){
                    fout_trace << "\t" << rtree.dup_rate << "\t" << rtree.del_rate << "\t" << rtree.chr_gain_rate << "\t" << rtree.chr_loss_rate << "\t" << rtree.wgd_rate;
                }
                if(cons==1){
                    fout_trace << "\t" << rtree.get_tree_height();
                }
                // print out the branch lengths

                if(cons){
                    for(int k=0; k<rtree.nintedge; ++k){
                        fout_trace << "\t" << rtree.intedges[k]->length;
                    }
                }
                else{
                    for(int k=0; k<rtree.nedge-1; k++){
                        fout_trace << "\t" << rtree.edges[k].length;
                    }
                }
                fout_trace << endl;

                string newick = rtree.make_newick(precision);
                fout_tree << "tree " << i << " = " << newick << ";" << endl;
            }
        }

        fout_tree << "end;" << endl;

        fout_trace.close();
        fout_tree.close();

        double n_keep = (n_draws - n_burnin)/n_gap;

        if(fix_topology==0){
            double accept_rate_topology = naccepts_toplogy / n_keep;
            cout << "Accept rate of topology proposal: " << accept_rate_topology << endl;
        }

        double accept_rate_blen = naccepts_blen / n_keep;
        cout << "Accept rate of branch length proposal: " << accept_rate_blen << endl;

        if(cons==1){
            double accept_rate_height = naccepts_height / n_keep;
            cout << "Accept rate of tree height proposal: " << accept_rate_height << endl;
        }

        if(model==0){
            double accept_rate_mrate = naccepts_mrate / n_keep;
            cout << "Accept rate of mutation rate proposal: " << accept_rate_mrate << endl;
        }

        if(model==1){
            double accept_rate_mrate = naccepts_dup / n_keep;
            cout << "Accept rate of duplication rate proposal: " << accept_rate_mrate << endl;
            accept_rate_mrate = naccepts_del / n_keep;
            cout << "Accept rate of deletion rate proposal: " << accept_rate_mrate << endl;

            if(only_seg==0){
                accept_rate_mrate = naccepts_gain / n_keep;
                cout << "Accept rate of chromosme gain rate proposal: " << accept_rate_mrate << endl;

                accept_rate_mrate = naccepts_loss / n_keep;
                cout << "Accept rate of chromosme loss rate proposal: " << accept_rate_mrate << endl;

                accept_rate_mrate = naccepts_wgd / n_keep;
                cout << "Accept rate of whole genome doubling rate proposal: " << accept_rate_mrate << endl;
            }
        }
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


int main (int argc, char ** const argv) {
    int n_draws, n_burnin, n_gap, model, cons, maxj, seed, clock, cn_max, correct_bias;
    double lambda, sigma_blen, mbactrian, lambda_topl, sigma_height;
    double dirichlet_param, tlen_shape, tlen_scale;
    double rate_shape, rate_scale;
    string datafile, timefile, trace_param_file, trace_tree_file, rtreefile, config_file;
    string dirichlet_alpha;
    int sample_prior, fix_topology, only_seg;
    double rmu, rdup_rate, rdel_rate, rgain_rate, rloss_rate, rwgd_rate;
    double mu, dup_rate, del_rate, chr_gain_rate, chr_loss_rate, wgd_rate;
    double sigma_mut, sigma_dup, sigma_del, sigma_gain, sigma_loss, sigma_wgd;
    double sigma_lmut, sigma_ldup, sigma_ldel, sigma_lgain, sigma_lloss, sigma_lwgd;

    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
    ;
    po::options_description required("Required parameters");
    required.add_options()
            ("cfile,c", po::value<string>(&datafile)->required(), "input copy number profile file")
            ("tfile,t", po::value<string>(&timefile)->required(), "input time information file")
    ;
    po::options_description optional("Optional parameters");
    optional.add_options()
            ("config_file", po::value<string>(&config_file)->default_value(""), "configuration file of input parameters")
            ("nsample,s", po::value<int>(&Ns)->default_value(5), "number of samples or regions")
            ("cn_max", po::value<int>(&cn_max)->default_value(4), "maximum copy number of a segment")

            ("trace_param_file", po::value<string>(&trace_param_file)->default_value("trace-mcmc-params.txt"), "output trace file of parameter values")
            ("trace_tree_file", po::value<string>(&trace_tree_file)->default_value("trace-mcmc-trees.txt"), "output trace file of trees")
            ("n_burnin,r", po::value<int>(&n_burnin)->default_value(9000), "number of burnin samples")
            ("n_draws,n", po::value<int>(&n_draws)->default_value(10000), "number of posterior draws to keep")
            ("n_gap", po::value<int>(&n_gap)->default_value(1), "sampling every kth samples ")
            ("sample_prior", po::value<int>(&sample_prior)->default_value(0), "whether or not to sample from the prior only")

            ("model,d", po::value<int>(&model)->default_value(0), "model of evolution (0: JC69, 1: 1-step bounded)")
            ("clock", po::value<int>(&clock)->default_value(0), "model of molecular clock (0: strict global, 1: random local clock)")
            ("fix_topology", po::value<int>(&fix_topology)->default_value(0), "whether or not to fix the topology of the tree")
            ("constrained", po::value<int>(&cons)->default_value(1), "constraints on branch length (0: none, 1: fixed total time)")
            ("estimate_mu", po::value<int>(&maxj)->default_value(1), "estimation of mutation rate (0: mutation rate fixed to be the given value, 1: estimating mutation rate)")
            ("correct_bias", po::value<int>(&correct_bias)->default_value(1), "correct ascertainment bias")

            ("rtreefile", po::value<string>(&rtreefile)->default_value(""), "reference tree file")
            ("rmu", po::value<double>(&rmu)->default_value(0.02), "mutation rate of the reference tree")
            ("rdup_rate", po::value<double>(&rdup_rate)->default_value(0.01), "segment duplication rate of the reference tree")
            ("rdel_rate", po::value<double>(&rdel_rate)->default_value(0.01), "segment deletion rate of the reference tree")
            ("rgain_rate", po::value<double>(&rgain_rate)->default_value(0.01), "chromosme gain rate of the reference tree")
            ("rloss_rate", po::value<double>(&rloss_rate)->default_value(0.01), "chromosme loss rate of the reference tree")
            ("rwgd_rate", po::value<double>(&rwgd_rate)->default_value(0.01), "whole genome doubling rate of the reference tree")

            ("sigma_blen,g", po::value<double>(&sigma_blen)->default_value(2.5), "sigma for proposal of branch length")
            ("sigma_height", po::value<double>(&sigma_height)->default_value(2.5), "sigma for prior of tree height")

            ("only_seg", po::value<int>(&only_seg)->default_value(0), "Whether or not to only consider segment-level mutations (0: include chromosome gain/loss and whole genome doubling, 1: only consider segment-level mutations)")
            ("mu,x", po::value<double>(&mu)->default_value(0.025), "mean of mutation rate")
            ("dup_rate", po::value<double>(&dup_rate)->default_value(0.01), "mean of segment duplication rate")
            ("del_rate", po::value<double>(&del_rate)->default_value(0.01), "mean of segment deletion rate")
            ("chr_gain_rate", po::value<double>(&chr_gain_rate)->default_value(0.001), "mean of chromosome gain rate")
            ("chr_loss_rate", po::value<double>(&chr_loss_rate)->default_value(0.001), "mean of chromosome loss rate")
            ("wgd_rate", po::value<double>(&wgd_rate)->default_value(0.0001), "mean of WGD (whole genome doubling) rate")

            ("sigma_mut", po::value<double>(&sigma_mut)->default_value(0.05), "sigma for proposal of mutation rate")
            ("sigma_dup", po::value<double>(&sigma_dup)->default_value(0.05), "sigma for proposal of segment duplication rate")
            ("sigma_del", po::value<double>(&sigma_del)->default_value(0.05), "sigma for proposal of segment deletion rate")
            ("sigma_gain", po::value<double>(&sigma_gain)->default_value(0.05), "sigma for proposal of chromosme gain rate")
            ("sigma_loss", po::value<double>(&sigma_loss)->default_value(0.05), "sigma for proposal of chromosme loss rate")
            ("sigma_wgd", po::value<double>(&sigma_wgd)->default_value(0.05), "sigma for proposal of whole genome doubling rate")
            ("lambda,l", po::value<double>(&lambda)->default_value(0.5), "lambda for multiplier proposal")
            ("lambda_topl", po::value<double>(&lambda_topl)->default_value(1), "lambda for topology proposal")
            ("mbactrian,m", po::value<double>(&mbactrian)->default_value(0.95), "m for Bactrian proposal")

            ("sigma_lmut", po::value<double>(&sigma_lmut)->default_value(0.05), "sigma for prior of log10 of mutation rate")
            ("sigma_ldup", po::value<double>(&sigma_ldup)->default_value(0.05), "sigma for prior of log10 of segment duplication rate")
            ("sigma_ldel", po::value<double>(&sigma_ldel)->default_value(0.05), "sigma for prior of log10 of segment deletion rate")
            ("sigma_lgain", po::value<double>(&sigma_lgain)->default_value(0.05), "sigma for prior of log10 of chromosme gain rate")
            ("sigma_lloss", po::value<double>(&sigma_lloss)->default_value(0.05), "sigma for prior of log10 of chromosme loss rate")
            ("sigma_lwgd", po::value<double>(&sigma_lwgd)->default_value(0.05), "sigma for prior of log10 of whole genome doubling rate")
            ("tlen_shape", po::value<double>(&tlen_shape)->default_value(10), "shape of gamma prior for total branch length")
            ("dirichlet_alpha", po::value<string>(&dirichlet_alpha)->default_value(""), "Dirichlet prior for branch length proportion")
            ("dirichlet_param", po::value<double>(&dirichlet_param)->default_value(0), "Symmetric Dirichlet prior for branch length proportion")
            ("tlen_scale", po::value<double>(&tlen_scale)->default_value(1), "scale of gamma prior for total branch length")
            ("rate_shape", po::value<double>(&rate_shape)->default_value(1), "shape of gamma prior for mutation rate")
            ("rate_scale", po::value<double>(&rate_scale)->default_value(0.1), "scale of gamma prior for mutation rate")

            // ("output_dir", po::value<string>(&output_dir)->default_value("./"), "output directory")
            ("verbose", po::value<int>(&debug)->default_value(0), "verbose level (0: default, 1: debug)")
            ("seed", po::value<int>(&seed)->default_value(0), "seed used for generating random numbers")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(required).add(optional);
    po::options_description config_file_options;
    config_file_options.add(generic).add(required).add(optional);
    po::variables_map vm;

    try {
            po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
            po::notify(vm);
            cout << "configuration file: " << config_file << endl;

            if(config_file!=""){
                ifstream ifs(config_file.c_str());
                po::store(po::parse_config_file(ifs, config_file_options), vm);
            }

            if(vm.count("help")) {
                    cout << cmdline_options << endl;
                    return 1;
            }

            if(vm.count("version")) {
                    cout << "svtreeml [version 0.1], a program to build a phylogenetic tree from copy number profile" << endl;
                    return 1;
            }
            po::notify(vm);
    } catch (const exception& e) {
            cerr << e.what() << endl;
            return 1;
    }

    setup_rng(seed);

    if(model==0){
        cout << "Assuming JC69 model " << endl;
    }
    if(model==1){
        cout << "Assuming Bounded model " << endl;
    }
    if (cons==0){
        cout << "Assuming the tree is unconstrained " << endl;
    }
    else{
        cout << "Assuming the tree is constrained by age " << endl;
    }

    if (maxj==0){
        cout << "Assuming mutation rate is fixed " << endl;
    }
    else{
        cout << "Estimating mutation rate" << endl;
    }

    if (correct_bias==0){
        cout << "Not correcting acquisition bias in likelihood computation " << endl;
    }
    else{
        cout << "Correcting acquisition bias in likelihood computation " << endl;
    }

    if (only_seg==0){
        cout << "Estimating mutation rates for segment duplication/deletion, chromosme gain/loss, and whole genome doubling " << endl;
    }
    else{
        cout << "Estimating mutation rates for segment duplication/deletion " << endl;
    }

    // vector<vector<int> > data = read_data_var_regions(datafile, Ns, cn_max);
    int num_invar_bins = 0;
    Nchar = 0;
    map<int, vector<vector<int>>> data = read_data_var_regions_by_chr(datafile, Ns, cn_max, num_invar_bins);
    // Nchar = data.size();

    // tobs already defined globally
    tobs = read_time_info(timefile, Ns, age);
    if(cons){
        cout << "The age of patient at the first sampling time: " << age << endl;
    }

    //vector<vector<int> > vobs; // already defined globally
    // for(int nc=0; nc<Nchar; ++nc) {
    //     vector<int> obs;
    //     for(int i=0; i<Ns; ++i) {
    //             obs.push_back(data[nc][i+3]);
    //     }
    //     vobs.push_back(obs);
    // }

    // Construct the CN matrix by chromosme
    int total_chr = data.rbegin()->first;
    for(int nchr=1; nchr<=total_chr; nchr++){
      vector<vector<int>> obs_chr;
      Nchar += data[nchr].size();
      for(int nc=0; nc<data[nchr].size(); ++nc){
          vector<int> obs;
          for(int i=0; i<Ns; ++i){
            obs.push_back( data[nchr][nc][i+3] );
          }
          obs_chr.push_back( obs );
      }
      vobs[nchr] = obs_chr;
    }

    // parameters for proposal distribution
    vector<double> proposal_parameters({lambda, mbactrian, sigma_blen});
    if(model==0){
        proposal_parameters.push_back(sigma_mut);
    }
    if(model==1){
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
    if(model==0){
        prior_parameters_mut.push_back(log10(mu));
        prior_parameters_mut.push_back(sigma_lmut);
    }
    if(model==1){
        prior_parameters_mut.push_back(log10(dup_rate));
        prior_parameters_mut.push_back(sigma_ldup);
        prior_parameters_mut.push_back(log10(del_rate));
        prior_parameters_mut.push_back(sigma_ldel);
        prior_parameters_mut.push_back(log10(chr_gain_rate));
        prior_parameters_mut.push_back(sigma_lgain);
        prior_parameters_mut.push_back(log10(chr_loss_rate));
        prior_parameters_mut.push_back(sigma_lloss);
        prior_parameters_mut.push_back(log10(wgd_rate));
        prior_parameters_mut.push_back(sigma_lwgd);
    }

    double min_height = *max_element(tobs.begin(), tobs.end());
    double max_height = age + min_height;
    vector<double> prior_parameters_height;
    prior_parameters_height.push_back(min_height);
    prior_parameters_height.push_back(max_height);

    if(rtreefile!="") {
            // MLE testing
            // read in true tree
            evo_tree test_tree = read_tree_info(rtreefile, Ns);
            test_tree.print();
            if(model==0){
                test_tree.mu = rmu;
            }
            if(model==1){
                test_tree.dup_rate = rdup_rate;
                test_tree.del_rate = rdel_rate;
                test_tree.chr_gain_rate = chr_gain_rate;
                test_tree.chr_loss_rate = chr_loss_rate;
                test_tree.wgd_rate = wgd_rate;
                test_tree.mu =dup_rate + del_rate + chr_gain_rate + chr_loss_rate + wgd_rate;
            }
            test_tree.tobs = tobs;
            test_tree.tree_height = test_tree.get_tree_height();
            test_tree.total_time = test_tree.get_total_time();
            cout << "Tree height " << test_tree.tree_height << endl;
            cout << "Total time " << test_tree.total_time << endl;
            double Ls = 0.0;
            Ls = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, test_tree, model, cons, cn_max, only_seg, correct_bias);
            // cout << "\nOriginal tree likelihood: " << Ls << endl;
            cout << 0 << "\t" << Ls  << "\t" << test_tree.mu << "\t" << test_tree.dup_rate << "\t" << test_tree.del_rate;
            // Output the original tree length
            double avg_blen = 0;
            for(int i=0; i<test_tree.nedge; ++i){
              avg_blen += test_tree.edges[i].length;
              cout << "\t" << test_tree.edges[i].length;
            }
            avg_blen = avg_blen/test_tree.nedge;
            cout << "\t" << avg_blen;
            cout << endl;
            cout << "Internal EDGES:" << endl;
            cout << "\tid\tstart\tend\tlength" << endl;
            for(int i=0; i<test_tree.nintedge; ++i){
              cout << "\t" << test_tree.intedges[i]->id+1 << "\t" << test_tree.intedges[i]->start+1 << "\t" << test_tree.intedges[i]->end+1 << "\t" << test_tree.intedges[i]->length << endl;
            }

            double Lf = 0;
            stringstream sstm;
            ofstream out_tree;

            cout << "Generate the start tree" << endl;
            // start with a random coalescence tree
            evo_tree rtree;
            // start with a random tree with the same topology
            if(fix_topology){
                gsl_vector* rblens = gsl_vector_alloc(test_tree.nedge);
                for(int i=0; i<test_tree.nedge; ++i){
                  gsl_vector_set(rblens, i, runiform(r, 1, age));
                }
                rtree = create_new_tree(rblens, test_tree, 0);
            }
            else{
                rtree = generate_coal_tree(Ns);
            }
            if(model==0){
                rtree.mu = mu;
            }
            if(model==1){
                rtree.dup_rate = dup_rate;
                rtree.del_rate = del_rate;
                rtree.mu = rtree.dup_rate + rtree.del_rate;
            }
            rtree.tobs = tobs;
            if(cons){
                adjust_tree_height(rtree);
                adjust_tree_tips(rtree);
                adjust_tree_blens(rtree);
            }
            rtree.tree_height = rtree.get_tree_height();
            rtree.total_time = rtree.get_total_time();
            cout << "Total time of random tree " << rtree.total_time << endl;


            sstm << "./test1/sim-data-" << cons << maxj << "-mcmc-tree-start.txt";
            out_tree.open(sstm.str());
            rtree.write(out_tree);
            out_tree.close();
            sstm.str("");

            // Ls = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
            Ls = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, 0, cn_max, only_seg, correct_bias);
            cout << "\nRandom tree likelihood: " << Ls << endl;

            // Estimate branch length with MCMC
            cout << "\n\n### Running MCMC" << endl;
            run_mcmc(rtree, model, n_draws, n_burnin, n_gap, proposal_parameters, prior_parameters_blen, prior_parameters_height, alphas, prior_parameters_mut, lambda_topl, trace_param_file, trace_tree_file, sample_prior, fix_topology, cons, maxj, cn_max, only_seg, correct_bias);
            // cout << "\nMinimised tree likelihood by MCMC / mu : " << Lf << "\t" << min_tree.mu*Nchar <<  endl;

            // cout << "\n\n### Running ML" << endl;
            // // evo_tree min_tree_ml = max_likelihood(rtree, model, Lf, 0.01, 0.01, 2000, cons, maxj);
            // evo_tree min_tree_ml = max_likelihood_BFGS(rtree, model, Lf, 0.01, 0.01, cons, maxj);
            // if(model==0){
            //     cout << "\nMinimised tree likelihood by ML / mu : " << Lf << "\t" << min_tree_ml.mu <<  endl;
            // }
            // if(model==1){
            //     cout << "\nMinimised tree likelihood by ML / dup_rate / del_rate : " << Lf << "\t" << min_tree_ml.dup_rate << "\t" << min_tree_ml.del_rate <<  endl;
            // }
            //
            // //min_tree.print();
            // sstm << "./test1/sim-data-" << cons << maxj << "-ML-tree.txt";
            // out_tree.open(sstm.str());
            // min_tree_ml.write(out_tree);
            // out_tree.close();
            // sstm.str("");
    }
    else{
        cout << "Generate the start tree" << endl;
        // start with a random coalescence tree
        evo_tree rtree = generate_coal_tree(Ns);

        if(model==0){
            rtree.mu = mu;
        }
        if(model==1){
            rtree.dup_rate = dup_rate;
            rtree.del_rate = del_rate;
            rtree.chr_gain_rate = chr_gain_rate;
            rtree.chr_loss_rate = chr_loss_rate;
            rtree.wgd_rate = wgd_rate;
            rtree.mu = rtree.dup_rate + rtree.del_rate + rtree.chr_gain_rate + rtree.chr_loss_rate + rtree.wgd_rate;
        }
        rtree.tobs = tobs;

        adjust_tree_height(rtree);
        adjust_tree_tips(rtree);
        adjust_tree_blens(rtree);

        rtree.tree_height = rtree.get_tree_height();
        rtree.total_time = rtree.get_total_time();
        cout << "Total hight of random tree " << rtree.tree_height << endl;

        // double Ls = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
        double Ls = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, only_seg, correct_bias);
        cout << "\nRandom tree likelihood: " << Ls << endl;

        // Estimate branch length with MCMC
        cout << "\n\n### Running MCMC" << endl;
        run_mcmc(rtree, model, n_draws, n_burnin, n_gap, proposal_parameters, prior_parameters_blen, prior_parameters_height, alphas, prior_parameters_mut, lambda_topl, trace_param_file, trace_tree_file, sample_prior, fix_topology, cons, maxj, cn_max, only_seg, correct_bias);
    }
}
