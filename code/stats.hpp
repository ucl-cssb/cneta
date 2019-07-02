#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>
#include <cstdio>

#include <algorithm>
#include <vector>
#include <map>

#include "lbfgsb_new.h"
// #include "expm.hpp"
#include "matexp/matrix_exponential.hpp"
#include "matexp/r8lib.hpp"
#include <boost/numeric/ublas/io.hpp>

extern "C" {
  gsl_rng * r;
  // vector<vector<int>> vobs;
  map<int, vector<vector<int>>> vobs;
  vector<double> tobs;
  // boost::numeric::ublas::matrix<double> qmat;
  int cn_max;
  // double qmat;
  int Ns;
  int Nchar;
  int num_invar_bins;
  int model;
  int correct_bias;
  int age;
  int has_wgd;
  int has_chr_gain;
  int has_chr_loss;
}

const double LARGE_LNL = -1e9;
const double SMALL_LNL = -1e20;
const double ERROR_X = 1.0e-4;
const double SMALL_PROB = 1.0e-20;

void setup_rng(int set_seed){
  gsl_rng_env_setup();

  const gsl_rng_type* T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  if( set_seed != 0 ){
    gsl_rng_set(r, set_seed);
  }else{
    int t = time(NULL);
    int pid = getpid();
    long s = t*pid;
    //cout << "pid:" << "\t" << getpid() << endl;
    cout << "seed:" << "\t" << t << "\t" << pid << "\t" << abs(s) << endl;
    gsl_rng_set (r, abs(s) );
  }
}

// unary function and pointer to unary function
// allows use of gsl rng for standard template algorithms
long unsigned myrng(long unsigned n)
{
  return gsl_rng_uniform_int(r, n);
}

long unsigned (*fp)(long unsigned) = myrng;

// factorial for choose(n,k)
int fact(int n){
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}

// wrapper for uniform
double runiform(gsl_rng* r, double a, double b){
  double myrandom = a + (b-a)*gsl_rng_uniform (r);

  while(myrandom==0){
    myrandom = a + (b-a)*gsl_rng_uniform (r);
  }
  return myrandom;
}


// sample an element according to a probability vector
// gsl_ran_multinomial?
int rchoose(gsl_rng* r, const vector<double>& rates){
  //cout << "rchoose, rates:";
  //for(int i=0; i<rates.size(); ++i) cout << "\t" << rates[i];
  //cout << endl;

  vector<double> p;
  double s = accumulate(rates.begin(), rates.end(), 0.0);

  for(int i=0; i<rates.size(); ++i){
    p.push_back(rates[i]/s);
  }

  vector<double> psum(p.size(),0.0);
  partial_sum(p.begin(),p.end(),psum.begin());

  double u = gsl_rng_uniform(r);
  int ret = -1;
  for(int i=0; i<rates.size(); ++i){
    if( u < psum[i] ){
      ret = i;
      break;
    }
  }

  //cout << "u=\t" << u << "\t" << ret << endl;
  return ret;
}


int get_num_params(){
    int num = 2;
    if (has_wgd == 0){
        num +=1;
    }
    if (has_chr_gain == 0){
        num +=1;
    }
    if (has_chr_loss == 0){
        num +=1;
    }
    return num;
}

// generate neutral coalescent trees
// here nsample is the number of cancer samples (not including germline node)
// here we directly calculate the edges in the tree
void generate_coal_tree(const int& nsample, vector<int>& edges, vector<double>& lengths, vector<double>& epoch_times, vector<double>& times){
  //cout << "GENERATING COAL TREE" << endl;
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
    double t = gsl_ran_exponential(r, 1/lambda);
    t_tot += t;

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
  double t = gsl_ran_exponential(r, 1/lambda);
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

  //cout << "total time of tree: " << t_tot << " : ";
  //for(int i=0; i<epoch_times.size(); ++i) cout << "\t" << epoch_times[i];
  //cout << endl;

  // invert the times
  //cout << "times of nodes:" << endl;
  //for(int i=0; i<times.size(); ++i){
    //cout << i+1 << "\t" << times[i] << endl;
  //}
}


// Scale the total time by given time
evo_tree generate_coal_tree(const int& nsample){
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
    double t = gsl_ran_exponential(r, 1/lambda);
    t_tot += t;

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
  double t = gsl_ran_exponential(r, 1/lambda);
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

//
// Likelihood and transition probabilities
//
double get_transition_prob(const double& mu, double blength, const int& sk, const int& sj){
  //if(debug) cout << "\tget_transition_prob" << endl;

  // assume a basic Markov model here cf JC
  // Qij = u/5 when i!=j
  // we have five states 0, 1, 2, 3, 4
  // Pij(t) = exp(-ut) k_ij + (1 - exp(-ut))*pi_j
  // pi_0 = pi_1 = pi_2 = pi_3 = pi_4 = 1/5

  double prob = 0;

  if( sk == sj ){
    prob = exp(-mu*blength) + (1 - exp(-mu*blength))*0.2;
  }else{
    prob = (1 - exp(-mu*blength))*0.2;
  }

  return prob;
}


// boost::numeric::ublas::matrix<double> get_rate_matrix_bounded(double dup_rate=0.01, double del_rate=0.01) {
//     boost::numeric::ublas::matrix<double> m(cn_max+1, cn_max+1);
//     for (unsigned i = 0; i < m.size1 (); ++ i){
//         for (unsigned j = 0; j < m.size2 (); ++ j){
//             m (i, j) = 0;
//         }
//     }
//     for( unsigned i = 1; i < cn_max; i++){
//         m (i,i-1) = 2 * i * del_rate;
//         m (i,i+1) = 2 * i * dup_rate;
//         m(i, i) = 0 - m (i,i-1) - m (i,i+1);
//     }
//     m (cn_max, cn_max - 1) = 2 * cn_max * del_rate;
//     m (cn_max, cn_max) = 0 - m (cn_max, cn_max - 1);
//
//     // std::cout << m << std::endl;
//     return m;
// }

void get_rate_matrix_bounded(double* m, double dup_rate, double del_rate, const int cn_max) {
    int debug = 0;
    int ncol = cn_max + 1;

    for (unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }
    for( unsigned i = 1; i < cn_max; i++){
        m[i+(i-1)*ncol] = 2 * i * del_rate;
        m[i+(i+1)*ncol] = 2 * i * dup_rate;
        m[i+i*ncol] = 0 - m[i+(i-1)*ncol] - m[i+(i+1)*ncol];
    }
    m[cn_max + (cn_max - 1)*ncol] = 2 * cn_max * del_rate;
    m[cn_max + cn_max*ncol] = 0 - m[cn_max + (cn_max - 1)*ncol];

    if(debug){
        std::cout << m << std::endl;
        r8mat_print ( ncol, ncol, m, "  A:" );
    }
}

// // http://www.guwi17.de/ublas/examples/
// // Use eigendecomposition to compute matrix exponential
// double get_transition_prob_bounded(const boost::numeric::ublas::matrix<double>& q, double t, const int& sk, const int& sj ){
//     //
//     boost::numeric::ublas::matrix<double> tmp = q * t;
//     boost::numeric::ublas::matrix<double> p = expm_pad(tmp);
//     cout << "t: " << t << endl;
//     cout << "Q matrix: " << q << endl;
//     cout << "tmp matrix: " << tmp << endl;
//     cout << "P matrix: " << p << endl;
//     return p(sk, sj);
// }
void get_transition_matrix_bounded(double* q, double* p, double t, int cn_max){
    int debug = 0;
    int n = cn_max + 1;
    double tmp[n*n]={0};

    for(int i = 0; i < n*n; i++){
        tmp[i] = q[i] * t;
    }
    double* res = r8mat_expm1(n, tmp);
    for(int i=0; i<n*n; i++){
        // if(res[i]<SMALL_PROB){
        //     p[i] = 0;
        // }
        // else{
        p[i] = res[i];
        // }
    }

    if(debug){
        cout << "t: " << t << endl;
        r8mat_print ( n, n, q, "  Q matrix:" );
        r8mat_print ( n, n, tmp, "  TMP matrix:" );
        r8mat_print ( n, n, p, "  P matrix:" );
    }
}

double get_transition_prob_bounded(double* p, const int& sk, const int& sj, int cn_max){
    int debug = 0;

    int n = cn_max + 1;
    int i = sk  + sj * n;
    double v = p[i];

    if(debug){
        r8mat_print ( n, n, p, "  P matrix before access:" );
        cout << "Prob at " << i << "("<< sk << ", " << sj << ") is: " << v << endl;
    }
    // if(v<SMALL_PROB){
    //     v = 0;
    // }
    return v;
}


// Check whether the set of branch lengths is valid under the time constraints
int is_tree_valid(evo_tree& rtree, int cons){
    for(int i=0; i<rtree.edges.size(); ++i){
        if (rtree.edges[i].length < 0 ){
            cout << "Negative branch length!" << endl;
            return 0;
        }
    }
    if(cons){
        double tot = rtree.get_total_time();
        if(tot > age + 1){
            cout << "Wrong total time (larger than age)!" << "\t" << tot << "\t" << age << endl;
            return 0;
        }
    }
    return 1;
}


// Compute the likelihood without grouping sites by chromosme
double get_likelihood(const int& Ns, const int& Nchar, const vector<vector<int>>& vobs, evo_tree& rtree, int model, int cons, const int cn_max){
  int debug = 0;
  if(debug) cout << "\tget_likelihood" << endl;
  int nstate = cn_max + 1;

  // return 0 if the tree is not valid
  if(!is_tree_valid(rtree, cons)){
      return SMALL_LNL;
  }

  double qmat[(nstate)*(nstate)]={0};
  // double pmat[(cn_max+1)*(cn_max+1)]={0};
  double pmati[(nstate)*(nstate)]={0};
  double pmatj[(nstate)*(nstate)]={0};

  if(model == 1){
      get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
  }

  double logL = 0;
  for(int nc=0; nc<Nchar; ++nc){
    vector<int> obs = vobs[nc];

    if(debug){
      cout << "char: " << nc << endl;
      for(int i=0; i<Ns; ++i) cout << "\t" << obs[i];
      cout << endl;
    }

    vector<vector<double> > L_sk_k( rtree.ntotn, vector<double>(nstate,0) );

    for(int i=0; i<Ns; ++i){
      for(int j=0; j<nstate; ++j){
	         if( j == obs[i] ) L_sk_k[i][j] = 1.0;
      }
    }
    // set unaltered
    L_sk_k[Ns][2] = 1.0;

    if(debug){
      cout << "\nLikelihood so far:\n";
      for(int i=0; i<rtree.ntotn; ++i){
    	for(int j=0; j<nstate; ++j){
    	  cout << "\t" << L_sk_k[i][j];
    	}
	    cout << endl;
      }
    }

    //create a list of nodes to loop over, making sure the root is last
    vector<int> knodes;
    for(int k=Ns+2; k<rtree.ntotn; ++k) knodes.push_back( k );
    knodes.push_back(Ns+1);

    for(int kn=0; kn<knodes.size(); ++kn){
      int k = knodes[kn];
      int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
      double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
      int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
      double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
      if(model==1){
           get_transition_matrix_bounded(qmat, pmati, bli, cn_max);
           get_transition_matrix_bounded(qmat, pmatj, blj, cn_max);
           if(debug){
               cout << "Get Pmatrix for branch length " << bli << endl;
               r8mat_print( nstate, nstate, pmati, "  P matrix after change:" );
               cout << "Get Pmatrix for branch length " << blj << endl;
               r8mat_print( nstate, nstate, pmatj, "  P matrix after change:" );
           }
      }

      if(debug) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni+1 << " , " << bli << "\t" <<  nj+1 << " , " << blj << endl;

      //loop over possible values of sk
      for(int sk=0; sk<nstate; ++sk){
    	double Li = 0;
    	// loop over possible si
    	for(int si=0; si<nstate; ++si){
            if (model == 0){
                Li += get_transition_prob(rtree.mu, bli, sk, si) * L_sk_k[ni][si];
            }
            if (model == 1){
                Li += get_transition_prob_bounded(pmati, sk, si, cn_max) * L_sk_k[ni][si];
            }
    	  //cout << "\tscoring: Li\t" << li << "\t" << get_transition_prob(mu, bli, sk, si ) << "\t" << L_sk_k[ni][si] << endl;
        }
    	double Lj = 0;
    	// loop over possible sj
    	for(int sj=0; sj<nstate; ++sj){
            if (model == 0){
    	         Lj += get_transition_prob(rtree.mu, blj, sk, sj) * L_sk_k[nj][sj];
            }
            if (model == 1){
                 Lj += get_transition_prob_bounded(pmatj, sk, sj, cn_max) * L_sk_k[nj][sj];
            }
    	}
	    //cout << "scoring: sk" << sk << "\t" <<  Li << "\t" << Lj << endl;
	    L_sk_k[k][sk] = Li*Lj;
      }

      if(debug){
    	cout << "\nLikelihood so far:\n";
    	for(int i=0; i<rtree.ntotn; ++i){
    	  for(int j=0; j<nstate; ++j){
    	    cout << "\t" << L_sk_k[i][j];
    	  }
    	  cout << endl;
    	}
      }
    }

    if(debug) cout << "Likelihood char : " << nc << "\t" << L_sk_k[Ns+1][2] << endl;
    if(L_sk_k[Ns+1][2] > 0) logL += log(L_sk_k[Ns+1][2]);
    else logL += LARGE_LNL;
  }

  if(debug) cout << "Final likelihood: " << logL << endl;
  return logL;
}


// Get the likelihood on one site of a chromosme
// z: possible changes in copy number caused by chromosme gain/loss
double get_likelihood_site(vector<vector<double>>& L_sk_k, evo_tree& rtree, const vector<int>& knodes, map<double, double*>& pmats, const int has_wgd, const int z, const int model, int cn_max){
  int debug = 0;
  int nstate = cn_max + 1;

  if(debug){
      cout << "Computing likelihood for one site" << endl;
  }

  for(int kn=0; kn<knodes.size(); ++kn){
    int k = knodes[kn];
    int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

    if(debug) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni+1 << " , " << bli << "\t" <<  nj+1 << " , " << blj << endl;

    //loop over possible values of sk
    for(int sk=0; sk<nstate; ++sk){
      int nsk = sk;  // state after changes by other large scale events
      if(has_wgd == 1) nsk = 2 * sk;
      nsk += z;
      // cout << "likelihood for state " << nsk << endl;
      if(nsk < 0 || nsk >= nstate) continue;
      double Li = 0;
      // loop over possible si
      for(int si=0; si<nstate; ++si){
          if (model == 0){
            Li += get_transition_prob(rtree.mu, bli, nsk, si) * L_sk_k[ni][si];
          }
          if (model == 1){
            // if(debug){
            //     cout << "1st branch length " << bli << endl;
            //     r8mat_print( nstate, nstate, pmats[bli], " Pmatrix for 1st branch length " );
            // }
            Li += get_transition_prob_bounded(pmats[bli], nsk, si, cn_max) * L_sk_k[ni][si];
          }
        //cout << "\tscoring: Li\t" << li << "\t" << get_transition_prob(mu, bli, sk, si ) << "\t" << L_sk_k[ni][si] << endl;
      }
      double Lj = 0;
      // loop over possible sj
      for(int sj=0; sj<nstate; ++sj){
          if (model == 0){
               Lj += get_transition_prob(rtree.mu, blj, nsk, sj) * L_sk_k[nj][sj];
          }
          if (model == 1){
              // if(debug){
              //     cout << "2nd branch length " << blj << endl;
              //     r8mat_print( nstate, nstate, pmats[blj], " Pmatrix for 2nd branch length ");
              // }
               Lj += get_transition_prob_bounded(pmats[blj], nsk, sj, cn_max) * L_sk_k[nj][sj];
          }
      }
      //cout << "scoring: sk" << sk << "\t" <<  Li << "\t" << Lj << endl;
      L_sk_k[k][nsk] = Li*Lj;
    }
    if(debug){
      cout << "\nLikelihood so far for one site:\n";
      for(int i=0; i<rtree.ntotn; ++i){
        for(int j=0; j<nstate; ++j){
          cout << "\t" << L_sk_k[i][j];
        }
        cout << endl;
      }
    }
  }
}


// Get the likelihood on a set of chromosmes
double get_likelihood_chr(map<int, vector<vector<int>>>& vobs, evo_tree& rtree, const vector<int>& knodes, map<double, double*>& pmats, const int has_wgd, const int model, const int cn_max){
    int debug = 0;
    int nstate = cn_max + 1;

    double logL = 0;    // for all chromosmes
    double chr_gain = rtree.chr_gain_rate;
    double chr_loss = rtree.chr_loss_rate;
    if(debug){
        cout << "chromosme gain rate " << chr_gain << endl;
        cout << "chromosme loss rate " << chr_loss << endl;
        cout << "Number of chr so far " << vobs.size() << endl;
    }

    for(int nchr=1; nchr<=vobs.size(); nchr++){
      // cout << "Chr " << nchr << endl;
      double chr_logL = 0;  // for one chromosme
      double chr_logL_normal = 0, chr_logL_gain = 0, chr_logL_loss = 0;
      double site_logL = 0;   // log likelihood for all sites on a chromosme
      int z = 0;    // no chr gain/loss
      // cout << " chromosme number change is " << 0 << endl;
      for(int nc=0; nc<vobs[nchr].size(); nc++){
          // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
          // for each site of the chromosme
          vector<int> obs = vobs[nchr][nc];
          // construct a table for each state of each node
          vector<vector<double>> L_sk_k( rtree.ntotn, vector<double>(nstate,0) );
          for(int i=0; i<Ns; ++i){
              for(int j=0; j<nstate; ++j){
                     if(j == obs[i])  L_sk_k[i][j] = 1.0;
              }
          }
          // set unaltered
          L_sk_k[Ns][2] = 1.0;

          if(debug){
            cout << "\nLikelihood for tips:\n";
            for(int i=0; i<rtree.ntotn; ++i){
            for(int j=0; j<nstate; ++j){
              cout << "\t" << L_sk_k[i][j];
            }
            cout << endl;
            }
          }

          get_likelihood_site(L_sk_k, rtree, knodes, pmats, has_wgd, z, model, cn_max);

          if(debug) cout << "Likelihood char : " << nc << "\t" << L_sk_k[Ns+1][2] << endl;
          if(L_sk_k[Ns+1][2] > 0) site_logL += log(L_sk_k[Ns+1][2]);
          else site_logL += LARGE_LNL;

          if(debug){
            cout << "\nLikelihood so far:\n";
            for(int i=0; i<rtree.ntotn; ++i){
                for(int j=0; j<nstate; ++j){
                  cout << "\t" << L_sk_k[i][j];
                }
                cout << endl;
            }
         }
      }

      double chr_normal = 1;
      if(abs(chr_loss-0) > SMALL_PROB){
         chr_normal -= chr_loss;
      }
      if(abs(chr_gain-0) > SMALL_PROB){
         chr_normal -= chr_gain;
      }

      chr_logL_normal = log(chr_normal) + site_logL;
      chr_logL += chr_logL_normal;
      if(debug){
          cout << "Likelihood without chr gain/loss for " << nchr << " is "  << chr_normal << endl;
          cout << "Site Likelihood for " << nchr << " is "  << site_logL << endl;
          cout << "Likelihood without chr gain/loss: " << chr_logL_normal << endl;
      }

      if(abs(chr_loss-0) > SMALL_PROB){
          z = -1;
          double site_logL = 0;   // log likelihood for all sites on a chromosme
          // cout << " chromosme number change is " << z << endl;
          for(int nc=0; nc<vobs[nchr].size(); nc++){
              // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
              // for each site of the chromosme
              vector<int> obs = vobs[nchr][nc];
              // construct a table for each state of each node
              vector<vector<double>> L_sk_k( rtree.ntotn, vector<double>(nstate,0) );
              for(int i=0; i<Ns; ++i){
                  for(int j=0; j<nstate; ++j){
                         if(j == obs[i])  L_sk_k[i][j] = 1.0;
                  }
              }
              // set unaltered
              L_sk_k[Ns][2] = 1.0;

              if(debug){
                cout << "\nLikelihood for tips:\n";
                for(int i=0; i<rtree.ntotn; ++i){
                for(int j=0; j<nstate; ++j){
                  cout << "\t" << L_sk_k[i][j];
                }
                cout << endl;
                }
              }

              get_likelihood_site(L_sk_k, rtree, knodes, pmats, has_wgd, z, model, cn_max);

              if(debug) cout << "Likelihood char : " << nc << "\t" << L_sk_k[Ns+1][2] << endl;
              if(L_sk_k[Ns+1][2] > 0) site_logL += log(L_sk_k[Ns+1][2]);
              else site_logL += LARGE_LNL;

              if(debug){
                cout << "\nLikelihood so far:\n";
                for(int i=0; i<rtree.ntotn; ++i){
                for(int j=0; j<nstate; ++j){
                  cout << "\t" << L_sk_k[i][j];
                }
                cout << endl;
                }
              }
          } // for all sites on a chromosme

          chr_logL_loss = log(chr_loss) + site_logL;
          chr_logL +=  log(1 + exp(chr_logL_loss-chr_logL_normal));
          if(debug){
              cout << "\nLikelihood before chr loss for " << nchr << " is " << site_logL << endl;
              cout << "\nLikelihood after chr loss: " << chr_logL_loss << endl;
          }
      } // for all chromosme loss

      if(abs(chr_gain-0) > SMALL_PROB){
          z = 1;
          double site_logL = 0;   // log likelihood for all sites on a chromosme
          // cout << " chromosme number change is " << z << endl;
          for(int nc=0; nc<vobs[nchr].size(); nc++){
              // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
              // for each site of the chromosme
              vector<int> obs = vobs[nchr][nc];
              // construct a table for each state of each node
              vector<vector<double>> L_sk_k( rtree.ntotn, vector<double>(nstate,0) );
              for(int i=0; i<Ns; ++i){
                  for(int j=0; j<nstate; ++j){
                         if(j == obs[i])  L_sk_k[i][j] = 1.0;
                  }
              }
              // set unaltered
              L_sk_k[Ns][2] = 1.0;

              if(debug){
                cout << "\nLikelihood for tips:\n";
                for(int i=0; i<rtree.ntotn; ++i){
                for(int j=0; j<nstate; ++j){
                  cout << "\t" << L_sk_k[i][j];
                }
                cout << endl;
                }
              }

              get_likelihood_site(L_sk_k, rtree, knodes, pmats, has_wgd, z, model, cn_max);

              if(debug) cout << "Likelihood char : " << nc << "\t" << L_sk_k[Ns+1][2] << endl;
              if(L_sk_k[Ns+1][2] > 0) site_logL += log(L_sk_k[Ns+1][2]);
              else site_logL += LARGE_LNL;

              if(debug){
                cout << "\nLikelihood so far:\n";
                for(int i=0; i<rtree.ntotn; ++i){
                for(int j=0; j<nstate; ++j){
                  cout << "\t" << L_sk_k[i][j];
                }
                cout << endl;
                }
              }
          } // for all sites on a chromosme

          chr_logL_gain = log(chr_gain) + site_logL;
          if(chr_logL_loss > 0){
              chr_logL += log(1 + 1 / (exp(chr_logL_normal-chr_logL_gain) + exp(chr_logL_loss-chr_logL_gain)));
          }
          else{
              chr_logL += log(1 + exp(chr_logL_gain-chr_logL_normal));
          }

          if(debug){
              cout << "\nLikelihood before chr gain for " << nchr << " is " << site_logL << endl;
              cout << "\nLikelihood after chr gain: " << chr_logL_gain << endl;
          }
      } // for all chromosme loss
      // chr_logL = chr_logL_normal + log(1 + exp(chr_logL_loss-chr_logL_normal)) + log(1 + 1 / (exp(chr_logL_normal-chr_logL_gain) + exp(chr_logL_loss-chr_logL_gain)));
      logL += chr_logL;
      if(debug){
          cout << "\nLikelihood with chr gain/loss for one chromosme: " << logL << endl;
      }
    } // for each chromosme
    if(debug){
        cout << "\nLikelihood with chr gain/loss for all chromosmes: " << logL << endl;
    }
    return logL;
}

// Compute the likelihood of dummy sites consisting entirely of 2s for the tree
double get_likelihood_bias(map<int, vector<vector<int>>>& vobs, evo_tree& rtree, const vector<int>& knodes, map<double, double*>& pmats, const int model, const int cn_max){
    int debug = 0;
    int nstate = cn_max + 1;
    double logL = 0;

    for(int nchr=1; nchr<=vobs.size(); nchr++){
        double site_logL = 0;   // log likelihood for all sites on a chromosme
        for(int nc=0; nc<vobs[nchr].size(); nc++){
            vector<vector<double>> L_sk_k( rtree.ntotn, vector<double>(nstate,0) );
            for(int i=0; i<Ns; ++i){
                L_sk_k[i][2] = 1.0;
            }
            // set unaltered
            L_sk_k[Ns][2] = 1.0;

          for(int kn=0; kn<knodes.size(); ++kn){
            int k = knodes[kn];
            int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
            double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
            int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
            double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

            // if(debug) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni+1 << " , " << bli << "\t" <<  nj+1 << " , " << blj << endl;

            //loop over possible values of sk
            for(int sk=0; sk<nstate; ++sk){
              double Li = 0;
              // loop over possible si
              for(int si=0; si<nstate; ++si){
                  if (model == 0){
                    Li += get_transition_prob(rtree.mu, bli, sk, si) * L_sk_k[ni][si];
                  }
                  if (model == 1){
                    Li += get_transition_prob_bounded(pmats[bli], sk, si, cn_max) * L_sk_k[ni][si];
                  }
              }
              double Lj = 0;
              // loop over possible sj
              for(int sj=0; sj<nstate; ++sj){
                  if (model == 0){
                       Lj += get_transition_prob(rtree.mu, blj, sk, sj) * L_sk_k[nj][sj];
                  }
                  if (model == 1){
                       Lj += get_transition_prob_bounded(pmats[blj], sk, sj, cn_max) * L_sk_k[nj][sj];
                  }
              }
              L_sk_k[k][sk] = Li*Lj;
              // if(debug){
              //     cout << "scoring: sk " << sk << "\t" << Li << "\t" << Lj << "\t" << L_sk_k[k][sk] << endl;
              // }
            }
        }

          if(L_sk_k[Ns+1][2] > 0) site_logL += log(L_sk_k[Ns+1][2]);
          else  site_logL += 0;
      }
      logL += site_logL;
    }
    // assert(exp(logL) < 1);
    double lkl = exp(logL);
    double bias = 0;
    if(lkl < 1)
        bias = log(1-lkl);
    if(debug){
        cout << "The log likelihood of invariant sites is " << logL << endl;
        cout << "The bias of invariant sites is " << bias << endl;
    }
    return bias;
}


// Incorporate chromosme gain/loss and WGD
double get_likelihood_revised(const int& Ns, const int& Nchar, const int& num_invar_bins, map<int, vector<vector<int>>>& vobs, evo_tree& rtree, int model, int cons, const int cn_max, int correct_bias = 0){
  int debug = 0;
  if(debug) cout << "\tget_likelihood" << endl;

  // return 0 if the tree is not valid
  if(!is_tree_valid(rtree, cons)){
      return SMALL_LNL;
  }

  int nstate = cn_max + 1;
  double qmat[(nstate)*(nstate)]={0};
  if(model == 1){
      if(debug){
          cout << "Getting rate matrix" << endl;
      }
      get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
  }

  //create a list of nodes to loop over, making sure the root is last
  vector<int> knodes;
  for(int k=Ns+2; k<rtree.ntotn; ++k) knodes.push_back( k );
  knodes.push_back(Ns+1);

  // Find the transition probability matrix for each branch
  map<double, double*> pmats;
  for(int kn=0; kn<knodes.size(); ++kn){
    int k = knodes[kn];
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
    if(model==1){
        double pmati[(nstate)*(nstate)]={0};
        double pmatj[(nstate)*(nstate)]={0};

        if(pmats.count(bli) == 0){
            get_transition_matrix_bounded(qmat, pmati, bli, cn_max);
            pmats[bli] = new double[(nstate)*(nstate)];
            for(int i=0; i<nstate*nstate; i++){
                pmats[bli][i] = pmati[i];
            }
        }
        if(pmats.count(blj) == 0){
            get_transition_matrix_bounded(qmat, pmatj, blj, cn_max);
            pmats[blj] = new double[(nstate)*(nstate)];
            for(int i=0; i<nstate*nstate; i++){
                pmats[blj][i] = pmatj[i];
            }
        }

        // if(debug){
        //      cout << "Get Pmatrix for 1st branch length " << bli << endl;
        //      r8mat_print( nstate, nstate, pmats[bli], "  P matrix after change in map:" );
        //      r8mat_print( nstate, nstate, pmati, "  P matrix after change:" );
        //      cout << "Get Pmatrix for 2nd branch length " << blj << endl;
        //      r8mat_print( nstate, nstate, pmats[blj], "  P matrix after change in map:" );
        //      r8mat_print( nstate, nstate, pmatj, "  P matrix after change:" );
        // }
     }
  }

  if(debug){
      for( auto it = pmats.begin(); it != pmats.end(); ++it )
      {
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print( nstate, nstate, it->second, "  P matrix:" );
      }
  }


  double logL = 0;

  double wgd = rtree.wgd_rate;
  // cout << "WGD rate is " << wgd << endl;
  if(abs(wgd-0) > SMALL_PROB){
       // cout << "Computing the likelihood with consideration of WGD" << endl;
       logL += (1-wgd) * get_likelihood_chr(vobs, rtree, knodes, pmats, 0, model, cn_max);
       logL += wgd * get_likelihood_chr(vobs, rtree, knodes, pmats, 1, model, cn_max);
  }
  else{
      logL += get_likelihood_chr(vobs, rtree, knodes, pmats, 0, model, cn_max);
  }


  if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;
  if(correct_bias == 1){
      double bias = get_likelihood_bias(vobs, rtree, knodes, pmats, model, cn_max);
      if(abs(bias-0) > SMALL_PROB)
        logL = logL - Nchar * bias;
  }

  if(debug) cout << "Final likelihood: " << logL << endl;

  return logL;
}


//
// Maximum likelihood and optimisation functions
//

double my_f (const gsl_vector *v, void *params){
  evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current branch lengths
  vector<edge> enew;
  for(int i=0; i<(tree->nedge); ++i){
    enew.push_back( tree->edges[i] );
  }
  for(int i=0; i<(tree->nedge)-1; ++i){
    enew[i].length = exp( gsl_vector_get(v, i) );
  }
  evo_tree new_tree(Ns+1, enew);
  if(model==0){
      new_tree.mu = tree->mu;
  }
  if(model==1){
      new_tree.dup_rate = tree->dup_rate;
      new_tree.del_rate = tree->del_rate;
      new_tree.chr_gain_rate = tree->chr_gain_rate;
      new_tree.chr_loss_rate = tree->chr_loss_rate;
      new_tree.wgd_rate = tree->wgd_rate;
  }

  // return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, model, 0, cn_max);
  return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, new_tree, model, 0, cn_max, correct_bias);
}

double my_f_mu (const gsl_vector *v, void *params){
  //cout << "my_f_mu" << endl;
  evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current branch lengths
  vector<edge> enew;
  for(int i=0; i<(tree->nedge); ++i){
    enew.push_back( tree->edges[i] );
  }
  for(int i=0; i<(tree->nedge)-1; ++i){
    enew[i].length = exp( gsl_vector_get(v, i) );
  }
  evo_tree new_tree(Ns+1, enew);
  if(model==0){
      new_tree.mu = exp( gsl_vector_get(v,tree->nedge-1) );  // 0 to nedge-2 are epars, nedge-1 is mu
  }
  if(model==1){
      new_tree.dup_rate = exp( gsl_vector_get(v,tree->nedge-1) );
      new_tree.del_rate = exp( gsl_vector_get(v,tree->nedge) );
      new_tree.chr_gain_rate = exp( gsl_vector_get(v,tree->nedge+1) );
      new_tree.chr_loss_rate = exp( gsl_vector_get(v,tree->nedge+2) );
      new_tree.wgd_rate = exp( gsl_vector_get(v,tree->nedge+3) );
  }

  // return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, model, 0, cn_max);
  return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, new_tree, model, 0, cn_max, correct_bias);
}

double my_f_cons (const gsl_vector *v, void *params){
  evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current parameters observing timing constraints
  vector<edge> enew;
  for(int i=0; i<(tree->nedge); ++i){
    enew.push_back( tree->edges[i] );
  }

  // parameters coming in are internal branch lengths followed by total time
  int count = 0;
  for(int i=0; i<(tree->nedge)-1; ++i){
    if(enew[i].end > Ns){
      enew[i].length = exp( gsl_vector_get(v, count) );
      count++;
    }else{
      enew[i].length = 0;
    }
  }

  evo_tree new_tree(tree->nleaf, enew, exp( gsl_vector_get(v, count)), tree->tobs);
  if(model==0){
      new_tree.mu = tree->mu;
  }
  if(model==1){
      new_tree.dup_rate = tree->dup_rate;
      new_tree.del_rate = tree->del_rate;
      new_tree.chr_gain_rate = tree->chr_gain_rate;
      new_tree.chr_loss_rate = tree->chr_loss_rate;
      new_tree.wgd_rate = tree->wgd_rate;
  }

  // return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, model, 1, cn_max);
  return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, new_tree, model, 1, cn_max, correct_bias);
}

double my_f_cons_mu (const gsl_vector *v, void *params){
  evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current parameters observing timing constraints
  vector<edge> enew;
  for(int i=0; i<(tree->nedge); ++i){
    enew.push_back( tree->edges[i] );
  }

  // parameters coming in are internal branch lengths followed by mutation rate
  int count = 0;
  for(int i=0; i<(tree->nedge)-1; ++i){
    if(enew[i].end > Ns){
      enew[i].length = exp( gsl_vector_get(v, count) );
      count++;
    }else{
      enew[i].length = 0;
    }
  }

  evo_tree new_tree(tree->nleaf, enew, exp( gsl_vector_get(v, count)), tree->tobs);
  if(model==0){
      new_tree.mu = exp( gsl_vector_get(v, count+1) );
  }
  if(model==1){
      new_tree.dup_rate = exp( gsl_vector_get(v,count+1) );
      new_tree.del_rate = exp( gsl_vector_get(v,count+2) );
      new_tree.chr_gain_rate = exp( gsl_vector_get(v,count+3) );
      new_tree.chr_loss_rate = exp( gsl_vector_get(v,count+4) );
      new_tree.wgd_rate = exp( gsl_vector_get(v,count+5) );
  }

  // return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, model, 1, cn_max);
  return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, new_tree, model, 1, cn_max, correct_bias);
}

// Output error informaton without aborting
void my_err_handler(const char * reason, const char * file, int line, int gsl_errno){
    fprintf (stderr, "failed, %s, reason: %s, file %s line %d \n", gsl_strerror(gsl_errno), reason, file, line);
    return;
}

void get_variables(evo_tree& rtree, int model, int cons, int maxj, double *x){
    int debug = 0;
    // create a new tree from current value of parameters
    vector<edge> enew;
    for(int i=0; i<rtree.nedge; ++i){
      enew.push_back(rtree.edges[i]);
    }

    if(cons == 0){
      // The first element of x is not used for optimization
      // The index of x is added by 1 compared with index for simplex method
      for(int i=0; i<rtree.nedge-1; ++i){
          enew[i].length = x[i + 1];
      }
      evo_tree new_tree(Ns+1, enew);
      if(maxj==0){
          if(model == 0){
              new_tree.mu = rtree.mu;
          }
          if(model == 1){
              new_tree.dup_rate = rtree.dup_rate;
              new_tree.del_rate = rtree.del_rate;
              new_tree.chr_gain_rate = rtree.chr_gain_rate;
              new_tree.chr_loss_rate = rtree.chr_loss_rate;
              new_tree.wgd_rate = rtree.wgd_rate;
              new_tree.mu = new_tree.dup_rate + new_tree.del_rate + new_tree.chr_gain_rate + new_tree.chr_loss_rate + new_tree.wgd_rate;
          }
      }else{
          if(model == 0){
              new_tree.mu = x[rtree.nedge];
              if(debug){
                  for(int i=0; i<rtree.nedge+1; i++){ cout << x[i] << '\n';}
                  cout << "mu value so far: " << new_tree.mu << endl;
              }
          }
          if(model == 1){
              new_tree.dup_rate = x[rtree.nedge];
              new_tree.del_rate = x[rtree.nedge+1];
              new_tree.chr_gain_rate = x[rtree.nedge+2];
              new_tree.chr_loss_rate = x[rtree.nedge+3];
              new_tree.wgd_rate = x[rtree.nedge+4];
              new_tree.mu = new_tree.dup_rate + new_tree.del_rate + new_tree.chr_gain_rate + new_tree.chr_loss_rate + new_tree.wgd_rate;
              if(debug){
                  for(int i=0; i<rtree.nedge+2; i++){ cout << x[i] << '\n';}
                  cout << "dup_rate value so far: " << new_tree.dup_rate << endl;
                  cout << "del_rate value so far: " << new_tree.del_rate << endl;
                  cout << "chr_gain_rate value so far: " << new_tree.chr_gain_rate << endl;
                  cout << "chr_loss_rate value so far: " << new_tree.chr_loss_rate << endl;
                  cout << "wgd_rate value so far: " << new_tree.wgd_rate << endl;
              }
          }
      }
      rtree = evo_tree(new_tree);
    }else{
      int count = 0;
      for(int i=0; i<rtree.nedge-1; ++i){
        if(enew[i].end > Ns){
          	enew[i].length = x[count + 1];
          	count++;
        }else{
  	        enew[i].length = 0;
        }
      }
      if(debug){
          cout << "total height so far: " << x[count + 1] << endl;
      }
      evo_tree new_tree(rtree.nleaf, enew, x[count + 1], rtree.tobs);
      if(maxj==0){
        if(model == 0){
            new_tree.mu = rtree.mu;
        }
        if(model == 1){
            new_tree.dup_rate = rtree.dup_rate;
            new_tree.del_rate = rtree.del_rate;
            new_tree.chr_gain_rate = rtree.chr_gain_rate;
            new_tree.chr_loss_rate = rtree.chr_loss_rate;
            new_tree.wgd_rate = rtree.wgd_rate;
            new_tree.mu = new_tree.dup_rate + new_tree.del_rate + new_tree.chr_gain_rate + new_tree.chr_loss_rate + new_tree.wgd_rate;
        }
      }else{
        if(model == 0){
            new_tree.mu = x[rtree.nintedge+2];
            if(debug){
                for(int i=0; i<=rtree.nintedge+2; i++){ cout << x[i] << '\n';}
                cout << "mu value so far: " << new_tree.mu << endl;
            }
        }
        if(model == 1){
            new_tree.dup_rate = x[rtree.nintedge+2];
            new_tree.del_rate = x[rtree.nintedge+3];
            new_tree.chr_gain_rate = x[rtree.nintedge+4];
            new_tree.chr_loss_rate = x[rtree.nintedge+5];
            new_tree.wgd_rate = x[rtree.nintedge+6];
            new_tree.mu = new_tree.dup_rate + new_tree.del_rate + new_tree.chr_gain_rate + new_tree.chr_loss_rate + new_tree.wgd_rate;
            if(debug){
                for(int i=0; i<=rtree.nintedge+6; i++){ cout << x[i] << '\n';}
                cout << "dup_rate value so far: " << new_tree.dup_rate << endl;
                cout << "del_rate value so far: " << new_tree.del_rate << endl;
                cout << "chr_gain_rate value so far: " << new_tree.chr_gain_rate << endl;
                cout << "chr_loss_rate value so far: " << new_tree.chr_loss_rate << endl;
                cout << "wgd_rate value so far: " << new_tree.wgd_rate << endl;
            }
         }
      }
      rtree = evo_tree(new_tree);
    }
}

/**
    the target function which needs to be optimized
    @param x the input vector x
    @return the function value at x
*/
double targetFunk(evo_tree& rtree, const int model, const int cons, const int maxj, const int cn_max, const int correct_bias, double x[]) {
    // negative log likelihood
    get_variables(rtree, model, cons, maxj, x);
    // return -1.0*get_likelihood(Ns, Nchar, vobs, rtree, model, cons);
    return -1.0*get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, rtree, model, cons, cn_max, correct_bias);
}

/**
	the approximated derivative function
	@param x the input vector x
	@param dfx the derivative at x
	@return the function value at x
*/
double derivativeFunk(evo_tree& rtree, const int model, const int cons, const int maxj, const int cn_max, const int correct_bias, int ndim, double x[], double dfx[]) {
    int debug = 0;
	double *h = new double[ndim+1];
    double temp;
    int dim;
	double fx = targetFunk(rtree, model, cons, maxj, cn_max, correct_bias, x);
	for (dim = 1; dim <= ndim; dim++ ){
		temp = x[dim];
		h[dim] = ERROR_X * fabs(temp);
		if (h[dim] == 0.0) h[dim] = ERROR_X;
		x[dim] = temp + h[dim];
		h[dim] = x[dim] - temp;
		dfx[dim] = (targetFunk(rtree, model, cons, maxj, cn_max, correct_bias, x));
		x[dim] = temp;
	}
	for (dim = 1; dim <= ndim; dim++ ){
        dfx[dim] = (dfx[dim] - fx) / h[dim];
        if(debug){
            cout << "dfx[dim] " << dfx[dim] << endl;
        }
    }
    delete [] h;
	return fx;
}

double optimFunc(evo_tree& rtree, const int model, const int cons, const int maxj, const int cn_max, const int correct_bias, int nvar, double *vars) {
    return targetFunk(rtree, model, cons, maxj, cn_max, correct_bias, vars-1);
}

double optimGradient(evo_tree& rtree, const int model, const int cons, const int maxj, const int cn_max, const int correct_bias, int nvar, double *x, double *dfx) {
    return derivativeFunk(rtree, model, cons, maxj, cn_max, correct_bias, nvar, x-1, dfx-1);
//    const double ERRORX = 1e-5;
//	double fx = optimFunc(nvar, x);
//	double h, temp;
//	for (int dim = 0; dim <= nvar; dim++ ){
//		temp = x[dim];
//		h = ERRORX * fabs(temp);
//		if (h == 0.0) h = ERRORX;
//		x[dim] = temp + h;
//		h = x[dim] - temp;
//		dfx[dim] = (optimFunc(nvar, x) - fx) / h;
//		x[dim] = temp;
//	}
//	return fx;
}


void lbfgsb(evo_tree& rtree, const int model, const int cons, const int maxj, const int cn_max, const int correct_bias, int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, int *fail,
		double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT)
{
	char task[60];
	double f, *g, dsave[29], *wa;
	int tr = -1, iter = 0, *iwa, isave[44], lsave[4];

	/* shut up gcc -Wall in 4.6.x */

	for(int i = 0; i < 4; i++) lsave[i] = 0;

	if(n == 0) { /* not handled in setulb */
		*fncount = 1;
		*grcount = 0;
		*Fmin = optimFunc(rtree, model, cons, maxj, cn_max, correct_bias, n, u);
		strcpy(msg, "NOTHING TO DO");
		*fail = 0;
		return;
	}
	if (nREPORT <= 0) {
		cerr << "REPORT must be > 0 (method = \"L-BFGS-B\")" << endl;
		exit(1);
	}
	switch(trace) {
	case 2: tr = 0; break;
	case 3: tr = nREPORT; break;
	case 4: tr = 99; break;
	case 5: tr = 100; break;
	case 6: tr = 101; break;
	default: tr = -1; break;
	}

	*fail = 0;
	g = (double*) malloc (n * sizeof(double));
	/* this needs to be zeroed for snd in mainlb to be zeroed */
	wa = (double *) malloc((2*m*n+4*n+11*m*m+8*m) * sizeof(double));
	iwa = (int *) malloc(3*n * sizeof(int));
	strcpy(task, "START");
	while(1) {
		/* Main workhorse setulb() from ../appl/lbfgsb.c : */
		setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task, tr, lsave, isave, dsave);
		/*    Rprintf("in lbfgsb - %s\n", task);*/
		if (strncmp(task, "FG", 2) == 0) {
			f = optimGradient(rtree, model, cons, maxj, cn_max, correct_bias, n, x, g);
			if (!isfinite(f)) {
				cerr << "L-BFGS-B needs finite values of 'fn'" << endl;
				exit(1);
			}

		} else if (strncmp(task, "NEW_X", 5) == 0) {
			iter++;
			if(trace == 1 && (iter % nREPORT == 0)) {
				cout << "iter " << iter << " value " << f << endl;
			}
			if (iter > maxit) {
				*fail = 1;
				break;
			}
		} else if (strncmp(task, "WARN", 4) == 0) {
			*fail = 51;
			break;
		} else if (strncmp(task, "CONV", 4) == 0) {
			break;
		} else if (strncmp(task, "ERROR", 5) == 0) {
			*fail = 52;
			break;
		} else { /* some other condition that is not supposed to happen */
			*fail = 52;
			break;
		}
	}
	*Fmin = f;
	*fncount = *grcount = isave[33];
	if (trace) {
		cout << "final value " << *Fmin << endl;
		if (iter < maxit && *fail == 0)
			cout << "converged" << endl;
		else
			cout << "stopped after " << iter << " iterations\n";
	}
	strcpy(msg, task);
	free(g);
	free(wa);
	free(iwa);
}


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
double L_BFGS_B(evo_tree& rtree, const int model, const int cons, const int maxj, const int cn_max, const int correct_bias, int n, double* x, double* l, double* u, double pgtol, int maxit) {
    int debug = 0;
    int i;
	double Fmin;
	int fail;
	int fncount;
	int grcount;
	char msg[100];

	int m = 10;          // number of BFGS updates retained in the "L-BFGS-B" method. It defaults to 5.

	int *nbd;           // 0: unbounded; 1: lower bounded; 2: both lower & upper; 3: upper bounded
	nbd = new int[n];
	for (i=0; i<n; i++)
		nbd[i] = 2;

	double factr = 1e+7; // control the convergence of the "L-BFGS-B" method.
	// Convergence occurs when the reduction in the object is within this factor
	// of the machine tolerance.
	// Default is 1e7, that is a tolerance of about 1e-8

//	double pgtol = 0;   // helps control the convergence of the "L-BFGS-B" method.
//    pgtol = 0.0;
	// It is a tolerance on the projected gradient in the current search direction.
	// Default is zero, when the check is suppressed

	int trace = 0;      // non-negative integer.
    if (debug)
        trace = 1;
	// If positive, tracing information on the progress of the optimization is produced.
	// Higher values may produce more tracing information.

	int nREPORT = 10;   // The frequency of reports for the "L-BFGS-B" methods if "trace" is positive.
	// Defaults to every 10 iterations.

/*#ifdef USE_OLD_PARAM
	lbfgsb(n, m, x, l, u, nbd, &Fmin, fn, gr1, &fail, ex,
			factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
#else*/

	lbfgsb(rtree, model, cons, maxj, cn_max, correct_bias, n, m, x, l, u, nbd, &Fmin, &fail,
			factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
//#endif

    if (fail == 51 || fail == 52) {
        cout << msg << endl;
    }

	delete[] nbd;

    return Fmin;
}


// Using BFGS method to get the maximum likelihood with lower and upper bounds
 evo_tree max_likelihood_BFGS(evo_tree& rtree, int model, double& minL, const double tolerance, const int miter,  int cons=0, int maxj=0, int cn_max=4, int correct_bias=1){
    int debug = 0;
    int npar_ne;
    int ndim = 0;
    double *variables, *upper_bound, *lower_bound;
    int i;

    // Set variables
    if(cons == 0){
      npar_ne = rtree.nedge-1;
      if(maxj==0){
        ndim = npar_ne;
      }else{
          if(model==0){
              ndim = npar_ne + 1;
          }
          if(model==1){
              ndim = npar_ne + 5;
          }
      }
      variables = new double[ndim+1];
      upper_bound = new double[ndim+1];
      lower_bound = new double[ndim+1];
      // initialise the best guess branch length and mu if required
      for(i=0; i<npar_ne; ++i){
        variables[i+1] = rtree.edges[i].length;
        lower_bound[i+1] = 0;
        upper_bound[i+1] = age + *max_element(rtree.tobs.begin(), rtree.tobs.end());
      }
      if(maxj==1){
          if(model == 0){
              i = npar_ne;
              variables[i+1] = rtree.mu;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;
          }
          if(model == 1){
              i = npar_ne;
              variables[i+1] = rtree.dup_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;

              i = npar_ne+1;
              variables[i+1] = rtree.del_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;

              i = npar_ne+2;
              variables[i+1] = rtree.chr_gain_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;

              i = npar_ne+3;
              variables[i+1] = rtree.chr_loss_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;

              i = npar_ne+4;
              variables[i+1] = rtree.wgd_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;
          }
      }
    }else{
      npar_ne = rtree.nintedge + 1;
      if(maxj==0){
        ndim = npar_ne;
      }else{
        if(model==0){
            ndim = npar_ne + 1;
        }
        if(model==1){
            ndim = npar_ne + 5;
        }
      }
      // initialise with internal edges
      variables = new double[ndim+1];
      upper_bound = new double[ndim+1];
      lower_bound = new double[ndim+1];

      // initialise the best guess branch length and mu if required
      // rtree.print();
      for(i=0; i< rtree.nintedge; ++i){
        variables[i+1] = rtree.intedges[i]->length;
        lower_bound[i+1] = 0;
        // upper_bound[i+1] = age + *max_element(rtree.tobs.begin(), rtree.tobs.end());
        upper_bound[i+1] = age;
      }

      // initialise with current total tree time until first sample
      i = rtree.nintedge;
      variables[i+1] = rtree.get_total_time();
      // cout << "total time: " << variables[i+1] << endl;
      vector<int> internal_lens = rtree.get_internal_lengths();
      double max_ilen =  *max_element(internal_lens.begin(), internal_lens.end());
      // lower_bound[i+1] = *max_element(rtree.tobs.begin(), rtree.tobs.end());
      lower_bound[i+1] = max_ilen;
      upper_bound[i+1] = age;

      if(maxj==1){
          if(model == 0){
              i = npar_ne;
              variables[i+1] = rtree.mu;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;
          }
          if(model == 1){
              i = npar_ne;
              variables[i+1] = rtree.dup_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;

              i = npar_ne+1;
              variables[i+1] = rtree.del_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;

              i = npar_ne+2;
              variables[i+1] = rtree.chr_gain_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;

              i = npar_ne+3;
              variables[i+1] = rtree.chr_loss_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;

              i = npar_ne+4;
              variables[i+1] = rtree.wgd_rate;
              lower_bound[i+1] = 0;
              upper_bound[i+1] = 1;
          }
       }
    }

    // variables contains the parameters to estimate (branch length, mutation rate)
    minL = L_BFGS_B(rtree, model, cons, maxj, cn_max, correct_bias, ndim, variables+1, lower_bound+1, upper_bound+1, tolerance, miter);
    if (debug){
        // cout << "mu of current ML tree: " << rtree.mu << endl;
        cout << "lnL of current ML tree: " << minL << endl;
        rtree.print();
    }

    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables;

    // rtree has been updated in the optimization process
    return rtree;
}

// given a tree, maximise the branch lengths (and optionally mu) assuming branch lengths are independent or constrained in time
evo_tree max_likelihood(evo_tree& rtree, int model, double& minL, const double ssize, const double tolerance, const int miter, int cons=0, int maxj=0, int cn_max=4, int correct_bias=1){
  int debug = 0;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  /* save original handler, install new handler */
  gsl_error_handler_t  *old_handler = gsl_set_error_handler (&my_err_handler);

  int npar_ne;
  int npar;
  gsl_vector *x;

  if(cons == 0){
    npar_ne = rtree.nedge-1;
    if(maxj==0){
      npar = npar_ne;
    }else{
      if(model==0){
          npar = npar_ne + 1;
      }
      if(model==1){
          npar = npar_ne + 5;
      }
    }

    // initialise the best guess branch length and mu if required
    x = gsl_vector_alloc (npar);
    for(int i=0; i<npar_ne; ++i){
      gsl_vector_set (x, i, log(rtree.edges[i].length));
    }
    if(maxj==1){
        if(model==0){
            gsl_vector_set (x, npar_ne, log(rtree.mu));
        }
        if(model==1){
            gsl_vector_set (x, npar_ne, log(rtree.dup_rate) );
            gsl_vector_set (x, npar_ne+1, log(rtree.del_rate) );
            gsl_vector_set (x, npar_ne+2, log(rtree.chr_gain_rate) );
            gsl_vector_set (x, npar_ne+3, log(rtree.chr_loss_rate) );
            gsl_vector_set (x, npar_ne+4, log(rtree.wgd_rate) );
        }
    }

    if(maxj==0){
      minex_func.f = my_f;
    }else{
      minex_func.f = my_f_mu;
    }

  }else{
    npar_ne = rtree.nintedge + 1;
    if(maxj==0){
      npar = npar_ne;
    }else{
        if(model==0){
            npar = npar_ne + 1;
        }
        if(model==1){
            npar = npar_ne + 5;
        }
    }

    x = gsl_vector_alloc (npar);
    // initialise with internal edges
    for(int i=0; i<rtree.nintedge; ++i){
      gsl_vector_set (x, i, log(rtree.intedges[i]->length));
    }

    // initialise with current total tree time
    gsl_vector_set (x, rtree.nintedge, log(rtree.get_total_time()) );
    if(maxj==1){
        if(model==0){
            gsl_vector_set (x, npar_ne, log(rtree.mu) );
        }
        if(model==1){
            gsl_vector_set (x, npar_ne, log(rtree.dup_rate) );
            gsl_vector_set (x, npar_ne+1, log(rtree.del_rate) );
            gsl_vector_set (x, npar_ne+2, log(rtree.chr_gain_rate) );
            gsl_vector_set (x, npar_ne+3, log(rtree.chr_loss_rate) );
            gsl_vector_set (x, npar_ne+4, log(rtree.wgd_rate) );
        }
    }

    if(maxj==0){
      minex_func.f = my_f_cons;
    }else{
      minex_func.f = my_f_cons_mu;
    }
  }
  evo_tree *p = &rtree;
  void * pv = p;

  // Set initial step sizes to 1
  //cout << "numbers: nedge / nintedge / npar / x->size: " << rtree.nedge << "\t" << rtree.nintedge << "\t" << npar << "\t" << x->size << endl;
  gsl_vector* ss = gsl_vector_alloc (npar);
  gsl_vector_set_all (ss, ssize);

  // Initialize method and iterate
  minex_func.n = npar;
  minex_func.params = pv;

  s = gsl_multimin_fminimizer_alloc (T, npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  size_t iter = 0;
  int status;
  double size;
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status) break;

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, tolerance);
    // gsl_set_error_handler_off();

    if(debug){
      printf ("%5lu, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n",
	      iter,
	      exp(gsl_vector_get (s->x, 0)),
	      exp(gsl_vector_get (s->x, s->x->size-1)),
	      s->fval, size);
    }

    if (status == GSL_SUCCESS){
      if(debug){
        	printf ("converged to minimum at\n");

        	printf ("%5lu, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n",
        		iter,
        		exp(gsl_vector_get (s->x, 0)),
        		exp(gsl_vector_get (s->x, s->x->size-1)),
        		s->fval, size);
        }
    }
  }
  while (status == GSL_CONTINUE && iter < miter);

  if( status == GSL_CONTINUE ){
    cout << "### WARNING: maximum likelihood did not converge" << endl;
  }

  // create a new tree
  vector<edge> enew;
  for(int i=0; i<rtree.nedge; ++i){
    enew.push_back( rtree.edges[i] );
  }

  if(cons == 0){
    for(int i=0; i<npar_ne; ++i){
      enew[i].length = exp(gsl_vector_get(s->x, i));
    }
    evo_tree new_tree(Ns+1, enew);
    if(maxj==0){
        if(model==0){
            new_tree.mu = rtree.mu;
        }
        if(model==1){
            new_tree.dup_rate = rtree.dup_rate;
            new_tree.del_rate = rtree.del_rate;
            new_tree.chr_gain_rate = rtree.chr_gain_rate;
            new_tree.chr_loss_rate = rtree.chr_loss_rate;
            new_tree.wgd_rate = rtree.wgd_rate;
            new_tree.mu = new_tree.dup_rate + new_tree.del_rate + new_tree.chr_gain_rate + new_tree.chr_loss_rate + new_tree.wgd_rate;
        }
    }else{
        if(model==0){
            new_tree.mu = exp(gsl_vector_get(s->x, npar_ne));
        }
        if(model==1){
            new_tree.dup_rate = exp(gsl_vector_get(s->x, npar_ne));
            new_tree.del_rate = exp(gsl_vector_get(s->x, npar_ne+1));
            new_tree.chr_gain_rate = exp(gsl_vector_get(s->x, npar_ne+2));
            new_tree.chr_loss_rate = exp(gsl_vector_get(s->x, npar_ne+3));
            new_tree.wgd_rate = exp(gsl_vector_get(s->x, npar_ne+4));
            new_tree.mu = new_tree.dup_rate + new_tree.del_rate + new_tree.chr_gain_rate + new_tree.chr_loss_rate + new_tree.wgd_rate;
        }
    }

    minL = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    /* restore original handler */
    gsl_set_error_handler (old_handler);

    return new_tree;

  }else{
    int count = 0;
    for(int i=0; i<rtree.nedge-1; ++i){
      if(enew[i].end > Ns){
        	enew[i].length = exp( gsl_vector_get(s->x, count) );
        	count++;
      }else{
	        enew[i].length = 0;
      }
    }
    // evo_tree new_tree(rtree.nleaf, enew, age, rtree.tobs);
    evo_tree new_tree(rtree.nleaf, enew, exp( gsl_vector_get(s->x, count)), rtree.tobs);
    if(maxj==0){
        if(model==0){
            new_tree.mu = rtree.mu;
        }
        if(model==1){
            new_tree.dup_rate = rtree.dup_rate;
            new_tree.del_rate = rtree.del_rate;
            new_tree.chr_gain_rate = rtree.chr_gain_rate;
            new_tree.chr_loss_rate = rtree.chr_loss_rate;
            new_tree.wgd_rate = rtree.wgd_rate;
            new_tree.mu = new_tree.dup_rate + new_tree.del_rate + new_tree.chr_gain_rate + new_tree.chr_loss_rate + new_tree.wgd_rate;
        }
    }else{
        if(model==0){
            new_tree.mu = exp(gsl_vector_get(s->x, npar_ne));
        }
        if(model==1){
            new_tree.dup_rate = exp(gsl_vector_get(s->x, npar_ne));
            new_tree.del_rate = exp(gsl_vector_get(s->x, npar_ne+1));
            new_tree.chr_gain_rate = exp(gsl_vector_get(s->x, npar_ne+2));
            new_tree.chr_loss_rate = exp(gsl_vector_get(s->x, npar_ne+3));
            new_tree.wgd_rate = exp(gsl_vector_get(s->x, npar_ne+4));
            new_tree.mu = new_tree.dup_rate + new_tree.del_rate + new_tree.chr_gain_rate + new_tree.chr_loss_rate + new_tree.wgd_rate;
        }
    }

    minL = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    /* restore original handler */
    gsl_set_error_handler (old_handler);

    // if (!is_tree_valid(new_tree)){
    //     // keep the original tree
    //     return rtree;
    // }

    return new_tree;
  }

}
