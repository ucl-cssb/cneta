#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>

#include <algorithm>
#include <vector>

extern "C" {
  gsl_rng * r;
  vector<vector<int> > vobs;
  vector<double> tobs;
  int Ns;
  int Nchar;
}
  
const double LARGE_LNL = -1e9;

void setup_rng(int set_seed){

  gsl_rng_env_setup();

  const gsl_rng_type* T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  if( set_seed != 0 ){
    gsl_rng_set(r, set_seed );
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
    double t = gsl_ran_exponential(r, lambda);
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
  double t = gsl_ran_exponential(r, lambda);
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
    double t = gsl_ran_exponential(r, lambda);
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
  double t = gsl_ran_exponential(r, lambda);
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

  evo_tree ret(nsample+1, edges, lengths);
  return ret;
}

//
// Likelihood and transition probabilities
// 
double get_transition_prob(const double& mu, double blength, const int& sk, const int& sj ){
  //if(debug) cout << "\tget_transition_prob" << endl;

  // assume a basic Markov model here cf JC
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

// Likelihood of a tree
//
double get_likelihood(const int& Ns, const int& Nchar, const vector<vector<int> >& vobs, evo_tree& rtree){
  int debug = 0;
  if(debug) cout << "\tget_likelihood" << endl;
  int nstate = 5;

  double logL = 0;
  for(int nc=0; nc<Nchar; ++nc){
    vector<int> obs = vobs[nc];

    if(0){
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

    if(0){
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

      if(0) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni+1 << " , " << bli << "\t" <<  nj+1 << " , " << blj << endl;

      //loop over possible values of sk
      for(int sk=0; sk<nstate; ++sk){

	double Li = 0;
	// loop over possible si
	for(int si=0; si<nstate; ++si){
	  Li += get_transition_prob(rtree.mu, bli, sk, si ) * L_sk_k[ni][si];
	  //cout << "\tscoring: Li\t" << li << "\t" << get_transition_prob(mu, bli, sk, si ) << "\t" << L_sk_k[ni][si] << endl;
	}
	double Lj = 0;
	// loop over possible sj
	for(int sj=0; sj<nstate; ++sj){
	  Lj += get_transition_prob(rtree.mu, blj, sk, sj ) * L_sk_k[nj][sj];
	}
	//cout << "scoring: sk" << sk << "\t" <<  Li << "\t" << Lj << endl;
	L_sk_k[k][sk] = Li*Lj;
      }

      if(0){
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
  new_tree.mu = tree->mu;

  return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree);
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
  new_tree.mu = exp( gsl_vector_get(v,tree->nedge-1) );  // 0 to nedge-2 are epars, nedge-1 is mu

  return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree);
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
  new_tree.mu = tree->mu;

  return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree);
}

double my_f_cons_mu (const gsl_vector *v, void *params){

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
  new_tree.mu = exp( gsl_vector_get(v, count+1) );

  return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree);
}



// given a tree, maximise the branch lengths (and optionally mu) assuming branch lengths are independent or constrained in time
evo_tree max_likelihood(evo_tree& rtree, double& minL, const double ssize, const double tolerance, const int miter, int cons=0, int maxj=0){

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;

  int npar_ne;
  int npar;
  gsl_vector *x;

  if(cons == 0){
    npar_ne = rtree.nedge-1;
    if(maxj==0){
      npar = npar_ne;
    }else{
      npar = npar_ne + 1;
    }

    // initialise the best guess branch length and mu if required
    x = gsl_vector_alloc (npar);
    for(int i=0; i<npar_ne; ++i){
      gsl_vector_set (x, i, log(rtree.edges[i].length));
    }
    if(maxj==1) gsl_vector_set (x, npar_ne, log(rtree.mu) );

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
      npar = npar_ne + 1;
    }

    x = gsl_vector_alloc (npar);

    // initialise with internal edges
    for(int i=0; i<rtree.nintedge; ++i){
      gsl_vector_set (x, i, log(rtree.intedges[i]->length));
    }

    // initialise with current total tree time
    gsl_vector_set (x, rtree.nintedge, log(rtree.get_total_time()) );

    if(maxj==1) gsl_vector_set (x, npar_ne, log(rtree.mu) );

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

    if(0){
      printf ("%5lu, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n",
	      iter,
	      exp(gsl_vector_get (s->x, 0)),
	      exp(gsl_vector_get (s->x, s->x->size-1)),
	      s->fval, size);
    }

    if (status == GSL_SUCCESS){
      if(0){
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
      new_tree.mu = rtree.mu;
    }else{
      new_tree.mu = exp(gsl_vector_get(s->x, npar_ne));
    }

    minL = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

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
    evo_tree new_tree(rtree.nleaf, enew, exp( gsl_vector_get(s->x, count)), rtree.tobs);
    if(maxj==0){
      new_tree.mu = rtree.mu;
    }else{
      new_tree.mu = exp(gsl_vector_get(s->x, npar_ne));
    }

    minL = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    return new_tree;
  }

}
