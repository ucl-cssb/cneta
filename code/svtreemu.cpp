// run maximum likelihood inference


#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <random>
#include <sstream>
#include <ctime>
#include <map>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>

#include <boost/algorithm/string.hpp>

#include "gzstream.h"

#include "evo_tree.hpp"
#include "genome.hpp"
#include "stats.hpp"

using namespace std;

gsl_rng * r;

int debug = 0;
double LARGE_LNL = -1e9;

// global values for gsl function minimization
vector<vector<int> > vobs;
vector<double> tobs;
int Ns;
int Nchar;

// global value for tree search
map<string,int> searched_trees;
bool exhausted_tree_search = false;


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

vector<double> read_time_info(const string& filename, const int& Ns){
  if(debug) cout << "\tread_time_info" << endl;
  vector<double> t_info;
  ifstream infile (filename.c_str());
  if (infile.is_open()){
    std::string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf; 
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

      int dt = atof(split[1].c_str());
      //cout << "read dt: " << dt << endl;
      t_info.push_back(dt);
    }
  }else{
    std::cerr << "Error: open of time data unsuccessful: " <<  filename << std::endl;
    exit(1);
  }

  if( t_info.size() != Ns ){
    std::cerr << "Error: timing information does not contain Ns entries: " <<  filename << std::endl;
    exit(1);
  }

  return t_info;
}

evo_tree read_tree_info(const string& filename, const int& Ns){
  if(debug) cout << "\tread_tree_info" << endl;

  vector<edge> edges;
  int id = 0;
  
  ifstream infile (filename.c_str());
  if (infile.is_open()){
    std::string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf; 
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

      if( id > 0 ){
	int start = atoi(split[0].c_str());
	int end   = atoi(split[1].c_str());
	double length = atof(split[2].c_str());
	//cout << "t: " << id << "\t" << start << "\t" << end << "\t" << length << endl;
	edges.push_back( edge( id-1, start-1, end-1, length) );
      }
      id++;
    }
  }else{
    std::cerr << "Error: open of tree data unsuccessful: " <<  filename << std::endl;
    exit(1);
  }

  evo_tree new_tree(Ns+1, edges);
  //new_tree.print();
  
  return new_tree;
}

vector<vector<int> > read_data_var_regions(const string& filename, const int& Ns, const int& max_cn){
  if(debug) cout << "\tread_data_var_regions" << endl;
  vector<vector<vector<int> > > s_info;

  // data indexed by [sample][data][ chr, bid, cn ]
  for(int i=0; i<Ns; ++i) s_info.push_back( vector< vector<int> >(4401, vector<int>(3,0) ) );
  
  //ifstream infile (filename.c_str());
  igzstream infile (filename.c_str());

  int counter = 0; 
  //if (infile.is_open()){
    std::string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;
      
      std::vector<std::string> split;
      std::string buf; 
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

      int sample = atoi(split[0].c_str());
      if( sample > Ns) break;

      //cout << sample-1 << "\t" << counter << endl;
      s_info[sample-1][counter][0] = atoi( split[1].c_str() );
      s_info[sample-1][counter][1] = atof( split[2].c_str() );
      s_info[sample-1][counter][2] = atof( split[3].c_str() );
      counter++;

      if(counter >= 4401) counter = 0;
    }
      
    //}else{
    // std::cerr << "Error: open of data file unsuccessful: " <<  filename << std::endl;
    //exit(1);
    //}
  cout << "\tSuccessfully read input file" << endl;
  
  // Loop over and output only the regions that have varied
  vector<int> var_bins(4401,0);
  for(int k=0; k<4401; ++k){ 
    int sum = 0;
    for(int i=0; i<Ns; ++i){
      sum += abs(s_info[i][k][2]);
    }
    if(sum != 2*Ns) var_bins[k] = 1;
  }

  if(debug){
    cout << "\tVariable bins found:" << endl;
    for(int k=0; k<4401; ++k){
      if(var_bins[k] == 1){
	cout << s_info[0][k][0] << "\t" << s_info[0][k][1];
	for(int i=0; i<Ns; ++i) cout << "\t" << s_info[i][k][2];
	cout << endl;
      }
    }
  }

  int nvar = accumulate(var_bins.begin(), var_bins.end(), 0);
  cout << "\tFound variable bins:\t" << nvar << endl;

  // We now need to convert runs of variable bins into segments of constant cn values
  
  vector<vector<int> > segs;
  for(int k=0; k<4401;){

    if(var_bins[k] == 1){
      //in_seg = true;
      int chr = s_info[0][k][0];
      int seg_start = s_info[0][k][1];
      int id_start = k;

      // hold the total copy number of the first bin. If this changes we need a new segment
      int seg_cn_tot = 0;
      for(int j=0; j<Ns; ++j) seg_cn_tot += s_info[j][k][2];

      //cout << "seg_start: " << chr << "\t" << seg_start << ", cn= " << seg_cn_tot << endl;

      bool const_cn = true;
      while( var_bins[k] == 1 && s_info[0][k][0] == chr && const_cn == true){

	// calculate new total cn of next bin
	int cn_tot = 0;
	for(int j=0; j<Ns; ++j) cn_tot += s_info[j][k][2];
	//cout << "\tbin:\t" << k+1 << "\t" << cn_tot << endl;
	if( cn_tot == seg_cn_tot){
	  const_cn = true;
	  ++k;
	}else{
	  const_cn = false;
	  //cout << "\tsplitting segment" << endl;
	}
      }
      int seg_end = s_info[0][k-1][1];
      int id_end = k-1;
      
      //cout << "seg_end:\t" << seg_end << "\t" << k << endl;
      //cout << endl;

      vector<int> seg;
      seg.push_back(chr);
      seg.push_back(id_start);
      seg.push_back(id_end);
      seg.push_back(seg_start);
      seg.push_back(seg_end); 
      segs.push_back(seg);

      // rewind k by one to get the split segment start correct
      if(const_cn == false) k--;
    }
    ++k;
  }
  cout << "\tFound segments:\t\t" << segs.size() << endl;  

  vector<vector<int> > ret;
  
  for(int i=0; i<segs.size(); ++i){
    vector<double> av_cn(Ns,0);
    bool valid = true;

    for(int j=0; j<Ns; ++j){
      for(int k=segs[i][1]; k<(segs[i][2]+1); ++k){
	av_cn[j] += s_info[j][k][2];
      }
      av_cn[j] = av_cn[j]/( segs[i][2] - segs[i][1] + 1 );

      // check all cns across the segment are integer valued
      if( ceil(av_cn[j]) != floor(av_cn[j]) ) valid = false;
    }

    if(debug){
      cout << "\t" << segs[i][0] << "\t" << segs[i][1] << "\t" << segs[i][2];
      for(int j=0; j<Ns; ++j) cout << "\t" << av_cn[j];
      cout << "\t" << valid << endl;
    }
      
    if( valid == true ){
      vector<int> vals;
      vals.push_back( segs[i][0] ); // chr
      vals.push_back( segs[i][1] ); // start
      vals.push_back( segs[i][2] ); // end
      for(int j=0; j<Ns; ++j){
	int cn = (int) av_cn[j];
	if( cn <= max_cn ) vals.push_back( cn );
	else vals.push_back( max_cn );
      }
      ret.push_back( vals );
    }
    
  }
  cout << "\tUsing segments:\t\t" << ret.size() << endl;
  for(int i=0; i<ret.size(); ++i){
    for(int j=0; j<ret[i].size(); ++j){
    cout << "\t" << ret[i][j];
    }
    cout << endl;
  }
  
  return ret;
  
}

double get_transition_prob(const double& mu, double blength, const int& sk, const int& sj ){
  //if(debug) cout << "\tget_transition_prob" << endl;

  int model = 0;

  if(model == 0){
    // assume a basic Markov model here cf JC
    // we have five states 0, 1, 2, 3, 4
    // Pij(t) = exp(-ut) k_ij + (1 - exp(-ut))*pi_j
    // pi_0 = pi_1 = pi_2 = pi_3 = pi_4 = 1/5

    if( sk == sj ){
      return exp(-mu*blength) + (1 - exp(-mu*blength))*0.2;
    }else{
      return (1 - exp(-mu*blength))*0.2;
    }
  }
  if(model == 1){
    // this model assumes a band structure on the Q matrix
    
    if( sk == sj ){
      if( sk == 0 || sk == 4 ){
	exp(-mu*blength) + (1 - exp(-mu*blength))/2.0;
      }
      else{
	exp(-mu*blength) + (1 - exp(-mu*blength))/3.0;
      }
    }else{
      if( (sk == 0 && sj == 1) || (sk == 1 && sj == 0) || (sk == 4 && sj == 3) || (sk == 3 && sj == 4) ){
	return (1 - exp(-mu*blength))/2.0;
      }
      else if( abs(sk - sj) == 1){
	return (1 - exp(-mu*blength))/3.0;
      }
      else{
	return 0;
      } 
    }  
  }
}


double get_likelihood(const int& Ns, const int& Nchar, const vector<vector<int> >& vobs, evo_tree& rtree){
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
evo_tree max_likelihood(evo_tree& rtree, double& minL, int cons=0, int maxj=0){

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
  //gsl_vector_set_all (ss, 0.01);
  gsl_vector_set_all (ss, 1.0);
  
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
    status = gsl_multimin_test_size (size, 1e-4);

    if(0){
      printf ("%5d, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n", 
	      iter,
	      exp(gsl_vector_get (s->x, 0)), 
	      exp(gsl_vector_get (s->x, s->x->size-1)), 
	      s->fval, size);
    }
      
    if (status == GSL_SUCCESS){
      if(0){
	printf ("converged to minimum at\n");
    
	printf ("%5d, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n", 
		iter,
		exp(gsl_vector_get (s->x, 0)), 
		exp(gsl_vector_get (s->x, s->x->size-1)), 
		s->fval, size);
      }
    }

  }
  while (status == GSL_CONTINUE && iter < 2000);

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

int main (int argc, char ** const argv) {
  setup_rng(0);

  // input cn file and temporal information
  string datafile(argv[1]);
  string timefile(argv[2]);
  string treefile(argv[3]);
  
  // number of regions
  Ns = atoi(argv[4]);

  double mu_0 = atoi(argv[5]);

  int cn_max = 4;
  vector<vector<int> > data = read_data_var_regions(datafile, Ns, cn_max);
  Nchar = data.size();

  // tobs already defined globally
  tobs = read_time_info(timefile,Ns);
  
  // read in mle tree
  evo_tree test_tree = read_tree_info(treefile,Ns);
  
  //vector<vector<int> > vobs; // already defined globally
  for(int nc=0; nc<Nchar; ++nc){
    vector<int> obs;
    for(int i=0; i<Ns; ++i){
      obs.push_back( data[nc][i+3] );
    }
    vobs.push_back( obs );
  }
  
  // estimate mutation rate
  test_tree.print();
  test_tree.tobs = tobs;
  test_tree.mu = 1.0;
  
  double Ls = get_likelihood(Ns, Nchar, vobs, test_tree);
  cout << "\nOriginal tree -ve likelihood: " << -Ls << endl;

  cout << "\n\n### Running optimisation: branches constrained, mu free" << endl;

  double minL = -1*LARGE_LNL;
  evo_tree min_tree_mu;
  
  for(int i=0; i<100; ++i){
    double Lf = 0;
    double mu_g = gsl_ran_lognormal(r, mu_0/Nchar, 1);
    cout << "testing mu_g = " << mu_g << endl;
    test_tree.mu = mu_g;
    evo_tree min_tree = max_likelihood(test_tree, Lf, 1, 1);

    if(Lf < minL){
      minL = Lf;
      cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree_mu.mu*Nchar <<endl;
      min_tree_mu = min_tree;
      min_tree_mu.print();
    }
    
  }
    
  stringstream sstm;
  ofstream out_tree;
  sstm << "results-maxL-mu-tree.txt";
  out_tree.open(sstm.str());
  min_tree_mu.write(out_tree);
  out_tree.close();
  sstm.str("");  
}
