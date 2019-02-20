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
      cout << "read dt: " << dt << endl;
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

string create_tree_string( evo_tree tree ){
  stringstream sstm;
  for(int i=0; i<tree.ntotn; ++i){
    sstm << tree.nodes[i].id+1;
    if(tree.nodes[i].daughters.size() == 2){
      sstm << ";" << tree.nodes[i].daughters[0]+1 << ";" << tree.nodes[i].daughters[1]+1;
    }
    sstm << ":";
  }  
  return sstm.str();
}

string order_tree_string( string tree ){
  stringstream sstm;

  vector<string> split1;
  boost::split(split1, tree, [](char c){return c == ':';});

  for(int i=0; i<split1.size()-1; ++ i){     // split creates an empty string at the end
    //sstm << split1[i];
    //cout << "\t" << split1[i] << endl;
    
    vector<string> split2;
    boost::split(split2, split1[i], [](char c){return c == ';';});

    if( split2.size() == 1){
      sstm << split1[i];
    }
    else{ 
      sstm << split2[0] << ";"; //  << split2[1] << ";" << split2[2];
      if( atoi(split2[1].c_str() ) < atoi(split2[2].c_str() ) ){
	sstm << split2[1] << ";" << split2[2];
      }else{
	sstm << split2[2] << ";" << split2[1];
      }
    }
    sstm << ":";
  }
  return sstm.str();
}


evo_tree perturb_tree( const int& Ns, const int& Nchar, const evo_tree& tree ){
  //if(debug) cout << "\tperturb_tree" << endl;
  
  // swap two leaf nodes: 0 to Ns-1 are the samples, Ns is germline, Ns+1 is root
  vector<int> n_0;
  for(int i=0; i<Ns; ++i) n_0.push_back( i );
  random_shuffle(n_0.begin(), n_0.end(),fp);
  int n1 = n_0[Ns-1];
  n_0.pop_back();
  random_shuffle(n_0.begin(), n_0.end(),fp);
  int n2 = n_0[0];

  //cout << "\ttree: ";
  //cout << create_tree_string( tree ) << endl;
  //cout << "\tswapping nodes:\t" << n1+1 << "\t" << n2+1 << endl;

  vector<edge> enew;
  for(int i=0; i<tree.nedge; ++i){
    enew.push_back( tree.edges[i] );
    bool changed = false;
    if(enew[i].end == n1 && changed == false){
      enew[i].end = n2;
      changed = true;
    }
    if(enew[i].end == n2 && changed == false){
      enew[i].end = n1;
      changed = true;
    }
  }

  evo_tree new_tree(Ns+1, enew);
  
  //cout << "\tntree: ";
  //cout << create_tree_string( new_tree ) << endl;

  //string ordered = order_tree_string( create_tree_string( new_tree ) );
  //cout << "\totree: ";
  //cout << ordered << endl;
  //cout << endl;
  
  return new_tree;
}

evo_tree perturb_tree( const int& Ns, const int& Nchar, vector<evo_tree> trees ){
  //if(debug) cout << "\tperturb_tree" << endl;

  int count = 0;
  while(true){
    // randomly sample the fit population
    int ind = gsl_rng_uniform_int(r, trees.size());

    // generate a new tree
    evo_tree ttree = perturb_tree( Ns, Nchar, trees[ind] );

    string tstring = order_tree_string( create_tree_string( ttree ) );

    if ( searched_trees.find(tstring) == searched_trees.end() ) {
      // add the tree
      searched_trees[ tstring ] = 1;

      return ttree;
      
    }
    else {
      //cout << "Tree already present" << endl;
      count++;
    }

    if(count > 100){
      //cout << "\tperturb_tree cannot find new topologies" << endl;
      exhausted_tree_search = true;
      return ttree;
    }
  }
  
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


double get_likelihood(const int& Ns, const int& Nchar, const vector<vector<int> >& vobs, evo_tree& rtree, const double& mu){
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
	  Li += get_transition_prob(mu, bli, sk, si ) * L_sk_k[ni][si];
	  //cout << "\tscoring: Li\t" << li << "\t" << get_transition_prob(mu, bli, sk, si ) << "\t" << L_sk_k[ni][si] << endl;
	}
	double Lj = 0;
	// loop over possible sj
	for(int sj=0; sj<nstate; ++sj){
	  Lj += get_transition_prob(mu, blj, sk, sj ) * L_sk_k[nj][sj];
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

/*
double prior(const evo_tree& tree){
  for(int i=0; i<tree.nedge; ++i){
    if( tree.edges[i].length < 0.0){
      return 0.0;
    }
  }
  return 1.0;
}
*/

double my_f (const gsl_vector *v, void *params){
  
  evo_tree *tree = (evo_tree*) params;

  double mu = 1.0/Nchar; // theta = 2*lambda

  //create a new tree with the current branch lengths
  vector<edge> enew;
  for(int i=0; i<(tree->nedge); ++i){
    enew.push_back( tree->edges[i] );
  }
  for(int i=0; i<(tree->nedge)-1; ++i){
    enew[i].length = exp( gsl_vector_get(v, i) );
  }
  evo_tree new_tree(Ns+1, enew);
    
  //return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, mu)*prior(new_tree);
  return -1.0*get_likelihood(Ns, Nchar, vobs, new_tree, mu);
}
  
  


// given a tree, maximise the mu and branch lengths assuming branch lengths are independent
evo_tree max_likelihood(evo_tree& rtree, double& minL){

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;

  int npar_ne = rtree.nedge-1;
  int npar = npar_ne;

  // initialise the best guess branch length
  gsl_vector *x = gsl_vector_alloc (npar);
  for(int i=0; i<npar_ne; ++i){
    gsl_vector_set (x, i, log(rtree.edges[i].length));
  }
  //gsl_vector_set (x, i, log(1.0));
  
  evo_tree *p = &rtree;
  void * pv = p;  
  //my_f (x, pv);
  
  // Set initial step sizes to 1
  gsl_vector* ss = gsl_vector_alloc (npar);
  gsl_vector_set_all (ss, 0.01);

  // Initialize method and iterate
  minex_func.n = npar;
  minex_func.f = my_f;
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
    status = gsl_multimin_test_size (size, 1e-3);

    if (status == GSL_SUCCESS){
      if(0){
	printf ("converged to minimum at\n");
    
	printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
		iter,
		gsl_vector_get (s->x, 0), 
		gsl_vector_get (s->x, 1), 
		s->fval, size);
      }
    }

  }
  while (status == GSL_CONTINUE && iter < 1000);

  // create a new tree
  vector<edge> enew;
  for(int i=0; i<rtree.nedge; ++i){
    enew.push_back( rtree.edges[i] );
  }
  for(int i=0; i<rtree.nedge-1; ++i){
    enew[i].length = exp(gsl_vector_get(s->x, i));
  }
  evo_tree new_tree(Ns+1, enew);

  minL = s->fval;
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return new_tree;
}

// map constrained and free branch lengths
vector<int> create_mapping( const evo_tree& rtree, const int& s0 ){
  vector<int> ret(rtree.nedge-1, 0);
  int count = 1;
  for(int i=0; i<rtree.nedge; ++i){
    if( rtree.edges[i].end < Ns ){
      if( rtree.edges[i].end == s0 ){
	ret[i] = 0;
      }
      else{
	ret[i] = -1;
      }
    }
    else{
      ret[i] = count;
      ++count;
    }
  }

  return ret;
}

/*
vector<double> extract_free_edges( const evo_tree& rtree, const int& s0 ){
  
  vector<double> ret;
  
  for(int i=0; i<rtree.nedge-1; ++i){
    if( rtree.edges[i].end == s0 ){
      ret.push_back(rtree.edges[i].length);
    }
  }
  // next free parameters are the branch lengths of the internal edges
  // internal edges are any edge that doesn't end on a leaf node
  for(int i=0; i<rtree.nedge-1; ++i){
    if( rtree.edges[i].end >= Ns+1 ){
      ret.push_back(rtree.edges[i].length);
    }
  }
  return ret;
}

vector<double> reconstruct_edges( const vector<double>& reduced_edges, const int& s0, const vector<double> t_info ){
  
  vector<double> ret;
  for(int i=0; i<rtree.nedge-1; ++i){
    if( rtree.edges[i].end < Ns ){
      ret.push_back( reduced_edges[0] + t_info[ rtree.edges[i].end - 1] );
    }
  }
  for(int i=0; i<rtree.nedge-1; ++i){
    if( rtree.edges[i].end >= Ns+1 ){
      ret.push_back(rtree.edges[i].length);
    }
  

  ret.push_back( reduced_edges[0]
    
  }
  
  return ret;
}
*/  


// given a tree, maximise the mu and branch lengths assuming branch lengths are related through tobs
evo_tree max_likelihood_cons(evo_tree& rtree, double& minL, const int& s0, const vector<double>& tobs){

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;

  // now there is essentially only one parameter for the leaf nodes corresponding to the length of t0
  int npar_ne = rtree.nedge - 1 - (Ns - 1);
  int npar = npar_ne;

  // initialise the best guess branch length
  vector<int> mapping = create_mapping( rtree, s0 );
  
  gsl_vector *x = gsl_vector_alloc (npar);
  for(int i=0; i<npar_ne; ++i){
    gsl_vector_set (x, i, log(rtree.edges[i].length) );
  }
 
  // next free parameter mutation rate
  //gsl_vector_set (x, i, log(1.0));
  
  evo_tree *p = &rtree;
  void * pv = p;  
  //my_f (x, pv);
  
  // Set initial step sizes to 1
  gsl_vector* ss = gsl_vector_alloc (npar);
  gsl_vector_set_all (ss, 0.01);

  // Initialize method and iterate
  minex_func.n = npar;
  minex_func.f = my_f;
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
    status = gsl_multimin_test_size (size, 1e-3);

    if (status == GSL_SUCCESS){
      if(0){
	printf ("converged to minimum at\n");
    
	printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
		iter,
		gsl_vector_get (s->x, 0), 
		gsl_vector_get (s->x, 1), 
		s->fval, size);
      }
    }

  }
  while (status == GSL_CONTINUE && iter < 1000);

  // create a new tree
  vector<edge> enew;
  for(int i=0; i<rtree.nedge; ++i){
    enew.push_back( rtree.edges[i] );
  }
  for(int i=0; i<rtree.nedge-1; ++i){
    enew[i].length = exp(gsl_vector_get(s->x, i));
  }
  evo_tree new_tree(Ns+1, enew);

  minL = s->fval;
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return new_tree;
}


evo_tree do_evolutionary_algorithm(const int& Npop, const int& Ngen, const int& max_static){

  cout << "Running evolutionary algorithm" << endl;
  
  // create initial population of trees. Sample from coalescent trees
  vector<evo_tree> trees;
  
  vector<double> lnLs(2*Npop,0);
  for(int i=0; i<Npop; ++i){
    evo_tree rtree = generate_coal_tree(Ns);
    trees.push_back( rtree );

    string tstring = order_tree_string( create_tree_string( rtree ) );
    searched_trees[ tstring ] = 1;
  }
  double min_lnL = 1e20;
  evo_tree min_lnL_tree;
  double Lf = 0;
  int count_static = 0;
  
  for(int g=0; g<Ngen; ++g){
    // Growth stage: create Npop copies + Npop copies with mutation
    // The top scoring trees have already been scored and optimised
    vector<evo_tree> new_trees;
    vector<evo_tree> opt_trees;

    if( g == 0 ){
      for(int i=0; i<Npop; ++i){
	new_trees.push_back( trees[i] );
      }
    
      for(int i=0; i<Npop; ++i){
	//new_trees.push_back( perturb_tree(Ns, Nchar, trees[i]) );
	new_trees.push_back( perturb_tree(Ns, Nchar, trees) );
      }
    
      // Selection: score the trees
      for(int i=0; i<2*Npop; ++i){
	evo_tree otree = max_likelihood( new_trees[i], Lf );
	otree.score = Lf;
	lnLs[i] = Lf;
	opt_trees.push_back( otree );
	cout << "g/i/lnL:\t" << g << "\t" << i << "\t" << lnLs[i] << endl;
      }
    }else{

      // Leave this subpopulation unchanged
      for(int i=0; i<Npop; ++i){
	new_trees.push_back( trees[i] );
	opt_trees.push_back( trees[i] );
	lnLs[i] = new_trees[i].score;
      }

      // Perturb this subpopulation
      for(int i=0; i<Npop; ++i){
	//new_trees.push_back( perturb_tree(Ns, Nchar, trees[i]) );       // trees is of size Npop

	new_trees.push_back( perturb_tree(Ns, Nchar, trees) );
		
	evo_tree otree = max_likelihood( new_trees[Npop + i], Lf );     // new_trees of size 2 Npop
	otree.score = Lf;
	lnLs[Npop + i] = Lf;
	opt_trees.push_back( otree );
	cout << "g/i/lnL:\t" << g << "\t" << Npop+i << "\t" << lnLs[i] << endl;
      }

      // Randomly sample again?

    }
    
    vector<int> index(2*Npop);
    int x=0;
    iota( index.begin(), index.end(), x++);
    sort( index.begin(), index.end(), [&](int i,int j){return lnLs[i]<lnLs[j];} );

    // Selection: calculate mean fitness of top half of population
    double meand = 0;
    for(int k=0; k<Npop; ++k){
      //cout << "\t\t" << lnLs[ index[k] ] << endl;
      meand += lnLs[ index[k] ];
    }
    meand = meand/Npop;
    if( g%1 == 0) cout << "g / av dist / top dist / trees searched \t" << g << "\t" << meand << "\t" << lnLs[ index[0] ] << "\t" << searched_trees.size() << endl;

    // Selection: select top half
    for(int i=0; i<Npop; ++i){
      trees[i] = opt_trees[ index[i] ];
    }

    // Selection: record the best (lowest) scoring tree
    if( lnLs[ index[0] ] < min_lnL ){
      min_lnL = lnLs[ index[0] ];
      min_lnL_tree = opt_trees[ index[0] ];
      count_static = 0;

      min_lnL_tree.print();
      
    }else{
      cout << "min static" << endl;
      count_static += 1;
    }

    if( max_static > 0 && count_static == max_static ){
      cout << "\t### static likelihood function. Finishing on ngen = " << g << endl;
      break;
    }
    
    fill(lnLs.begin(), lnLs.end(), 0);
    
    if( exhausted_tree_search == true ){
      cout << "\tperturb_tree struggling to find new topologies. Either exhausted possible trees or local minimum" << endl;
      break;
    }
  }

  cout << "FINISHED. MIN -ve logL = " << min_lnL << endl;

  return min_lnL_tree;
}

int main (int argc, char ** const argv) {
  setup_rng(0);

  // input cn file and temporal information
  string datafile(argv[1]);
  string timefile(argv[2]);
  
  // number of regions
  Ns = atoi(argv[3]);
  
  // control parameters
  //int Npop = 10;
  //int Ngen = 5000;
  //int max_static = 10;

  int Npop = atoi(argv[4]);
  int Ngen = atoi(argv[5]);
  int max_static = atoi(argv[6]);

  int cn_max = 4;
  vector<vector<int> > data = read_data_var_regions(datafile, Ns, cn_max);
  Nchar = data.size();

  vector<double> del_t = read_time_info(timefile,Ns);
  int s0 = 0;
  for(int i=0; i<del_t.size();++i){
    if(del_t[i] == 0){
      s0 = i;
      break;
    }
  }
  cout << "assigned ref branch to leaf: " << s0 << endl;

  //exit(0);
  
  //vector<vector<int> > vobs; // already defined globally
  for(int nc=0; nc<Nchar; ++nc){
    vector<int> obs;
    for(int i=0; i<Ns; ++i){
      obs.push_back( data[nc][i+3] );
    }
    vobs.push_back( obs );
  }

  if(0){
    static const int arr1[] = {8,5, 8,1, 9,2, 9,3, 10,9, 10,8, 11,4, 11,10, 7,11, 7,6 };
    vector<int> e (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
    for(int i=0; i<e.size();++i) e[i] = e[i] - 1;
  
    static const double arr2[] = {18.49, 38.49, 51.71, 31.71, 0.51, 3.73, 22.2, 0.013, 0.99, 0};
    vector<double> l (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
    //for(int i=0; i<l.size();++i) l[i] = l[i]/10.0;
  
    evo_tree test_tree(Ns+1, e, l);
    test_tree.print();

    double mu = 1.0/Nchar;
    double Ls = 0.0;

    Ls = get_likelihood(Ns, Nchar, vobs, test_tree, mu);
    cout << "\nOriginal tree -ve likelihood: " << -Ls << endl;

    //for(int i=0; i<10; ++i){
    //double mu_p = mu + (i-5)*0.01*mu;
    //Ls = get_likelihood(Ns, Nchar, vobs, test_tree, mu_p );
    //cout << "\n -ve likelihood: " << mu_p << "\t" << -Ls << endl;
    //}

    test_tree.print();

    double Lf = 0;
    cout << "\n\nRunning optimisation free" << endl;
    evo_tree min_tree = max_likelihood(test_tree, Lf);
    cout << "\nMinimised tree likelihood: " << Lf << endl;
    min_tree.print();

    vector<int> mapping = create_mapping( test_tree, s0 );
    cout << "##### created mapping #####" << endl;
    for(int i=0; i< mapping.size(); ++i){
      cout << "\t" << i+1 << "\t" << mapping[i] << endl;
    }
    
    //cout << "\n\nRunning optimisation constrained" << endl;
    //evo_tree min_tree2 = max_likelihood_cons(test_tree, Lf, s0, del_t);
    //cout << "\nMinimised tree likelihood: " << Lf << endl;
    //min_tree2.print();
    
    
  }
  
  
  if(1){
    evo_tree min_lnL_tree = do_evolutionary_algorithm(Npop, Ngen, max_static);
  
    // Write out the top tree
    stringstream sstm;
    //cout << "Best fitting tree, -ve lnL = " << global_min << endl;
    min_lnL_tree.print();
    sstm << "sim-data-" << "1" << "-tree.txt";
    ofstream out_tree(sstm.str());
    min_lnL_tree.write(out_tree);
    out_tree.close();
    sstm.str("");
  }
  
  //evo_tree min_lnL_tree = do_greedy_algorithm();
}
