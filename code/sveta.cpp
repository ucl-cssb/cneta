// to build the program:
//   g++ sveta.cpp -o sveta -L/usr/local/lib/ -lgsl -I/usr/local/include
//
// to run the model:
//   ./sveta


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

#include "gzstream.h"

#include "sveta.hpp"

#include "evo_tree.hpp"
#include "genome.hpp"
#include "stats.hpp"

using namespace std;

vector<mutation> generate_mutation_times(const int& edge_id, const double& blength, const double& node_time, const vector<double>& rate_constants ){
  bool print = false;

  if(print) cout << "\tgenerate_mutations, blength:" << "\t" << blength << endl;

  vector<mutation> ret;
  
  // model this as a poisson process
  double time = 0.0;
  while( time < blength){
    vector<double> rates;
    for(int i=0; i<rate_constants.size(); ++i){
      rates.push_back(rate_constants[i]);
      //cout << "##:\t" << rates[i] << endl;
    }

    double rate = accumulate(rates.begin(), rates.end(), 0.0);
    double tevent = gsl_ran_exponential (r, 1/rate);
    time += tevent;

    // choose type of event
    int e = rchoose(r,rates);
    if(print) cout << "mut times, tevent, total time, time/branch len, event\t" << tevent << "\t" << time << "\t" << blength << "\t" << e << endl;

    ret.push_back( mutation( edge_id, e, time/blength, node_time+time ) );
  }
  ret.pop_back();
  
  return ret;
}

void apply_mutations( const vector<mutation>& muts, genome& g, bool print = false ){
  //print = true;
  
  // mean variant size in bins
  double mdup = g.mean_dup_size;
  double mdel = g.mean_del_size;
  
  for(int i=0; i<muts.size(); ++i){

    g.nmuts[ muts[i].type ]++;
    g.mutations.push_back( muts[i] );
    
    if(muts[i].type == 0){   //duplications

      if(print == true){
	cout << "GENOME BEFORE duplication " << endl;
	//g.print();
	//g.print_cn();
      }
	
      // make sure a duplication somewhere is possible
      bool possible = false;
      for(int j=0; j < g.chrs.size(); ++j){
	if( g.chrs[j].size() > 1 ) possible = true; 
      }

      if( possible == true ){
	bool done = false;
	while( done == false ){
	  int c = gsl_rng_uniform_int(r, g.chrs.size() );
	  int len = gsl_ran_exponential(r, mdup);
	  //int len =  1 + gsl_ran_poisson (r, ldup);
	  //cout << "dup len:" << len << endl;
	  
	  if( len < g.chrs[c].size() ){
	    // choose a random location up to size-len
	    int loc = gsl_rng_uniform_int(r,g.chrs[c].size()-len);
      
	    if(print == true) cout << "\tSV: duplicating segment, chr, start, len: " << c << "\t" << loc << "\t" << len+1 << endl;
	    //cout << "SV: insert: " << loc+len+1 << "\t" << loc << "\t" << loc+len+1 << endl;
	    g.chrs[c].insert(g.chrs[c].begin()+loc+len+1, g.chrs[c].begin()+loc, g.chrs[c].begin()+loc+len+1);
	    done = true;
	  }
	}
      }else{
	if(print == true) cout << "\tSV: duplication failed:" << endl;
      }

      if(print == true){
	cout << "GENOME AFTER duplication " << endl;
	//g.print();
	//g.print_cn();
      }
    }
    else if( muts[i].type == 1){   //deletions

      if(print == true){
	cout << "GENOME BEFORE deletion " << endl;
	//g.print();
	//g.print_cn();
      }
	
      // make sure a deletion somewhere is possible
      bool possible = false;
      for(int j=0; j < g.chrs.size(); ++j){
	if( g.chrs[j].size() > 1 ) possible = true; 
      }

      if( possible == true ){
	bool done = false;
	while( done == false ){
	  int c = gsl_rng_uniform_int(r, g.chrs.size() );
	  int len = (int) gsl_ran_exponential(r, mdel);
	  //int len =  1 + gsl_ran_poisson (r, ldel);
	  //cout << "del len:" << len << endl;
	
	  if( len < g.chrs[c].size() ){
	    int loc = gsl_rng_uniform_int(r,g.chrs[c].size()-len);
	    if(print == true) cout << "\tSV: deleting segment, chrs, start, len: " << c << "\t" << loc << "\t" << len+1 << endl;
	    g.chrs[c].erase(g.chrs[c].begin()+loc, g.chrs[c].begin()+loc+len+1);
	    done = true;
	  }
	}
      }else{
	if(print == true) cout << "\tSV: deletion failed:" << endl;
      }

      if(print == true){
	cout << "GENOME AFTER deletion " << endl;
	//g.print();
	//g.print_cn();
      }
    }
    else if( muts[i].type == 2){   //chr loss
      if( g.chrs.size() > 0 ){
	int c = gsl_rng_uniform_int(r, g.chrs.size() );
	g.chrs.erase(g.chrs.begin()+c);
      }
	
    }
    else if( muts[i].type == 3){   //chr gain
      if( g.chrs.size() > 0 ){
	int c = gsl_rng_uniform_int(r, g.chrs.size() );
	g.chrs.insert(g.chrs.end(), g.chrs.begin()+c, g.chrs.begin()+c+1);
      }	
    }
    else if( muts[i].type == 4){   // whole genome duplication
      if( g.chrs.size() > 0 ){
	g.chrs.insert(g.chrs.end(), g.chrs.begin(), g.chrs.end());
      }	
    }
    else{
      cerr << "Unknown mutation type" << endl;
      exit(1);
    }
  }
}

void traverse_tree_mutating(const int& node_id, const evo_tree& tree, const vector<double>& rate_constants,
			    map<int, vector<mutation> >& all_muts, vector<genome>& genomes ){
  //cout << "\ttraverse_tree: " << node_id+1 << endl;

  if( !tree.nodes[node_id].isRoot ){
    // copy the parent
    genomes[ node_id ] = genomes[ tree.nodes[node_id].parent ];
    genomes[ node_id ].node_id = node_id;

    // apply the mutations from parent -> daughter
    //cout << "MUTATING genome: " << tree.nodes[node_id].parent+1 << " -> " << node_id+1 << "\t edge id: " << tree.nodes[node_id].e_in+1 << endl;
    int edge_id = tree.nodes[node_id].e_in;
    vector<mutation> muts;
    if( tree.edges[edge_id].length > 0 ){
      muts = generate_mutation_times(edge_id, tree.edges[edge_id].length, tree.node_times[ tree.nodes[node_id].parent ], rate_constants);
      apply_mutations( muts, genomes[ node_id ] );
    }
    all_muts[edge_id] = muts;
  }
  
  if( tree.nodes[ node_id ].daughters.size() == 0 ){
    // we are done
    return;      
  }else{
    traverse_tree_mutating(tree.nodes[ node_id ].daughters[0], tree, rate_constants, all_muts, genomes );
    traverse_tree_mutating(tree.nodes[ node_id ].daughters[1], tree, rate_constants, all_muts, genomes );
  }
  return;
}

void simulate_samples(vector<genome>& genomes, const evo_tree& tree, genome& germline, const vector<double>& rate_constants, bool print = true ){
  
  // assign the germline to the root of the tree
  germline.node_id = tree.root_node_id;
  
  for(int i=0; i<(tree.nnode+tree.nleaf); ++i){
    if(i == tree.root_node_id){
      genomes.push_back( germline );
    }
    else{
      genomes.push_back( genome() );
    }
  }

  // hold all mutations by edge id
  map<int, vector<mutation> > muts;
  
  // move through the evolutionary tree mutating genomes
  traverse_tree_mutating( tree.root_node_id, tree, rate_constants, muts, genomes );
  
  // final samples returned, print out leaf nodes
  if(print == true){
    cout << "MUTATIONS:" << endl;
    for(map<int, vector<mutation> >::iterator it = muts.begin(); it != muts.end(); it++){
      cout << "EDGE, id: " << it->first+1
	   << "\t" << tree.edges[it->first].start+1 << " -> " << tree.edges[it->first].end+1
	   << "\t" << it->second.size() << endl;
      vector<mutation> v = it->second;
      for(int i=0; i<v.size(); ++i){
	v[i].print();
      }
    }
    cout << endl;
    cout << "LEAF GENOMES:" << endl;
    for(int i=0; i<tree.nleaf; ++i){
      //genomes[i].print();
      genomes[i].print_muts();
      //genomes[i].print_cn();
      tree.print_ancestral_edges( genomes[i].node_id );
      cout << endl;
    }
  }
}


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

/*
// args:
// Ns: number of samples (not including germline)
// prc: 5 rate constants for duplication, deletion, chromosome gain, chromosome loss, wgd
// pvs: mean dup size, mean del size
// ptree: parameters for the tree toplogy
// ret: 2d array of size nbins x Ns
evo_tree construct_tree( const int& Ns, const vector<double>& epars, const vector<double>& lengths, vector<double>& node_times ){

  vector<int> edges;
  edges.push_back( 6 );
  edges.push_back( 5 );
  edges.push_back( 6 );
  edges.push_back( 9 );
  edges.push_back( 7 );
  edges.push_back( epars[0] );
  edges.push_back( 7 );
  edges.push_back( epars[1] );
  edges.push_back( 7 );
  edges.push_back( epars[2] );
  edges.push_back( 8 );
  edges.push_back( epars[3] );
  edges.push_back( 8 );
  edges.push_back( epars[4] );
  edges.push_back( 9 );
  edges.push_back( epars[5] );
  edges.push_back( 9 );
  edges.push_back( epars[6] );

  node_times[5] = 0;
  node_times[4] = 0;
  node_times[8] = lengths[1]; // root -> 9
  node_times[7] = lengths[1];
  
  evo_tree test_tree(Ns+1, edges, lengths);
  return test_tree;
}
*/
//void run_sample_set(int Ns, double* prc, double* pvs, double* ptree, double* pl, int* ret){
void run_sample_set(int Ns, double* prc, double* pvs, int* ret){

  static const int arr[] = {367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};
  vector<int> chr_lengths (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    
  vector<double> rate_consts (prc, prc + sizeof(prc) / sizeof(prc[0]) );    
  genome germline(chr_lengths,2);
  germline.mean_dup_size = pvs[0];
  germline.mean_del_size = pvs[1];
    
  vector<int> edges;
  
  vector<double> lengths;
  vector<double> epoch_times;
  vector<double> node_times;
  stringstream sstm;
  vector<genome> results;  
  //cout << "\n\n###### New sample collection ######" << endl;
  //cout << "###### Ns+1= " << Ns+1 << endl;
      
  generate_coal_tree(Ns, edges, lengths, epoch_times, node_times);
  evo_tree test_tree(Ns+1, edges, lengths);

  //for(int i=0; i<6; ++i) epars.push_back( ptree[i] );
  //for(int i=0; i<8; ++i) lengths.push_back( pl[i] );
  //evo_tree test_tree = construct_tree(Ns, epars, lengths, node_times );
  
  test_tree.node_times = node_times;
  simulate_samples(results, test_tree, germline, rate_consts, false);

  int nbins = 4401;
  for(int i=0; i<(test_tree.nleaf-1); ++i){
    vector<int> cns = results[i].get_cn_vector();
    for(int j=0; j<nbins; ++j){
      ret[nbins*i + j] = cns[j];
    }
  }

  if(0){
    sstm << "test-data-cn.txt";
    ofstream out_cn(sstm.str());
    for(int j=0; j<test_tree.nleaf; ++j){
      results[j].write(out_cn);
    }
    out_cn.close();
    sstm.str("");
  }
  return;
}	    

//////////////////////////////////////////////////////////
///                                                    ///
///   MAIN                                             ///
///                                                    ///
//////////////////////////////////////////////////////////

// Tree representation
// if we have n samples these trees always have the germline variation as the n+1 samples
// the final edge length is always zero
// edges are always directed from -> to going from the node/MRCA
// therefore leaf nodes are always in the second column of edges and the root is only in the first

// Internally nodes and edges have ids starting from 0
// to match up with ape in R we always print id+1


/*
//order of parameters expected:  Npop, ngen, g_d, mu_i, p_tran, svp2, trp2, SV_min, svgradient, svtransgrad, fitness, SV_mean, mu_k, gnmdou, maxchr, minchr
void read_params( const string& filename, vector<double>& prior_rates_l, vector<double>& prior_rates_u){

  ifstream infile (filename.c_str());

  int counter = 0; 
  if (infile.is_open()){

    std::string line;
    
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf; 
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

      prior_rates_l.push_back( atof( split[0].c_str() ) );
      prior_rates_u.push_back( atof( split[1].c_str() ) );	      
    
      counter++;
    }
      
  }else{
    std::cerr << "Error: open of input parameter file unsuccessful: " <<  filename << std::endl;
  }
  
  return counter;
}
*/

int main (int argc, char ** const argv) {
  setup_rng(0);

  int mode = 0;

  // output directory
  string dir(argv[1]);
  dir = dir + "/";

  // number of regions
  int Ns = atoi(argv[2]);
  
  // number of multi-region samples
  int Nsims = atoi(argv[3]);

  // priors on rates
  //vector<double> pars;
  //read_params(argv[4],pars);

  // four event types: duplication, deletion, chromosome gain, chromosome loss, wgd
  // rates are 1/mean
  vector<double> rate_consts; 
  rate_consts.push_back( atof(argv[4]) );
  rate_consts.push_back( atof(argv[5]) );
  rate_consts.push_back( atof(argv[6]) );
  rate_consts.push_back( atof(argv[7]) );
  rate_consts.push_back( atof(argv[8]) );

  // parameters for mean of dup/del size distributions
  vector<double> var_size; 
  var_size.push_back( atof(argv[9]) );
  var_size.push_back( atof(argv[10]) );

  cout << "rates:\t" << rate_consts[0] << "\t" << rate_consts[1]  << "\t" << rate_consts[2]  << "\t" << rate_consts[3]  << "\t" << rate_consts[4] << endl;
  cout << "sizes:\t" << var_size[0] << "\t" << var_size[1] << endl;
  
  // simulate and coalescent tree and apply SVs
  if(mode == 0){
    
    static const int arr[] = {367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};
    vector<int> chr_lengths (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    
    // can specify the germline in a number of different ways
    //genome germline(1,10);
    //genome germline(42,1000);
    //genome germline(chr_lengths);
    
    genome germline(chr_lengths,2);
    germline.mean_dup_size = var_size[0];
    germline.mean_del_size = var_size[1];
    
    vector<int> edges;
    vector<double> lengths;
    vector<double> epoch_times;
    vector<double> node_times;
    stringstream sstm;
    vector<genome> results;
    for(int i=0; i<Nsims; ++i){
      //cout << "\n\n###### New sample collection ######" << endl;
      //int Ns = 2 + gsl_rng_uniform_int(r, 8);
      //int Ns = 4;
      //cout << "###### Ns+1= " << Ns+1 << endl;

      generate_coal_tree(Ns, edges, lengths, epoch_times, node_times);
      //cout << "EDGES:" << endl;
      //for(int j=0; j<edges.size(); j=j+2) cout << edges[j]+1 << " -> " << edges[j+1]+1 << endl;

      evo_tree test_tree(Ns+1, edges, lengths);
      test_tree.node_times = node_times;
      simulate_samples(results, test_tree, germline, rate_consts, false);

      sstm << dir << "sim-data-" << i+1 << "-info.txt";
      ofstream out_info(sstm.str());
      out_info << "NODE TIMES:";
      out_info << "\tnid\ttime" << endl;
      for(int j=0; j< node_times.size(); ++j){
	out_info << "\t" << j+1 << "\t" << node_times[j] << endl;
      }
      out_info << endl;
     
      for(int i=0; i<test_tree.nleaf; ++i){
	test_tree.print_ancestral_edges( results[i].node_id, out_info );
      }
      out_info << endl;
      for(int i=0; i<test_tree.nleaf; ++i){
	results[i].print_muts(out_info);
      }
      out_info.close();
      sstm.str("");
      
      sstm << dir << "sim-data-" << i+1 << "-tree.txt";
      ofstream out_tree(sstm.str());
      test_tree.write(out_tree);
      out_tree.close();
      sstm.str("");
      
      sstm << dir << "sim-data-" << i+1 << "-cn.txt.gz";
      //ofstream out_cn(sstm.str());
      ogzstream out_cn(sstm.str().c_str());
      for(int j=0; j<test_tree.nleaf; ++j){
	results[j].write(out_cn);
      }
      out_cn.close();
      sstm.str("");

      sstm << dir << "sim-data-" << i+1 << "-mut.txt";
      ofstream out_mut(sstm.str());
      for(int j=0; j<test_tree.nleaf; ++j){
	for(int i=0; i<results[j].mutations.size(); ++i){
	  out_mut << j+1 << "\t" << results[j].mutations[i].edge_id+1
		 << "\t" << results[j].mutations[i].type << "\t" << results[j].mutations[i].btime << "\t" << results[j].mutations[i].gtime << endl;
	}
      }
      out_mut.close();
      sstm.str("");
      
      edges.clear();
      lengths.clear();
      results.clear();
      epoch_times.clear();
      node_times.clear();
    }
   
  }

  /*
  //create test tree
  if(mode == 1){
    //static const int arr1[] = {7,8, 6,7, 8,1, 8,2, 7,9, 9,3, 9,4, 6,5 };
    static const int arr1[] = {6,7, 5,6, 7,0, 7,1, 6,8, 8,2, 8,3, 5,4 };

    vector<int> e (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );

    static const double arr2[] = {5.1,6.3,10.2,9.5,5.2,3.2,5.4,0};
    vector<double> l (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  
    evo_tree test_tree(5, e, l);

    test_tree.print();
  }
  
  // test out the genome classes and functionality
  if(mode == 2){

    if(0){
      vector<mutation> dels;
      dels.push_back( mutation( 0, 1, 0.5, 0 ) );
      dels.push_back( mutation( 0, 1, 0.5, 0 ) );
      dels.push_back( mutation( 0, 1, 0.5, 0 ) );
      dels.push_back( mutation( 0, 1, 0.5, 0 ) );

      vector<mutation> dups;
      dups.push_back( mutation( 0, 0, 0.5, 0 ) );
      dups.push_back( mutation( 0, 0, 0.5, 0 ) );
      dups.push_back( mutation( 0, 0, 0.5, 0 ) );
      dups.push_back( mutation( 0, 0, 0.5, 0 ) );

      vector<mutation> wgds;
      wgds.push_back( mutation( 0, 4, 0.5, 0 ) );
    
      genome g1(2,10);
      g1.node_id = 0;
      g1.print();
      g1.print_cn();
      apply_mutations( dups, g1 );
      apply_mutations( dels, g1 );
      apply_mutations( wgds, g1 );
      g1.print();
      g1.print_cn();
    }

    if(1){
      static const int arr3[] = {10,4,3};
      vector<int> chr_lengths (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]) );

      genome g2(chr_lengths);
      g2.node_id = 0;
      g2.print();
      g2.print_cn();

      genome g2d(chr_lengths,2);
      g2d.node_id = 0;
      g2d.print();
      g2d.print_cn();

    }
    
    //genome g3 = g2;
    //g3.print();
    //g3.print_cn();
  }
  if(mode == 3){
    stringstream sstm;
    
    int Ns = 4;
    double prc[] = {0.1, 0.1, 0.1, 0.1, 0.05};
    double pvs[] = {30.0, 30.0};

    for(int i=0; i<10; ++i){
    
      int* ret = new int[Ns*4401];
      run_sample_set(Ns, prc, pvs, &(ret[0]) );

      sstm << dir << "sim-data-" << i+1 << "-cn.txt";
      ofstream out_cn(sstm.str());
      for(int i=0; i<(Ns*4401); ++i){
	out_cn << ret[i] << endl;
      }
      out_cn.close();
      sstm.str("");
      delete[] ret;
    }
    
  }
  */
}
