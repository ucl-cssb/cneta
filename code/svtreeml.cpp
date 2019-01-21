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

#include "gzstream.h"

#include "evo_tree.hpp"
#include "genome.hpp"
#include "stats.hpp"

using namespace std;

gsl_rng * r;

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
   
vector<vector<int> > read_data_var_regions(const string& filename, const int& Ns){

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

  if(0){
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

  vector<vector<int> > segs;
  for(int k=0; k<4401;){

    if(var_bins[k] == 1){
      //in_seg = true;
      int chr = s_info[0][k][0];
      int seg_start = s_info[0][k][1];
      int id_start = k;
      //cout << "seg_start: " << chr << "\t" << seg_start << endl;

      while( var_bins[k] == 1 && s_info[0][k][0] == chr){
	++k;
      }
      int seg_end = s_info[0][k-1][1];
      int id_end = k-1;
      
      //cout << "seg_end: " << s_info[0][k][0] << "\t" << seg_end << endl;
      //cout << endl;
      vector<int> seg;
      seg.push_back(chr);
      seg.push_back(id_start);
      seg.push_back(id_end);
      seg.push_back(seg_start);
      seg.push_back(seg_end); 
      segs.push_back(seg);
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

    //cout << "\t" << segs[i][0] << "\t" << segs[i][1] << "\t" << segs[i][2];
    //for(int j=0; j<Ns; ++j) cout << "\t" << av_cn[j];
    //cout << "\t" << valid << endl;

    if( valid == true ){
      vector<int> vals;
      vals.push_back( segs[i][0] );
      vals.push_back( segs[i][1] );
      vals.push_back( segs[i][2] );
      for(int j=0; j<Ns; ++j) vals.push_back( (int) av_cn[j] );
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

evo_tree perturb_mutation_tree( const int& Ns, const int& Nchar, const evo_tree& tree, const int& max_cn ){
  
  double u = runiform(r, 0, 1);

  if(u < 0.3){
    // swap two leaf nodes: 0 to Ns-1 are the samples, Ns is germline, Ns+1 is root
    vector<int> n_0;
    for(int i=0; i<Ns; ++i) n_0.push_back( i );
    random_shuffle(n_0.begin(), n_0.end(),fp);
    int n1 = n_0[Ns-1];
    n_0.pop_back();
    random_shuffle(n_0.begin(), n_0.end(),fp);
    int n2 = n_0[0];
    
    //cout << "swapping nodes:\t" << n1+1 << "\t" << n2+1 << endl;

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
    new_tree.chars = tree.chars;
    return new_tree;

  }else{
    evo_tree new_tree = tree;

    // change one of the internal node characters
    int nintn = new_tree.nnode - 1;
    int np = gsl_rng_uniform_int(r, nintn);
    int nc = gsl_rng_uniform_int(r, Nchar);

    //cout << "switching chars:\t" << np << "\t" << nc << endl;
    
    new_tree.chars[tree.nleaf+1 + np][nc] =  gsl_rng_uniform_int(r, max_cn);
    
    return new_tree;
  }
  
}

int dist( const int& Nchar, const vector<int>& ai, const vector<int>& af ){
  int ret = 0;
  for(int i=0; i<Nchar; ++i) ret += abs(af[i] - ai[i]);
  return ret;
}

void print_assignment(const int& Nchar, const int& ntotalnodes, const vector<vector<int> >& node_ass){
  for(int i=0; i<ntotalnodes; ++i){
    for(int j=0; j<Nchar; ++j){
      cout << "\t" << node_ass[i][j];
    }
    cout << endl;
  }
}

double get_parsimony_score(const int& Nchar, evo_tree& tree){
  //cout << "\tget_parsimony score" << endl;
  //print_assignment(Nchar, 9, tree.chars);

  // Loop over edges and calculate distances
  int score = 0;
  for(int i=0; i<tree.nedge; ++i){
    int score_e = dist( Nchar, tree.chars[ tree.edges[i].start ], tree.chars[ tree.edges[i].end ] );

    // set branch lengths
    tree.edges[i].length = score_e;
    score += score_e;
    
    //cout << "edge score" << endl;
    //cout << "\t" << tree.edges[i].id+1 << "\t" << tree.edges[i].start+1 << "\t" << tree.edges[i].end+1 << endl;
    //for(int j=0; j<Nchar; ++j) cout << "\t" << node_ass[tree.edges[i].start][j];
    //cout << endl;
    //for(int j=0; j<Nchar; ++j) cout << "\t" << node_ass[tree.edges[i].end][j];
    //cout << endl;
  }
  //cout << "total score: " << score << endl;
  return score;
}


int main (int argc, char ** const argv) {
  setup_rng(0);

  // output directory
  string datafile(argv[1]);

  // number of regions
  int Ns = atoi(argv[2]);
  
  // control parameters
  //int Npop = 1000;
  //int Ngen = 5000;
  //int max_cn = 12;

  int Npop = atoi(argv[3]);
  int Ngen = atoi(argv[4]);
  int max_cn = atoi(argv[5]);
  
  vector<vector<int> > data = read_data_var_regions(datafile, Ns);
  int Nchar = data.size();

  // node_ass is the initial assignments of the characters
  int ntotn = 2*(Ns+1) - 1;
  vector<vector<int> > node_ass(ntotn, vector<int>(Nchar, 0));

  for(int i=0; i<Ns; ++i){
    for(int j=0; j<Nchar; ++j){
      // recall that the first three values are info
      node_ass[i][j] = data[j][i+3];
    }
  }
  for(int i=Ns; i<Ns+2; ++i){
    for(int j=0; j<Nchar; ++j){
      node_ass[i][j] = 2;
    }
  }
  
  for(int i=Ns+2; i<ntotn; ++i){
    for(int j=0; j<Nchar; ++j){
      node_ass[i][j] = 2 + gsl_rng_uniform_int(r, 2);
    }
  }
  //tree.chars = node_ass;
  
  //double score = get_parsimony_score(Nchar, tree);
  //cout << "Initial assignment score:\t" << score << endl;
  //print_assignment( Nchar, ntotn, node_ass );
  
  //print_assignment(Nchar, ntotn, tree.chars);
  //cout << "tree 2" << endl;
  //evo_tree tree2 = tree;
  //print_assignment(Nchar, ntotn, tree2.chars);
  //evo_tree tree3 = perturb_tree(Ns, Nchar, tree2);
  //cout << "tree 3" << endl;
  //print_assignment(Nchar, ntotn, tree3.chars);

  
  // create initial population of trees. Sample from coalescent trees
  vector<evo_tree> trees;
  vector<double> distances(2*Npop,0);
  for(int i=0; i<Npop; ++i){
    //trees.push_back( evo_tree(Ns+1, e, l) );

    evo_tree rtree = generate_coal_tree(Ns);
    trees.push_back( rtree );
    trees[i].chars = node_ass;
    //cout << "Initial assignment score:\t" << get_parsimony_score(Nchar, trees[i]) << endl;
  }

  double min_dist = 1e20;
  evo_tree min_dist_tree;

  int count_static = 0;
  
  for(int g=0; g<Ngen; ++g){
    // Growth stage: create Npop copies + Npop copies with mutation
    vector<evo_tree> new_trees;
    for(int i=0; i<Npop; ++i){
      new_trees.push_back( trees[i] );
    }
    
    for(int i=0; i<Npop; ++i){
      new_trees.push_back( perturb_mutation_tree(Ns, Nchar, trees[i], max_cn) );
    }

    // Selection: score and keep the top Npop
    for(int i=0; i<2*Npop; ++i){
      distances[i] = get_parsimony_score(Nchar, new_trees[i]);
      //cout << "g/i/distance:\t" << g << "\t" << i << "\t" << distances[i] << endl;
    }
    
    vector<int> index(2*Npop);
    int x=0;
    iota( index.begin(), index.end(), x++);
    sort( index.begin(), index.end(), [&](int i,int j){return distances[i]<distances[j];} );

    // Selection: calculate mean fitness of top half of population
    double meand = 0;
    for(int k=0; k<Npop; ++k){
      meand += distances[ index[k] ];
    }
    meand = meand/Npop;
    if( g%100 == 0) cout << "g / av dist / top dist \t" << g << "\t" << meand << "\t" << distances[ index[0] ] << endl;

    // Selection: select top half
    for(int i=0; i<Npop; ++i){
      trees[i] = new_trees[ index[i] ];
    }

    // Selection: record the best (lowest) scoring tree
    if( distances[ index[0] ] < min_dist ){
      min_dist = distances[ index[0] ];
      min_dist_tree = new_trees[ index[0] ];
      count_static = 0;
    }else{
      //cout << "incrementing\t" << min_dist << endl;
      count_static += 1;
    }

    if( count_static == 100 ){
      cout << "\t### static distance function. Finishing on ngen = " << g << endl;
      break;
    }
    
    fill(distances.begin(), distances.end(), 0);
  }

  // Write out the top tree
  stringstream sstm;
  min_dist_tree.print();
  sstm << "sim-data-" << "1" << "-tree.txt";
  ofstream out_tree(sstm.str());
  min_dist_tree.write(out_tree);
  out_tree.close();
  sstm.str("");

  //cout << "Final population" << endl;
  //for(int i=0; i<Npop; ++i){
  //  double d = get_parsimony_score(Nchar, trees[i]);
  //  cout << i << "\t" << d << endl;
  //}
  
  //for(int i=0; i<10; ++i){
  //  trees[i].print();
  //  sstm << "sim-data-" << i+1 << "-tree.txt";
  //  ofstream out_tree(sstm.str());
  //  trees[i].write(out_tree);
  //  out_tree.close();
  //  sstm.str("");
  //}
  
  //vector<int> sim(4401*Ns,0);
  //run_sample_set(Ns, prc, pvs, e, l, sim);
  //double dst = edistance( data, sim );
  //cout << "distance:\t" << dst << endl;

  /*
  if(0){
  evo_tree t1(Ns+1, e, l);
  cout << "printing t1" << endl;
  t1.print();
  evo_tree t2 = t1;
  cout << "printing t2" << endl;
  t2.print();
  cout << "done" << endl;

  stringstream sstm;
  for(int i=0; i<10; ++i){
    evo_tree new_tree = perturb_tree( Ns, t2 );
    cout << "perturbed" << endl;
    new_tree.print();

    sstm << "sim-data-" << i+1 << "-tree.txt";
    ofstream out_tree(sstm.str());
    new_tree.write(out_tree);
    out_tree.close();
    sstm.str("");
  }
  }
  */

}
