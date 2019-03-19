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
//#include <gsl/gsl_multimin.h>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "gzstream.h"

#include "evo_tree.hpp"
#include "genome.hpp"
#include "stats.hpp"
#include "utilities.hpp"

using namespace std;

int debug = 0;
const int CN_MAX = 4;

// global values for gsl function minimization
//vector<vector<int> > vobs;
//vector<double> tobs;
//int Ns;
//int Nchar;

// global value for tree search
map<string,int> searched_trees;
bool exhausted_tree_search = false;

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
  new_tree.mu = tree.mu;
  new_tree.tobs = tree.tobs;

  //cout << "\tntree: ";
  //cout << create_tree_string( new_tree ) << endl;

  //string ordered = order_tree_string( create_tree_string( new_tree ) );
  //cout << "\totree: ";
  //cout << ordered << endl;
  //cout << endl;

  return new_tree;
}

evo_tree perturb_tree( const int& Ns, const int& Nchar, vector<evo_tree> trees ){
  if(debug) cout << "\tperturb_tree" << endl;

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

evo_tree do_evolutionary_algorithm(const int& Npop, const int& Ngen, const int& max_static, const double& mu, const double ssize, const double tolerance, const int miter, int cons=0, int maxj=0){

  //cout << "Running evolutionary algorithm" << endl;

  // create initial population of trees. Sample from coalescent trees
  vector<evo_tree> trees;

  vector<double> lnLs(2*Npop,0);
  for(int i=0; i<Npop; ++i){
    evo_tree rtree = generate_coal_tree(Ns);
    rtree.mu = mu;
    rtree.tobs = tobs;
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
	evo_tree otree = max_likelihood(new_trees[i], Lf, ssize, tolerance, miter, cons, maxj);
	otree.score = Lf;
	lnLs[i] = Lf;
	opt_trees.push_back( otree );
	//cout << "g/i/lnL:\t" << g << "\t" << i << "\t" << lnLs[i] << endl;
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

	evo_tree otree = max_likelihood(new_trees[Npop + i], Lf, ssize, tolerance, miter, cons, maxj);     // new_trees of size 2 Npop
	otree.score = Lf;
	lnLs[Npop + i] = Lf;
	opt_trees.push_back( otree );
	//cout << "g/i/lnL:\t" << g << "\t" << Npop+i << "\t" << lnLs[i] << endl;
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

/*
vector<edge> create_edges_from_nodes( const vector<node>& nodes, const vector<double>& node_times ){
  vector<edge> enew;

  // Add leaf and internal nodes
  int root_id = 0;
  int id = 0;
  for(int i=0; i<nodes.size(); ++i){
    if(nodes[i].parent == -1) root_id = nodes[i].id;
  }

  for(int i=0; i<nodes.size(); ++i){
    if(nodes[i].id != root_id && nodes[i].id != root_id-1 && nodes[i].parent != root_id){
      enew.push_back( edge(id, nodes[i].parent, nodes[i].id, node_times[nodes[i].id] - node_times[nodes[i].parent]) );
      id++;
    }
  }
  // Add root nodes
  for(int i=0; i<nodes.size(); ++i){
  if( nodes[i].parent == root_id && nodes[i].id != root_id-1){
    enew.push_back( edge(id, nodes[i].parent, nodes[i].id, node_times[nodes[i].id] - node_times[nodes[i].parent]) );
    id++;
  }
  }
  enew.push_back( edge(id, root_id, Ns, 0) );

  return enew;
}
*/

int main (int argc, char ** const argv) {
  int Npop, Ngen, max_static, miter, bootstrap;
  double tolerance, ssize, mu;
  string datafile, timefile, ofile;

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
    ("nsample,s", po::value<int>(&Ns)->default_value(5), "number of samples or regions")
    ("npop,p", po::value<int>(&Npop)->default_value(100), "number of population")
    ("ngen,g", po::value<int>(&Ngen)->default_value(50), "number of generation")
    ("nstop,e", po::value<int>(&max_static)->default_value(3), "stop condition")
    ("tolerance,r", po::value<double>(&tolerance)->default_value(1e-2), "tolerance value")
    ("miter,m", po::value<int>(&miter)->default_value(2000), "maximum number of iterations in maximization")
    ("ssize,z", po::value<double>(&ssize)->default_value(0.01), "initial step size")
    ("bootstrap,b", po::value<int>(&bootstrap)->default_value(0), "doing bootstrap or not")
    ("ofile,o", po::value<string>(&ofile)->default_value("results-maxL-tree.txt"), "output tree file")
    ("mu,x", po::value<double>(&mu)->default_value(0.025), "mutation rate (SCA/locus/time)")
    ;

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(required).add(optional);
  po::variables_map vm;

  try {
      po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);

      if(vm.count("help")){
          cout << cmdline_options << endl;
          return 1;
      }

      if(vm.count("version")){
          cout << "svtreeml [version 0.1], a program to build a phylogenetic tree from copy number profile" << endl;
          return 1;
      }

      po::notify(vm);

      datafile = vm["cfile"].as<string>();
      timefile = vm["tfile"].as<string>();
      Ns = vm["nsample"].as<int>();
      Npop = vm["npop"].as<int>();
      Ngen = vm["ngen"].as<int>();
      max_static = vm["nstop"].as<int>();
      tolerance = vm["tolerance"].as<double>();
      miter = vm["miter"].as<int>();
      ssize = vm["ssize"].as<double>();
      // cout << "Input: " << endl;
      // cout << " Data file: " << datafile << endl;
      // cout << " Time file: " << timefile << endl;
      // cout << " Number of samples: " << Ns << endl;
      // cout << " Number of population: " << Npop << endl;
      // cout << " Number of generation: " << Ngen << endl;
      // cout << " Tolerance value: " << tolerance << endl;
  } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
  }

  setup_rng(0);

  vector<vector<int> > data = read_data_var_regions(datafile, Ns, CN_MAX);
  Nchar = data.size();

  // tobs already defined globally
  tobs = read_time_info(timefile,Ns);

  //vector<vector<int> > vobs; // already defined globally
  for(int nc=0; nc<Nchar; ++nc){
    vector<int> obs;
    for(int i=0; i<Ns; ++i){
      obs.push_back( data[nc][i+3] );
    }
    vobs.push_back( obs );
  }

  // Do a bootstrap replicate
  if(bootstrap){
      // create a copy of vobs to resample
      vector<vector<int> > vobs_copy = vobs;
      vobs.clear();
      // cout << "rows selected: ";
      for(int nc=0; nc<Nchar; ++nc){
          // randomly select a row
         int i = gsl_rng_uniform_int(r, Nchar);
         // cout << i << "\t";
         vobs.push_back(vobs_copy[i]);
      }
      // cout << endl;
  }

  if(0){
    // MLE testing

    //static const int arr1[] = {8,5, 8,1, 9,2, 9,3, 10,9, 10,8, 11,4, 11,10, 7,11, 7,6 };
    //vector<int> e (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
    //for(int i=0; i<e.size();++i) e[i] = e[i] - 1;

    //static const double arr2[] = {18.49, 38.49, 51.71, 31.71, 0.51, 3.73, 22.2, 0.013, 0.99, 0};
    //vector<double> l (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

    //evo_tree test_tree(Ns+1, e, l, 1);
    //test_tree.tobs = tobs;
    //test_tree.print();

    // read in true tree
    evo_tree test_tree = read_tree_info("./test/sim-data-1-tree.txt",Ns);
    test_tree.print();

    test_tree.mu = 1.0/Nchar;
    double Ls = 0.0;

    Ls = get_likelihood(Ns, Nchar, vobs, test_tree);
    cout << "\nOriginal tree -ve likelihood: " << -Ls << endl;

    //for(int i=0; i<10; ++i){
    //double mu_p = mu + (i-5)*0.01*mu;
    //Ls = get_likelihood(Ns, Nchar, vobs, test_tree, mu_p );
    //cout << "\n -ve likelihood: " << mu_p << "\t" << -Ls << endl;
    //}

    double Lf = 0;
    stringstream sstm;
    ofstream out_tree;

    cout << "\n\n### Running optimisation: branches free, mu fixed" << endl;
    evo_tree min_tree = max_likelihood(test_tree, Lf, ssize, tolerance, miter, 0, 0);
    cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree.mu*Nchar <<  endl;
    min_tree.print();

    sstm << "sim-data-" << "00" << "-tree.txt";
    out_tree.open(sstm.str());
    min_tree.write(out_tree);
    out_tree.close();
    sstm.str("");

    cout << "\n\n### Running optimisation: branches free, mu free" << endl;
    //reinitialize estimate of mu
    test_tree.mu = 1.0/Nchar;
    //test_tree.mu = 1;

    evo_tree min_tree_mu = max_likelihood(test_tree, Lf, ssize, tolerance, miter, 0, 1);
    cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree_mu.mu*Nchar << endl;
    min_tree_mu.print();

    sstm << "sim-data-" << "01" << "-tree.txt";
    out_tree.open(sstm.str());
    min_tree_mu.write(out_tree);
    out_tree.close();
    sstm.str("");

    cout << "\n\n### Running optimisation: branches constrained, mu fixed" << endl;
    test_tree.tobs = tobs;

    evo_tree min_tree2 = max_likelihood(test_tree, Lf, ssize, tolerance, miter, 1, 0);
    cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree2.mu*Nchar <<endl;
    min_tree2.print();

    sstm << "sim-data-" << "10" << "-tree.txt";
    out_tree.open(sstm.str());
    min_tree2.write(out_tree);
    out_tree.close();
    sstm.str("");

    cout << "\n\n### Running optimisation: branches constrained, mu free" << endl;
    test_tree.tobs = tobs;

    evo_tree min_tree2_mu = max_likelihood(test_tree, Lf, ssize, tolerance, miter, 1, 1);
    cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree2_mu.mu*Nchar <<endl;
    min_tree2_mu.print();

    sstm << "sim-data-" << "11" << "-tree.txt";
    out_tree.open(sstm.str());
    min_tree2_mu.write(out_tree);
    out_tree.close();
    sstm.str("");
  }

  if(1){
    cout << "Assuming mu (SCA/locus/time):  " << mu << endl;
    evo_tree min_lnL_tree = do_evolutionary_algorithm(Npop, Ngen, max_static,mu, ssize, tolerance, miter, 1, 0);

    // Write out the top tree
    //cout << "Best fitting tree, -ve lnL = " << global_min << endl;
    min_lnL_tree.print();
    ofstream out_tree(ofile);
    min_lnL_tree.write(out_tree);
    out_tree.close();
  }

}
