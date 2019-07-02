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

#include <boost/program_options.hpp>

#include "gzstream.h"

#include "evo_tree.hpp"
#include "genome.hpp"
#include "stats.hpp"
#include "utilities.hpp"

using namespace std;

// The number of trees to search before terminating
const int MAX_TREE = 100;
const double MAX_NLNL = 1e20;
int debug = 0;

// global values for gsl function minimization
//vector<vector<int> > vobs;
//vector<double> tobs;
//int Ns;
//int Nchar;

// global value for tree search
map<string,int> searched_trees;
bool exhausted_tree_search = false;


evo_tree perturb_tree( const int& Ns, const evo_tree& tree ){
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
  new_tree.dup_rate = tree.dup_rate;
  new_tree.del_rate = tree.del_rate;
  new_tree.tobs = tree.tobs;

  //cout << "\tntree: ";
  //cout << create_tree_string( new_tree ) << endl;

  //string ordered = order_tree_string( create_tree_string( new_tree ) );
  //cout << "\totree: ";
  //cout << ordered << endl;
  //cout << endl;

  return new_tree;
}

evo_tree perturb_tree( const int& Ns, vector<evo_tree> trees ){
  if(debug) cout << "\tperturb_tree" << endl;

  int count = 0;
  while(true){
    // randomly sample the fit population
    int ind = gsl_rng_uniform_int(r, trees.size());

    // generate a new tree
    evo_tree ttree = perturb_tree( Ns, trees[ind] );
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

    if(count > MAX_TREE){
      //cout << "\tperturb_tree cannot find new topologies" << endl;
      exhausted_tree_search = true;
      return ttree;
    }
  }
}

evo_tree do_evolutionary_algorithm(const int& Npop, const int& Ngen, const int& max_static, const vector<double>& rates, const double ssize, const double tolerance, const int miter, const int optim, const int model, const int cons, const int maxj, const int cn_max, const int correct_bias){
  //cout << "Running evolutionary algorithm" << endl;
  // create initial population of trees. Sample from coalescent trees
  vector<evo_tree> trees;

  vector<double> lnLs(2*Npop,0);
  for(int i=0; i<Npop; ++i){
    evo_tree rtree = generate_coal_tree(Ns);
    rtree.tobs = tobs;
    double tree_height = rtree.get_tree_height();
    double min_height = *max_element(tobs.begin(), tobs.end());
    if(cons){
        double max_height = age + min_height;
        double new_tree_height = runiform(r, min_height, max_height);
        double ratio = new_tree_height/tree_height;
        rtree.scale_time_internal(ratio);
        tree_height = new_tree_height;
    }

    if(rates.size()>1){
        rtree.dup_rate = rates[0];
        rtree.del_rate = rates[1];
        rtree.chr_gain_rate = rates[2];
        rtree.chr_loss_rate = rates[3];
        rtree.wgd_rate = rates[4];
        rtree.mu = accumulate(rates.begin(), rates.end(), 0.0);
    }
    else{
        rtree.mu = rates[0];
    }

    // double curr_ttime = rtree.get_total_time();
    // evo_tree rtree_adjusted(rtree.nleaf, rtree.edges, curr_ttime, tobs);
    trees.push_back( rtree );

    string tstring = order_tree_string( create_tree_string( rtree ) );
    searched_trees[ tstring ] = 1;
  }
  double min_lnL = MAX_NLNL;
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
    	new_trees.push_back( perturb_tree(Ns, trees) );
      }

      // Selection: score the trees
      for(int i=0; i<2*Npop; ++i){
        evo_tree otree;
        if(optim == 0){
          otree = max_likelihood(new_trees[i], model, Lf, ssize, tolerance, miter, cons, maxj, cn_max, correct_bias);
        }
        if(optim == 1){
          otree = max_likelihood_BFGS(new_trees[i], model, Lf, tolerance, miter, cons, maxj, cn_max, correct_bias);
        }
        otree.score = Lf;
        // cout << "otree tobs " << otree.tobs[0] << endl;
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
        	new_trees.push_back( perturb_tree(Ns, trees) );
            // new_trees of size 2 Npop
            evo_tree otree;
            if(optim == 0){
        	    otree = max_likelihood(new_trees[Npop + i], model, Lf, ssize, tolerance, miter, cons, maxj, cn_max, correct_bias);
            }
            if(optim == 1){
                otree = max_likelihood_BFGS(new_trees[Npop + i], model, Lf, tolerance, miter, cons, maxj, cn_max, correct_bias);
            }
        	otree.score = Lf;
            // cout << "otree tobs " << otree.tobs[0] << endl;
        	lnLs[Npop + i] = Lf;
        	opt_trees.push_back( otree );
        	//cout << "g/i/lnL:\t" << g << "\t" << Npop+i << "\t" << lnLs[i] << endl;
        }
      // Randomly sample again?
    }

    vector<int> index(2*Npop);
    int x=0;
    iota( index.begin(), index.end(), x++);
    sort( index.begin(), index.end(), [&](int i,int j){ return lnLs[i]<lnLs[j];} );

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


void run_test(const string& tree_file, int Ns, int Nchar, int model, int cn_max, vector<double> tobs, vector<double> rates, double ssize, double tolerance, double miter){
    // MLE testing
    //static const int arr1[] = {8,5, 8,1, 9,2, 9,3, 10,9, 10,8, 11,4, 11,10, 7,11, 7,6 };
    //vector<int> e (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
    //for(int i=0; i<e.size();++i) e[i] = e[i] - 1;

    //static const double arr2[] = {18.49, 38.49, 51.71, 31.71, 0.51, 3.73, 22.2, 0.013, 0.99, 0};
    //vector<double> l (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

    //evo_tree test_tree(Ns+1, e, l, 1);
    //test_tree.tobs = tobs;
    //test_tree.print();

    // string tree_file = "./test/sim-data-1-tree.txt";
    double dup_rate = rates[0];
    double del_rate = rates[1];
    double chr_gain_rate = rates[2];
    double chr_loss_rate = rates[3];
    double wgd_rate = rates[4];

    // read in true tree
    evo_tree test_tree = read_tree_info(tree_file, Ns);
    test_tree.print();
    if(model==0){
        test_tree.mu = 1.0/Nchar;
    }
    if(model==1){
        test_tree.dup_rate = dup_rate;
        test_tree.del_rate = del_rate;
        test_tree.mu = dup_rate + del_rate;
    }
    test_tree.tobs = tobs;

    double Ls = 0.0;

    // Ls = get_likelihood(Ns, Nchar, vobs0, test_tree, model, 0);
    // cout << "\nOriginal tree -ve likelihood without incorporating chr gain/loss and WGD: " << -Ls << endl;

    Ls = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, test_tree, model, 0, cn_max, correct_bias);
    cout << "\nOriginal tree -ve likelihood with revised method: " << -Ls << endl;

    //for(int i=0; i<10; ++i){
    //double mu_p = mu + (i-5)*0.01*mu;
    //Ls = get_likelihood(Ns, Nchar, vobs, test_tree, mu_p, model, 0);
    //cout << "\n -ve likelihood: " << mu_p << "\t" << -Ls << endl;
    //}

    double Lf = 0;
    stringstream sstm;
    ofstream out_tree;

    cout << "\n\n### Running optimisation: branches free, mu fixed" << endl;
    evo_tree min_tree = max_likelihood(test_tree, model, Lf, ssize, tolerance, miter, 0, 0, cn_max, correct_bias);
    cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree.mu <<  endl;
    min_tree.print();

    sstm << "./test/sim-data-" << "00" << "-tree.txt";
    out_tree.open(sstm.str());
    min_tree.write(out_tree);
    out_tree.close();
    sstm.str("");

    test_tree = read_tree_info(tree_file, Ns);
    test_tree.dup_rate = dup_rate;
    test_tree.del_rate = del_rate;
    test_tree.mu = dup_rate + del_rate;
    evo_tree min_tree_bfgs = max_likelihood_BFGS(test_tree, model, Lf, tolerance, miter, 0, 0, cn_max, correct_bias);
    cout << "\nMinimised tree likelihood / mu by BFGS : " << Lf << "\t" << min_tree_bfgs.dup_rate << "\t" << min_tree_bfgs.del_rate << endl;
    min_tree_bfgs.print();

    sstm << "./test/sim-data-" << "00" << "-BFGS-tree.txt";
    out_tree.open(sstm.str());
    min_tree_bfgs.write(out_tree);
    out_tree.close();
    sstm.str("");

    cout << "\n\n### Running optimisation: branches free, mu free" << endl;

    test_tree = read_tree_info(tree_file, Ns);
    evo_tree min_tree_mu = max_likelihood(test_tree, model, Lf, ssize, tolerance, miter, 0, 1, cn_max, correct_bias);
    cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree_mu.mu << endl;
    min_tree_mu.print();

    sstm << "./test/sim-data-" << "01" << "-tree.txt";
    out_tree.open(sstm.str());
    min_tree_mu.write(out_tree);
    out_tree.close();
    sstm.str("");

    test_tree = read_tree_info(tree_file, Ns);
    min_tree_mu = max_likelihood_BFGS(test_tree, model, Lf, tolerance, miter, 0, 1, cn_max, correct_bias);
    cout << "\nMinimised tree likelihood / mu by BFGS : " << Lf << "\t" << min_tree_mu.dup_rate << "\t" << min_tree_mu.del_rate << endl;
    min_tree_mu.print();

    sstm << "./test/sim-data-" << "01" << "-BFGS-tree.txt";
    out_tree.open(sstm.str());
    min_tree_mu.write(out_tree);
    out_tree.close();
    sstm.str("");


    cout << "\n\n### Running optimisation: branches constrained, mu fixed" << endl;
    test_tree = read_tree_info(tree_file, Ns);
    test_tree.tobs = tobs;
    test_tree.mu = 1.0/Nchar;
    evo_tree min_tree2 = max_likelihood(test_tree, model, Lf, ssize, tolerance, miter, 1, 0, cn_max, correct_bias);
    cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree2.mu <<endl;
    min_tree2.print();

    sstm << "./test/sim-data-" << "10" << "-tree.txt";
    out_tree.open(sstm.str());
    min_tree2.write(out_tree);
    out_tree.close();
    sstm.str("");

    test_tree = read_tree_info(tree_file, Ns);
    test_tree.tobs = tobs;
    test_tree.dup_rate = dup_rate;
    test_tree.del_rate = del_rate;
    test_tree.mu = dup_rate + del_rate;
    min_tree2 = max_likelihood_BFGS(test_tree, model, Lf, tolerance, miter, 1, 0, cn_max, correct_bias);
    cout << "\nMinimised tree likelihood / mu by BFGS: " << Lf << "\t" << min_tree2.dup_rate << "\t" << min_tree2.del_rate << endl;
    min_tree2.print();

    sstm << "./test/sim-data-" << "10" << "-BFGS-tree.txt";
    out_tree.open(sstm.str());
    min_tree2.write(out_tree);
    out_tree.close();
    sstm.str("");


    cout << "\n\n### Running optimisation: branches constrained, mu free" << endl;
    test_tree = read_tree_info(tree_file, Ns);
    test_tree.tobs = tobs;
    evo_tree min_tree2_mu = max_likelihood(test_tree, model, Lf, ssize, tolerance, miter, 1, 1, cn_max, correct_bias);
    cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree2_mu.mu*Nchar <<endl;
    min_tree2_mu.print();

    sstm << "./test/sim-data-" << "11" << "-tree.txt";
    out_tree.open(sstm.str());
    min_tree2_mu.write(out_tree);
    out_tree.close();
    sstm.str("");

    test_tree = read_tree_info(tree_file, Ns);
    test_tree.tobs = tobs;
    min_tree2_mu = max_likelihood_BFGS(test_tree, model, Lf, tolerance, miter, 1, 1, cn_max, correct_bias);
    cout << "\nMinimised tree likelihood / mu by BFGS : " << Lf << "\t" << min_tree2_mu.dup_rate << "\t" << min_tree2_mu.del_rate << endl;
    min_tree2_mu.print();

    sstm << "./test/sim-data-" << "11" << "-BFGS-tree.txt";
    out_tree.open(sstm.str());
    min_tree2_mu.write(out_tree);
    out_tree.close();
    sstm.str("");
}


// Compute the likelihood of a tree given the observed copy number profile
double compute_tree_likelihood(const string& tree_file, int Ns, int Nchar, int num_invar_bins, map<int, vector<vector<int>>>& vobs, vector<double>& tobs, vector<double>& rates, int model, int cons, int cn_max, int correct_bias){
    evo_tree tree;
    // string format = tree_file.substr(tree_file.find_last_of(".") + 1);
    // cout << "Format of the input tree file is " << format << endl;

    // if (format == "txt"){
        tree = read_tree_info(tree_file, Ns);
        tree.print();
        tree.tobs = tobs;
    // }
    // if (format == "nex"){
    //     tree = read_nexus(tree_file);
    // }

    if(model==0){
        tree.mu = 1.0/Nchar;
    }
    if(model==1){
        tree.dup_rate = rates[0];
        tree.del_rate = rates[1];
        tree.chr_gain_rate = rates[2];
        tree.chr_loss_rate = rates[3];
        tree.wgd_rate = rates[4];
        tree.mu = accumulate(rates.begin(), rates.end(), 0.0);
    }

    double Ls = 0.0;
    Ls = get_likelihood_revised(Ns, Nchar, num_invar_bins, vobs, tree, model, cons, cn_max, correct_bias);

    return Ls;
}

double maximize_tree_likelihood(const string& tree_file, const string& ofile, int Ns, int Nchar, vector<double>& tobs, vector<double>& rates, double ssize, double tolerance, int miter, int model, int optim, int cons, int maxj, int cn_max, int correct_bias){
    evo_tree tree;
    // string format = tree_file.substr(tree_file.find_last_of(".") + 1);
    // cout << "Format of the input tree file is " << format << endl;

    // if (format == "txt"){
        tree = read_tree_info(tree_file, Ns);
        tree.print();
        tree.tobs = tobs;
    // }
    // if (format == "nex"){
    //     tree = read_nexus(tree_file);
    // }

    if(model==0){
        tree.mu = 1.0/Nchar;
    }
    if(model==1){
        tree.dup_rate = rates[0];
        tree.del_rate = rates[1];
        tree.chr_gain_rate = rates[2];
        tree.chr_loss_rate = rates[3];
        tree.wgd_rate = rates[4];
        tree.mu = accumulate(rates.begin(), rates.end(), 0.0);
    }

    double Lf = 0;
    evo_tree min_tree;
    if(optim == 1){
        min_tree = max_likelihood_BFGS(tree, model, Lf, tolerance, miter, cons, maxj, cn_max, correct_bias);
    }
    if(optim == 0){
        min_tree = max_likelihood(tree, model, Lf, ssize, tolerance, miter, cons, maxj, cn_max, correct_bias);
    }
    cout << "\nMinimised tree likelihood: " << Lf << endl;
    if(maxj == 1){
        cout << "Estimated mutation rates: " << endl;
        cout << "\t" << min_tree.dup_rate << "\t" << min_tree.del_rate << "\t" << min_tree.chr_gain_rate << "\t" << min_tree.chr_loss_rate << "\t" << min_tree.wgd_rate << endl;
    }
    min_tree.print();

    ofstream out_tree(ofile);
    min_tree.write(out_tree);
    out_tree.close();
}


int main (int argc, char ** const argv) {
  int Npop, Ngen, max_static, miter, bootstrap, seed, cons, maxj, optim, correct_bias, mode, cn_max;
  double tolerance, ssize, mu, dup_rate, del_rate, chr_gain_rate, chr_loss_rate, wgd_rate;
  string datafile, timefile, ofile, tree_file;

  int test = 0; // Whether or not to run test mode

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
    ("ofile,o", po::value<string>(&ofile)->default_value("results-maxL-tree.txt"), "output tree file with maximum likelihood")
    ("tree_file", po::value<string>(&tree_file)->default_value(""), "input tree file ")
    ("cn_max", po::value<int>(&cn_max)->default_value(4), "maximum copy number of a segment")

    ("npop,p", po::value<int>(&Npop)->default_value(100), "number of population")
    ("ngen,g", po::value<int>(&Ngen)->default_value(50), "number of generation")
    ("nstop,e", po::value<int>(&max_static)->default_value(3), "stop condition")
    ("tolerance,r", po::value<double>(&tolerance)->default_value(1e-2), "tolerance value")
    ("miter,m", po::value<int>(&miter)->default_value(2000), "maximum number of iterations in maximization")
    ("ssize,z", po::value<double>(&ssize)->default_value(0.01), "initial step size")

    ("mu,x", po::value<double>(&mu)->default_value(0.025), "mutation rate (allele/locus/year)")
    ("dup_rate", po::value<double>(&dup_rate)->default_value(0.01), "duplication rate (allele/locus/year)")
    ("del_rate", po::value<double>(&del_rate)->default_value(0.01), "deletion rate (allele/locus/year)")
    ("chr_gain_rate", po::value<double>(&chr_gain_rate)->default_value(0), "chromosome gain rate")
    ("chr_loss_rate", po::value<double>(&chr_loss_rate)->default_value(0), "chromosome loss rate")
    ("wgd_rate", po::value<double>(&wgd_rate)->default_value(0), "WGD (whole genome doubling) rate")

    ("model,d", po::value<int>(&model)->default_value(1), "model of evolution (0: JC69, 1: 1-step bounded)")
    ("constrained", po::value<int>(&cons)->default_value(1), "constraints on branch length (0: none, 1: fixed total time)")
    ("fixm", po::value<int>(&maxj)->default_value(0), "estimation of mutation rate (0: mutation rate fixed to be the given value, 1: estimating mutation rate)")
    ("optim", po::value<int>(&optim)->default_value(1), "method of optimization (0: Simplex, 1: L-BFGS-B)")
    ("correct_bias", po::value<int>(&correct_bias)->default_value(1), "correct ascertainment bias")

    ("bootstrap,b", po::value<int>(&bootstrap)->default_value(0), "doing bootstrap or not")
    ("mode", po::value<int>(&mode)->default_value(1), "running mode of the program (0: Test on example data; 1: Compute maximum likelihood tree from copy number profile; 2: Compute the likelihood of a given tree with branch length; 3: Compute the maximum likelihood of a given tree)")
    ("verbose", po::value<int>(&debug)->default_value(0), "verbose level (0: default, 1: debug)")
    ("seed", po::value<int>(&seed)->default_value(0), "seed used for generating random numbers")
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
          cout << "svtreeml [version 0.1], a program to build a maximum likelihood phylogenetic tree from copy number profile" << endl;
          cout << "This program can also be used to compute the likelihood of a tree with fixed branch lengths or maximum likelihood of a tree with fixed topology (and/or mutation rates), given the copy number profile." << endl;
          return 1;
      }
      po::notify(vm);
  } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
  }

  setup_rng(seed);

  int num_invar_bins = 0;
  Nchar = 0;
  map<int, vector<vector<int>>> data = read_data_var_regions_by_chr(datafile, Ns, cn_max, num_invar_bins);
  // vector<vector<int>> data0 = read_data_var_regions(datafile, Ns, CN_MAX);
  // Nchar = data.size();

  // tobs already defined globally
  tobs = read_time_info(timefile, Ns, age);
  if(cons){
      cout << "The age of patient at the first sampling time: " << age << endl;
  }

  vector<double> rates;
  if(model==0){
      rates.push_back(mu);
      cout << "Assuming JC69 model " << endl;
  }
  if(model==1){
      rates.push_back(dup_rate);
      rates.push_back(del_rate);
      rates.push_back(chr_gain_rate);
      rates.push_back(chr_loss_rate);
      rates.push_back(wgd_rate);
      cout << "Assuming Bounded model " << endl;
  }

  if (cons==0){
      cout << "Assuming the tree is unconstrained " << endl;
  }
  else{
      cout << "Assuming the tree is constrained by age " << endl;
  }

  if (maxj==0){
      cout << "Assuming mutation rate (allele/locus/year) is fixed " << endl;
  }
  else{
      cout << "Estimating mutation rate (allele/locus/year)" << endl;
  }

  if (optim==0){
      cout << "Using Simplex method for optimization " << endl;
  }
  else{
      cout << "Using L-BFGS-B method for optimization" << endl;
  }

  if (correct_bias==0){
      cout << "Not correcting acquisition bias in likelihood computation " << endl;
  }
  else{
      cout << "Correcting acquisition bias in likelihood computation " << endl;
  }

  // Construct the CN matrix
  // vector<vector<int>> vobs0;
  // for(int nc=0; nc<data0.size(); ++nc){
  //   vector<int> obs;
  //   for(int i=0; i<Ns; ++i){
  //     obs.push_back( data0[nc][i+3] );
  //   }
  //   vobs0.push_back( obs );
  // }
  // cout << "The number of sites used: " << data0.size() << endl;
  // Construct the CN matrix by chromosme
  for(int nchr=1; nchr<=data.size(); nchr++){
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


  if(mode == 0){
      cout << "Running test " << endl;
      run_test(tree_file, Ns, Nchar, model, cn_max, tobs, rates, ssize, tolerance, miter);
  }

  if (mode == 1){
      // cout << "The number of sites used after grouping by chromosme: " << Nchar << endl;
      // Do a bootstrap replicate
      // if(bootstrap){
      //     // create a copy of vobs to resample
      //     vector<vector<int> > vobs_copy = vobs;
      //     vobs.clear();
      //     // cout << "rows selected: ";
      //     for(int nc=0; nc<Nchar; ++nc){
      //         // randomly select a site
      //        int i = gsl_rng_uniform_int(r, Nchar);
      //        // cout << i << "\t";
      //        vobs.push_back(vobs_copy[i]);
      //     }
      //     // cout << endl;
      // }
      cout << "Building maximum likelihood tree from copy number profile " << endl;
      if(bootstrap){
          // create a copy of vobs to resample
          map<int, vector<vector<int>>> vobs_copy = vobs;
          vobs.clear();
          for(int nchr=1; nchr<=data.size(); nchr++){
            vector<vector<int>> obs_chr;
            for(int nc=0; nc<data[nchr].size(); ++nc){
                  // randomly select a site
                 int i = gsl_rng_uniform_int(r, data[nchr].size());
                 // cout << i << "\t";
                 obs_chr.push_back(vobs_copy[nchr][i]);
            }
            vobs[nchr] = obs_chr;
         }
          // cout << endl;
      }

      evo_tree min_lnL_tree = do_evolutionary_algorithm(Npop, Ngen, max_static, rates, ssize, tolerance, miter, optim, model, cons, maxj, cn_max, correct_bias);
      if(maxj==1){
          if(model==0){
              cout << "Estimated mutation rate (allele/locus/year):  " << min_lnL_tree.mu << endl;
          }
          if(model==1){
              cout << "Estimated duplication rate (allele/locus/year):  " << min_lnL_tree.dup_rate << endl;
              cout << "Estimated deletion rate (allele/locus/year):  " << min_lnL_tree.del_rate << endl;
              cout << "Estimated chromosme gain rate (year):  " << min_lnL_tree.chr_gain_rate << endl;
              cout << "Estimated chromosme loss rate (year):  " << min_lnL_tree.chr_loss_rate << endl;
              cout << "Estimated whole genome doubling rate (year):  " << min_lnL_tree.wgd_rate << endl;
          }
      }
      // Write out the top tree
      //cout << "Best fitting tree, -ve lnL = " << global_min << endl;
      min_lnL_tree.print();

      // Check the validity of the tree
      if(cons==1){
          if(!is_tree_valid(min_lnL_tree, cons)){
              cout << "The final tree is not valid!" << endl;
          }
      }
      ofstream out_tree(ofile);
      min_lnL_tree.write(out_tree);
      out_tree.close();
  }

  if(mode == 2){
      cout << "Computing the likelihood of a given tree from copy number profile " << endl;
      double Lf = compute_tree_likelihood(tree_file, Ns, Nchar, num_invar_bins, vobs, tobs, rates, model, cons, cn_max, correct_bias);
      cout << "The negative log likelihood of the input tree is " << -Lf << endl;
  }

  if(mode == 3){
      cout << "Computing maximum likelihood of a given tree from copy number profile " << endl;
      maximize_tree_likelihood(tree_file, ofile, Ns, Nchar, tobs, rates, ssize, tolerance, miter, model, optim, cons, maxj, cn_max, correct_bias);
  }
}
