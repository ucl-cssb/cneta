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

#include <boost/program_options.hpp>

#include "gzstream.h"

#include "evo_tree.hpp"
#include "genome.hpp"
#include "stats.hpp"
#include "utilities.hpp"

using namespace std;

int debug = 0;

// global values for gsl function minimization
//vector<vector<int> > vobs;
//vector<double> tobs;
//int Ns;
//int Nchar;

int main (int argc, char ** const argv) {
  int miter, nmax, seed, cn_max;
  double tolerance, ssize, mu_0, vlnorm;
  string datafile, timefile, treefile, ofile;

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
     ("pfile,p", po::value<string>(&treefile)->required(), "input tree file")
     ;
  po::options_description optional("Optional parameters");
  optional.add_options()
    ("nsample,s", po::value<int>(&Ns)->default_value(5), "number of samples or regions")
    ("cn_max", po::value<int>(&cn_max)->default_value(4), "maximum copy number of a segment")
    ("tolerance,r", po::value<double>(&tolerance)->default_value(1e-2), "tolerance value")
    ("miter,m", po::value<int>(&miter)->default_value(2000), "maximum number of iterations in maximization")
    ("nmax,n", po::value<int>(&nmax)->default_value(100), "number of maximizations to attempt")
    ("ssize,z", po::value<double>(&ssize)->default_value(0.01), "initial step size")
    ("ofile,o", po::value<string>(&ofile)->default_value("results-maxL-mu-tree.txt"), "output tree file")
    ("mu,x", po::value<double>(&mu_0)->default_value(1.0), "initial mutation rate estimate (SCA/locus/time)")
    ("vlnorm,l", po::value<double>(&vlnorm)->default_value(1.0), "scale of lognorm for initial value sampling")
    ("model,d", po::value<int>(&model)->default_value(0), "model of evolution (0: JC69, 1: 1-step bounded)")
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
          cout << "svtreeml [version 0.1], a program to build a phylogenetic tree from copy number profile" << endl;
          return 1;
      }
      po::notify(vm);
  } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
  }

  setup_rng(seed);

  // vector<vector<int>> data = read_data_var_regions(datafile, Ns, CN_MAX);
  // Nchar = data.size();
  int num_invar_bins = 0;
  Nchar = 0;
  map<int, vector<vector<int>>> data = read_data_var_regions_by_chr(datafile, Ns, cn_max, num_invar_bins);

  // tobs already defined globally
  tobs = read_time_info(timefile, Ns, age);
  cout << "The age of patient at the first sampling time: " << age << endl;

  // read in mle tree
  evo_tree test_tree = read_tree_info(treefile, Ns);

  //vector<vector<int> > vobs; // already defined globally
  // for(int nc=0; nc<Nchar; ++nc){
  //   vector<int> obs;
  //   for(int i=0; i<Ns; ++i){
  //     obs.push_back( data[nc][i+3] );
  //   }
  //   vobs.push_back( obs );
  // }

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

  // estimate mutation rate
  //test_tree.print();
  test_tree.tobs = tobs;
  test_tree.mu = mu_0;

  //double Ls = get_likelihood(Ns, Nchar, vobs, test_tree);
  //cout << "\nOriginal tree -ve likelihood: " << -Ls << endl;

  cout << "\n\n### Running optimisation: branches constrained, mu free" << endl;
  cout << "vlnorm: " << vlnorm << endl;

  double minL = -1*LARGE_LNL;
  evo_tree min_tree_mu(test_tree);

  for(int i=0; i<nmax; ++i){
    double Lf = 0;
    double mu_g;
    if(i==0){
      mu_g = mu_0;
    }else{
      mu_g = gsl_ran_lognormal(r, log(mu_0), vlnorm);
    }
    test_tree.mu = mu_g;

    //evo_tree min_tree = max_likelihood(test_tree, Lf, 1, 1);
    evo_tree min_tree = max_likelihood(test_tree, Lf, ssize, tolerance, miter, 1, 1);
    cout << "Testing mu_g / -lnL / mu: " << mu_g << " / " << Lf << " / " << min_tree.mu << endl;

    if(Lf < minL){
      minL = Lf;
      min_tree_mu = min_tree;
      //cout << "\nMinimised tree likelihood / mu : " << Lf << "\t" << min_tree_mu.mu*Nchar <<endl;
      //min_tree_mu.print();
    }

  }

  cout << "\nMinimised tree likelihood / Nchar / mu (SCA/locus/time): " << minL << " / " << Nchar << " / " <<  min_tree_mu.mu << endl;
  min_tree_mu.print();

  stringstream sstm;
  ofstream out_tree;
  out_tree.open(ofile);
  min_tree_mu.write(out_tree);
  out_tree.close();
  sstm.str("");
}
