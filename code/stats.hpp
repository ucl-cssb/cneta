#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include <algorithm>
#include <vector>

extern gsl_rng * r;


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
