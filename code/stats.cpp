#include "stats.hpp"

// long unsigned (*myrng)(long unsigned){ &myrng };


void setup_rng(gsl_rng* r, unsigned set_seed){
  if(set_seed){
    gsl_rng_set(r, set_seed);
  }else{
    int t = time(NULL);
    int pid = getpid();
    long s = t * pid;
    //cout << "pid:" << "\t" << getpid() << endl;
    // cout << "seed:" << "\t" << t << "\t" << pid << "\t" << abs(s) << endl;
    cout << "seed:" << "\t" << abs(s) << endl;
    gsl_rng_set(r, abs(s));
  }
}


// wrapper for uniform
double runiform(gsl_rng* r, double a, double b){
  double myrandom = a + (b - a) * gsl_rng_uniform(r);

  while(myrandom == 0){
    myrandom = a + (b - a) * gsl_rng_uniform(r);
  }

  return myrandom;
}



// sample an element according to a probability vector
// gsl_ran_multinomial?
int rchoose(gsl_rng* r, const vector<double>& rates){
  //cout << "rchoose, rates:";
  //for(int i = 0; i < rates.size(); ++i) cout << "\t" << rates[i];
  //cout << endl;

  vector<double> p;
  double s = accumulate(rates.begin(), rates.end(), 0.0);

  for(int i = 0; i < rates.size(); ++i){
    p.push_back(rates[i]/s);
  }

  vector<double> psum(p.size(),0.0);
  partial_sum(p.begin(),p.end(),psum.begin());

  double u = gsl_rng_uniform(r);
  int ret = -1;
  for(int i = 0; i < rates.size(); ++i){
    if(u < psum[i]){
      ret = i;
      break;
    }
  }

  //cout << "u=\t" << u << "\t" << ret << endl;
  return ret;
}
