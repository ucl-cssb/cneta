extern "C" {

  gsl_rng * r;

  void setup_rng(int set_seed);
  
  //void run_sample_set(int Ns, double* prc, double* pvs, double* ptree, int* ret);

  void run_sample_set(int Ns, double* prc, double* pvs, int* ret);
}
