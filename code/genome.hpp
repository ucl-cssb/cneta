#ifndef GENOME_HPP
#define GENOME_HPP

#include "common.hpp"
#include "stats.hpp"

// classes used to simulate mutations across genome along a phylogenetic tree

const int N_MUT_TYPE = 5;


class mutation{
public:
  int type;
  double btime;
  double gtime;
  int edge_id;
  int chr;  // chromosome on which the mutation occurs
  int seg; // segment on which the mutation occurs

  mutation(const int& _eid, const int& _type, const double& _btime, const double& _gtime, int _chr, int _seg);

  void print() const;
};


class segment{
public:
  int chr;
  int seg_id;

  segment();
  segment(const int& _chr, const int& _seg_id);

  segment(const segment& _s2);
};



class genome{
public:
  vector<vector<segment>> chrs;
  int node_id;    // associated tree node

  // map [chr][seg_id][count]
  copy_number cn_profile;
  copy_number allele_cn_profile;

  // count number of mutation for each type (5 for now)
  vector<int> nmuts;
  // store all the mutations generated
  vector<mutation> mutations;

  genome();
  genome(const genome& _g2);
  // create a genome with nchr x nseg
  genome(const int& _nchr, const int& _nseg);
  // create a genome with nchr x nseg_i, allowing different germline ploidy
  // allow different number of segs on different chromosmes
  genome(const int& _nchr, const vector<int>& _nsegs, const int& ploidy);
  // create a genome with varying chromosome sizes
  genome(const vector<int>& _chr_lens);
  // allow different germline ploidy
  genome(const vector<int>& _chr_lens, const int& ploidy);

  // Compute total copy number
  void calculate_cn();
 // Compute allele-specific copy number
  void calculate_allele_cn();

  void print_muts() const;
  void print_muts(ostream& stream) const;
  void print() const;

  void print_cn();
  void write(ogzstream& of);
  void write_rcn(ogzstream& of, gsl_rng* r);
  void write_rcn_baseline(ogzstream& of, const vector<int>& chr_lengths, gsl_rng* r);
  void write_allele_cn(ogzstream& of, const vector<int>& chr_lengths);

  vector<int> get_cn_vector();   // get the CNP of the genome
};


int get_ploidy_at_chr(const genome& g, int c);
int get_max_cn_chr(genome& g, int c);
int get_max_seg_num_chr(genome& g, int c);
int get_max_cn_chr_seg(genome& g, int c, int loc, int len);
int get_max_seg_num_seg(genome& g, int c, int ins_start, int len);
int get_max_cn_genome(genome& g);
int get_num_chr(const genome& g);
void get_available_chr(const genome& g, vector<int>& available_chrs);

// The WGD rate of the genome changes as the copy numbers change, due to the upper limit on total copy number
// WGD is impossible when it is approaching cn_max
void update_wgd_rate(genome& g, double& wgd_rate, int cn_max);
// The chromosome gain rate of the genome changes as the copy numbers change, due to the upper limit on total copy number
void update_chr_gain_rate(genome& g, int i, double& chr_gain_rate, int cn_max);
void update_chr_loss_rate(genome& g, int i, double& chr_loss_rate);


// rates: The rate of all events (segment duplication, segment deletion) for finite sites model
double get_site_rates(const vector<double>& rate_consts, int model, int stateA, int stateB, int cn_max, double& site_dup_rate, double& site_del_rate);
// rates: The total mutation rate at each site in the genome, used for randomly selecting a site
// The site is pre-specified at the beginning, denoted by each segment ID
// site_*_rates: duplication/deletion probabilities for all sites on the genome, used for picking a site for duplication/deletion
// chr_*_rates: chromosome gain/loss probabilities for all sites on the genome, used for picking a site for chromosome gain/loss
// type_rates: used for selecting different mutational types
double get_total_rates_allele_specific(genome& g, vector<double>& site_dup_rates, vector<double>& site_del_rates, vector<double>& chr_gain_rates, vector<double>& chr_loss_rates, vector<double>& type_rates, const vector<double>& rate_consts, int model, int cn_max, int debug);


void select_haplotype(genome& g, int& c, long unsigned (*fp_myrng)(long unsigned));
void select_haplotype_by_seg(genome& g, int& c, int seg_id, long unsigned (*fp_myrng)(long unsigned));

int generate_duplication(genome& g, int c, int seg_id, int mean_dup_size, int model, int cn_max, gsl_rng* r, int debug);
int generate_deletion(genome& g, int c, int seg_id, int mean_del_size, int model, int cn_max, gsl_rng* r, int debug);
int generate_chr_gain(genome& g, int c, int debug);
int generate_chr_loss(genome& g, int c, int debug);
int generate_wgd(genome& g, int cn_max, int debug);

// Get relative copy number from absolute copy number based on the baseline, using baseline strategy to pre-process data for data with WGD
int get_rcn_baseline(int cn, int baseline, gsl_rng* r);

#endif
