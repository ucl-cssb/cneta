#ifndef GENOME_HPP
#define GENOME_HPP

#include "common.hpp"
#include "stats.hpp"

// classes used to simulate mutations across genome along a phylogenetic tree

const int BIN_SIZE[] = {367, 385, 335, 316, 299, 277, 251, 243, 184, 210, 215, 213, 166, 150, 134, 118, 121, 127, 79, 106, 51, 54};
const int N_MUT_TYPE = 5;
const vector<string> MUT_TYPES {"duplication", "deletion", "gain", "loss", "WGD"};

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


// a region in the genome, which may be a fixed-size bin or a segment of unknown size
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
  /*
   chrs stores the segments that are currently in the genome, assuming sequential order;
   Once a chr is appended to the vector, its position is assumed to be fixed, according to NUM_CHR and PLOIDY, although there may be no segments on this chr due to deletion or loss.
   This order assumption is important to compute CNP with index of chrs.
   */
  vector<vector<segment>> chrs;

  int node_id;    // associated tree node
  int num_site;   // a genome is divideded into a number of sites for computing CNP

  // the cn profile is calculated wrt to the original unmutated genome
  // map [chr][seg_id] = copy_number
  // in a Map, every key is mapped with default value zero when the map is declared.
  copy_number cn_profile;
  copy_number allele_cn_profile;  // another allele is stored NUM_CHR after

  // count number of mutation for each type (5 for now)
  vector<int> nmuts;
  // store all the mutations generated
  vector<mutation> mutations;

  /* 
  can specify the germline in a number of different ways: genome germline(1,10); genome germline(chr_lengths);
  A chromosome (haplotype) is composed of nseg sites.
  Each segment in the chromosome is specified by chrID (from 0 to nchr - 1) and segID (from 0 to nseg - 1)
  Different haplotypes of the same chr are distinguished by position. 
  e.g. Hapolotype 0 is from 0 to nchr - 1; Hapolotype 1 is from nchr to 2* nchr - 1
  */
  genome();
  // create a genome with nchr x nseg
  genome(const int& _nchr, const int& _nseg);
  // create a genome with nchr x nseg_i, allowing different germline ploidy and different number of segs on different chromosomes
  genome(const int& _nchr, const vector<int>& _nsegs, const int& ploidy);
  // create a genome with varying chromosome sizes and allow different germline ploidy
  genome(const vector<int>& _chr_lens, const int& ploidy = 1);
  genome(const genome& _g2);

  void initialize_cnp();    // initialize CNP at the beginning so that all positions are recorded
  int reset_cn(copy_number& cnp);   // reset counter to 0 for recounting after CN changes without losing positions
  void calculate_cn();   // Compute total copy number
  void calculate_allele_cn();   // Compute haplotype-specific copy number

  void print_muts(ostream& stream) const;
  void print() const;
  void print_cn() const;
  void write(ogzstream& of);
  void write_rcn(ogzstream& of, gsl_rng* r);
  // write haplotype-specific relative CNs
  void write_allele_rcn(ogzstream& of, gsl_rng* r);
  void write_allele_cn(ogzstream& of);

  vector<int> get_cn_vector();   // get the CNP of the genome
};

void print_cnp(copy_number& cn_profile);

int get_ploidy_at_chr(const genome& g, int c);


/***** functions to get maximum CNs at different scales (used in updating mutation rate) **********/
// Check the maximum copy number of the interval to duplicate, [start, end), when the interval has multiple different segments. Assume the copy number profile is updated (for checking computation of max CN on chr c seg [start, end))
int get_max_cn_dup_seg(genome& g, int c, int start, int end);
// Find the current maximum copy number for one interval (loc, loc + len) on a chromosome, assuming the interval contains several segments (not used)
int get_max_cn_chr_seg(genome& g, int c, int loc, int len);
// Find the current maximum copy number for one chromosome
int get_max_cn_chr(genome& g, int c);
// Find the current maximum copy number for the whole genome
int get_max_cn_genome(genome& g);

int get_num_available_chr(const genome& g);
void get_available_chr(const genome& g, vector<int>& available_chrs);

// The WGD rate of the genome changes as the copy numbers change, due to the upper limit on total copy number
// WGD is impossible when it is approaching cn_max
void update_wgd_rate(genome& g, double& wgd_rate, int cn_max);
// The chromosome gain rate of the genome changes as the copy numbers change, due to the upper limit on total copy number
void update_chr_gain_rate(genome& g, int c, double& chr_gain_rate, int cn_max);
void update_chr_loss_rate(genome& g, int c, double& chr_loss_rate);

/*
Return the rate, q_i = sum(q_{ij}), where i!=j, of segment duplication and deletion at a site
stateA, stateB: CNs for allele A and B
*/
double get_site_rates(const vector<double>& rate_consts, int model, int stateA, int stateB, int cn_max, double& site_dup_rate, double& site_del_rate);
/*
Return the total mutation rate across all sites in the genome, used for randomly selecting a site
The site is pre-specified at the beginning, denoted by each segment ID
site_*_rates: duplication/deletion probabilities for all sites on the genome, used for picking a site for duplication/deletion
chr_*_rates: chromosome gain/loss probabilities for all sites on the genome, used for picking a site for chromosome gain/loss
type_rates: used for selecting different mutational types
*/
double get_total_rates_allele_specific(genome& g, vector<double>& site_dup_rates, vector<double>& site_del_rates, vector<double>& chr_gain_rates, vector<double>& chr_loss_rates, vector<double>& type_rates, const vector<double>& rate_consts, int model, int cn_max, int debug);

/*
 Note: the sites are counted for all the haplotypes. 
 The positions on different haplotypes are counted as one site.
 The starting point of sites is 0
 Used when segment ID starting from 0 for each chromosome
*/
int site2chr(int site, int& chr, const vector<int>& chr_lengths, int debug);
// Randomly choose a haplotype of a chromosome, used in generating chromosome gain/loss
void select_haplotype(genome& g, int& c, long unsigned (*fp_myrng)(long unsigned));
// Randomly select a haplotype of a chromosome with the chosen segments, used in generating segment duplication/deletion
void select_haplotype_by_seg(genome& g, int& c, int seg_id, long unsigned (*fp_myrng)(long unsigned));

int generate_duplication(genome& g, int c, int seg_id, int mean_dup_size, int model, int cn_max, gsl_rng* r, int debug);
int generate_deletion(genome& g, int c, int seg_id, int mean_del_size, int model, int cn_max, gsl_rng* r, int debug);
// c is the chromosome ID, which reflects its haplotype, e.g. c = 23 means it is chr 1 haplotype 1
int generate_chr_gain(genome& g, int c, int cn_max, int debug);
int generate_chr_loss(genome& g, int c, int debug);
int generate_wgd(genome& g, int cn_max, int debug);

// Get relative copy number from absolute copy number based on the baseline, using baseline strategy to pre-process data for data with WGD
int get_rcn_baseline(int cn, int baseline, gsl_rng* r);


#endif
