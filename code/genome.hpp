#ifndef GENOME_HPP
#define GENOME_HPP

#include "common.hpp"

// classes used to simulate mutations across genome along a phylogenetic tree

const int N_MUT_TYPE = 5;


class mutation {
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

  segment(const int& _chr, const int& _seg_id);

  segment(const segment& _s2);
};



class genome {
public:
  vector<int> chr_lengths;
  vector<vector<segment>> chrs;
  int node_id;

  // map [chr][seg_id][count]
  map<int, map<int,int>> cn_profile;
  map<int, map<int,int>> allele_cn_profile;

  // count types of mutation
  vector<int> nmuts;
  // store mutations
  vector<mutation> mutations;

  // mean of exp size distributions
  double mean_dup_size;
  double mean_del_size;

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
  void write_allele_cn(ogzstream& of);

  vector<int> get_cn_vector();
};


int get_ploidy_at_chr(const genome& g, int c);
int get_max_cn_chr(genome& g, int c);
int get_max_seg_num_chr(genome& g, int c);
int get_max_cn_chr_seg(genome& g, int c, int loc, int len);
int get_max_seg_num_seg(genome& g, int c, int ins_start, int len);
int get_max_cn_genome(genome& g);
int get_num_chr(const genome& g);
void get_available_chr(const genome& g, vector<int>& available_chrs);

#endif
