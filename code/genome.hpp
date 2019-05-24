// genome.hpp

#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "gzstream.h"

using namespace std;

//enum mutation_type {duplication, deletion, chr_loss, chr_gain, wgd };
int n_mut_type = 5;
const int CN_MAX = 4;

class mutation {
public:
  int type;
  double btime;
  double gtime;
  int edge_id;

  mutation(const int& _eid, const int& _type, const double& _btime, const double& _gtime){
    type = _type;
    btime = _btime;
    gtime = _gtime;
    edge_id = _eid;
  }

  void print() const{
    //cout << "\t" << type << "\t" << chr << "\t" << loc << "\t" << size << fixed << setprecision(3) << "\t" << time << endl;
    cout << "\t" << type << fixed << setprecision(3) << "\t" << btime << "\t" << gtime << endl;
  }

};

class segment{
public:
  int chr;
  int seg_id;

  segment(const int& _chr, const int& _seg_id){
    chr    = _chr;
    seg_id = _seg_id;
  }

  segment(const segment& _s2){
    chr    = _s2.chr;
    seg_id = _s2.seg_id;
  }
};

class genome {
public:
  vector<int> chr_lengths;
  vector<vector<segment>> chrs;
  int node_id;

  // map [chr][seg_id][count]
  map<int, map<int,int> > cn_profile;

  // count types of mutation
  vector<int> nmuts;
  // store mutations
  vector<mutation> mutations;

  // mean of exp size distributions
  double mean_dup_size;
  double mean_del_size;

  genome(){};

  genome(const genome& _g2) {
    chrs.clear();
    chrs.insert(chrs.end(), _g2.chrs.begin(), _g2.chrs.end() );

    chr_lengths.clear();
    chr_lengths.insert(chr_lengths.end(), _g2.chr_lengths.begin(), _g2.chr_lengths.end() );
    node_id = _g2.node_id;

    cn_profile.clear();
    cn_profile = _g2.cn_profile;

    map<int, map<int,int> >::iterator nit1;
    map<int,int>::iterator nit2;
    for(nit1 = cn_profile.begin(); nit1 != cn_profile.end(); ++nit1){
      for(nit2 = nit1->second.begin(); nit2 != nit1->second.end(); ++nit2){
          nit2->second = 0;
      }
    }

    nmuts.clear();
    for(int i=0; i<_g2.nmuts.size(); ++i) nmuts.push_back(_g2.nmuts[i]);

    mutations.clear();
    for(int i=0; i<_g2.mutations.size(); ++i) mutations.push_back(_g2.mutations[i]);

    mean_dup_size = _g2.mean_dup_size;
    mean_del_size = _g2.mean_del_size;
  }

  // create a genome with nchr x nseg
  genome(const int& _nchr, const int& _nseg){
    for(int i=0; i<_nchr; ++i){
      chr_lengths.push_back(_nseg);

      vector<segment> chr;
      for(int j=0; j < _nseg; ++j){
        chr.push_back( segment(i,j) );
        cn_profile[i][j] = 0;
      }
      //fill_n (std::back_inserter(chrs), 1, chr);
      chrs.push_back(chr);
    }

    for(int i=0; i<n_mut_type; ++i) nmuts.push_back(0);
  }

  // create a genome with varying chromosome sizes
  genome(const vector<int>& _chr_lens){
    for(int i=0; i<_chr_lens.size(); ++i){
      chr_lengths.push_back(_chr_lens[i]);

      vector<segment> chr;
      for(int j=0; j < _chr_lens[i]; ++j){
        chr.push_back( segment(i,j) );
        cn_profile[i][j] = 0;
      }
      chrs.push_back(chr);
    }

    for(int i=0; i<n_mut_type; ++i) nmuts.push_back(0);
  }

  // allow different germline ploidy
  genome(const vector<int>& _chr_lens, const int& ploidy){

    for( int p=0; p<ploidy; ++p){
      for(int i=0; i<_chr_lens.size(); ++i){
        chr_lengths.push_back(_chr_lens[i]);

        vector<segment> chr;
        for(int j=0; j < _chr_lens[i]; ++j){
          chr.push_back( segment(i,j) );
          cn_profile[i][j] = 0;
        }
        chrs.push_back(chr);
      }
    }

    for(int i=0; i<n_mut_type; ++i) nmuts.push_back(0);
  }


  void calculate_cn(){

    map<int, map<int,int> >::iterator nit1;
    map<int,int>::iterator nit2;
    for(nit1 = cn_profile.begin(); nit1 != cn_profile.end(); ++nit1){
      for(nit2 = nit1->second.begin(); nit2 != nit1->second.end(); ++nit2){
          nit2->second = 0;
      }
    }

    /* // check the zeroing
    map<int, map<int,int> >::const_iterator it1;
    map<int,int>::const_iterator it2;
    cout << "\tCOPY NUMBER (zeros):" << endl;
    for(it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
      for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        cout << "\t" << it1->first << "_" << it2->first;
      }
      cout << endl;
      for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        cout << "\t" << it2->second;
      }
      cout << endl;
    }
    */

    // here chrs stores the segments that are currently in the genome
    // the cn profile is calculated wrt to the original unmutated genome
    for(int i=0; i<chrs.size(); ++i){
      for(int j=0; j<chrs[i].size(); ++j){
        segment s = chrs[i][j];
        //cout << "#####:" << "\t" << s.chr << "\t" << s.seg_id << endl;
        cn_profile[s.chr][s.seg_id]++;
      }
    }
  }

  void print_muts(){
    cout << "\tMUTATIONS   (" << node_id+1 << ")  dup, deln, crl, crg, wgd:";
    for(int i=0; i< nmuts.size(); ++i) cout << "\t" << nmuts[i];
    cout << endl;
    //cout << "\t";
    //for(int i=0; i< mutations.size(); ++i) cout << "\t" << mutations[i].edge_id+1 << "," << mutations[i].type << "," << mutations[i].btime << "," << mutations[i].gtime;
    //cout << endl;
  }

  void print_muts(ostream& stream){
    vector<string> stype;
    stype.push_back("dup");
    stype.push_back("del");
    stype.push_back("crl");
    stype.push_back("crg");
    stype.push_back("wgd");

    stream << "MUTATIONS   (" << node_id+1 << ")  dup, deln, crl, crg, wgd:";
    for(int i=0; i< nmuts.size(); ++i) stream << "\t" << nmuts[i];
    stream << endl;
    stream << "\teid\ttype\ttime" << endl;
    for(int i=0; i< mutations.size(); ++i){
      stream << "\t" << mutations[i].edge_id+1 << "\t" << stype[mutations[i].type] << "\t" << mutations[i].btime << "\t" << mutations[i].gtime << endl;
    }
      stream << endl;
  }

  void print_cn(){
    calculate_cn();

    // cn is stored as [chr][seg_id][count]

    map<int, map<int,int> >::const_iterator it1;
    map<int,int>::const_iterator it2;
    cout << "\tCOPY NUMBER (" << node_id+1 << "):" << endl;
    for(it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
      for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        cout << "\t" << it1->first << "_" << it2->first;
      }
      cout << endl;
      for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        cout << "\t" << it2->second;
      }
      cout << endl;
    }
  }

  void print(){
    cout << "\tGENOME (node): " << node_id+1 << endl;
    for(int i=0; i<chrs.size(); ++i){
      cout << "\t\tchr / nsegs: " << i << "\t" << chrs[i].size() << endl;
      cout << "\t";
      for(int j=0; j<chrs[i].size(); ++j){
    cout << "\t" << chrs[i][j].chr << "_" << chrs[i][j].seg_id;
      }
      cout << endl;
    }
  }

  void write(ofstream& of){
    calculate_cn();

    map<int, map<int,int> >::const_iterator it1;
    map<int,int>::const_iterator it2;
    for(it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
      for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        of << node_id+1 << "\t" << it1->first+1 << "\t" << it2->first << "\t" << it2->second << endl;
      }
    }
  }

  void write(ogzstream& of){
    calculate_cn();

    map<int, map<int,int> >::const_iterator it1;
    map<int,int>::const_iterator it2;
    for(it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
      for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        of << node_id+1 << "\t" << it1->first+1 << "\t" << it2->first << "\t" << it2->second << endl;
      }
    }
  }

  vector<int> get_cn_vector(){
    calculate_cn();
    vector<int> cns;

    map<int, map<int,int> >::const_iterator it1;
    map<int,int>::const_iterator it2;
    for(it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
      for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        cns.push_back(it2->second);
      }
    }
    return cns;
  }
};
