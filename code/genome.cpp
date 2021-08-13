#include "genome.hpp"



/****************** mutation *********************/
mutation::mutation(const int& _eid, const int& _type, const double& _btime, const double& _gtime, int _chr, int _seg):
type(_type), btime(_btime), gtime(_gtime), edge_id(_eid), chr(_chr), seg(_seg){
}

void mutation::print() const{
  //cout << "\t" << type << "\t" << chr << "\t" << loc << "\t" << size << fixed << setprecision(3) << "\t" << time << endl;
  cout << "\t" << type << fixed << setprecision(3) << "\t" << btime << "\t" << gtime << "\t" << chr << "\t" << seg << endl;
}


/****************** segment *********************/
segment::segment(){};

segment::segment(const int& _chr, const int& _seg_id):
chr(_chr), seg_id(_seg_id){
  chr    = _chr;
  seg_id = _seg_id;
}

segment::segment(const segment& _s2){
  chr    = _s2.chr;
  seg_id = _s2.seg_id;
}


/****************** genome *********************/
genome::genome(){};

genome::genome(const genome& _g2) {
  chrs.assign(_g2.chrs.begin(), _g2.chrs.end() );

  // chr_lengths.assign(_g2.chr_lengths.begin(), _g2.chr_lengths.end() );
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

  nmuts.assign(_g2.nmuts.begin(), _g2.nmuts.end());
  mutations.assign(_g2.mutations.begin(), _g2.mutations.end());
}

// create a genome with nchr x nseg
genome::genome(const int& _nchr, const int& _nseg){
  for(int i = 0; i < _nchr; ++i){
    // chr_lengths.push_back(_nseg);

    vector<segment> chr;
    chr.resize(_nseg);
    for(int j = 0; j < _nseg; ++j){
      chr.push_back(segment(i,j));
      cn_profile[i][j] = 0;
      allele_cn_profile[i][j] = 0;
    }
    //fill_n (std::back_inserter(chrs), 1, chr);
    chrs.push_back(chr);
  }

  for(int i = 0; i < N_MUT_TYPE; ++i) nmuts.push_back(0);
}

// create a genome with nchr x nseg_i, allowing different germline ploidy
// allow different number of segs on different chromosmes
genome::genome(const int& _nchr, const vector<int>& _nsegs, const int& ploidy){
  for(int p = 0; p < ploidy; ++p){
      for(int i = 0; i < _nchr; ++i){
        // chr_lengths.push_back(_nsegs[i]);

        vector<segment> chr;
        for(int j = 0; j < _nsegs[i]; ++j){
          chr.push_back( segment(i,j) );
          cn_profile[i][j] = 0;
          allele_cn_profile[i+_nchr*ploidy][j] = 0;
          cout << i+_nchr*ploidy << "\t" << j << endl;
        }
        //fill_n (std::back_inserter(chrs), 1, chr);
        chrs.push_back(chr);
      }
  }
  assert(allele_cn_profile.size() == ploidy * _nchr);
  for(int i = 0; i < N_MUT_TYPE; ++i) nmuts.push_back(0);
}

// create a genome with varying chromosome sizes
genome::genome(const vector<int>& _chr_lens){
  for(int i = 0; i < _chr_lens.size(); ++i){
    // chr_lengths.push_back(_chr_lens[i]);

    vector<segment> chr;
    for(int j = 0; j < _chr_lens[i]; ++j){
      chr.push_back( segment(i,j) );
      cn_profile[i][j] = 0;
      allele_cn_profile[i][j] = 0;
    }
    chrs.push_back(chr);
  }

  for(int i = 0; i < N_MUT_TYPE; ++i) nmuts.push_back(0);
}

// allow different germline ploidy
genome::genome(const vector<int>& _chr_lens, const int& ploidy){
  for(int p=0; p<ploidy; ++p){
    for(int i = 0; i < _chr_lens.size(); ++i){
      // chr_lengths.push_back(_chr_lens[i]);

      vector<segment> chr;
      for(int j = 0; j < _chr_lens[i]; ++j){
        chr.push_back( segment(i,j) );
        cn_profile[i][j] = 0;
        allele_cn_profile[i+_chr_lens.size()*ploidy][j] = 0;
      }
      chrs.push_back(chr);
    }
  }

  for(int i = 0; i < N_MUT_TYPE; ++i) nmuts.push_back(0);
}


// Compute total copy number
void genome::calculate_cn(){
  for(map<int, map<int,int>>::iterator nit1 = cn_profile.begin(); nit1 != cn_profile.end(); ++nit1){
    for(map<int,int>::iterator nit2 = nit1->second.begin(); nit2 != nit1->second.end(); ++nit2){
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
  for(int i = 0; i < chrs.size(); ++i){
    for(int j = 0; j < chrs[i].size(); ++j){
      segment s = chrs[i][j];
      //cout << "#####:" << "\t" << s.chr << "\t" << s.seg_id << endl;
      cn_profile[s.chr][s.seg_id]++;
    }
  }
}


// Compute allele-specific copy number
void genome::calculate_allele_cn(){
  for(map<int, map<int,int> >::iterator nit1 = allele_cn_profile.begin(); nit1 != allele_cn_profile.end(); ++nit1){
      for(map<int,int>::iterator nit2 = nit1->second.begin(); nit2 != nit1->second.end(); ++nit2){
          nit2->second = 0;
      }
  }

  // cout << "Computing allele-specific copy number " << endl;
  // cout << chrs.size() << endl;
  // here chrs stores the segments that are currently in the genome
  // the cn profile is calculated wrt to the original unmutated genome
  for(int i = 0; i < chrs.size(); ++i){
    for(int j = 0; j < chrs[i].size(); ++j){
      segment s = chrs[i][j];
      // cout << "#####:" << "\t" << i%44 << "\t" << s.seg_id << endl;
      allele_cn_profile[i % (NUM_CHR * NORM_PLOIDY)][s.seg_id]++;
    }
  }
}

void genome::print_muts() const{
  cout << "\tMUTATIONS   (" << node_id+1 << ") duplication, deletion, gain, loss, WGD:";
  for(int i = 0; i < nmuts.size(); ++i) cout << "\t" << nmuts[i];
  cout << endl;
  //cout << "\t";
  //for(int i = 0; i <  mutations.size(); ++i) cout << "\t" << mutations[i].edge_id+1 << "," << mutations[i].type << "," << mutations[i].btime << "," << mutations[i].gtime;
  //cout << endl;
}


void genome::print_muts(ostream& stream) const{
  vector<string> stype{"duplication", "deletion", "gain", "loss", "WGD"};

  stream << "MUTATIONS   (" << node_id+1 << ") duplication, deletion, gain, loss, WGD:";
  for(int i = 0; i <  nmuts.size(); ++i) stream << "\t" << nmuts[i];
  stream << endl;
  stream << "\teid\ttype\ttime" << endl;
  for(int i = 0; i <  mutations.size(); ++i){
    stream << "\t" << mutations[i].edge_id+1 << "\t" << stype[mutations[i].type] << "\t" << mutations[i].btime << "\t" << mutations[i].gtime  << "\t" << mutations[i].chr + 1 << "\t" << mutations[i].seg << endl;
  }
    stream << endl;
}

void genome::print_cn(){
  calculate_cn();

  // cn is stored as [chr][seg_id][count]
  map<int,int>::const_iterator it2;
  cout << "\tCOPY NUMBER (" << node_id+1 << "):" << endl;
  for(map<int, map<int,int> >::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      cout << "\t" << it1->first+1 << "_" << it2->first;
    }
    cout << endl;
    for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      cout << "\t" << it2->second;
    }
    cout << endl;
  }
}

void genome::print() const{
  cout << "\tGENOME (node): " << node_id+1 << endl;
  for(int i = 0; i < chrs.size(); ++i){
    cout << "\t\tchr / nsegs: " << i << "\t" << chrs[i].size() << endl;
    cout << "\t";
    for(int j = 0; j < chrs[i].size(); ++j){
        cout << "\t" << chrs[i][j].chr << "_" << chrs[i][j].seg_id;
    }
    cout << endl;
  }
}


void genome::write(ogzstream& of){
  calculate_cn();

  for(map<int, map<int,int> >::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(map<int,int>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      // sample and chr starts with 1; segment starts with 0
      of << node_id+1 << "\t" << it1->first+1 << "\t" << it2->first << "\t" << it2->second << endl;
    }
  }
}

void genome::write_allele_cn(ogzstream& of, const vector<int>& chr_lengths){
  int debug = 0;
  calculate_allele_cn();

  // cout << "Writing allele-specific copy number " << endl;
  for(int i = 0; i < allele_cn_profile.size()/2; i++){
     map<int,int> segs1 = allele_cn_profile[i];
     map<int,int> segs2 = allele_cn_profile[i + NUM_CHR];
     // cout << i << "\t" << segs1.size() << "\t" << segs2.size() << endl;
     // Find the id of the last segment since some segment may get lost
     int max_size = chr_lengths[i];
     if(debug) cout << i << "\t" << max_size << endl;
     for(int j = 0; j < max_size; j++){
         int cn1, cn2;
         if (segs1.find(j) == segs1.end()) {
             cn1 = 0;
         }else{
             cn1 = segs1[j];
         }
         if (segs2.find(j) == segs2.end()) {
             cn2 = 0;
         }else{
             cn2 = segs2[j];
         }
         of << node_id+1 << "\t" << i+1 << "\t" << j << "\t" << cn1 << "\t" << cn2 << endl;
    }
  }
}

vector<int> genome::get_cn_vector(){
  calculate_cn();

  vector<int> cns;
  for(map<int, map<int,int> >::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(map<int,int>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      cns.push_back(it2->second);
    }
  }
  return cns;
}


// Find the number of a specific chromosome
int get_ploidy_at_chr(const genome& g, int c){
    int ploidy_at_c = 0;
    for(int i = 0; i < g.chrs.size(); i++){
        if(i % NUM_CHR == c){
            ploidy_at_c += 1;
        }
    }
    return ploidy_at_c;
}


// Find the current maximum copy number for one chromosome
int get_max_cn_chr(genome& g, int c){
    int debug = 0;
    g.calculate_cn();
    vector<segment> segs = g.chrs[c];
    int max_cn = 0;
    for(int i = 0; i < segs.size(); i++){
        int cn = g.cn_profile[c%NUM_CHR][g.chrs[c][i].seg_id];
        if(cn > max_cn){
            max_cn = cn;
        }
    }
    if(debug){
        cout << "maximum copy number for chr " << c+1 << " is " << max_cn << endl;
    }
    return max_cn;
}

// Find the maximum copy number of a segment on a specified chromosome
int get_max_seg_num_chr(genome& g, int c){
    g.calculate_cn();
    map<int, int> seg_cn;     // record seg_id and current copy number
    vector<int> segs;     // record seg_id and current multiplicity
    for(int i = 0; i < g.chrs[c].size(); i++){
        int sid = g.chrs[c][i].seg_id;
        int cn = g.cn_profile[c%NUM_CHR][sid];
        seg_cn.insert(pair<int, int>(sid, cn));
        segs.push_back(sid);
    }
    map<int,int> dup;
    for_each(segs.begin(), segs.end(), [&dup]( int val ){ dup[val]++; } );

    int num_mseg = 0;
    for(auto p : dup){
        int total_cn = p.second + seg_cn[p.first];
        // if(debug){
        //     cout << "Number of segments: " << p.first << ' ' << p.second << ' ' << total_cn << endl;
        // }
        if(total_cn > num_mseg){
            num_mseg = total_cn;
        }
    }
    return num_mseg;
}

// Find the current maximum copy number for one interval on a chromosome, assuming the interval contains several segments
int get_max_cn_chr_seg(genome& g, int c, int loc, int len){
    int debug = 0;
    if(g.chrs[c].size() <= 0){
        return 0;
    }
    g.calculate_cn();
    int max_state = g.cn_profile[c%NUM_CHR][g.chrs[c][loc].seg_id];
    // Check the state of all the regions
    for(int i = loc + 1; i < loc + len + 1; i++){
        int state = g.cn_profile[c%NUM_CHR][g.chrs[c][i].seg_id];
        // cout << "copy number for " << i << " is " << state << endl;
        if(state > max_state){
            max_state = state;
        }
    }
    if(debug){
        cout << "maximum copy number for chr " << c+1 << ", loc " << loc << ", length " << len << " is " << max_state << endl;
    }
    return max_state;
}


// Check the maximum copy number of the interval to duplicate when the interval has multiple different segments
// Assume the copy number profile is updated
int get_max_seg_num_seg(genome& g, int c, int ins_start, int len){
    // g.calculate_cn();
    // Some segments may not have the largest copy number for now. But there current copy number plus the multiplicity will make the new copy number larger than specified.
    map<int, int> seg_cn;     // record seg_id and current copy number
    vector<int> segs;     // record seg_id and current multiplicity
    for(int i=ins_start; i < ins_start + len + 1; i++){
        int sid = g.chrs[c][i].seg_id;
        int cn = g.cn_profile[c%NUM_CHR][sid];
        seg_cn.insert(pair<int, int>(sid, cn));
        segs.push_back(sid);
    }
    map<int,int> dup;
    for_each(segs.begin(), segs.end(), [&dup](int val){ dup[val]++; } );

    int num_mseg = 0;
    for(auto p : dup){
        int total_cn = p.second + seg_cn[p.first];
        // if(debug){
        //     cout << "Number of segments: " << p.first << ' ' << p.second << ' ' << total_cn << endl;
        // }
        if(total_cn > num_mseg){
            num_mseg = total_cn;
        }
    }
    return num_mseg;
}


// Find the current maximum copy number for the whole genome
int get_max_cn_genome(genome& g){
    int debug = 0;
    g.calculate_cn();
    int max_cn = 0;
    for(int i = 0; i < g.cn_profile.size(); i++){
        int cn = get_max_cn_chr(g, i);
        if(cn > max_cn){
            max_cn = cn;
        }
    }
    if(debug){
        cout << "maximum copy number for the whole genome is " << max_cn << endl;
    }
    return max_cn;
}


// Find the number of chromosomes with segments
int get_num_chr(const genome& g){
    int n = 0;
    for(int i = 0; i < g.chrs.size(); i++){
        if(g.chrs[i].size() > 0){
            n++;
        }
    }
    return n;
}

// Find the ID of chromosomes with segments
void get_available_chr(const genome& g, vector<int>& available_chrs){
    available_chrs.clear();
    for(int i = 0; i < g.chrs.size(); i++){
        // cout << "size of chr " << i << " is " << g.chrs[i].size() << endl;
        if(g.chrs[i].size() > 0){
            available_chrs.push_back(i);
        }
    }
}
