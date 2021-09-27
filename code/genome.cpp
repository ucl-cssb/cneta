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

  node_id = _g2.node_id;

  cn_profile.clear();
  cn_profile = _g2.cn_profile;

  copy_number::iterator nit1;
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
        vector<segment> chr;
        for(int j = 0; j < _nsegs[i]; ++j){
          chr.push_back(segment(i,j));
          cn_profile[i][j] = 0;
          allele_cn_profile[i + _nchr * ploidy][j] = 0;
          // cout << i + _nchr * ploidy << "\t" << j << endl;
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
    vector<segment> chr;
    for(int j = 0; j < _chr_lens[i]; ++j){
      chr.push_back(segment(i,j));
      cn_profile[i][j] = 0;
      allele_cn_profile[i][j] = 0;
    }
    chrs.push_back(chr);
  }

  for(int i = 0; i < N_MUT_TYPE; ++i) nmuts.push_back(0);
}

// allow different germline ploidy
genome::genome(const vector<int>& _chr_lens, const int& ploidy){
  for(int p = 0; p < ploidy; ++p){
    for(int i = 0; i < _chr_lens.size(); ++i){
      vector<segment> chr;
      for(int j = 0; j < _chr_lens[i]; ++j){
        chr.push_back(segment(i,j));
        cn_profile[i][j] = 0;
        allele_cn_profile[i + _chr_lens.size() * ploidy][j] = 0;
      }
      chrs.push_back(chr);
    }
  }

  for(int i = 0; i < N_MUT_TYPE; ++i) nmuts.push_back(0);
}


// Compute total copy number
void genome::calculate_cn(){
  for(copy_number::iterator nit1 = cn_profile.begin(); nit1 != cn_profile.end(); ++nit1){
    for(map<int,int>::iterator nit2 = nit1->second.begin(); nit2 != nit1->second.end(); ++nit2){
        nit2->second = 0;
    }
  }

  /* // check the zeroing
  copy_number::const_iterator it1;
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
  for(copy_number::iterator nit1 = allele_cn_profile.begin(); nit1 != allele_cn_profile.end(); ++nit1){
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
  cout << "\tMUTATIONS   (" << node_id + 1 << ") duplication, deletion, gain, loss, WGD:";
  for(int i = 0; i < nmuts.size(); ++i) cout << "\t" << nmuts[i];
  cout << endl;
  //cout << "\t";
  //for(int i = 0; i <  mutations.size(); ++i) cout << "\t" << mutations[i].edge_id + 1 << "," << mutations[i].type << "," << mutations[i].btime << "," << mutations[i].gtime;
  //cout << endl;
}


void genome::print_muts(ostream& stream) const{
  vector<string> stype{"duplication", "deletion", "gain", "loss", "WGD"};

  stream << "MUTATIONS   (" << node_id + 1 << ") duplication, deletion, gain, loss, WGD:";
  for(int i = 0; i < nmuts.size(); ++i) stream << "\t" << nmuts[i];
  stream << endl;
  stream << "\teid\ttype\ttime\trelative time (at the branch)\tabsolute time\tchr\tseg" << endl;
  for(int i = 0; i < mutations.size(); ++i){
    stream << "\t" << mutations[i].edge_id + 1 << "\t" << stype[mutations[i].type] << "\t" << mutations[i].btime << "\t" << mutations[i].gtime << "\t" << mutations[i].chr + 1 << "\t" << mutations[i].seg + 1 << endl;
  }
  stream << endl;
}


void genome::print_cn(){
  calculate_cn();

  // cn is stored as [chr][seg_id][count]
  map<int,int>::const_iterator it2;
  cout << "\tCOPY NUMBER (" << node_id + 1 << "):" << endl;
  for(copy_number::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      cout << "\t" << it1->first + 1 << "_" << it2->first;
    }
    cout << endl;
    for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      cout << "\t" << it2->second;
    }
    cout << endl;
  }
}

void genome::print() const{
  cout << "\tGENOME (node): " << node_id + 1 << endl;
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

  for(copy_number::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(map<int,int>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      // sample and chr starts with 1; segment starts with 0
      of << node_id + 1 << "\t" << it1->first + 1 << "\t" << it2->first + 1 << "\t" << it2->second << endl;
    }
  }
}


void genome::write_rcn(ogzstream& of, gsl_rng* r){
  calculate_cn();

  int nwgd = nmuts[N_MUT_TYPE - 1];
  int ploidy = NORM_PLOIDY;
  if(nwgd > 0) ploidy = pow(NORM_PLOIDY, 1 + nwgd);

  for(copy_number::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(map<int,int>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        // sample and chr starts with 1; segment starts with 0
        int rcn = it2->second;
        if(nwgd > 0){
          // normalize CN so that current ploidy is normal
          int incs = pow(NORM_PLOIDY, nwgd);
          // randomly round up when there are WGDs
          rcn = get_rcn_baseline(rcn, incs, r);
        }
        rcn = rcn - NORM_PLOIDY;
        if(rcn < -2) rcn = -2;
        if(rcn > 2) rcn = 2;
        of << node_id + 1 << "\t" << it1->first + 1 << "\t" << it2->first + 1 << "\t" << rcn << endl;
    }
  }
}


void genome::write_rcn_baseline(ogzstream& of, const vector<int>& chr_lengths, gsl_rng* r){
    calculate_allele_cn();

    int nwgd = nmuts[N_MUT_TYPE - 1];
    int ploidy = NORM_PLOIDY;
    if(nwgd > 0) ploidy = pow(NORM_PLOIDY, 1 + nwgd);
    // std::cout << "Number of WGD events is " << nwgd << ", so ploidy is " << ploidy << endl;
    int allele_ploidy = ploidy / 2;

  // cout << "Writing allele-specific copy number " << endl;
  for(int i = 0; i < allele_cn_profile.size() / 2; i++){
     map<int,int> segs1 = allele_cn_profile[i];
     map<int,int> segs2 = allele_cn_profile[i + NUM_CHR];
     // cout << i << "\t" << segs1.size() << "\t" << segs2.size() << endl;
     // Find the id of the last segment since some segment may get lost
     int max_size = chr_lengths[i];
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
        int rcnA = get_rcn_baseline(cn1, allele_ploidy, r);
        int rcnB = get_rcn_baseline(cn2, allele_ploidy, r);
        // int rcn = rcnA + rcnB;
         of << node_id + 1 << "\t" << i + 1 << "\t" << j + 1 << "\t" << rcnA << "\t" << rcnB << endl;
    }
  }
}


void genome::write_allele_cn(ogzstream& of, const vector<int>& chr_lengths){
  int debug = 0;
  calculate_allele_cn();

  // cout << "Writing allele-specific copy number " << endl;
  for(int i = 0; i < allele_cn_profile.size() / 2; i++){
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
         of << node_id + 1 << "\t" << i + 1 << "\t" << j + 1 << "\t" << cn1 << "\t" << cn2 << endl;
    }
  }
}

vector<int> genome::get_cn_vector(){
  calculate_cn();

  vector<int> cns;
  for(copy_number::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
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
        cout << "maximum copy number for chr " << c + 1 << " is " << max_cn << endl;
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
    if(g.chrs[c].empty()){
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
        cout << "maximum copy number for chr " << c + 1 << ", loc " << loc << ", length " << len << " is " << max_state << endl;
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
        if(!g.chrs[i].empty()){
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
        if(!g.chrs[i].empty()){
            available_chrs.push_back(i);
        }
    }
}


void update_wgd_rate(genome& g, double& wgd_rate, int cn_max){
    int max_cn = get_max_cn_genome(g);
    if(max_cn <= 0 || 2 * max_cn > cn_max){
        // cout << "Allowed MAX CN is " << cn_max << endl;
        // cout << "Maximum copy number on genome " << g.node_id << " is " << max_cn << endl;
        wgd_rate = 0.0;
    }
}


void update_chr_gain_rate(genome& g, int i, double& chr_gain_rate, int cn_max){
    double max_cn_c = get_max_cn_chr(g, i);
    double max_snum = get_max_seg_num_chr(g, i);
    if(max_cn_c <= 0 || max_snum >= cn_max){
        chr_gain_rate = 0.0;
    }
}


void update_chr_loss_rate(genome& g, int i, double& chr_loss_rate){
    double max_cn_c = get_max_cn_chr(g, i);
    if(max_cn_c <= 0){
        chr_loss_rate = 0.0;
    }
}


double get_site_rates(const vector<double>& rate_consts, int model, int stateA, int stateB, int cn_max, double& site_dup_rate, double& site_del_rate){
    double rate = 0.0;

    if(model == MK){ // Mk
        // assume duplication rate equals deletion rate
        assert(rate_consts[0] == rate_consts[1]);
        rate = 4 * rate_consts[0] / 5;
    }else if(model == BOUNDT){  // bounded total
        double dup_rate = rate_consts[0];
        double del_rate = rate_consts[1];
        int state = stateA + stateB;

        if(state == 0){
            site_dup_rate = 0.0;
            site_del_rate = 0.0;
        }else if(state == cn_max){
            site_dup_rate = 0.0;
            site_del_rate = 2 * cn_max * del_rate;
        }else{
            site_dup_rate = 2 * state * dup_rate;
            site_del_rate = 2 * state * del_rate;
        }
        rate = site_dup_rate + site_del_rate;
    }else if(model == BOUNDA){   // bounded allele specific
        double dup_rate = rate_consts[0];
        double del_rate = rate_consts[1];

        if(stateA == 0 && stateB == 0){
            site_dup_rate = 0.0;
            site_del_rate = 0.0;
        }else if(stateA + stateB == cn_max){   // only deletion can occur
            site_dup_rate = 0.0;
            if(stateA == 0 || stateB == 0){   // deletion can only occur at one allele
                site_del_rate = del_rate;
            }else{
                site_del_rate = 2 * del_rate;
            }
        }else{
            if(stateA == 0 || stateB == 0){
                site_dup_rate = dup_rate;
                site_del_rate = del_rate;
            }else{
                site_dup_rate = 2 * dup_rate;
                site_del_rate = 2 * del_rate;
                // cout << "assign rate " << site_dup_rate << "\t" << site_del_rate << endl;
            }
        }
        rate = site_dup_rate + site_del_rate;
    }else{

    }

    // cout << stateA << "\t" << stateB << "\t" << site_dup_rate << "\t" << site_del_rate << endl;

    return rate;
}


double get_total_rates_allele_specific(genome& g, vector<double>& site_dup_rates, vector<double>& site_del_rates, vector<double>& chr_gain_rates, vector<double>& chr_loss_rates, vector<double>& type_rates, const vector<double>& rate_consts, int model, int cn_max, int debug){
    assert(rate_consts.size() == 5);
    double chr_gain_rate = rate_consts[2];
    double chr_loss_rate = rate_consts[3];
    double wgd_rate = rate_consts[4];
    double dup_rate_all = 0.0;
    double del_rate_all = 0.0;
    double gain_rate_all = 0.0;
    double loss_rate_all = 0.0;

    g.calculate_cn(); // Find the current state of the genome before computing rates (affected by copy number)
    g.calculate_allele_cn();

    vector<double> rates_chrs;  // mutation rates on each chromosome
    int s = 0;
    // Tranverse all the sites in the genome to get their rates
    for(int i = 0; i < g.cn_profile.size(); ++i){   // For each chromosome
        vector<double> rates_chr_sites; // mutation rates on all sites of one chromosome
        for(int j = 0; j < g.cn_profile[i].size(); ++j){    // For each segment in the chromosome
            // Get the copy number of this segment
            int stateA = g.allele_cn_profile[i][j];
            int stateB = g.allele_cn_profile[i + NUM_CHR][j];
            double site_dup_rate = 0.0;
            double site_del_rate = 0.0;
            double rate = get_site_rates(rate_consts, model, stateA, stateB, cn_max, site_dup_rate, site_del_rate);
            // if(debug){
            //     cout << "Rate on chr " << i + 1 << " seg " << j + 1 << " site " << s++ << " stateA " << stateA << " stateB " << stateB << " is " << rate << endl;
            //     cout << "Site duplication rate " << site_dup_rate << "\tSite deletion rate " << site_del_rate << endl;
            // }

            dup_rate_all += site_dup_rate;
            del_rate_all += site_del_rate;

            site_dup_rates.push_back(site_dup_rate);     // mutation rates on all sites of all chromosomes
            site_del_rates.push_back(site_del_rate);     // mutation rates on all sites of all chromosomes

            rates_chr_sites.push_back(rate);
        }
        // Get the rate on one chromosome
        double rate_c = accumulate(rates_chr_sites.begin(), rates_chr_sites.end(), 0.0);    // rates for duplication/deletion

        update_chr_gain_rate(g, i, chr_gain_rate, cn_max);
        update_chr_loss_rate(g, i, chr_loss_rate);
        rate_c += chr_loss_rate;
        rate_c += chr_gain_rate;

        if(debug){
            cout << "Rate on chr " << i + 1 << " is " << rate_c << endl;
            cout << "Chromosome gain rate " << chr_gain_rate << "\tChromosome loss rate " << chr_loss_rate << endl;
        }
        rates_chrs.push_back(rate_c);

        chr_gain_rates.push_back(chr_gain_rate);
        chr_loss_rates.push_back(chr_loss_rate);

        gain_rate_all += chr_gain_rate;
        loss_rate_all += chr_loss_rate;
    }
    // Get the total rate of mutation for the current genome (the sum of rates accross sites, chromosomes)
    double rate = accumulate(rates_chrs.begin(), rates_chrs.end(), 0.0);
    // cout << "WGD rate before " << wgd_rate << endl;
    update_wgd_rate(g,  wgd_rate, cn_max);
    // cout << "WGD rate after " << wgd_rate << endl;
    rate += wgd_rate;

    type_rates.push_back(dup_rate_all);
    type_rates.push_back(del_rate_all);
    type_rates.push_back(gain_rate_all);
    type_rates.push_back(loss_rate_all);
    type_rates.push_back(wgd_rate);

    if(debug){
        cout << "The total rate at all sites in the genome: " << rate << endl;
        // cout << site_dup_rates.size() << endl;
        // cout << site_del_rates.size() << endl;
        // cout << chr_gain_rates.size() << endl;
        // cout << chr_loss_rates.size() << endl;
        // cout << type_rates.size() << endl;

        cout << "Site rates:";
        for(int i = 0; i <  site_dup_rates.size(); i++){
            cout << "\t" << site_dup_rates[i] << "\t" << site_del_rates[i];
        }
        cout << endl;

        cout << "Chromosome rates:";
        for(int i = 0; i <  chr_gain_rates.size(); i++){
            cout << "\t" << chr_gain_rates[i] << "\t" << chr_loss_rates[i];
        }
        cout << endl;

        cout << "Rates of different mutational types on genome " << g.node_id;
        for(int i = 0; i < type_rates.size(); i++){
            cout << "\t" << type_rates[i];
        }
        cout << endl;
    }

    return rate;
}

// Randomly choose a haplotype of a chromosome
void select_haplotype(genome& g, int& c, long unsigned (*fp_myrng)(long unsigned)){
    // Find the number of possible haplotype of this chromosome
    int ploidy_at_c = get_ploidy_at_chr(g, c);
    assert(ploidy_at_c > 0);

    int hap_c = c;
    vector<int> haps;   // possible haplotypes
    for(int h = 0; h < ploidy_at_c; h++){
        haps.push_back(h);
    }

    // shuffle(haps.begin(), haps.end(), default_random_engine(seed));
    random_shuffle(haps.begin(), haps.end(), fp_myrng);

    for(int h = 0; h < ploidy_at_c; h++){   // randomly choose a haplotype
        hap_c = c + NUM_CHR * haps[h];
        // cout << "Size of chr " << hap_c << " is " << g.chrs[hap_c].size() << endl;
        if(g.chrs[hap_c].size() > 0)    break;
    }
    c = hap_c;
}


// Randomly choose a haplotype of a chromosome based on available segments
void select_haplotype_by_seg(genome& g, int& c, int seg_id, long unsigned (*fp_myrng)(long unsigned)){
    // Find the number of possible haplotype of this chromosome
    int ploidy_at_c = get_ploidy_at_chr(g, c);
    assert(ploidy_at_c > 0);

    int hap_c = c;
    vector<int> haps;   // possible haplotypes
    for(int h = 0; h < ploidy_at_c; h++){
        haps.push_back(h);
    }

    // shuffle(haps.begin(), haps.end(), default_random_engine(seed));
    random_shuffle(haps.begin(), haps.end(), fp_myrng);

    for(int h = 0; h < ploidy_at_c; h++){   // randomly choose a haplotype
        hap_c = c + NUM_CHR * haps[h];
        // cout << "Size of chr " << hap_c << " is " << g.chrs[hap_c].size() << endl;
        // check if haplotype has the chosen segment
        int has_seg = 0;
        if(g.chrs[hap_c].size() > 0){
            for(int i = 0; i < g.chrs[hap_c].size(); i++){
                if(g.chrs[hap_c][i].seg_id == seg_id){
                    has_seg = 1;
                    break;
                }
            }
            if(has_seg) break;
        }
    }
    c = hap_c;
}

// input: c -- chromosome number (for each haplotype); seg_id -- location on a chromosome (should be intepreted as segment ID)
// if return value > 0,  mutations are generated successfully
int generate_duplication(genome& g, int c, int seg_id, int mean_dup_size, int model, int cn_max, gsl_rng* r, int debug){
    if(g.chrs[c].size() <= 0){
        cout << "All segments on chr " << c + 1 << " has lost" << endl;
        return 0;
    }
    if(debug){
      cout << "GENOME BEFORE duplication " << endl;
      // g.print();
      g.print_cn();
    }
    // mean variant size in bins
    double mdup = mean_dup_size;
    int len = 0;
    if(mdup > 1) len = gsl_ran_exponential(r, mdup);

    if(debug){
        cout << "before duplication, Chr " << c + 1 << " has " << g.chrs[c].size() << " segments" << endl;
        // for(int i = 0; i < g.chrs[c].size(); i++){
        //     cout << "\t" << g.chrs[c][i].seg_id;
        // }
        // cout << endl;
        cout << "dup duplication:" << len << endl;
        cout << "\tSV: duplicating segment, chr, seg_id, len: " << c + 1 << "\t" << seg_id + 1 << "\t" << len + 1 << endl;
        cout << "Previous state: ";
        cout << "\t " << g.cn_profile[c % NUM_CHR][seg_id];
        cout << endl;
    }
    // Find the number of copies
    // int state = g.cn_profile[c % NUM_CHR][g.chrs[c][loc].seg_id];
    int state = g.cn_profile[c % NUM_CHR][seg_id];
    if(state >= cn_max){
        return 0;
    }

    int ins_start = 0;
    double u = runiform(r, 0, 1);
    if(u < 0.5){
        // randomly select a position
        if(debug) cout << "nontandem duplication" << endl;
        ins_start = gsl_rng_uniform_int(r, g.chrs[c].size());
    }else{ // tandem duplication
        // Find the location of seg_id in the genome
        for(int i = 0; i < g.chrs[c].size(); i++){
            if(g.chrs[c][i].seg_id == seg_id){
                ins_start = i;
                break;
            }
        }
    }
    int ins_end = ins_start + len + 1; //assume tandem duplication
    if(debug) cout << "   start " << ins_start << ", end " << ins_end << endl;

    if(model == MK){
        vector<int> possible_states;
        for (int i = state + 1; i <= cn_max; i++){
            possible_states.push_back(i);
        }
        // Randomly pick a possible state
        int ssize = possible_states.size();
        int pstate;
        if(debug){
            cout << "There are " << ssize << " possible state(s) for duplication" << endl;
        }
        if(ssize > 1){
            int k = gsl_rng_uniform_int(r, ssize);
            pstate = possible_states[k];
        }else{
            pstate = state + 1;
        }
        // Ensure there are "state" copy of the region
        // if(debug){
        //     cout << "copy start " << g.chrs[c][loc].seg_id << endl;
        //     cout << "copy end " << g.chrs[c][loc + len].seg_id << endl;
        // }
        if(pstate - state <= 0){
            cout << "no possible states" << '\n';
            return 0;
        }
        // assume tandem duplication
        for(int j = 0; j < pstate - state; j++){
            // g.chrs[c].insert(g.chrs[c].begin() + loc+(len + 1)*(j + 1), g.chrs[c].begin() + loc, g.chrs[c].begin() + loc + len + 1);
            int ins_idx = ins_start + (len + 1) * j;
            g.chrs[c].insert(g.chrs[c].begin() + ins_idx, len + 1, segment(c % NUM_CHR, seg_id));

            if(debug){
                cout << j + 1 << "th dup" << endl;
                int it;
                it = g.chrs[c][ins_idx].seg_id;
                cout << "dup start " << it << endl;
                it = g.chrs[c][ins_idx + len].seg_id;
                cout << "dup end " << it << endl;
            }
        }
    }else if(model == BOUNDT || model == BOUNDA){
        // int max_state = get_max_cn_chr_seg(g, c, seg_id, len);
        // int max_snum = get_max_seg_num_seg(g, c, seg_id, len);
        // // cout << "cn_max is " << cn_max << endl;
        // if(max_snum > cn_max){
        //     return 0;
        // }
        // g.chrs[c].insert(g.chrs[c].begin()+ins_end, g.chrs[c].begin() + loc, g.chrs[c].begin() + loc + len + 1);
        g.chrs[c].insert(g.chrs[c].begin() + ins_start, len + 1, segment(c % NUM_CHR, seg_id));

        if(debug){
            cout << "End position " << ins_start+len + 1 << endl;
            cout << "Insertion position " << ins_start << endl;
            cout << "Current state: ";
            for(int i = 0; i < len + 1; i++){
                cout << "\t " << g.cn_profile[c % NUM_CHR][g.chrs[c][ins_start + i].seg_id];
            }
            cout << endl;
            set<double> segs;
            for(int i = 0; i < g.chrs[c].size(); i++){
                segs.insert(g.chrs[c][i].seg_id);
            }
            cout << "After duplication, Chr " << c + 1 << " has " << g.chrs[c].size() << " segments (unique: " << segs.size() << ")" << endl;
        }
    }else{
        cerr << "Model not supported!" << endl;
        exit(EXIT_FAILURE);
    }

    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
        cout << "GENOME AFTER duplication " << endl;
        // g.print();
        g.print_cn();
    }
    return 1;
}


// Assume the chosen segment is in the specified chromosome
int generate_deletion(genome& g, int c, int seg_id, int mean_del_size, int model, int cn_max, gsl_rng* r, int debug){
    if(g.chrs[c].size() <= 0){
        cout << "All segments on chr " << c + 1 << " has lost" << endl;
        return 0;
    }
    if(debug){
      cout << "GENOME BEFORE deletion " << endl;
      // g.print();
      g.print_cn();
    }

    double mdel = mean_del_size;
    int len = 0;
    if(mdel > 1)   len = (int) gsl_ran_exponential(r, mdel);

    if(debug){
      cout << "before deletion, Chr " << c + 1 << " has " << g.chrs[c].size() << " segments" << endl;
      // for(int i = 0; i < g.chrs[c].size(); i++){
      //     cout << "\t" << g.chrs[c][i].seg_id;
      // }
      // cout << endl;
      cout << "deletion length:" << len << endl;
    }

    // Find the regions to delete
    int del_start = 0;
    int del_msize = 0;
    for(int i = 0; i < g.chrs[c].size(); i++){
        if(g.chrs[c][i].seg_id == seg_id){
            if(g.chrs[c][i-1].seg_id != seg_id){
                del_start = i;
                // cout << i << "\t" << g.chrs[c][i].seg_id << "\t" << seg_id + 1 << endl;
                del_msize = 0;
            }
            del_msize++;
        }
    }
    // if(len + del_start >= g.chrs[c].size()){
    if(del_msize < len){
        cout << "Not enough size to delete " << len << ", delete " << del_msize << endl;
        len = del_msize;
        // return 0;
    }
    int del_end = del_start + len + 1;
    if(debug) cout << "deletion start " << del_start << ", end " << del_end << endl;

    if(model == MK){
        // int state = g.cn_profile[c % NUM_CHR][g.chrs[c][loc].seg_id];
        int state = g.cn_profile[c % NUM_CHR][seg_id];
        // Find the number of segments in other haplotypes
        int hap_state = 0;
        int ploidy_at_c = get_ploidy_at_chr(g, c);
        // cout << "Ploidy at " << c << " is " << ploidy_at_c << endl;
        for(int j = 0; j < ploidy_at_c; j++){
            if(j == c / NUM_CHR) continue;
            int hap = c % NUM_CHR + j * NUM_CHR;
            // cout << "Alternative haplotype " << hap << endl;
            for(int i = 0; i < g.chrs[hap].size();i++){
                if(g.chrs[hap][i].seg_id == seg_id)
                {
                    hap_state++;
                }
            }
        }
        if(debug){
            cout << "Copies in other haplotypes: " << hap_state << endl;
            cout << "\tSV: deleting segment, chrs, seg_id, len: " << c << "\t" << seg_id + 1 << "\t" << len + 1 << endl;
        }

        vector<int> possible_states;
        int pstate;
        for (int i = 0; i <= state - hap_state - 1; i++){
            possible_states.push_back(i);
        }
        // Randomly pick a possible state
        int ssize = possible_states.size();
        if(debug){
            cout << "There are " << ssize << " possible state(s) for deletion" << endl;
        }
        if(ssize <= 0){
            if(debug) cout << "Impossible to do deletion on Chr " << c << endl;
            return 0;
        }else{
            if(ssize > 1){
                int k = gsl_rng_uniform_int(r, ssize);
                pstate = possible_states[k];
            }else{
                pstate = possible_states[0];
            }
            // Ensure there are "state" copy of the region
            int start = g.chrs[c][del_start].seg_id;
            int end = g.chrs[c][del_start + len].seg_id;
            if(debug){
                cout << "del start " << start << endl;
                cout << "del end " << end << endl;
            }

            if(pstate < 0){
                return 0;
            }

            for(int j = 0; j <= pstate; j++){
                // erase all the segment with seg_id in the specified region
                for(int i = 0; i < g.chrs[c].size();i++){
                    int seg = g.chrs[c][i].seg_id;
                    if(seg >= start && seg <= end){
                        g.chrs[c].erase(g.chrs[c].begin() + i);
                        if(debug){
                            cout << "deleting " << seg << " at position " << i << endl;
                        }
                        break;
                    }
                }
            }
        }
    }else if(model == BOUNDT || model == BOUNDA){
        // Find the minimum copy number of the region to delete when the interval to delete has different segments
        // int min_state = g.cn_profile[c % NUM_CHR][g.chrs[c][loc].seg_id];
        // for (int i = loc; i < loc + len + 1; i++){
        //     int state = g.cn_profile[c % NUM_CHR][g.chrs[c][i].seg_id];
        //     if(state < min_state){
        //         min_state = state;
        //     }
        // }
        // if(min_state <= 0){
        //     return 0;
        // }
        g.chrs[c].erase(g.chrs[c].begin() + del_start, g.chrs[c].begin() + del_end);
    }else{
        cerr << "Model not supported!" << endl;
        exit(EXIT_FAILURE);
    }

    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
        cout << "Current state: " << g.cn_profile[c % NUM_CHR][seg_id] << endl;
        set<double> segs;
        for(int i = 0; i < g.chrs[c].size(); i++){
            segs.insert(g.chrs[c][i].seg_id);
        }
        cout << "After deletion, Chr " << c + 1 << " has " << g.chrs[c].size() << " segments (unique: " << segs.size() << ")" << endl;

        cout << "GENOME AFTER deletion " << endl;
        // g.print();
        g.print_cn();
    }
    return 1;
}


int generate_chr_gain(genome& g, int c, int debug){
    // vector<int> available_chrs;
    // get_available_chr(g, available_chrs);
    // if(available_chrs.size()<=0){
    //     // cout << "No available chromosomes!" << endl;
    //     return 0;
    // }
    // int ci = gsl_rng_uniform_int(r, available_chrs.size());
    // int c = available_chrs[ci];
    // int max_cn_c = get_max_cn_chr(g, c);
    // Some segments with maximum copy number may have several copies and hence have copy number larger than specified value when having a chromosome gain (with overlapping events)
    // int max_snum = get_max_seg_num_chr(g, c);
    // if(g.chrs[c].size() <= 0 || max_snum > cn_max){
    //     return 0;
    // }

    // Make sure that the maximum copy number is smaller than the specified value
    int orig_num_chr = get_num_chr(g);
    int orig_size = g.chrs.size();
    assert(orig_size % NUM_CHR == 0);
    // cout << "Size of chrs " << orig_size << endl;
    // assert(orig_size % NUM_CHR == 0);
    int is_in = 0;
    int new_chr = (orig_size / NUM_CHR) * NUM_CHR + c % NUM_CHR; // If c is not in the genome
    // To make sure the chromosome IDs match in different set of copies, add all the haploid chromosomes with others being empty
    for(int i = 0; i < orig_size % NUM_CHR; i++){  // If c is in the genome, but has lost all segments
      int idx = c % NUM_CHR + i * NUM_CHR;
      if(g.chrs[idx].size() <= 0){
          new_chr = idx;
          is_in = 1;
          break;
      }
    }
    if(!is_in){
        for(int i = orig_size; i < orig_size+NUM_CHR; i++){ // Fill the position for other chromosomes before the one to duplicate so that the chromosome ID follows the pattern
          vector<segment> e;
          g.chrs.insert(g.chrs.end(), e);
          // cout << "Size of chr " <<  i << " is " << g.chrs[i].size() << endl;
        }
    }

    g.chrs[new_chr].insert(g.chrs[new_chr].begin(), g.chrs[c].begin(), g.chrs[c].end());
    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
      cout << "Chromosome gain in Chromosome ID " << c << endl;
      cout << "ID for gained Chromosome " << new_chr << endl;
      cout << "Copies to make in chr " << c << endl;
      for(int i = 0; i < g.chrs[c].size(); i++){
          cout << "\t" << g.chrs[c][i].chr << "," << g.chrs[c][i].seg_id;
      }
      cout << endl;
      cout << "There are " << orig_num_chr << " non-empty chromosome IDs before chr gain" << endl;
      cout << "There are " << get_num_chr(g) << " non-empty chromosome IDs now" << endl;
      g.print_cn();
    }

    return 1;
}


int generate_chr_loss(genome& g, int c, int debug){
    // vector<int> available_chrs;
    // get_available_chr(g, available_chrs);
    // if(available_chrs.size()<=0){
    //     // cout << "No available chromosomes!" << endl;
    //     return 0;
    // }
    // int ci = gsl_rng_uniform_int(r, available_chrs.size());
    // int c = available_chrs[ci];
    // int c = gsl_rng_uniform_int(r, g.chrs.size());
    int orig_num_chr = get_num_chr(g);
    // Delete the segments in the chromosome, but keep the chromosome ID so that it can be remapped by NUM_CHR
    g.chrs[c].erase(g.chrs[c].begin(), g.chrs[c].end());
    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
        cout << "Chromosome loss in " << c + 1 << endl;
        cout << "There are " << orig_num_chr << " non-empty chromosomes before chr loss" << endl;
        cout << "There are " << get_num_chr(g) << " non-empty chromosomes now" << endl;
        g.print_cn();
    }
    return 1;
}


int generate_wgd(genome& g, int cn_max, int debug){
    int max_cn_g = get_max_cn_genome(g);
    if(g.chrs.size() <= 0 || 2 * max_cn_g > cn_max){
        return 0;
    }

    int orig_num_chr = get_num_chr(g);
    int orig_size = g.chrs.size();
    // g.chrs.insert(g.chrs.end(), g.chrs.begin(), g.chrs.end());   // cause duplicated insertions on Mac
    vector<segment> s;
    g.chrs.insert(g.chrs.end(), orig_size, s);
    // replicate the segments for each chromosome
    for(int i = 0; i < orig_size; i++){
         // cout << "Insert chr " << i << " for chr " << i + orig_size << endl;
         g.chrs[i + orig_size].insert(g.chrs[i + orig_size].begin(), g.chrs[i].begin(), g.chrs[i].end());
    }
    g.calculate_cn();   // compute copy number after each mutation event
    g.calculate_allele_cn();

    if(debug){
        cout << "Whole genome doubling" << endl;
        cout << "There are " << orig_num_chr << " non-empty chromosomes before WGD" << endl;
        cout << "There are " << get_num_chr(g) << " non-empty chromosomes now" << endl;
        g.print_cn();
        g.print();
    }
    return 1;
}


int get_rcn_baseline(int cn, int baseline, gsl_rng* r){
    int rcn = cn / baseline;   // quotient

    // remainder is always 0 when ploidy = 1
    // randomly round up when remainder is not 0
    if(cn % baseline != 0){
        float rand = runiform(r, 0, 1);
        float val = cn / float(baseline);
        // check it is really float
        // cout << val << endl;
        if(rand < 0.5){
            rcn = ceil(val);
        }else{            
            rcn = floor(val);
            if(rcn == 0) rcn = 1;           
        }
    }

    return rcn;
}
