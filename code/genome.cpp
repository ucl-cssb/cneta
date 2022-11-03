#include "genome.hpp"



/****************** mutation *********************/
mutation::mutation(const int& _eid, const int& _type, const double& _btime, const double& _gtime, int _chr, int _seg):
type(_type), btime(_btime), gtime(_gtime), edge_id(_eid), chr(_chr), seg(_seg){
}

void mutation::print() const{
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
genome::genome(): node_id(-1), num_site(0){};


genome::genome(const int& _nchr, const int& _nseg): node_id(-1){
  int k = 0; // make segment ID unique across the whole genome
  num_site = 0;
  for(int i = 0; i < _nchr; ++i){
    vector<segment> chr;
    chr.resize(_nseg);
    for(int j = 0; j < _nseg; ++j){
      chr[j] = segment(i, k++);
      num_site++;
    }
    //fill_n (std::back_inserter(chrs), 1, chr);
    chrs.push_back(chr);
  }

  initialize_cnp();
  for(int i = 0; i < N_MUT_TYPE; ++i) nmuts.push_back(0);
}


genome::genome(const int& _nchr, const vector<int>& _nsegs, const int& ploidy): node_id(-1){ 
  for(int p = 0; p < ploidy; ++p){
      int k = 0;     
      for(int i = 0; i < _nchr; ++i){
        vector<segment> chr;
        chr.resize(_nsegs[i]);
        for(int j = 0; j < _nsegs[i]; ++j){
          chr[j] = segment(i + _nchr * p, k++);
        }
        chrs.push_back(chr);
      }
  }

  initialize_cnp();
  num_site = accumulate(_nsegs.begin(), _nsegs.end(), 0);
  for(int i = 0; i < N_MUT_TYPE; ++i) nmuts.push_back(0);
}


genome::genome(const vector<int>& _chr_lens, const int& ploidy): node_id(-1){ 
  int _nchr = _chr_lens.size();
  for(int p = 0; p < ploidy; ++p){
    int k = 0;
    for(int i = 0; i < _nchr; ++i){
      vector<segment> chr;
      chr.resize(_chr_lens[i]);
      for(int j = 0; j < _chr_lens[i]; ++j){
        chr[j] = segment(i + _nchr * p, k++);
      }
      chrs.push_back(chr);
    }
  }

  initialize_cnp();
  num_site = accumulate(_chr_lens.begin(), _chr_lens.end(), 0);
  for(int i = 0; i < N_MUT_TYPE; ++i) nmuts.push_back(0);
}


genome::genome(const genome& _g2) {
  chrs = _g2.chrs;
  node_id = _g2.node_id;
  num_site = _g2.num_site;

  cn_profile = _g2.cn_profile;
  allele_cn_profile = _g2.allele_cn_profile;

  nmuts = _g2.nmuts;
  mutations = _g2.mutations;
}


// TODO: to change CNP key from chr to haplotype to allow duplication across chromosomes
void genome::initialize_cnp(){
   int total_chr = NUM_CHR * NORM_PLOIDY;

   for(int i = 0; i < chrs.size(); ++i){
    for(int j = 0; j < chrs[i].size(); ++j){
      segment s = chrs[i][j];
    //   cout << "#####:" << "\t" << s.chr << "\t" << s.seg_id << endl;
      // s.chr is from 0 to NUM_CHR
      cn_profile[s.chr % NUM_CHR][s.seg_id]++;
      allele_cn_profile[s.chr % total_chr][s.seg_id]++;
    }
  }  

  if(2 * cn_profile.size() != allele_cn_profile.size()){
      cout << "CNP size not match! " << cn_profile.size() << ", " << allele_cn_profile.size() << endl;
      exit(EXIT_FAILURE);
  } 

  int size = 0;
  for(auto cnp: cn_profile){
    size += cnp.second.size();
  }
  int asize = 0;
  for(auto cnp: allele_cn_profile){
    asize += cnp.second.size();
  }
  assert(asize == size * 2);
}


int genome::reset_cn(copy_number& cnp){
  int size = 0;
  for(copy_number::iterator it1 = cnp.begin(); it1 != cnp.end(); ++it1){
    for(map<int, int>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
         it2->second = 0;
         size++;
    }
  }
//   cout << "size of CNP " << size << endl;
  return size;
}


// Compute total copy number
void genome::calculate_cn(){
  int size = reset_cn(cn_profile);
  if(size != num_site){
      cout << "Wrong CNP!" << endl;
      for(int i = 0; i < NUM_CHR; i++){
        if(cn_profile[i].size() != BIN_SIZE[i]){
            cout << "chr " << i + 1 << " (should have " << BIN_SIZE[i] << " sites) has " << cn_profile[i].size() << " sites:" << endl;  
            for(auto cnp: cn_profile[i]){
                cout << "\t" << cnp.first + 1;
            } 
            cout << endl;          
        }         
      }
      exit(EXIT_FAILURE);
  }

  // a completely lost segment will disappear from "chrs[i]"
  for(int i = 0; i < chrs.size(); ++i){
    for(int j = 0; j < chrs[i].size(); ++j){
      segment s = chrs[i][j];
      cn_profile[s.chr % NUM_CHR][s.seg_id]++;
    }
  }
}


// Compute haplotype-specific copy number
void genome::calculate_allele_cn(){
  int size = reset_cn(allele_cn_profile);
  if(size != num_site * 2){
      cout << "Wrong haplotype-specific CNP!" << endl;
      print_cnp(allele_cn_profile);
      for(auto cnp: allele_cn_profile){
          cout << cnp.second.size() << endl;
      }      
  }
  
  // cout << "Computing haplotype-specific copy number " << endl;
  int total_chr = NUM_CHR * NORM_PLOIDY;
  // i is from 0 to total_chr - 1 given that chr gain/loss or WGD will only lead to empty chr vector
  for(int i = 0; i < chrs.size(); ++i){
    for(int j = 0; j < chrs[i].size(); ++j){
      segment s = chrs[i][j];
      allele_cn_profile[s.chr % total_chr][s.seg_id]++;
    }
  }
}


// converting absolute copy number relative to the baseline ploidy
int genome::get_rcn_baseline(int cn, int baseline, gsl_rng* r, int random_round){
    int rcn = cn / baseline;   // quotient

    // remainder is always 0 when ploidy = 1
    // randomly round up when remainder is not 0
    if(cn % baseline != 0){
        double rand = runiform(r, 0, 1);
        double val = cn / double(baseline);       
        // cout << val << endl;    // check it is really float
        if(random_round){
            if(rand < 0.5){
                rcn = ceil(val);
            }else{
                rcn = floor(val);
                if(rcn == 0) rcn = 1;    
            }
        }else{
            rcn = round(val);
        }
    }

    return rcn;
}


void genome::print_muts(ostream& stream) const{
  stream << "MUTATIONS   (" << node_id + 1 << ") " << MUT_TYPES[0] << ", "  << MUT_TYPES[1] << ", " << MUT_TYPES[2] << ", "  << MUT_TYPES[3] << ", "  << MUT_TYPES[4] << ":";

  for(int i = 0; i < nmuts.size(); ++i) stream << "\t" << nmuts[i];
  stream << endl;

  stream << "\teid\ttype\ttime\trelative time (at the branch)\tabsolute time\tchr\tseg" << endl;
  for(int i = 0; i < mutations.size(); ++i){
    stream << "\t" << mutations[i].edge_id + 1 << "\t" << MUT_TYPES[mutations[i].type] << "\t" << mutations[i].btime << "\t" << mutations[i].gtime << "\t" << mutations[i].chr + 1 << "\t" << mutations[i].seg + 1 << endl;
  }
  stream << endl;
}


void genome::print_cn() const{
  // cn is stored as [chr][seg_id][count]
  map<int, int>::const_iterator it2;
  cout << "\tCOPY NUMBER (" << node_id + 1 << "):" << endl;
  for(copy_number::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      cout << "\t" << it1->first + 1 << "_" << it2->first + 1;
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
    cout << "\t\tchr / nsegs: " << i + 1 << "\t" << chrs[i].size() << endl;
    cout << "\t";
    for(int j = 0; j < chrs[i].size(); ++j){
        cout << "\t" << chrs[i][j].chr + 1 << "_" << chrs[i][j].seg_id + 1;
    }
    cout << endl;
  }
}


void genome::write(ogzstream& of){
  for(copy_number::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(map<int, int>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      // sample and chr starts with 1; segment starts with 0      
      of << node_id + 1 << "\t" << it1->first + 1 << "\t" << it2->first + 1 << "\t" << it2->second << endl;
    }
  }
}


 // normalize CN so that current ploidy is normal
void genome::write_rcn(ogzstream& of, gsl_rng* r, int use_nwgd, int random_round){
  int nwgd = nmuts[N_MUT_TYPE - 1];
  int baseline = 0;
  if(use_nwgd){
    assert(nwgd >= 0);
    // ploidy = pow(NORM_PLOIDY, 1 + nwgd);
    baseline = pow(NORM_PLOIDY, nwgd);
  }else{
    // use rounded mean copy number as base line
    vector<int> cns = get_cn_vector();
    baseline = round(get_mean_cn(cns) / 2);
  }

  for(copy_number::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(map<int, int>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        // sample and chr starts with 1; segment starts with 0
        int rcn = it2->second;
        rcn = get_rcn_baseline(rcn, baseline, r, random_round);
        rcn = rcn - NORM_PLOIDY;
        if(rcn < -2) rcn = -2;
        if(rcn > 2) rcn = 2;
        of << node_id + 1 << "\t" << it1->first + 1 << "\t" << it2->first + 1 << "\t" << rcn << endl;
    }
  }
}


void genome::write_allele_rcn(ogzstream& of, gsl_rng* r, int use_nwgd, int random_round){
    int nwgd = nmuts[N_MUT_TYPE - 1];
    int alleleA_baseline = 0;
    int alleleB_baseline = alleleA_baseline;

    // get the baseline
    if(use_nwgd){
        assert(nwgd >= 0);
        // std::cout << "Number of WGD events is " << nwgd << ", so ploidy is " << ploidy << endl;
        alleleA_baseline = pow(NORM_PLOIDY, nwgd);
        alleleB_baseline = alleleA_baseline;
    }else{
        // use rounded mean copy number as base line
        vector<int> cnsA;
        vector<int> cnsB;
        get_allele_cn_vector(cnsA, cnsB);
        alleleA_baseline = round(get_mean_cn(cnsA));
        alleleB_baseline = round(get_mean_cn(cnsB));
    }

    // cout << "Writing haplotype-specific copy number " << endl;
    for(int i = 0; i < allele_cn_profile.size() / 2; i++){
        map<int, int> segs1 = allele_cn_profile[i];
        map<int, int> segs2 = allele_cn_profile[i + NUM_CHR];
        // cout << i << "\t" << segs1.size() << "\t" << segs2.size() << endl;
        // Find the id of the last segment since some segment may get lost
        assert(segs1.size() == segs2.size());
        for(auto seg: segs1){
            int seg_id = seg.first;
            int cn1 = seg.second;        
            int cn2 = segs2[seg_id];
            int rcnA = get_rcn_baseline(cn1, alleleA_baseline, r, random_round);
            int rcnB = get_rcn_baseline(cn2, alleleB_baseline, r, random_round);
            // int rcn = rcnA + rcnB;
            of << node_id + 1 << "\t" << i + 1 << "\t" << seg_id + 1 << "\t" << rcnA << "\t" << rcnB << endl;
        }
    }
}


void genome::write_allele_cn(ogzstream& of){
  assert(allele_cn_profile.size() == 2 * NUM_CHR);
  // cout << "Writing haplotype-specific copy number " << endl;
  for(int i = 0; i < NUM_CHR; i++){
     map<int, int> segs1 = allele_cn_profile[i];
     map<int, int> segs2 = allele_cn_profile[i + NUM_CHR];
     // cout << i << "\t" << segs1.size() << "\t" << segs2.size() << endl;
     // Find the id of the last segment since some segment may get lost
     assert(segs1.size() == segs2.size());
    for(auto seg: segs1){
        int seg_id = seg.first;
        int cn1 = seg.second;        
        int cn2 = segs2[seg_id];
        of << node_id + 1 << "\t" << i + 1 << "\t" << seg_id + 1 << "\t" << cn1 << "\t" << cn2 << endl;
    }
  }
}


vector<int> genome::get_cn_vector(){
  vector<int> cns;
  for(copy_number::const_iterator it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(map<int, int>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      cns.push_back(it2->second);
    }
  }
  return cns;
}


void genome::get_allele_cn_vector(vector<int>& cnsA, vector<int>& cnsB){
    assert(cnsA.size() == 0 && cnsB.size() == 0);

    for(int i = 0; i < allele_cn_profile.size() / 2; i++){
        map<int, int> segs1 = allele_cn_profile[i];
        map<int, int> segs2 = allele_cn_profile[i + NUM_CHR];
        assert(segs1.size() == segs2.size());
        for(auto seg: segs1){
            int seg_id = seg.first;
            int cn1 = seg.second;        
            int cn2 = segs2[seg_id];
            cnsA.push_back(cn1);
            cnsB.push_back(cn2);
        }
    }
}


// get average total copy number by sample
double genome::get_mean_cn(const vector<int>& cns){
    int tcn = accumulate(cns.begin(), cns.end(), 0);
    int nele = cns.size();
    double avg_tcn = (double) tcn / nele;
    // int avg_tcn = round(favg_tcn);
    // cout << "Average total copy number " << favg_tcn << "\t" << avg_tcn << endl;
    
    return avg_tcn;
}



void print_cnp(copy_number& cn_profile){
  copy_number::const_iterator it1;
  map<int, int>::const_iterator it2;
  cout << "\tCOPY NUMBER PROFILE:" << endl;
  for(it1 = cn_profile.begin(); it1 != cn_profile.end(); ++it1){
    for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      cout << " " << it1->first + 1 << "_" << it2->first + 1 << ":" << it2->second;
    }
     cout << endl;
  } 
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


int get_max_cn_dup_seg(genome& g, int c, int start, int end){
    // Some segments may not have the largest copy number for now. But the current copy number plus the multiplicity will make the new copy number larger than specified after duplication.
    map<int, int> seg_cn;     // record seg_id and current copy number
    vector<int> segs;     

    for(int i = start; i < end; i++){
        int chr = g.chrs[c][i].chr;
        int sid = g.chrs[c][i].seg_id;
        int cn = g.cn_profile[chr % NUM_CHR][sid];
        if(seg_cn.find(sid) == seg_cn.end()){
            seg_cn.insert(pair<int, int>(sid, cn));
        }else{  // same segment should have the same CN in cn_profile
            assert(seg_cn[sid] == cn);
        }        
        segs.push_back(sid);
    }
    
    map<int, int> dup;  // record seg_id and current multiplicity
    for_each(segs.begin(), segs.end(), [&dup](int val){ dup[val]++; } );

    int mcn_dup = 0;   // maximum total CN for the segments after duplication
    for(auto p : dup){
      int sid = p.first;
      int multiplicity = p.second;
      int total_cn = p.second + seg_cn[p.first];
      if(total_cn > mcn_dup){
          mcn_dup = total_cn;
      }
      // if(debug){
      //     cout << "Number of segments: " << p.first << ' ' << p.second << ' ' << total_cn << endl;
      // }
    }

    return mcn_dup;
}


int get_max_cn_chr_seg(genome& g, int c, int loc, int len){
    if(g.chrs[c].empty()){
        return 0;
    }

    int debug = 0;
    int max_cn = g.cn_profile[g.chrs[c][loc].chr % NUM_CHR][g.chrs[c][loc].seg_id];

    // Check the state of all the regions
    for(int i = loc + 1; i < loc + len; i++){
        int cn = g.cn_profile[g.chrs[c][i].chr % NUM_CHR][g.chrs[c][i].seg_id];
        // cout << "copy number for " << i << " is " << state << endl;
        if(cn > max_cn){
            max_cn = cn;
        }
    }

    if(debug){
        cout << "maximum copy number for chr " << c + 1 << ", loc " << loc << ", length " << len << " is " << max_cn << endl;
    }

    return max_cn;
}


int get_max_cn_chr(genome& g, int c){
    int debug = 0;
    int max_cn = 0;

    for(int i = 0; i < g.chrs[c].size(); i++){
        int cn = g.cn_profile[g.chrs[c][i].chr % NUM_CHR][g.chrs[c][i].seg_id];
        if(cn > max_cn){
            max_cn = cn;
        }
    }

    if(debug){
        cout << "maximum copy number for chr " << c + 1 << " is " << max_cn << endl;
    }

    return max_cn;
}


int get_max_cn_genome(genome& g){
    int debug = 0;
    int max_cn = 0;

    for(int chr = 0; chr < NUM_CHR; chr++){
        int cn = get_max_cn_chr(g, chr);
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
int get_num_available_chr(const genome& g){
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
    // max_cn <= 0 -- all segments have been lost
    if(max_cn <= 0 || 2 * max_cn > cn_max){
        // cout << "Allowed MAX CN is " << cn_max << endl;
        // cout << "Maximum copy number on genome " << g.node_id << " is " << max_cn << endl;
        wgd_rate = 0.0;
    }
}


void update_chr_gain_rate(genome& g, int c, double& chr_gain_rate, int cn_max){
    double max_cn_c = get_max_cn_chr(g, c);
    // double max_snum = get_max_seg_num_chr(g, c);
    // max_cn_c <= 0 -- all segments on chr c have been lost
    if(max_cn_c <= 0 || 2 * max_cn_c > cn_max){
        chr_gain_rate = 0.0;
    }
}


void update_chr_loss_rate(genome& g, int c, double& chr_loss_rate){
    double max_cn_c = get_max_cn_chr(g, c);
    if(max_cn_c <= 0){
        chr_loss_rate = 0.0;
    }
}


double get_site_rates(const vector<double>& rate_consts, int model, int stateA, int stateB, int cn_max, double& site_dup_rate, double& site_del_rate){
    double rate = 0.0;

    if(model == MK){
        // assume duplication rate equals deletion rate
        assert(rate_consts[0] == rate_consts[1]);
        rate = (cn_max - 1) * rate_consts[0];
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
    assert(g.cn_profile.size() == NUM_CHR);

    double chr_gain_rate = rate_consts[2];
    double chr_loss_rate = rate_consts[3];
    double wgd_rate = rate_consts[4];

    double dup_rate_all = 0.0;
    double del_rate_all = 0.0;
    double gain_rate_all = 0.0;
    double loss_rate_all = 0.0;

    if(debug){
        cout << "CNP before updating rates" << endl;
        g.print_cn();
    }

    vector<double> rates_chrs(g.cn_profile.size(), 0.0);  // mutation rates on each chromosome
    int i = 0;  // count sites on all chromosomes
    // Tranverse all the sites in the genome to get their rates
    for(int chr = 0; chr < NUM_CHR; ++chr){   // For each chromosome
        vector<double> rates_chr_sites(g.cn_profile[chr].size(), 0.0); // mutation rates on all sites of one chromosome 
        int j = 0;  // count sites on each chr       
        for(auto seg : g.cn_profile[chr]){    // For each segment in the chromosome
            int seg_id = seg.first;  // should be in sequencial order
            // Get the copy number of this segment
            int stateA = g.allele_cn_profile[chr][seg_id];
            int stateB = g.allele_cn_profile[chr + NUM_CHR][seg_id];
            double site_dup_rate = 0.0;
            double site_del_rate = 0.0;
            double rate = get_site_rates(rate_consts, model, stateA, stateB, cn_max, site_dup_rate, site_del_rate);
            // if(debug){
            //     cout << "Rate on chr " << chr + 1 << " seg " << seg_id + 1 << " site on genome " << i << " site on chr " << j << " stateA " << stateA << " stateB " << stateB << " rate " << rate << endl;
            //     cout << "Site duplication rate " << site_dup_rate << "\tSite deletion rate " << site_del_rate << endl;
            // }

            dup_rate_all += site_dup_rate;
            del_rate_all += site_del_rate;

            site_dup_rates[i] = site_dup_rate;     
            site_del_rates[i] = site_del_rate;     
            rates_chr_sites[j] = rate;

            i++;
            j++;
        }

        // Get the rate on one chromosome
        double rate_c = accumulate(rates_chr_sites.begin(), rates_chr_sites.end(), 0.0);    // rates for duplication/deletion
        update_chr_gain_rate(g, chr, chr_gain_rate, cn_max);
        update_chr_loss_rate(g, chr, chr_loss_rate);
        rate_c += chr_loss_rate;
        rate_c += chr_gain_rate;
        rates_chrs[chr] = rate_c;

        if(debug){
            cout << "Rate on chr " << chr + 1 << " is " << rate_c << endl;
            cout << "Chromosome gain rate " << chr_gain_rate << "\tChromosome loss rate " << chr_loss_rate << endl;
        }

        chr_gain_rates[chr] = chr_gain_rate;
        chr_loss_rates[chr] = chr_loss_rate;
        gain_rate_all += chr_gain_rate;
        loss_rate_all += chr_loss_rate;
    } 
    // assert(i == 4401);
    // Get the total rate of mutation for the current genome (the sum of rates accross sites, chromosomes)
    double rate = accumulate(rates_chrs.begin(), rates_chrs.end(), 0.0);
    // cout << "WGD rate before " << wgd_rate << endl;
    update_wgd_rate(g,  wgd_rate, cn_max);
    // cout << "WGD rate after " << wgd_rate << endl;
    rate += wgd_rate;

    type_rates.insert(type_rates.end(), {dup_rate_all, del_rate_all, gain_rate_all, loss_rate_all, wgd_rate});

    if(debug){
        cout << "The total rate at all sites in the genome: " << rate << endl;

        // cout << "Non-zero Site rates:";
        // for(int i = 0; i <  site_dup_rates.size(); i++){
        //     if(site_dup_rates[i] > 0 ||  site_del_rates[i] > 0)
        //         cout << "\t" << i << "_" << site_dup_rates[i] << "_" << site_del_rates[i];
        // }
        // cout << endl;

        // cout << "Chromosome rates:";
        // for(int i = 0; i <  chr_gain_rates.size(); i++){
        //     cout << "\t" << chr_gain_rates[i] << "\t" << chr_loss_rates[i];
        // }
        // cout << endl;

        // cout << "Rates of different mutational types on genome " << g.node_id;
        // for(int i = 0; i < type_rates.size(); i++){
        //     cout << "\t" << type_rates[i];
        // }
        // cout << endl;

        cout << "CNP after updating rates" << endl;
        g.print_cn();
    }

    return rate;
}


void select_haplotype(genome& g, int& c, long unsigned (*fp_myrng)(long unsigned)){
    // Find the number of possible haplotype of this chromosome
    int ploidy_at_c = get_ploidy_at_chr(g, c);
    assert(ploidy_at_c > 0);

    int hap_c = c;
    vector<int> haps(ploidy_at_c);   // possible haplotypes
    iota(haps.begin(), haps.end(), 0);
    // shuffle(haps.begin(), haps.end(), default_random_engine(seed));
    random_shuffle(haps.begin(), haps.end(), fp_myrng);

    for(int h = 0; h < ploidy_at_c; h++){   // randomly choose a haplotype
        hap_c = c + NUM_CHR * haps[h];
        // cout << "Size of chr " << hap_c << " is " << g.chrs[hap_c].size() << endl;
        if(g.chrs[hap_c].size() > 0)    break;
    }
    c = hap_c;
}


void select_haplotype_by_seg(genome& g, int& c, int seg_id, long unsigned (*fp_myrng)(long unsigned)){
    // Find the number of possible haplotype of this chromosome
    int ploidy_at_c = get_ploidy_at_chr(g, c);
    assert(ploidy_at_c > 0);

    int hap_c = c;
    vector<int> haps(ploidy_at_c);   // possible haplotypes
    iota(haps.begin(), haps.end(), 0);
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
// allow duplication of size > 1 to across segs with different IDs
int generate_duplication(genome& g, int c, int seg_id, int mean_dup_size, int model, int cn_max, gsl_rng* r, int debug){
    if(g.chrs[c].size() <= 0){
        if(debug) cout << "All sites on chr " << c + 1 << " has lost, cannot duplicate!" << endl;
        return 0;
    }

    // mean variant size in bins
    int len = 1;
    if(mean_dup_size > 1) len = gsl_ran_exponential(r, mean_dup_size);
    // Find the original position of the segment to duplicate
    int orig_start = 0;
    // Find the location of first seg_id in the genome
    for(int i = 0; i < g.chrs[c].size(); i++){
        if(g.chrs[c][i].seg_id == seg_id){
            orig_start = i;
            break;
        }
    }
    int orig_end = orig_start + len;       
    if(orig_end > g.chrs[c].size()){
        orig_end = g.chrs[c].size();
    }

    if(debug){
        cout << "before duplication, Chr " << c + 1 << " has " << g.chrs[c].size() << " segments" << endl;
        for(int i = 0; i < g.chrs[c].size(); i++){
            cout << "\t" << g.chrs[c][i].seg_id;
        }
        cout << endl;
        cout << "\tSV: duplicating site, chr, start, size: " << c + 1 << "\t" << seg_id + 1 << "\t" << len << endl;
        cout << "Previous state of selected sites: ";
        for(int i = orig_start; i < orig_end; i++){
            cout << "\t " << g.cn_profile[g.chrs[c][i].chr % NUM_CHR][g.chrs[c][i].seg_id];
        }
        cout << endl;       
    }
   
    // Find the number of copies. There may be >1 copies of the same segment
    int mcn_dup = get_max_cn_dup_seg(g, c, orig_start, orig_end);
    if(mcn_dup > cn_max){
        if(debug) cout << "Reaching maximum CN on chr " << c + 1 << " site [" << orig_start << "," << orig_end << "), cannot duplicate!" << endl;
        return 0;
    }

    // Find the position to insert the new copy
    int ins_start = 0;  
    int ins_chr = c; 
    double u = runiform(r, 0, 1);
    if(u < 0.5){
        // randomly select a position at the same chr?
        if(debug) cout << "Interspersed duplication" << endl;
        while(g.chrs[ins_chr].size() <= 0){
            ins_chr = gsl_rng_uniform_int(r, g.chrs.size());
            ins_start = gsl_rng_uniform_int(r, g.chrs[ins_chr].size());
        }
    }else{ 
        if(debug) cout << "Tandem duplication" << endl;
        ins_start = orig_end;
    }

    if(debug){
        cout << "   inserting following sites at chr : " << ins_chr + 1 << " from position " << ins_start << " to " << ins_start + len << endl;
        for(int i = orig_start; i < orig_end; i++){
            cout << g.chrs[c][i].chr + 1 << "\t" << g.chrs[c][i].seg_id + 1 << endl;
        }        
    }

    if(model == MK){
        int state = g.cn_profile[c % NUM_CHR][seg_id];
        if(state >= cn_max){
            return 0;
        }       
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
            int ins_idx = ins_start + (len) * j;
            g.chrs[c].insert(g.chrs[c].begin() + ins_idx, g.chrs[c].begin() + orig_start, g.chrs[c].begin() + orig_end);
            // g.chrs[c].insert(g.chrs[c].begin() + ins_idx, len, segment(c % NUM_CHR, seg_id));

            if(debug){
                cout << j + 1 << "th dup" << endl;
                int it = g.chrs[c][ins_idx].seg_id;
                cout << "dup start " << it << endl;
                it = g.chrs[c][ins_idx + len].seg_id;
                cout << "dup end " << it << endl;
            }
        }
    }else if(model == BOUNDT || model == BOUNDA){
        g.chrs[ins_chr].insert(g.chrs[ins_chr].begin() + ins_start, g.chrs[c].begin() + orig_start, g.chrs[c].begin() + orig_end);
        // g.chrs[c].insert(g.chrs[c].begin() + ins_start, len, segment(c % NUM_CHR, seg_id));
    }else{
        cerr << "Model not supported!" << endl;
        exit(EXIT_FAILURE);
    }

    if(debug){
        set<double> segs;
        for(int i = 0; i < g.chrs[ins_chr].size(); i++){
            segs.insert(g.chrs[ins_chr][i].seg_id);
        }       
        cout << "after duplication, Chr " << ins_chr + 1 << " has " << g.chrs[ins_chr].size() << " segments (unique: " << segs.size() << ")" << endl;  
        for(int i = 0; i < g.chrs[ins_chr].size(); i++){
            cout << "\t" << g.chrs[ins_chr][i].seg_id;
        }
        cout << endl;                        
    }

    return 1;
}


// Assume the chosen segment is in the specified chromosome
int generate_deletion(genome& g, int c, int seg_id, int mean_del_size, int model, int cn_max, gsl_rng* r, int debug){
    if(g.chrs[c].size() <= 0){              // make sure a deletion somewhere is possible
        if(debug) cout << "All sites on chr " << c + 1 << " has lost, cannot delete!" << endl;
        return 0;
    }

    int len = 1;
    if(mean_del_size > 1)  len = (int) gsl_ran_exponential(r, mean_del_size);

    if(debug){
      cout << "before deletion, Chr " << c + 1 << " has " << g.chrs[c].size() << " segments" << endl;
      for(int i = 0; i < g.chrs[c].size(); i++){
          cout << "\t" << g.chrs[c][i].seg_id;
      }
      cout << endl;
      cout << "deletion length:" << len << endl;
    }

    // Find the regions to delete
    int del_start = 0;
    // int del_msize = 0;
    for(int i = 0; i < g.chrs[c].size(); i++){
        if(g.chrs[c][i].seg_id == seg_id){    // only delete regions/bins with the same ID
            // if(g.chrs[c][i-1].seg_id != seg_id){
            //     del_start = i;
            //     // cout << i << "\t" << g.chrs[c][i].seg_id << "\t" << seg_id + 1 << endl;
            //     del_msize = 0;
            // }
            // del_msize++;
            del_start = i;
            break;
        }
    }
    // if(len + del_start >= g.chrs[c].size()){
    // if(del_msize < len){
    //     cout << "Not enough size to delete " << len << ", delete " << del_msize << endl;
    //     len = del_msize;
    //     // return 0;
    // }
    int del_end = del_start + len;
    if(del_end > g.chrs[c].size()){
      del_end = g.chrs[c].size();
    }
    if(debug) cout << "deletion at " << c + 1 << ", start " << del_start << ", end " << del_end << endl;

    if(model == MK){
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
                if(g.chrs[hap][i].seg_id == seg_id){
                    hap_state++;
                }
            }
        }
        if(debug){
            cout << "Copies in other haplotypes: " << hap_state << endl;
            cout << "\tSV: deleting segment, chrs, seg_id, len: " << c << "\t" << seg_id + 1 << "\t" << len << endl;
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
        g.chrs[c].erase(g.chrs[c].begin() + del_start, g.chrs[c].begin() + del_end);
    }else{
        cerr << "Model not supported!" << endl;
        exit(EXIT_FAILURE);
    }

    return 1;
}


int generate_chr_gain(genome& g, int c, int cn_max, int debug){
    double max_cn_c = get_max_cn_chr(g, c);
    if(max_cn_c <= 0 || 2 * max_cn_c > cn_max){
        return 0;
    }

    // Make sure that the maximum copy number is smaller than the specified value
    int orig_size = g.chrs.size();
    assert(orig_size % NUM_CHR == 0);
    // cout << "Size of chrs " << orig_size << endl;
    
    int new_chr = orig_size + c % NUM_CHR;   // insert new chr at the end by default
    int is_in = 0;
    // If some haplotype of c is in the genome, but has lost all segments, just fill the segments
    for(int i = 0; i < orig_size / NUM_CHR; i++){  
      int idx = c % NUM_CHR + i * NUM_CHR;
      if(g.chrs[idx].size() <= 0){
          new_chr = idx;
          is_in = 1;
          break;
      }
    }
    // If c is not in the genome, add all the haploid chromosomes with others being empty to make sure the chromosome IDs match in different set of copies
    if(!is_in){
        for(int i = orig_size; i < orig_size + NUM_CHR; i++){ // Fill the position for other chromosomes before the one to duplicate so that the chromosome ID follows the pattern
          vector<segment> e;
          g.chrs.insert(g.chrs.end(), e);
          // cout << "Size of chr " <<  i << " is " << g.chrs[i].size() << endl;
        }
    }

    g.chrs[new_chr].insert(g.chrs[new_chr].begin(), g.chrs[c].begin(), g.chrs[c].end());

    if(debug){
      cout << "ID for gained Chromosome " << new_chr << endl;
      cout << "Copies to make for gain in chr " << c + 1 << endl;
      for(int i = 0; i < g.chrs[c].size(); i++){
          cout << "\t" << g.chrs[c][i].chr + 1 << "," << g.chrs[c][i].seg_id + 1;
      }
      cout << endl;     
    }

    return 1;
}


int generate_chr_loss(genome& g, int c, int debug){
  double max_cn_c = get_max_cn_chr(g, c);
  if(max_cn_c <= 0){
      return 0;
  }

  // Delete the segments in the chromosome, but keep the chromosome ID so that it can be remapped by NUM_CHR
  g.chrs[c].erase(g.chrs[c].begin(), g.chrs[c].end());

  return 1;
}


int generate_wgd(genome& g, int cn_max, int debug){
    int max_cn_g = get_max_cn_genome(g);
    if(g.chrs.size() <= 0 || 2 * max_cn_g > cn_max){
        return 0;
    }

    int orig_size = g.chrs.size();
    // g.chrs.insert(g.chrs.end(), g.chrs.begin(), g.chrs.end());   // cause duplicated insertions on Mac
    vector<segment> s;
    g.chrs.insert(g.chrs.end(), orig_size, s);
    // replicate the segments for each chromosome
    for(int i = 0; i < orig_size; i++){
         // cout << "Insert chr " << i << " for chr " << i + orig_size << endl;
         g.chrs[i + orig_size].insert(g.chrs[i + orig_size].begin(), g.chrs[i].begin(), g.chrs[i].end());
    }

    return 1;
}


int site2chr(int site, int& chr, const vector<int>& chr_lengths, int debug){
    int seg = -1;
    if(debug) cout << "Site to convert is " << site << endl;
    if(site >= 0 && site < chr_lengths[0]){
        chr = 0;
        seg = site;
    }else{
        if(debug) cout << "There are " << chr_lengths.size() << " chromosomes" << endl;
        for(int i = 0; i < chr_lengths.size(); i++){
            int sum1 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 1, 0);
            int sum2 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 2, 0);
            // if(debug){
            //     cout << "Size for chr " << i + 1 << " is " << chr_lengths[i] << endl;
            //     cout << "sum until chromosome " << i + 1 << " is " << sum1 << endl;
            //     cout << "sum until chromosome " << i + 2 << " is " << sum2 << endl;
            // }
            if(site >= sum1 && site < sum2){
                chr = i + 1;
                seg = site - sum1;
                break;
            }
        }
    }
    assert(seg < chr_lengths[chr]);

    return seg;
}
