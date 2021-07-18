#include "parse_cn.hpp"


// collection of functions for reading/writing


// Convert the allele-specific state (ordered by 0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2 2/1	3/0	 0/4 1/3 2/2 3/1	4/0) to total copy number
int state_to_total_cn(int state, int cn_max){
    int cn = 0;
    vector<int> sums;
    // cout << "Sums of state: ";
    for(int i = 1; i <= (cn_max + 1); i++){
        int s = i * (i+1) / 2;
        sums.push_back(s);
        // cout << "\t" << s;
    }
    // cout << endl;
    if(state < sums[0]) return 0;
    int i = 1;
    do{
        if(state >= sums[i-1] && state < sums[i]) return i;
        i++;
    }while(i <= cn_max);
    // cout << state << "\t" <<
    return 0;
}


// Convert the allele-specific state (ordered by 0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2 2/1	3/0	 0/4 1/3 2/2 3/1	4/0) to total copy number
int state_to_allele_cn(int state, int cn_max, int& cnA, int& cnB){
    int cn = 0;
    vector<int> sums;
    // cout << "Sums of state: ";
    for(int i = 1; i <= (cn_max + 1); i++){
        int s = i * (i+1) / 2;
        sums.push_back(s);
        // cout << "\t" << s;
    }
    // cout << endl;
    if(state < sums[0]) return 0;
    int i = 1;
    do{
        if(state >= sums[i-1] && state < sums[i]){
             // total copy number is i;
             int diff = state - sums[i-1];
             cnA = diff;
             cnB = i - cnA;
        }
        i++;
    }while(i <= cn_max);
    // cout << state << "\t" <<
    return 0;
}


// Read the samping time and patient age of each sample
// Ns: number of samples
// age: age at earliest sample
vector<double> read_time_info(const string& filename, const int& Ns, int& age, int debug){
  if(debug) cout << "\tread_time_info" << endl;

  vector<double> t_info;
  vector<int> ages;
  ifstream infile (filename.c_str());
  if(infile.is_open()){
    std::string line;
    while(!getline(infile, line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf;
      stringstream ss(line);
      while(ss >> buf) split.push_back(buf);

      double dt = atof(split[1].c_str());
      //cout << "read dt: " << dt << endl;
      t_info.push_back(dt);

      if(split.size() > 2){
          int a = atoi(split[2].c_str());
          ages.push_back(a);
      }
    }
    if(ages.size() > 0){
        age = *min_element(ages.begin(), ages.end());
    }
  }else{
    std::cerr << "Error: open of time data unsuccessful: " << filename << std::endl;
    exit(1);
  }

  if(t_info.size() != Ns){
    std::cerr << "Error: timing information does not contain Ns entries: " << filename << std::endl;
    exit(1);
  }

  return t_info;
}


// Change the allele specific copy number to the state used in substitution rate matrix
int allele_cn_to_state(int cnA, int cnB){
    int tcn = cnA + cnB;
    int s = 0;

    int nprev = 0;
    // There are i+1 combinations for a total copy number of i
    for(int i = 0; i < tcn; i++){
        nprev += i + 1;
    }
    // cout << nprev << " cases before " << cnA << "," << cnB << endl;
    s = nprev;

    for(int j = 0; j < cnA; j++){
        // cout << j << endl;
        s += 1;
    }
    // cout << "State is " << s << endl;
    return s;
}

// Format of input total copy number: sample, chr, seg, copy_number
// Format of input allele specific copy number: sample, chr, seg, copy_number A, copy_number B -> converted into a specific number by position
vector<vector<vector<int>>> read_cn(const string& filename, int Ns, int &num_total_bins, int cn_max, int is_total, int debug){
    vector<vector<vector<int>>> s_info;
    num_total_bins = 0;
    // data indexed by [sample][data][ chr, bid, cn ]
    for(int i = 0; i < Ns; ++i) s_info.push_back(vector<vector<int>>());

    igzstream infile(filename.c_str());
    int counter = 0;
    std::string line;
    int prev_sample = 1;

    while(!getline(infile, line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf;
      stringstream ss(line);
      while(ss >> buf) split.push_back(buf);

      int sample = atoi(split[0].c_str());
      // Read next sample
      if(sample > Ns){
          break;
      }
      if(prev_sample != sample){
          num_total_bins = counter;
          counter = 0;
      }

      //cout << sample-1 << "\t" << counter << endl;
      int chr = atoi( split[1].c_str());  // chr
      int sid = atoi( split[2].c_str());  // segment ID
      int cn = -1;
      if(is_total){
          cn = atoi( split[3].c_str());  // copy number
          if(cn > cn_max){
              cout << "copy number " << cn << " is decreased to " << cn_max << endl;
              cn = cn_max;
          }
      }else{
          int cn1 = atoi( split[3].c_str());  // copy number
          int cn2 = atoi( split[4].c_str());  // copy number
          cn = allele_cn_to_state(cn1, cn2);
          if(cn > cn_max){
              cout << "copy number " << cn << " is larger than " << cn_max << "! Please decrease the copy number or increase maximum copy number allowed!" << endl;
          }
      }
      vector<int> vcn{chr, sid, cn};
      s_info[sample-1].push_back(vcn);

      counter++;

      // if(counter >= num_total_bins) counter = 0;
      prev_sample = sample;
    }

    if(debug) cout << "\tSuccessfully read input file" << endl;

    return s_info;
}

// Find the potential number of WGDs for each sample
void get_num_wgd(const vector<vector<vector<int>>>& s_info, int cn_max, vector<int>& obs_num_wgd, int is_total, int debug){
    cout << "Getting the potential number of WGDs for each sample" << endl;
    for(int i = 0; i < s_info.size(); i++){
        vector<vector<int>> s_cn = s_info[i];
        // count the presence of each copy number
        map<int, int> cn_count;
        for(int k = 0; k < cn_max; k++){
            cn_count[k] = 0;
        }
        for (int j = 0; j < s_cn.size(); j++){
            vector<int> cp = s_cn[j];
            int cn = cp[2];
            if(!is_total) cn = state_to_total_cn(cn, cn_max);
            cn_count[cn] += 1;
        }
        // Find the most frequenct copy number
        int most_freq_cn = 0;
        int max_count = 0;
        for(auto cnc : cn_count){
            if(debug) cout << cnc.first << "\t" << cnc.second << endl;
            if(cnc.second > max_count){
                max_count = cnc.second;
                most_freq_cn = cnc.first;
            }
        }

        int nwgd = 0;
        int mode_logcn = ceil(log2(most_freq_cn));
        if(mode_logcn > 1) nwgd = mode_logcn - 1;
        obs_num_wgd.push_back(nwgd);
        if(debug) cout << "Sample " << i+1 << " probably has " << nwgd << " WGD event(s)" << endl;
    }
}


// Find the potential number of chromosome changes for each sample
void get_change_chr(const vector<vector<vector<int>>>& s_info, vector<vector<int>>& obs_change_chr, int cn_max, int is_total, int debug){
    cout << "Getting the potential number of chromosome changes for each sample" << endl;
    for(int i = 0; i < s_info.size(); i++){
        // chr, seg, CN
        vector<vector<int>> s_cn = s_info[i];
        int num_seg = 0;
        double avg_cn = 0;
        map<int, vector<int>> chr_cn;
        vector<int> chr_change;

        // iterate through all records for a sample
        for(int j = 0; j < s_cn.size(); j++){
            vector<int> cp = s_cn[j];
            int cn = cp[2];
            if(!is_total) cn = state_to_total_cn(cn, cn_max);
            chr_cn[cp[0]].push_back(cn);
            avg_cn += cn;
            num_seg++;
        }

        avg_cn = avg_cn / num_seg;

        for(auto c : chr_cn){
            vector<int> cp = c.second;
            int chr_sum_cn = accumulate(cp.begin(), cp.end(), 0);
            double avg_chr_cn = (double) chr_sum_cn / cp.size();
            double num_change = avg_chr_cn - avg_cn;
            int round_num_change = (int) (num_change + 0.5 - (num_change<0));
            // cout << "Number of segments in chromosome " << c.first << " is " << cp.size() << "; avg cn: " << avg_chr_cn << "; exact num changes: " << num_change  << "; num changes: " << round_num_change << endl;
            chr_change.push_back(round_num_change);
        }

        if(debug){
            cout << "Sample " << i+1 << endl;
            // cout << "Number of segments " << num_seg << endl;
            int num_gain = count_if(chr_change.begin(), chr_change.end(), [](int c){return c >= 1;});
            int num_loss = count_if(chr_change.begin(), chr_change.end(), [](int c){return c <= -1;});
            cout << "There are probably " << num_gain << " chromosome gain events and " << num_loss << " chromosome loss events" << endl;
            cout << "Average copy number in the genome " << avg_cn << endl << endl;
        }
        obs_change_chr.push_back(chr_change);
    }
}



// Find the largest copy number in a sample
void get_sample_mcn(const vector<vector<vector<int>>>& s_info, vector<int>& sample_max_cn, int cn_max, int is_total, int debug){
    cout << "Getting the largest copy number for each sample" << endl;
    for(int i = 0; i < s_info.size(); i++){
        vector<vector<int>> s_cn = s_info[i];
        vector<int> cns;
        // iterate through all records for a sample
        for(int j = 0; j < s_cn.size(); j++){
            int cn = s_cn[j][2];
            if(!is_total) cn = state_to_total_cn(cn, cn_max);
            cns.push_back(cn);
        }
        int mcn = *max_element(cns.begin(), cns.end());
        sample_max_cn.push_back(mcn);
        if(debug) cout << "Largest copy number in sample " << i + 1 << " is " << mcn << endl;
    }
}


// Distinguish invariable and variable sites;
// Combine adjacent invariable sites
// TODO: Distinguish different types of invariant sites: 2 vs 4, which have different likelihoods
vector<vector<int>> get_invar_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int is_total, int debug){
    num_invar_bins = 0;
    // Find the number of invariable sites for each character (state)
    // Loop over and output only the regions that have varied, which provide information for tree building
    vector<int> var_bins(num_total_bins, 0);
    for(int k = 0; k < num_total_bins; ++k){
        // // using sum to detect variant bins -- not work for aneuploid genomes
        // int sum = 0;
        // for(int i = 0; i < Ns; ++i){
        //     sum += abs(s_info[i][k][2]);
        // }
        // if((is_total == 1 && sum != 2 * Ns) || (!is_total && sum != 4 * Ns)){    // each site has number 2 when it is total CN or 4 when it is allele-specific CN
        //     var_bins[k] = 1;
        // }
        // else{
        //     num_invar_bins += 1;
        // }
        int has_diff = 0;
        int cn1 = s_info[0][k][2];
        for(int i = 1; i < Ns; ++i){
            int cn2 = s_info[i][k][2];
            if(cn2 != cn1){
              has_diff = 1;
              var_bins[k] = 1;
              break;
            }
        }
        if(has_diff == 0){
          num_invar_bins += 1;
        }
    }

    if(debug){
        cout << "\tVariable bins found:" << endl;
        for(int k = 0; k < num_total_bins; ++k){
            if(var_bins[k] == 1){
            cout << s_info[0][k][0] << "\t" << s_info[0][k][1];
            for(int i = 0; i < Ns; ++i) cout << "\t" << s_info[i][k][2];
            cout << endl;
            }
        }
    }

    int nvar = accumulate(var_bins.begin(), var_bins.end(), 0);
    cout << "\tTotal number of bins:\t" << num_total_bins << endl;
    cout << "\tFound variable bins:\t" << nvar << endl;
    cout << "\tFound invariable bins:\t" << num_invar_bins << endl;

    vector<vector<int>> segs;
    for(int k = 0; k < num_total_bins;){
        if(var_bins[k] == 1){
          //in_seg = true;
          int chr = s_info[0][k][0];
          int seg_start = s_info[0][k][1];
          int id_start = k;

          // hold all the sites in a bin
          vector<int> prev_bin;
          for(int j = 0; j < Ns; ++j) prev_bin.push_back(s_info[j][k][2]);
          //cout << "seg_start: " << chr << "\t" << seg_start << ", cn =  " << seg_cn_tot << endl;

          // Check the subsequent bins
          bool const_cn = true;
          k++;
          while(var_bins[k] == 1 && const_cn){
              vector<int> curr_bin;
              for(int j = 0; j < Ns; ++j) curr_bin.push_back(s_info[j][k][2]);
              if(is_equal_vector(prev_bin, curr_bin) && s_info[0][k][0] == chr){
            	  const_cn = true;
            	  ++k;
              }else{
            	  const_cn = false;
            	  //cout << "\tsplitting segment" << endl;
              }
          }
          int seg_end = s_info[0][k-1][1];
          int id_end = k-1;
          //cout << "seg_end:\t" << seg_end << "\t" << k << endl;
          //cout << endl;

          vector<int> seg{chr, id_start, id_end, seg_start, seg_end};
          segs.push_back(seg);

          // rewind k by one to get the split segment start correct
          if(!const_cn) k--;
        }
        ++k;
    }
    cout << "\tFound segments:\t\t" << segs.size() << endl;

    return segs;
}


vector<vector<int>> get_all_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int incl_all, int is_total, int debug){
    num_invar_bins = 0;
    // Find the number of invariable sites for each character (state)
    // Loop over and output only the regions that have varied
    vector<int> var_bins(num_total_bins, 0);
    for(int k = 0; k < num_total_bins; ++k){
        int sum = 0;
        for(int i = 0; i < Ns; ++i){
          sum += abs(s_info[i][k][2]);
        }
        if((is_total==1 && sum != 2*Ns) || (is_total==0 && sum != 4*Ns)){    // each site has number 2 when it is total CN or 4 when it is allele-specific CN
            var_bins[k] = 1;
        }
        else{
            num_invar_bins += 1;
        }
    }

    if(debug){
        cout << "\tVariable bins found:" << endl;
        for(int k = 0; k < num_total_bins; ++k){
          if(var_bins[k] == 1){
            cout << s_info[0][k][0] << "\t" << s_info[0][k][1];
            for(int i = 0; i < Ns; ++i) cout << "\t" << s_info[i][k][2];
            cout << endl;
          }
        }
    }

    int nvar = accumulate(var_bins.begin(), var_bins.end(), 0);
    cout << "\tTotal number of bins:\t" << num_total_bins << endl;
    cout << "\tFound variable bins:\t" << nvar << endl;
    cout << "\tFound invariable bins:\t" << num_invar_bins << endl;

    vector<vector<int>> segs;
    if(incl_all){
        for(int k = 0; k < num_total_bins; k++){
              int chr = s_info[0][k][0];
              int seg_start = s_info[0][k][1];
              int id_start = k;
              int seg_end = s_info[0][k][1];
              int id_end = k;
              //cout << "seg_end:\t" << seg_end << "\t" << k << endl;
              //cout << endl;
              vector<int> seg{chr, id_start, id_end, seg_start, seg_end};
              segs.push_back(seg);
        }
    }
    else{
        for(int k = 0; k < num_total_bins; k++){
            if(var_bins[k] == 1){
                int chr = s_info[0][k][0];
                int seg_start = s_info[0][k][1];
                int id_start = k;
                int seg_end = s_info[0][k][1];
                int id_end = k;
                //cout << "seg_end:\t" << seg_end << "\t" << k << endl;
                //cout << endl;
                vector<int> seg{chr, id_start, id_end, seg_start, seg_end};
                segs.push_back(seg);
            }
        }
    }

    cout << "\tFound segments:\t\t" << segs.size() << endl;

    return segs;
}


map<int, vector<vector<int>>>  group_segs_by_chr(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug){
    map<int, vector<vector<int>>> ret;
    int Nchar = 0;

    for(int i = 0; i < segs.size(); ++i){
        vector<double> av_cn(Ns,0);     // Compute the average copy number of a segment
        bool valid = true;

        for(int j = 0; j < Ns; ++j){
          for(int k=segs[i][1]; k < (segs[i][2]+1); ++k){
                 av_cn[j] += s_info[j][k][2];
          }
          av_cn[j] = av_cn[j] / ( segs[i][2] - segs[i][1] + 1);
          // The average should be the same as the value of each bin
          assert(av_cn[j] == s_info[j][segs[i][1]][2]);
          // check all cns across the segment are integer valued
          if(ceil(av_cn[j]) != floor(av_cn[j]) ) valid = false;
        }

        if(debug){
          cout << "\t" << segs[i][0] << "\t" << segs[i][1] << "\t" << segs[i][2];
          for(int j = 0; j < Ns; ++j) cout << "\t" << av_cn[j];
          cout << "\t" << valid << endl;
        }

        if(valid){
          vector<int> vals;
          vals.push_back(segs[i][0]); // chr
          vals.push_back(segs[i][1]); // start
          vals.push_back(segs[i][2]); // end
          for(int j = 0; j < Ns; ++j){
          	int cn = (int) av_cn[j];
          	if(cn <= cn_max ) vals.push_back(cn);
          	else vals.push_back(cn_max);
          }
          ret[segs[i][0]].push_back(vals);
          Nchar += 1;
        }
    }

    cout << "\tUsing segments:\t\t" << Nchar << endl;
    if(debug){
        for(auto it : ret){
            vector<vector<int>> sites = it.second;
            for(int j = 0; j < sites.size(); ++j){
                 cout << "\t" << sites[j][0];
                 for(int k = 0; k < Ns; ++k){
                     cout << "\t" << sites[j][k+3];
                 }
                 cout << endl;
            }
        }
    }

    return ret;
}

vector<vector<int>> group_segs(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug){
    vector<vector<int>> ret;

    for(int i = 0; i < segs.size(); ++i){
        vector<double> av_cn(Ns,0);
        bool valid = true;

        for(int j = 0; j < Ns; ++j){
          for(int k=segs[i][1]; k < (segs[i][2]+1); ++k){
                 av_cn[j] += s_info[j][k][2];
          }
          av_cn[j] = av_cn[j]/( segs[i][2] - segs[i][1] + 1);

          // check all cns across the segment are integer valued
          if(ceil(av_cn[j]) != floor(av_cn[j]) ) valid = false;
        }

        if(debug){
          cout << "\t" << segs[i][0] << "\t" << segs[i][1] << "\t" << segs[i][2];
          for(int j = 0; j < Ns; ++j) cout << "\t" << av_cn[j];
          cout << "\t" << valid << endl;
        }

        if(valid){
          // chr, start, end
          vector<int> vals{segs[i][0], segs[i][1], segs[i][2]};
          for(int j = 0; j < Ns; ++j){
          	int cn = (int) av_cn[j];
          	if(cn <= cn_max ) vals.push_back(cn);
          	else vals.push_back(cn_max);
          }
          ret.push_back(vals);
        }
    }

    cout << "\tUsing segments:\t\t" << ret.size() << endl;
    if(debug){
        for(int j = 0; j < ret.size(); ++j){
             for(int k = 0; k < Ns; ++k){
                 cout << "\t" << ret[j][k+3];
             }
             cout << endl;
        }
    }

    return ret;
}

// Read the input copy numbers
vector<vector<int>> read_data_var_regions(const string& filename, const int& Ns, const int& cn_max, int& num_invar_bins, int& num_total_bins, int& seg_size, vector<int>&  obs_num_wgd, vector<vector<int>>& obs_change_chr, vector<int>& sample_max_cn, int model, int is_total, int debug){
    cout << "reading data and calculating CNA regions" << endl;
    vector<vector<vector<int>>> s_info = read_cn(filename, Ns, num_total_bins, cn_max, is_total, debug);
    if(model == DECOMP){
        get_num_wgd(s_info, cn_max, obs_num_wgd, is_total, debug);
        get_change_chr(s_info, obs_change_chr, cn_max, is_total, debug);
    }
    get_sample_mcn(s_info, sample_max_cn, cn_max, is_total);
    // We now need to convert runs of variable bins into segments of constant cn values, grouped by chromosome
    vector<vector<int>> segs = get_invar_segs(s_info, Ns, num_total_bins, num_invar_bins, debug);
    seg_size = segs.size();
    int max_cn_val = cn_max;
    if(!is_total){
        max_cn_val = (cn_max + 1) * (cn_max + 2) / 2 - 1;
    }
    vector<vector<int>> ret = group_segs(segs, s_info, Ns, max_cn_val, debug);

    return ret;
}



// Read the input copy numbers while converting runs of variable bins into segments of constant cn values and group them by chromosome
map<int, vector<vector<int>>> read_data_var_regions_by_chr(const string& filename, const int& Ns, const int& cn_max, int& num_invar_bins, int& num_total_bins, int &seg_size, vector<int>&  obs_num_wgd, vector<vector<int>>& obs_change_chr, vector<int>& sample_max_cn, int model, int is_total, int debug){
    cout << "reading data and calculating CNA regions by chromosome" << endl;
    vector<vector<vector<int>>> s_info = read_cn(filename, Ns, num_total_bins, cn_max, is_total, debug);

    if(model == DECOMP){
        get_num_wgd(s_info, cn_max, obs_num_wgd, is_total);
        get_change_chr(s_info, obs_change_chr, cn_max, is_total);
    }
    get_sample_mcn(s_info, sample_max_cn, cn_max, is_total, debug);

    vector<vector<int>> segs = get_invar_segs(s_info, Ns, num_total_bins, num_invar_bins, is_total, debug);
    seg_size = segs.size();
    int max_cn_val = cn_max;
    if(!is_total){
        max_cn_val = (cn_max + 1) * (cn_max + 2) / 2 - 1;
    }

    map<int, vector<vector<int>>> ret = group_segs_by_chr(segs, s_info, Ns, max_cn_val, debug);

    return ret;
}


// Read the input copy numbers as they are and group them by chromosome
map<int, vector<vector<int>>> read_data_regions_by_chr(const string& filename, const int& Ns, const int& cn_max, int& num_invar_bins, int& num_total_bins, int& seg_size, vector<int>&  obs_num_wgd, vector<vector<int>>& obs_change_chr, vector<int>& sample_max_cn, int model, int incl_all, int is_total, int debug){
    cout << "reading data and calculating CNA regions by chromosome" << endl;
    vector<vector<vector<int>>> s_info = read_cn(filename, Ns, num_total_bins, cn_max, is_total, debug);

    if(model == DECOMP){
        get_num_wgd(s_info, cn_max, obs_num_wgd, is_total, debug);
        get_change_chr(s_info, obs_change_chr, cn_max, is_total, debug);
    }
    get_sample_mcn(s_info, sample_max_cn, cn_max, is_total, debug);

    // We now need to convert runs of variable bins into segments of constant cn values, grouped by chromosome
    vector<vector<int>> segs = get_all_segs(s_info, Ns, num_total_bins, num_invar_bins, incl_all, is_total, debug);
    seg_size = segs.size();
    int max_cn_val = cn_max;
    if(!is_total){
        max_cn_val = (cn_max + 1) * (cn_max + 2) / 2 - 1;
    }

    map<int, vector<vector<int>>> ret = group_segs_by_chr(segs, s_info, Ns, max_cn_val, debug);

    return ret;
}


// Get the input matrix of copy numbers by chromosome
map<int, vector<vector<int>>> get_obs_vector_by_chr(map<int, vector<vector<int>>>& data, const int& Ns){
    map<int, vector<vector<int>>> vobs;
    // Construct the CN matrix by chromosome
    int total_chr = data.rbegin()->first;
    // int total_chr = data.size();   // Some chromosomes got lost in the segment merging
    for(int nchr = 1; nchr <= total_chr; nchr++){
        vector<vector<int>> obs_chr;
        // Nchar += data[nchr].size();
        for(int nc = 0; nc < data[nchr].size(); ++nc){
            vector<int> obs;
            for(int i = 0; i < Ns; ++i){
              obs.push_back(data[nchr][nc][i+3]);
            }
            obs_chr.push_back(obs);
        }
        vobs[nchr] = obs_chr;
    }
    return vobs;
}


void get_bootstrap_vector_by_chr(map<int, vector<vector<int>>>& data, map<int, vector<vector<int>>>& vobs, gsl_rng* r){
    // create a copy of vobs to resample
    map<int, vector<vector<int>>> vobs_copy = vobs;
    vobs.clear();
    int total_chr = data.rbegin()->first;
    // cout << "Total number of chromosomes " << total_chr << endl;
    for(int nchr = 1; nchr <= total_chr; nchr++){
      // cout << "Chr " << nchr << "\t";
      vector<vector<int>> obs_chr;
      for(int nc = 0; nc < data[nchr].size(); ++nc){
            // randomly select a site
           int i = gsl_rng_uniform_int(r, data[nchr].size());
           // cout << i << ":" << vobs_copy[nchr][i].size() << "\t" ;
           obs_chr.push_back(vobs_copy[nchr][i]);
      }
      // cout << obs_chr.size() << endl;
      vobs[nchr] = obs_chr;
    }
}
