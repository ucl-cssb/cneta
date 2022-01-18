#include "parse_cn.hpp"


// collection of functions for reading/writing copy number data


int state_to_total_cn(int state, int cn_max){
    int cn = 0;
    vector<int> sums;
    // cout << "Sums of state: ";
    for(int i = 1; i <= (cn_max + 1); i++){
        int s = i * (i + 1) / 2;
        sums.push_back(s);
        // cout << "\t" << s;
    }
    // cout << endl;

    if(state < sums[0]) return 0;

    int i = 1;
    do{
        if(state >= sums[i - 1] && state < sums[i]) return i;
        i++;
    }while(i <= cn_max);

    return 0;
}


int state_to_allele_cn(int state, int cn_max, int& cnA, int& cnB){
    int cn = 0;
    vector<int> sums;
    // cout << "Sums of state: ";
    for(int i = 1; i <= (cn_max + 1); i++){
        int s = i * (i + 1) / 2;
        sums.push_back(s);
        // cout << "\t" << s;
    }
    // cout << endl;

    if(state < sums[0]) return 0;

    int i = 1;
    do{
        if(state >= sums[i - 1] && state < sums[i]){
             // total copy number is i;
             int diff = state - sums[i-1];
             cnA = diff;
             cnB = i - cnA;
        }
        i++;
    }while(i <= cn_max);

    return 0;
}


vector<double> read_time_info(const string& filename, const int& Ns, int& age, int debug){
  if(debug) cout << "\tread_time_info" << endl;

  vector<double> t_info;
  vector<int> ages;
  ifstream infile(filename.c_str());
  if(infile.is_open()){
    string line;
    while(!getline(infile, line).eof()){
      if(line.empty()) continue;

      vector<string> split;
      string buf;
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
    cerr << "Error: open of time data unsuccessful: " << filename << endl;
    exit(EXIT_FAILURE);
  }

  if(t_info.size() != Ns){
    cerr << "Error: timing information does not contain " << Ns << " entries: " << filename << endl;
    exit(EXIT_FAILURE);
  }

  return t_info;
}


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


vector<vector<vector<int>>> read_cn(const string& filename, int Ns, int& num_total_bins, int cn_max, int is_total, int is_rcn, int debug){
    vector<vector<vector<int>>> s_info;
    num_total_bins = 0;
    // data indexed by [sample][data][ chr, bid, cn ]
    for(int i = 0; i < Ns; ++i) s_info.push_back(vector<vector<int>>());

    igzstream infile(filename.c_str());
    int counter = 0;
    string line;
    int prev_sample = 1;

    while(!getline(infile, line).eof()){
      if(line.empty()) continue;

      vector<string> split;
      string buf;
      stringstream ss(line);
      while(ss >> buf) split.push_back(buf);

      const char* sstr = split[0].c_str();
      if(!isdigit(*sstr)){
          cout << "Sample ID must be an integer, ordered from 1 to the specified number of patient samples!" << endl;
          exit(EXIT_FAILURE);
      }
      int sample = atoi(sstr);
      // Read next sample
      if(sample > Ns){
          cout << "Skipping sample with ID larger than " << Ns << endl;
          break;
      }
      if(prev_sample != sample){
          num_total_bins = counter;
          counter = 0;
      }

      int chr = atoi(split[1].c_str());  // chr
      int sid = atoi(split[2].c_str());  // segment ID
      int cn = -1;

      if(is_rcn){
          cn = atoi(split[3].c_str());  // copy number
          cn = cn + NORM_PLOIDY;
      }else{
        if(is_total){
            if(split.size() != 4){
                cout << "There should be 4 columns in the file for total copy numbers" << endl;
                exit(EXIT_FAILURE);
            }
            cn = atoi(split[3].c_str());  // copy number
            if(cn < 0){
                cout << "Negative copy numbers in line " << line << "!" << endl;
                exit(EXIT_FAILURE);
            }
            if(cn > cn_max){
                if(debug) cout << "copy number " << cn << " is decreased to " << cn_max << endl;
                cn = cn_max;
            }
        }else{
            if(split.size() != 5){
                cout << "Current file has " << split.size() << " columns " << endl;
                cout << "There should be 5 columns in the file for haplotype-specific copy numbers" << endl;
                exit(EXIT_FAILURE);
            }
            int cn1 = atoi(split[3].c_str());  // copy number
            int cn2 = atoi(split[4].c_str());  // copy number
            if(cn1 < 0 || cn2 < 0){
                cout << "Negative copy numbers in line " << line << "!" << endl;
                exit(EXIT_FAILURE);
            }
            int tcn = cn1 + cn2;
            if(tcn > cn_max){   // a bit hard to decease haplotype-specific CNs
                cout << "total copy number " << tcn << " is larger than " << cn_max << "! Please adjust input accordingly!";
              //   int larger_cn = (cn1 > cn2) ? cn1 : cn2;
              //   int diff = tcn - cn_max;
                exit(EXIT_FAILURE);
            }
            cn = allele_cn_to_state(cn1, cn2);
        }
      }

      vector<int> vcn{chr, sid, cn};
      s_info[sample - 1].push_back(vcn);

      counter++;

      prev_sample = sample;
    }

    if(debug) cout << "\tSuccessfully read input file with " << num_total_bins << " sites" << endl;

    return s_info;
}




void get_num_wgd(const vector<vector<vector<int>>>& s_info, int cn_max, vector<int>& obs_num_wgd, int is_total, int debug){
    cout << "Getting the potential number of WGDs for each sample" << endl;

    for(int i = 0; i < s_info.size(); i++){
        int nwgd = 0;
        int sum_cn = 0;
        vector<vector<int>> s_cn = s_info[i];
        // count the presence of each copy number
        // map<int, int> cn_count;
        // for(int k = 0; k < cn_max; k++){
        //     cn_count[k] = 0;
        // }
        for (int j = 0; j < s_cn.size(); j++){
            vector<int> cp = s_cn[j];
            int cn = cp[2];
            if(!is_total) cn = state_to_total_cn(cn, cn_max);
            // cn_count[cn] += 1;
            sum_cn += cn;
        }
        // // Find the most frequenct copy number
        // int most_freq_cn = 0;
        // int max_count = 0;
        // for(auto cnc : cn_count){
        //     if(debug) cout << cnc.first << "\t" << cnc.second << endl;
        //     if(cnc.second > max_count){
        //         max_count = cnc.second;
        //         most_freq_cn = cnc.first;
        //     }
        // }
        //
        // int mode_logcn = ceil(log2(most_freq_cn));
        // if(mode_logcn > 1) nwgd = mode_logcn - 1;
        // use ploidy, to be applicable to real data, assume at most 1 WGD events
        float ploidy = (float) sum_cn / s_cn.size();
        if(ploidy > WGD_CUTOFF) nwgd = 1;

        obs_num_wgd.push_back(nwgd);
        // if(debug)
          cout << "Sample " << i + 1 << " probably has " << nwgd << " WGD event, with ploidy " << ploidy << endl;
    }
}


void get_change_chr(const vector<vector<vector<int>>>& s_info, vector<vector<int>>& obs_change_chr, int cn_max, int is_total, int debug){
    cout << "Getting the potential number of chromosome changes for each sample" << endl;
    for(int i = 0; i < s_info.size(); i++){
        // chr, seg, CN
        vector<vector<int>> s_cn = s_info[i];
        int num_seg = 0;
        double avg_cn = 0.0;
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
            int round_num_change = (int) (num_change + 0.5 - (num_change < 0));
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



void get_sample_mcn(const vector<vector<vector<int>>>& s_info, vector<int>& sample_max_cn, int cn_max, int is_total, int debug){
    cout << "Getting the largest copy number for each sample" << endl;
    for(int i = 0; i < s_info.size(); i++){
        vector<vector<int>> s_cn = s_info[i];
        vector<int> cns;
        // iterate through all records for a sample
        for(int j = 0; j < s_cn.size(); j++){
            int cn = s_cn[j][2];
            // if(debug) cout << s_cn[j][0] << "\t" << s_cn[j][1] << "\t" << s_cn[j][2] << "\n";
            if(!is_total) cn = state_to_total_cn(cn, cn_max);
            cns.push_back(cn);
        }
        int mcn = *max_element(cns.begin(), cns.end());
        sample_max_cn.push_back(mcn);
        if(debug) cout << "Largest copy number in sample " << i + 1 << " is " << mcn << endl;
    }
}


void get_var_bins(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, vector<int>& var_bins, int is_total, int debug){
    for(int k = 0; k < num_total_bins; ++k){
        // using sum of CNs across samples to detect variant bins does not work for special cases
        vector<int> cns;
        for(int i = 0; i < Ns; ++i){
            int cn = s_info[i][k][2];
            cns.push_back(cn);
        }

        bool is_invar = true;

        if(is_total){
            is_invar = all_of(cns.begin(), cns.end(), [&] (int i) {return i == NORM_PLOIDY;});
        }else{ // For haplotype-specific CN, normal state is 1/1, with ID 4
            is_invar = all_of(cns.begin(), cns.end(), [&] (int i) {return i == NORM_ALLElE_STATE;});
        }

        if(is_invar){
          num_invar_bins += 1;
        }else{
          var_bins[k] = 1;
        }
    }

    if(debug){
        cout << "\tVariable bins found:" << endl;
        for(int k = 0; k < num_total_bins; ++k){
            if(var_bins[k]){
              cout << s_info[0][k][0] << "\t" << s_info[0][k][1];
              for(int i = 0; i < Ns; ++i) cout << "\t" << s_info[i][k][2];
              cout << endl;
            }
        }
    }

    int nvar = accumulate(var_bins.begin(), var_bins.end(), 0);
    cout << "\tTotal number of bins:\t" << num_total_bins << endl;
    cout << "\tNumber of variable bins:\t" << nvar << endl;
    cout << "\tNumber of invariable bins:\t" << num_invar_bins << endl;
}


vector<vector<int>> get_invar_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int is_total, int debug){
    num_invar_bins = 0;
    vector<int> var_bins(num_total_bins, 0);
    if(debug) cout << "\tGetting all the variable bins" << endl;
    get_var_bins(s_info, Ns, num_total_bins, num_invar_bins, var_bins, is_total, debug);

    vector<vector<int>> segs;
    for(int k = 0; k < num_total_bins;){  
        if(var_bins[k]){ // only starting from variable bins
          int chr = s_info[0][k][0];
          int seg_start = s_info[0][k][1];
          int id_start = k;

          // hold all the sites in a bin
          vector<int> prev_bin;
          for(int j = 0; j < Ns; ++j) prev_bin.push_back(s_info[j][k][2]);
        //   if(debug) cout << "seg_start: " << chr << "\t" << seg_start << ", cn =  " << s_info[0][k][2] << endl;

          // Check the subsequent bins
          bool const_cn = true;
          k++;
          while(k < num_total_bins && var_bins[k] && const_cn){ // break if next bin is invarible or has different CN
              vector<int> curr_bin;
              for(int j = 0; j < Ns; ++j) curr_bin.push_back(s_info[j][k][2]);
              if(is_equal_vector(prev_bin, curr_bin) && s_info[0][k][0] == chr){
            	  const_cn = true;
            	  ++k;
              }else{
            	  const_cn = false;
            	//   if(debug) cout << "\tsplitting segment" << endl;
              }
          }
          int seg_end = s_info[0][k-1][1];
          int id_end = k - 1;
        //   if(debug) cout << "seg_end:\t" << seg_end << "\t" << k << endl;

          // id_*: start from 0 until #segments; seg_*: original segment ID
          vector<int> seg{chr, id_start, id_end, seg_start, seg_end};
          segs.push_back(seg);

          // rewind k by one to get the split segment start correct
          if(!const_cn) k--;
        }
        ++k;
    }
    cout << "\tFound segments (after merging consecutive bins):\t" << segs.size() << endl;

    return segs;
}


vector<vector<int>> get_all_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int incl_all, int is_total, int debug){
    num_invar_bins = 0;
    vector<int> var_bins(num_total_bins, 0);
    get_var_bins(s_info, Ns, num_total_bins, num_invar_bins, var_bins, is_total, debug);

    vector<vector<int>> segs;
    if(incl_all){   // bins with normal CNs across all samples are also included in the CN matrix
        for(int k = 0; k < num_total_bins; k++){
            int chr = s_info[0][k][0];
            int seg_start = s_info[0][k][1];
            int id_start = k;
            int seg_end = s_info[0][k][1];
            int id_end = k;
            vector<int> seg{chr, id_start, id_end, seg_start, seg_end};
            segs.push_back(seg);
        }
    }else{  // normal sites will be accounted for based on the count
        for(int k = 0; k < num_total_bins; k++){
            if(var_bins[k]){
                int chr = s_info[0][k][0];
                int seg_start = s_info[0][k][1];
                int id_start = k;
                int seg_end = s_info[0][k][1];
                int id_end = k;
                vector<int> seg{chr, id_start, id_end, seg_start, seg_end};
                segs.push_back(seg);
            }
        }
    }

    cout << "\tFound segments:\t\t" << segs.size() << endl;

    return segs;
}


vector<double> compute_segment_cn(int i, const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug){
  vector<double> av_cn(Ns, 0.0);

  for(int j = 0; j < Ns; ++j){
    for(int k = segs[i][1]; k < (segs[i][2] + 1); ++k){
      av_cn[j] += s_info[j][k][2];
    }
    av_cn[j] = av_cn[j] / (segs[i][2] - segs[i][1] + 1);
    // The average should be the same as the value of each bin
    assert(av_cn[j] == s_info[j][segs[i][1]][2]);

    // check all CNs across the segment are integers
    if(ceil(av_cn[j]) != floor(av_cn[j])){
      cout << "Fractional copy number << " << av_cn[j] << " at " << segs[i][1] << ", " << segs[i][2] << endl;
      exit(EXIT_FAILURE);
    }

    if(av_cn[j] > cn_max){
      cout << "INVALID (larger than maximum allowed) copy number << " << av_cn[j] << " at " << segs[i][1] << ", " << segs[i][2] << endl;
      exit(EXIT_FAILURE);
    }
  }

  if(debug){
    cout << "\n" << segs[i][0] << "\t" << segs[i][1] << "\t" << segs[i][2];
    for(int j = 0; j < Ns; ++j) cout << "\t" << av_cn[j];
    cout << endl;
  }

  return av_cn;
}


map<int, vector<vector<int>>> group_segs_by_chr(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, const string& seg_file, int debug){
    map<int, vector<vector<int>>> ret;
    int Nchar = 0;

    ofstream fcn;
    if(seg_file != "")  fcn.open(seg_file);

    for(int i = 0; i < segs.size(); ++i){
        vector<double> av_cn = compute_segment_cn(i, segs, s_info, Ns, cn_max, debug);

        if(seg_file != "") fcn << segs[i][0] << "\t" << segs[i][1] + 1 << "\t" << segs[i][2] + 1 << "\t" << segs[i][3] << "\t" << segs[i][4];

        vector<int> vals{segs[i][0], segs[i][1], segs[i][2]};   // chr, start, end
        for(int j = 0; j < Ns; ++j){
            int cn = (int) av_cn[j];
            vals.push_back(cn);

            if(seg_file != "") fcn << "\t" << cn;
        }
        if(seg_file != "") fcn << endl;

        ret[segs[i][0]].push_back(vals);
        Nchar += 1;
    }

    if(seg_file != "") fcn.close();
    // output segments in a file for reference by site, which can be converted to same format as original input
    cout << "\tUsing segments:\t\t" << Nchar << endl;
    if(debug){
        // row: sites, column: samples, value: CN
        for(auto it : ret){
            vector<vector<int>> sites = it.second;
            for(int j = 0; j < sites.size(); ++j){
                 cout << sites[j][0] << "\t" << sites[j][1] + 1 << "\t" << sites[j][2] + 1;
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
        vector<double> av_cn = compute_segment_cn(i, segs, s_info, Ns, cn_max, debug);

        // chr, start, end
        vector<int> vals{segs[i][0], segs[i][1], segs[i][2]};
        for(int j = 0; j < Ns; ++j){
        	int cn = (int) av_cn[j];
        	vals.push_back(cn);
        }

        ret.push_back(vals);
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


int get_segs_cn(vector<vector<vector<int>>>& s_info, vector<vector<int>>& segs, const string& filename, const INPUT_PROPERTY& input_prop, INPUT_DATA& input_data, int debug){
    int cn_max = input_prop.cn_max;

    // get CNs for each bins or sites
    s_info = read_cn(filename, input_prop.Ns, input_data.num_total_bins, cn_max, input_prop.is_total, input_prop.is_rcn, debug);

    if(input_prop.model == DECOMP){
        get_num_wgd(s_info, cn_max, input_data.obs_num_wgd, input_prop.is_total, debug);
        get_change_chr(s_info, input_data.obs_change_chr, cn_max, input_prop.is_total, debug);
    }

    get_sample_mcn(s_info, input_data.sample_max_cn, cn_max, input_prop.is_total, debug);

    // find location of segments (consecutive sites) first
    if(input_prop.is_bin){
        // Read the input copy numbers while converting runs of variable bins into segments of constant cn values and group them by chromosome
        segs = get_invar_segs(s_info, input_prop.Ns, input_data.num_total_bins, input_data.num_invar_bins, input_prop.is_total, debug);
    }else{
        // Read the input copy numbers as they are and group them by chromosome
        segs = get_all_segs(s_info, input_prop.Ns, input_data.num_total_bins, input_data.num_invar_bins, input_prop.incl_all, input_prop.is_total, debug);
    }

    input_data.seg_size = segs.size();
    int max_cn_val = cn_max;
    if(!input_prop.is_total){
        max_cn_val = (cn_max + 1) * (cn_max + 2) / 2 - 1;
    }

  return max_cn_val;
}


vector<vector<int>> read_data_var_regions(const string& filename, const INPUT_PROPERTY& input_prop, INPUT_DATA& input_data, int debug){
  cout << "\nReading data and calculating CNA regions" << endl;

  vector<vector<vector<int>>> s_info;
  vector<vector<int>> segs;
  int max_cn_val = get_segs_cn(s_info, segs, filename, input_prop, input_data, debug);

  vector<vector<int>> ret = group_segs(segs, s_info, input_prop.Ns, max_cn_val, debug);

  return ret;
}


map<int, vector<vector<int>>> read_data_var_regions_by_chr(const string& filename, const INPUT_PROPERTY& input_prop, INPUT_DATA& input_data, const string& seg_file, int debug){
    cout << "\nReading data and group regions by chromosome" << endl;

    vector<vector<vector<int>>> s_info;
    vector<vector<int>> segs;
    int max_cn_val = get_segs_cn(s_info, segs, filename, input_prop, input_data, debug);

    // combine segment locations and site-level CNs to get the final CN matrix
    map<int, vector<vector<int>>> ret = group_segs_by_chr(segs, s_info, input_prop.Ns, max_cn_val, seg_file, debug);

    if(debug){
        cout << "Reading input finished" << endl;
    }

    return ret;
}


map<int, vector<vector<int>>> get_obs_vector_by_chr(map<int, vector<vector<int>>>& data, const int& Ns){
    map<int, vector<vector<int>>> vobs;
    // Construct the CN matrix by chromosome
    // Assume chromosomes in data are ordered numberically
    // int total_chr = data.rbegin()->first;  // Some chromosomes got lost in the segment merging, so total_chr may not equal to data.size()
    // int nchr = data.begin()->first;
    // for(; nchr <= total_chr; nchr++){
    for(auto dcn : data){
        int nchr = dcn.first;
        vector<vector<int>> obs_chr;
        for(int nc = 0; nc < data[nchr].size(); ++nc){
            vector<int> obs;
            for(int i = 0; i < Ns; ++i){
              obs.push_back(data[nchr][nc][i + 3]);
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

    // cout << "Total number of chromosomes " << total_chr << endl;
    for(auto dcn : data){
        int nchr = dcn.first;
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
