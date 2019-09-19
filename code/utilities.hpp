// collection of functions for reading/writing
extern int debug;
// const int num_total_bins = 4401;

// Read the samping time and patient age of each sample
vector<double> read_time_info(const string& filename, const int& Ns, int& age){
  if(debug) cout << "\tread_time_info" << endl;
  vector<double> t_info;
  vector<int> ages;
  ifstream infile (filename.c_str());
  if (infile.is_open()){
    std::string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf;
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

      double dt = atof(split[1].c_str());
      //cout << "read dt: " << dt << endl;
      t_info.push_back(dt);

      if(split.size()>2){
          int a = atoi(split[2].c_str());
          ages.push_back(a);
      }
    }
    if(ages.size()>0){
        age = *min_element(ages.begin(), ages.end());
    }
  }else{
    std::cerr << "Error: open of time data unsuccessful: " <<  filename << std::endl;
    exit(1);
  }

  if( t_info.size() != Ns ){
    std::cerr << "Error: timing information does not contain Ns entries: " <<  filename << std::endl;
    exit(1);
  }

  return t_info;
}

evo_tree read_tree_info(const string& filename, const int& Ns){
  if(debug) cout << "\tread_tree_info" << endl;

  vector<edge> edges;
  int id = 0;

  ifstream infile (filename.c_str());
  if (infile.is_open()){
    std::string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf;
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

      if( id > 0 ){
    	int start = atoi(split[0].c_str());
    	int end   = atoi(split[1].c_str());
    	double length = atof(split[2].c_str());
        if(end == Ns+1) length = 0;
        if( !(length>0) && end != Ns+1 ) length = 1;
    	// cout << "t: " << id << "\t" << start << "\t" << end << "\t" << length << endl;
    	edges.push_back( edge( id-1, start-1, end-1, length) );
      }
      id++;
    }
  }else{
    std::cerr << "Error: open of tree data unsuccessful: " <<  filename << std::endl;
    exit(1);
  }

  evo_tree new_tree(Ns+1, edges);
  //new_tree.print();

  return new_tree;
}


// Read a newick tree
// evo_tree read_newick(const string& filename){
//     ifstream infile (filename.c_str());
//     if (infile.is_open()){
//       std::string line;
//       istringstream newickstream(incomingNewick);
//     }
//     else{
//         std::cerr << "Error: open of tree data unsuccessful: " <<  filename << std::endl;
//         exit(1);
//     }
//
//     evo_tree new_tree(Ns+1, edges);
//     //new_tree.print();
//
//     return new_tree;
// }


bool is_equal_vector(const vector<int>& bin1, const vector<int>& bin2){
    assert(bin1.size() == bin2.size());
    for (int i=0; i < bin1.size(); i++){
        if(bin1[i] != bin2[i]){
            return false;
        }
    }
    return true;
}


vector<vector<vector<int>>> read_cn(const string& filename, int Ns, int &num_total_bins){
    vector<vector<vector<int>>> s_info;
    num_total_bins = 0;
    // data indexed by [sample][data][ chr, bid, cn ]
    for(int i=0; i<Ns; ++i) s_info.push_back(vector<vector<int>>());

    igzstream infile (filename.c_str());
    int counter = 0;
    std::string line;
    int prev_sample = 1;

    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf;
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

      int sample = atoi(split[0].c_str());
      // Read next sample
      if( sample > Ns) {
          break;
      }
      if(prev_sample != sample){
          num_total_bins = counter;
          counter = 0;
      }

      //cout << sample-1 << "\t" << counter << endl;
      int chr = atoi( split[1].c_str() );  // chr
      int sid = atoi( split[2].c_str() );  // segment ID
      int cn = atoi( split[3].c_str() );  // copy number
      vector<int> vcn{chr, sid, cn};
      s_info[sample-1].push_back(vcn);
      counter++;

      // if(counter >= num_total_bins) counter = 0;
      prev_sample = sample;
    }

    // cout << "\tSuccessfully read input file" << endl;

    return s_info;
}


vector<vector<int>> get_invar_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins){
    num_invar_bins = 0;
    // Find the number of invariable sites for each character (state)
    // Loop over and output only the regions that have varied
    vector<int> var_bins(num_total_bins, 0);
    for(int k=0; k<num_total_bins; ++k){
        int sum = 0;
        for(int i=0; i<Ns; ++i){
          sum += abs(s_info[i][k][2]);
        }
        if(sum != 2*Ns){
            var_bins[k] = 1;
        }
        else{
            num_invar_bins += 1;
        }
    }

    if(debug){
        cout << "\tVariable bins found:" << endl;
        for(int k=0; k<num_total_bins; ++k){
          if(var_bins[k] == 1){
            cout << s_info[0][k][0] << "\t" << s_info[0][k][1];
            for(int i=0; i<Ns; ++i) cout << "\t" << s_info[i][k][2];
            cout << endl;
          }
        }
    }

    int nvar = accumulate(var_bins.begin(), var_bins.end(), 0);
    cout << "\tTotal number of bins:\t" << num_total_bins << endl;
    cout << "\tFound variable bins:\t" << nvar << endl;
    cout << "\tFound invariable bins:\t" << num_invar_bins << endl;

    vector<vector<int>> segs;
    for(int k=0; k<num_total_bins;){
        if(var_bins[k] == 1){
          //in_seg = true;
          int chr = s_info[0][k][0];
          int seg_start = s_info[0][k][1];
          int id_start = k;

          // hold all the sites in a bin
          vector<int> prev_bin;
          for(int j=0; j<Ns; ++j) prev_bin.push_back(s_info[j][k][2]);
          //cout << "seg_start: " << chr << "\t" << seg_start << ", cn= " << seg_cn_tot << endl;

          // Check the subsequent bins
          bool const_cn = true;
          k++;
          while( var_bins[k] == 1 && const_cn == true){
              vector<int> curr_bin;
              for(int j=0; j<Ns; ++j) curr_bin.push_back(s_info[j][k][2]);
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

          vector<int> seg;
          seg.push_back(chr);
          seg.push_back(id_start);
          seg.push_back(id_end);
          seg.push_back(seg_start);
          seg.push_back(seg_end);
          segs.push_back(seg);

          // rewind k by one to get the split segment start correct
          if(const_cn == false) k--;
        }
        ++k;
    }
    cout << "\tFound segments:\t\t" << segs.size() << endl;

    return segs;
}


vector<vector<int>> get_all_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int incl_all){
    num_invar_bins = 0;
    // Find the number of invariable sites for each character (state)
    // Loop over and output only the regions that have varied
    vector<int> var_bins(num_total_bins, 0);
    for(int k=0; k<num_total_bins; ++k){
        int sum = 0;
        for(int i=0; i<Ns; ++i){
          sum += abs(s_info[i][k][2]);
        }
        if(sum != 2*Ns){
            var_bins[k] = 1;
        }
        else{
            num_invar_bins += 1;
        }
    }

    if(debug){
        cout << "\tVariable bins found:" << endl;
        for(int k=0; k<num_total_bins; ++k){
          if(var_bins[k] == 1){
            cout << s_info[0][k][0] << "\t" << s_info[0][k][1];
            for(int i=0; i<Ns; ++i) cout << "\t" << s_info[i][k][2];
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
        for(int k=0; k<num_total_bins; k++){
              int chr = s_info[0][k][0];
              int seg_start = s_info[0][k][1];
              int id_start = k;
              int seg_end = s_info[0][k][1];
              int id_end = k;
              //cout << "seg_end:\t" << seg_end << "\t" << k << endl;
              //cout << endl;
              vector<int> seg;
              seg.push_back(chr);
              seg.push_back(id_start);
              seg.push_back(id_end);
              seg.push_back(seg_start);
              seg.push_back(seg_end);
              segs.push_back(seg);
        }
    }
    else{
        for(int k=0; k<num_total_bins; k++){
            if(var_bins[k] == 1){
                int chr = s_info[0][k][0];
                int seg_start = s_info[0][k][1];
                int id_start = k;
                int seg_end = s_info[0][k][1];
                int id_end = k;
                //cout << "seg_end:\t" << seg_end << "\t" << k << endl;
                //cout << endl;
                vector<int> seg;
                seg.push_back(chr);
                seg.push_back(id_start);
                seg.push_back(id_end);
                seg.push_back(seg_start);
                seg.push_back(seg_end);
                segs.push_back(seg);
            }
        }
    }

    cout << "\tFound segments:\t\t" << segs.size() << endl;

    return segs;
}


map<int, vector<vector<int>>>  group_segs_by_chr(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int max_cn){
    map<int, vector<vector<int>>> ret;
    int Nchar = 0;

    for(int i=0; i<segs.size(); ++i){
        vector<double> av_cn(Ns,0);     // Compute the average copy number of a segment
        bool valid = true;

        for(int j=0; j<Ns; ++j){
          for(int k=segs[i][1]; k<(segs[i][2]+1); ++k){
                 av_cn[j] += s_info[j][k][2];
          }
          av_cn[j] = av_cn[j]/( segs[i][2] - segs[i][1] + 1 );
          // The average should be the same as the value of each bin
          assert(av_cn[j] == s_info[j][segs[i][1]][2]);
          // check all cns across the segment are integer valued
          if( ceil(av_cn[j]) != floor(av_cn[j]) ) valid = false;
        }

        if(debug){
          cout << "\t" << segs[i][0] << "\t" << segs[i][1] << "\t" << segs[i][2];
          for(int j=0; j<Ns; ++j) cout << "\t" << av_cn[j];
          cout << "\t" << valid << endl;
        }

        if( valid == true ){
          vector<int> vals;
          vals.push_back( segs[i][0] ); // chr
          vals.push_back( segs[i][1] ); // start
          vals.push_back( segs[i][2] ); // end
          for(int j=0; j<Ns; ++j){
        	int cn = (int) av_cn[j];
        	if( cn <= max_cn ) vals.push_back( cn );
        	else vals.push_back( max_cn );
          }
          ret[segs[i][0]].push_back( vals );
          Nchar += 1;
        }

    }

    cout << "\tUsing segments:\t\t" << Nchar << endl;
    if(debug){
        for(auto it : ret){
            vector<vector<int>> sites = it.second;
            for(int j=0; j<sites.size(); ++j){
                 cout << "\t" << sites[j][0];
                 for(int k=0; k<Ns; ++k){
                     cout << "\t" << sites[j][k+3];
                 }
                 cout << endl;
            }
        }
    }

    return ret;
}

vector<vector<int>> group_segs(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int max_cn){
    vector<vector<int>> ret;

    for(int i=0; i<segs.size(); ++i){
        vector<double> av_cn(Ns,0);
        bool valid = true;

        for(int j=0; j<Ns; ++j){
          for(int k=segs[i][1]; k<(segs[i][2]+1); ++k){
                 av_cn[j] += s_info[j][k][2];
          }
          av_cn[j] = av_cn[j]/( segs[i][2] - segs[i][1] + 1 );

          // check all cns across the segment are integer valued
          if( ceil(av_cn[j]) != floor(av_cn[j]) ) valid = false;
        }

        if(debug){
          cout << "\t" << segs[i][0] << "\t" << segs[i][1] << "\t" << segs[i][2];
          for(int j=0; j<Ns; ++j) cout << "\t" << av_cn[j];
          cout << "\t" << valid << endl;
        }

        if( valid == true ){
          vector<int> vals;
          vals.push_back( segs[i][0] ); // chr
          vals.push_back( segs[i][1] ); // start
          vals.push_back( segs[i][2] ); // end
          for(int j=0; j<Ns; ++j){
        	int cn = (int) av_cn[j];
        	if( cn <= max_cn ) vals.push_back( cn );
        	else vals.push_back( max_cn );
          }
          ret.push_back( vals );
        }
    }

    cout << "\tUsing segments:\t\t" << ret.size() << endl;
    if(debug){
        for(int j=0; j<ret.size(); ++j){
             for(int k=0; k<Ns; ++k){
                 cout << "\t" << ret[j][k+3];
             }
             cout << endl;
        }
    }

    return ret;
}

// Read the input copy numbers
vector<vector<int>> read_data_var_regions(const string& filename, const int& Ns, const int& max_cn, int &num_invar_bins, int &num_total_bins, int &seg_size){
    cout << "reading data and calculating CNA regions" << endl;
    vector<vector<vector<int>>> s_info = read_cn(filename, Ns, num_total_bins);
    // We now need to convert runs of variable bins into segments of constant cn values, grouped by chromosme
    vector<vector<int>> segs = get_invar_segs(s_info, Ns, num_total_bins, num_invar_bins);
    seg_size = segs.size();

    vector<vector<int>> ret = group_segs(segs, s_info, Ns, max_cn);

    return ret;
}



// Read the input copy numbers and group them by chromosme
map<int, vector<vector<int>>> read_data_var_regions_by_chr(const string& filename, const int& Ns, const int& max_cn, int &num_invar_bins, int &num_total_bins, int &seg_size){
    cout << "reading data and calculating CNA regions by chromosme" << endl;
    vector<vector<vector<int>>> s_info = read_cn(filename, Ns, num_total_bins);
    // We now need to convert runs of variable bins into segments of constant cn values, grouped by chromosme
    vector<vector<int>> segs = get_invar_segs(s_info, Ns, num_total_bins, num_invar_bins);
    seg_size = segs.size();

    map<int, vector<vector<int>>> ret = group_segs_by_chr(segs, s_info, Ns, max_cn);

    return ret;
}


// Read the input copy numbers as they are and group them by chromosme
map<int, vector<vector<int>>> read_data_regions_by_chr(const string& filename, const int& Ns, const int& max_cn, int &num_invar_bins, int &num_total_bins, int &seg_size, int incl_all=1){
    cout << "reading data and calculating CNA regions by chromosme" << endl;
    vector<vector<vector<int>>> s_info = read_cn(filename, Ns, num_total_bins);
    // We now need to convert runs of variable bins into segments of constant cn values, grouped by chromosme
    vector<vector<int>> segs = get_all_segs(s_info, Ns, num_total_bins, num_invar_bins, incl_all);
    seg_size = segs.size();

    map<int, vector<vector<int>>> ret = group_segs_by_chr(segs, s_info, Ns, max_cn);

    return ret;
}



//vector<vector<int> > vobs; // already defined globally
// for(int nc=0; nc<Nchar; ++nc) {
//     vector<int> obs;
//     for(int i=0; i<Ns; ++i) {
//             obs.push_back(data[nc][i+3]);
//     }
//     vobs.push_back(obs);
// }


// Get the input matrix of copy numbers by chromosme
map<int, vector<vector<int>>> get_obs_vector_by_chr(map<int, vector<vector<int>>>& data){
    map<int, vector<vector<int>>> vobs;
    // Construct the CN matrix by chromosome
    int total_chr = data.rbegin()->first;
    // int total_chr = data.size();   // Some chromosmes got lost in the segment merging
    for(int nchr=1; nchr <= total_chr; nchr++){
        vector<vector<int>> obs_chr;
        // Nchar += data[nchr].size();
        for(int nc=0; nc<data[nchr].size(); ++nc){
            vector<int> obs;
            for(int i=0; i<Ns; ++i){
              obs.push_back( data[nchr][nc][i+3] );
            }
            obs_chr.push_back( obs );
        }
        vobs[nchr] = obs_chr;
    }
    return vobs;
}


void get_bootstrap_vector_by_chr(map<int, vector<vector<int>>>& data, map<int, vector<vector<int>>>& vobs){
    // create a copy of vobs to resample
    map<int, vector<vector<int>>> vobs_copy = vobs;
    vobs.clear();
    int total_chr = data.rbegin()->first;
    // cout << "Total number of chromosmes " << total_chr << endl;
    for(int nchr=1; nchr<=total_chr; nchr++){
      // cout << "Chr " << nchr << "\t";
      vector<vector<int>> obs_chr;
      for(int nc=0; nc<data[nchr].size(); ++nc){
            // randomly select a site
           int i = gsl_rng_uniform_int(r, data[nchr].size());
           // cout << i << ":" << vobs_copy[nchr][i].size() << "\t" ;
           obs_chr.push_back(vobs_copy[nchr][i]);
      }
      // cout << obs_chr.size() << endl;
      vobs[nchr] = obs_chr;
    }
}
