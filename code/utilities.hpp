// collection of functions for reading/writing
extern int debug;

vector<double> read_time_info(const string& filename, const int& Ns){
  if(debug) cout << "\tread_time_info" << endl;
  vector<double> t_info;
  ifstream infile (filename.c_str());
  if (infile.is_open()){
    std::string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf;
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

      int dt = atof(split[1].c_str());
      //cout << "read dt: " << dt << endl;
      t_info.push_back(dt);
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
	//cout << "t: " << id << "\t" << start << "\t" << end << "\t" << length << endl;
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

vector<vector<int> > read_data_var_regions(const string& filename, const int& Ns, const int& max_cn){
  cout << "reading data and calculating CNA regions" << endl;
  vector<vector<vector<int> > > s_info;

  // data indexed by [sample][data][ chr, bid, cn ]
  for(int i=0; i<Ns; ++i) s_info.push_back( vector< vector<int> >(4401, vector<int>(3,0) ) );

  //ifstream infile (filename.c_str());
  igzstream infile (filename.c_str());

  int counter = 0;
  //if (infile.is_open()){
    std::string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf;
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

      int sample = atoi(split[0].c_str());
      if( sample > Ns) break;

      //cout << sample-1 << "\t" << counter << endl;
      s_info[sample-1][counter][0] = atoi( split[1].c_str() );
      s_info[sample-1][counter][1] = atof( split[2].c_str() );
      s_info[sample-1][counter][2] = atof( split[3].c_str() );
      counter++;

      if(counter >= 4401) counter = 0;
    }

    //}else{
    // std::cerr << "Error: open of data file unsuccessful: " <<  filename << std::endl;
    //exit(1);
    //}
  cout << "\tSuccessfully read input file" << endl;

  // Loop over and output only the regions that have varied
  vector<int> var_bins(4401,0);
  for(int k=0; k<4401; ++k){
    int sum = 0;
    for(int i=0; i<Ns; ++i){
      sum += abs(s_info[i][k][2]);
    }
    if(sum != 2*Ns) var_bins[k] = 1;
  }

  if(debug){
    cout << "\tVariable bins found:" << endl;
    for(int k=0; k<4401; ++k){
      if(var_bins[k] == 1){
	cout << s_info[0][k][0] << "\t" << s_info[0][k][1];
	for(int i=0; i<Ns; ++i) cout << "\t" << s_info[i][k][2];
	cout << endl;
      }
    }
  }

  int nvar = accumulate(var_bins.begin(), var_bins.end(), 0);
  cout << "\tFound variable bins:\t" << nvar << endl;

  // We now need to convert runs of variable bins into segments of constant cn values

  vector<vector<int> > segs;
  for(int k=0; k<4401;){

    if(var_bins[k] == 1){
      //in_seg = true;
      int chr = s_info[0][k][0];
      int seg_start = s_info[0][k][1];
      int id_start = k;

      // hold the total copy number of the first bin. If this changes we need a new segment
      int seg_cn_tot = 0;
      for(int j=0; j<Ns; ++j) seg_cn_tot += s_info[j][k][2];

      //cout << "seg_start: " << chr << "\t" << seg_start << ", cn= " << seg_cn_tot << endl;

      bool const_cn = true;
      while( var_bins[k] == 1 && s_info[0][k][0] == chr && const_cn == true){

	// calculate new total cn of next bin
	int cn_tot = 0;
	for(int j=0; j<Ns; ++j) cn_tot += s_info[j][k][2];
	//cout << "\tbin:\t" << k+1 << "\t" << cn_tot << endl;
	if( cn_tot == seg_cn_tot){
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

  vector<vector<int> > ret;

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
  //for(int i=0; i<ret.size(); ++i){
  //for(int j=0; j<ret[i].size(); ++j){
  //cout << "\t" << ret[i][j];
  //}
  //cout << endl;
  //}

  return ret;

}