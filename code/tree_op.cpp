#include "tree_op.hpp"


// check if the age of tip node is maintained
bool is_tip_age_valid(const vector<double>& node_ages, const vector<double>& tobs){
  // cout << "tobs:";
  // for(int i = 0; i < tobs.size(); i++){
  //   cout << "\t" << tobs[i];
  // }
  // cout << endl;
  // cout << "node_ages:";
  // for(int i = 0; i < node_ages.size(); i++){
  //   cout << "\t" << node_ages[i];
  // }
  // cout << endl;

  double max_time = *max_element(tobs.begin(), tobs.end());
  for(int i = 0; i < tobs.size(); i++){
    if(fabs(node_ages[i] - (max_time - tobs[i])) > SMALL_DIFF_BRANCH){
    //   cout << i << "\t" << node_ages[i] << "\t" << tobs[i] << endl;
      return false;
    }
  }
  return true;
}

// check if node_ages and node_times are consistent
bool is_age_time_consistent(const vector<double>& node_times, const vector<double>& node_ages){
  double max_age = *max_element(node_ages.begin(), node_ages.end());
  double max_time = *max_element(node_times.begin(), node_times.end());

  if(fabs(max_age - max_time) > SMALL_DIFF_BRANCH){
    return false;
  }else{
    return true;
  }
}


void check_node_age_ratio(evo_tree& tree, const vector<int>& knodes){
  // To confirm the conversion of edge to ratio is correct
  vector<double> ratios = tree.get_ratio_from_age();
  cout << "ratios of node times: ";
  for(int i = 0; i < ratios.size(); i++){
      cout << i + 1 << "\t" << "\t" << ratios[i] << endl;
  }
  tree.update_edges_from_ratios(ratios, knodes);
  cout << "New branch lengths: " << endl;
  for(int i = 0; i < tree.edges.size(); i++){
      cout << i + 1 << "\t" << "\t" << tree.edges[i].length << endl;
  }
}

void check_nodes_transverse(evo_tree& tree){
  cout << "Postorder traversal of internal nodes in original tree";
  vector<int> knodes;
  for(int k = (tree.nleaf + 1); k < (2 * tree.nleaf - 1); ++k){
    knodes.push_back(k);
  }
  knodes.push_back(tree.nleaf);
  for(auto k: knodes){
    cout << "\t" << k;
  }
  cout << endl;

  vector<int> knodes2;
  tree.get_inodes_postorder(&tree.nodes[tree.root_node_id], knodes2);
  cout << "Postorder traversal of internal nodes";
  for(auto k: knodes2){
    cout << "\t" << k;
  }
  cout << endl;
}


// The string will be different for different labeled histories (topologies may be the same but timings are different)
string create_tree_string(const evo_tree& tree){
  stringstream sstm;
  int ntotn = 2 * tree.nleaf - 1;

  for(int i = 0; i < ntotn; ++i){
    sstm << tree.nodes[i].id + 1;
    if(tree.nodes[i].daughters.size() == 2){
      sstm << ";" << tree.nodes[i].daughters[0] + 1 << ";" << tree.nodes[i].daughters[1] + 1;
    }
    sstm << ":";
  }

  return sstm.str();
}


// The string will be unique for different topologies
// The key is to name the internal node based on tips
string create_tree_string_uniq(const evo_tree& tree){
  int debug = 0;

  stringstream sstm;
  // Use a dictory to store the names of internal nodes;
  map<int, string> names;
  int ntotn = 2 * tree.nleaf - 1;

  for(int i = 0; i < ntotn; ++i){
    int nid = tree.nodes[i].id + 1;
    if(nid <= tree.nleaf){
        sstm << nid;
    }
    else if(nid == tree.nleaf + 1){  // root
        sstm << nid << ";" << tree.nodes[i].daughters[0] + 1 << ";" << tree.nodes[i].daughters[1] + 1;
    }else{    // internal nodes
        assert((tree.nodes[i].daughters.size() == 2));
        // cout << "Converting internal nodes" << endl;
        string pnode = "n_";
        int d1 = tree.nodes[i].daughters[0] + 1;
        int d2 = tree.nodes[i].daughters[1] + 1;
        if(d1 <= tree.nleaf && d2 <= tree.nleaf){   // The first time seeing this internal node
            // Sort the children so that the name is unique
            int min_d = d1;
            int max_d = d2;
            if(d1 > d2){
                 min_d = d2;
                 max_d = d1;
            }
            pnode = pnode + to_string(min_d) + "-" + to_string(max_d);
            // cout << nid << "\t" << pnode << endl;
            names[nid] = pnode;
            sstm << pnode << ";" << min_d << ";" << max_d;
        }else{
            // The internal nodes are numbered increasingly
            // cout << nid << ";" << d1 << ";" << d2 << endl;
            if(d1 > tree.nleaf && d2 <= tree.nleaf){
                pnode = pnode + to_string(d2) + "-[" + names[d1] + "]";
                // cout << nid << "\t" << pnode << "\t" <<  names[d1] << endl;
                names[nid] = pnode;
                sstm << pnode << ";" << d2 << ";" << names[d1];
            }else if(d2 > tree.nleaf && d1 <= tree.nleaf){  // d2 > rtree.nleaf - 1
                pnode = pnode + to_string(d1) + "-[" + names[d2] + "]";
                // cout << nid << "\t" << pnode << "\t" <<  names[d2] << endl;
                names[nid] = pnode;
                sstm << pnode << ";" << d1 << ";" << names[d2];
            }else{  // both nodes are not tips
                // Ensure pnode is unique by sorting names of daughters
                if(names[d1].compare(names[d2]) < 0){
                    pnode = pnode + "[" + names[d1] + "]-[" + names[d2] + "]";
                }else{
                    pnode = pnode + "[" + names[d2] + "]-[" + names[d1] + "]";
                }
                // cout << nid << "\t" << pnode << "\t" <<  names[d2] << endl;
                names[nid] = pnode;
                sstm << pnode << ";" << names[d1] << ";" << names[d2];
            }
        }
    }
    sstm << ":";
  }

  if(debug){
      for(auto it: names){
          cout << it.first << "\t" << it.second << endl;
      }
  }

  return sstm.str();
}


// for different labeled histories
string order_tree_string(const string& tree){
  stringstream sstm;
  vector<string> split1;
  boost::split(split1, tree, [](char c){ return c == ':'; });

  for(int i = 0; i < split1.size() - 1; ++ i){     // split creates an empty string at the end
    //sstm << split1[i];
    //cout << "\t" << split1[i] << endl;
    vector<string> split2;
    boost::split(split2, split1[i], [](char c){ return c == ';'; });

    if(split2.size() == 1){
      sstm << split1[i];
    }
    else{
      sstm << split2[0] << ";"; //  << split2[1] << ";" << split2[2];
      string s1 = split2[1];
      string s2 = split2[2];
      // if( atoi(split2[1].c_str() ) < atoi(split2[2].c_str() ) ){
      if(s1.compare(s2) < 0){
	         sstm << split2[1] << ";" << split2[2];
      }else{
	         sstm << split2[2] << ";" << split2[1];
      }
    }
    sstm << ":";
  }
  return sstm.str();
}


// for different topologies
// Sort the parent nodes at first
string order_tree_string_uniq(const string& tree){
  stringstream sstm;
  vector<string> split1;

  boost::split(split1, tree, [](char c){ return c == ':'; });
  // Sort split1
  sort(split1.begin(), split1.end());

  for(int i = 0; i < split1.size()-1; ++ i){     // split creates an empty string at the end
    //sstm << split1[i];
    //cout << "\t" << split1[i] << endl;
    vector<string> split2;
    boost::split(split2, split1[i], [](char c){ return c == ';'; });

    if(split2.size() == 1){
      sstm << split1[i];
    }else{
      sstm << split2[0] << ";"; //  << split2[1] << ";" << split2[2];
      string s1 = split2[1];
      string s2 = split2[2];
      // if( atoi(split2[1].c_str() ) < atoi(split2[2].c_str() ) ){
      if( s1.compare(s2) < 0 ){
	         sstm << split2[1] << ";" << split2[2];
      }else{
	         sstm << split2[2] << ";" << split2[1];
      }
    }
    sstm << ":";
  }

  return sstm.str();
}




// Randomly swap two leaves
evo_tree perturb_tree(evo_tree& tree, long unsigned (*fp_myrng)(long unsigned)){
    int debug = 0;
    if(debug) cout << "\tperturb one tree" << endl;

    // Ns = nleaf - 1
    // swap two leaf nodes: 0 to Ns-1 are the samples, Ns is germline, Ns + 1 is root
    vector<int> n_0;
    for(int i = 0; i < tree.nleaf - 1; ++i) n_0.push_back(i);

    // shuffle(n_0.begin(), n_0.end(), default_random_engine(seed));
    random_shuffle(n_0.begin(), n_0.end(), fp_myrng);

    int n1 = n_0[tree.nleaf - 2];
    n_0.pop_back();

    // shuffle(n_0.begin(), n_0.end(), default_random_engine(seed));
    random_shuffle(n_0.begin(), n_0.end(), fp_myrng);

    int n2 = n_0[0];

    if(debug){
        cout << "\toriginal tree: ";
        cout << order_tree_string_uniq(create_tree_string_uniq(tree)) << endl;
        cout << tree.make_newick(5) << endl;
        cout << "\tswapping nodes:\t" << n1 + 1 << "\t" << n2 + 1 << endl;
    }

    vector<edge> enew;
    int e1, e2;
    int end1, end2;
    double len1, len2;
    for(int i = 0; i < tree.edges.size(); ++i){
        enew.push_back(tree.edges[i]);
        bool changed = false;
        if(enew[i].end == n1 && !changed){
          enew[i].end = n2;
          e1 = i;
          end1 = n1;
          len1 = enew[i].length;
          changed = true;
        }
        if(enew[i].end == n2 && !changed){
          enew[i].end = n1;
          e2 = i;
          end2 = n2;
          len2 = enew[i].length;
          changed = true;
        }
    }

    // cout << "edge 1 " << "\t" << e1  << "\t" << enew[e1].start + 1 << "\t" << end1 + 1 << "\t" <<
    // enew[e1].length << endl;
    // cout << "edge 2 " << "\t" << e2  << "\t" << enew[e2].start + 1 << "\t" << end2 + 1 << "\t" <<
    // enew[e2].length << endl;

    enew[e1].length = len2;
    enew[e2].length = len1;

    // cout << "edge 1 " << "\t" << e1  << "\t" << enew[e1].start + 1 << "\t" << enew[e1].end + 1 << "\t" <<
    // enew[e1].length << endl;
    // cout << "edge 2 " << "\t" << e2  << "\t" << enew[e2].start + 1 << "\t" << enew[e2].end + 1 << "\t" <<
    // enew[e2].length << endl;

    evo_tree new_tree(tree.nleaf, enew, 1);

    new_tree.mu = tree.mu;
    new_tree.dup_rate = tree.dup_rate;
    new_tree.del_rate = tree.del_rate;
    new_tree.chr_gain_rate = tree.chr_gain_rate;
    new_tree.chr_loss_rate = tree.chr_loss_rate;
    new_tree.wgd_rate = tree.wgd_rate;


    if(debug){
        // cout << "\tnew tree: ";
        // cout << create_tree_string_uniq( new_tree ) << endl;
        string ordered = order_tree_string_uniq(create_tree_string_uniq(new_tree));
        cout << "\tordered new tree: ";
        cout << ordered << endl;
        cout << new_tree.make_newick(5) << endl;
        cout << endl;
    }

    return new_tree;
}


// test if the tree is correctly built
void test_evo_tree(const evo_tree& tree){
  // generate the list of ancestral nodes belonging to leaf nodes
  cout << "leaf ancestral nodes:" << endl;
  for(int i = 0; i < tree.nleaf-1; ++i){
    vector<int> anodes = tree.get_ancestral_nodes( tree.nodes[i].id);
    cout << "\tnode " << tree.nodes[i].id+1;
    for(int j = 0; j < anodes.size(); ++j) cout << "\t" << tree.nodes[anodes[j]].id+1;
    cout << endl;
  }

  // generate the list of ancestral edges belonging to leaf nodes
  cout << "leaf ancestral edges:" << endl;
  for(int i = 0; i < tree.nleaf-1; ++i){
    vector<int> aedges = tree.get_ancestral_edges( tree.nodes[i].id);
    reverse(aedges.begin(),aedges.end());
    cout << "\tnode " << tree.nodes[i].id+1;
    for(int j = 0; j < aedges.size(); ++j) cout << "\t" << tree.edges[aedges[j]].id+1
					    << " (" << tree.edges[aedges[j]].start+1 << " -> " << tree.edges[aedges[j]].end+1 << ") ";
    cout << endl;
  }
}


// generate neutral coalescent trees
// here nsample is the number of cancer samples (not including germline node)
// here we directly calculate the edges in the tree
void generate_coal_tree(const int& nsample, gsl_rng* r, long unsigned (*fp_myrng)(long unsigned), vector<int>& edges, vector<double>& lengths, vector<double>& epoch_times, vector<double>& times, const int& Ne, const double& beta, const double& gtime){
  vector<int> nodes;

  // For compatibility with ape in R, root node must be labelled
  // 0 to nsample-1 are leaf, nsample is germline, nsample + 1 is root
  // all other nodes ids start from here
  int node_count = nsample + 2;

  // create leaf nodes
  for(int i = 0; i < nsample; ++i) nodes.push_back(i);

  // create vector of event times
  // total nodes = 2*(nsample + 1) -1
  for(int i = 0; i < (2*nsample + 1); ++i) times.push_back(0.0);

  double t_tot = 0;
  double curr_time = 0;
  int nlin = nsample;
  while(nlin > 1){
    // sample a time from Exp( combinations(k,2) )
    // double lambda = fact(nlin)/( 2*fact(nlin-2));
    double lambda = nlin * (nlin - 1) / 2;
    double tn = gsl_ran_exponential(r, 1/lambda) * Ne;
    // cout << "Normal coalescence time  " << tn << endl;
    double t = tn;
    if(beta > 0){  // simulate exponential growth
        // t = (log(1 + beta * tn * exp(- beta * t_tot))) / beta;  // Formula in book "Gene Genealogies, Variation and Evolution", P99
        t = ((log(exp(beta*t_tot) + beta * tn)) / beta)- t_tot;    // Formula in CoalEvol7.3.5.c, equivalent to the one above
    }
    // cout << "Exponential coalescence time  " << t << endl;

    t_tot += t;
    // choose two random nodes from available list
    random_shuffle(nodes.begin(), nodes.end(), fp_myrng);

    // edge node_count -> node
    edges.push_back(node_count);
    edges.push_back(nodes[nodes.size()-1]);
    lengths.push_back(t_tot - times[ nodes[nodes.size()-1] ] );

    edges.push_back(node_count);
    edges.push_back(nodes[nodes.size()-2]);
    lengths.push_back(t_tot - times[ nodes[nodes.size()-2] ] );

    // update time for this node
    //cout << "\t" << node_count + 1 << "\t" << t_tot << endl;
    epoch_times.push_back(t_tot);
    times[ node_count ] = t_tot;

    nodes.pop_back();
    nodes.pop_back();

    nodes.push_back(node_count);
    node_count++;
    nlin--;
  }

  cout << "TMRCA of tumour samples (scaled by N): " << t_tot << endl;
  // create the root and germline nodes and edges
  double lambda = runiform(r, 0, 1);
  double tn = gsl_ran_exponential(r, 1/lambda) * Ne;
  double t = tn;
  if(beta > 0){  // simulate exponential growth
      // t = (log(1 + beta * tn * exp(- beta * t_tot))) / beta;
      t = ((log(exp(beta*t_tot) + beta * tn)) / beta)- t_tot;    // Formula in CoalEvol7.3.5.c
  }
  t_tot += t;
  epoch_times.push_back(t_tot);

  // add in the time for the root and germline nodes
  times[nsample + 1] = t_tot;
  times[nsample] = t_tot;

  edges.push_back(nsample + 1);
  edges.push_back(node_count-1);
  lengths.push_back(t);

  edges.push_back(nsample + 1);
  edges.push_back(nsample);
  lengths.push_back(0);

  // invert the times
  for(int i = 0; i < times.size(); ++i){
    times[i] = t_tot - times[i];
  }

  for(int l = 0; l < lengths.size(); ++l){
      lengths[l] = lengths[l] * gtime;
  }
  //cout << "total time of tree: " << t_tot << " : ";
  //for(int i = 0; i < epoch_times.size(); ++i) cout << "\t" << epoch_times[i];
  //cout << endl;

  // invert the times
  //cout << "times of nodes:" << endl;
  //for(int i = 0; i < times.size(); ++i){
    //cout << i + 1 << "\t" << times[i] << endl;
  //}
}



// Scale the total time by given time
evo_tree generate_coal_tree(const int& nsample, gsl_rng* r, long unsigned (*fp_myrng)(long unsigned), int Ne, double beta, double gtime){
   int debug = 0;
   if(debug) cout << "GENERATING COALESCENCE TREE" << endl;
   vector<int> edges;
   vector<double> lengths;
   vector<double> epoch_times;
   vector<double> times;
   vector<int> nodes;

   // For compatibility with ape in R, root node must be labelled
   // 0 to nsample-1 are leaf, nsample is germline, nsample + 1 is root
   // all other nodes ids start from here
   int node_count = nsample + 2;

   // create leaf nodes
   for(int i = 0; i < nsample; ++i) nodes.push_back(i);

   // create vector of event times
   // total nodes = 2*(nsample + 1) -1
   for(int i = 0; i < (2 * nsample + 1); ++i) times.push_back(0.0);

   double t_tot = 0;
   int nlin = nsample;
   while(nlin > 1){
     // sample a time from Exp( combinations(k,2) )
     //  double lambda = fact(nlin) / ( 2 * fact(nlin - 2) );    // may become overflow using fact
     double lambda = nlin * (nlin - 1) / 2;
     double tn = gsl_ran_exponential(r, 1 / lambda) * Ne;
     if(debug) cout << "Normal coalescence time  " << tn << endl;
     double t = tn;
     if(beta > 0){  // simulate exponential growth
         // t = (log(1 + beta * tn * exp(- beta * t_tot))) / beta;  // Formula in book "Gene Genealogies, Variation and Evolution", P99
         t = ((log(exp(beta*t_tot) + beta * tn)) / beta)- t_tot;    // Formula in CoalEvol7.3.5.c, equivalent to the one above
     }
     if(debug) cout << "exponential coalescence time  " << t << endl;
     t_tot += t;

     // choose two random nodes from available list
     random_shuffle(nodes.begin(), nodes.end(), fp_myrng);

     // edge node_count -> node
     edges.push_back(node_count);
     edges.push_back(nodes[nodes.size() - 1]);
     lengths.push_back(t_tot - times[nodes[nodes.size() - 1]]);

     edges.push_back(node_count);
     edges.push_back(nodes[nodes.size()-2]);
     lengths.push_back(t_tot - times[nodes[nodes.size() - 2]]);

     // update time for this node
     //cout << "\t" << node_count + 1 << "\t" << t_tot << endl;
     epoch_times.push_back(t_tot);
     times[node_count] = t_tot;

     nodes.pop_back();
     nodes.pop_back();

     nodes.push_back(node_count);
     node_count++;
     nlin--;
   }

   // create the root and germline nodes and edges
   // double lambda = 1;
   double lambda = runiform(r, 0, 1);
   double tn = gsl_ran_exponential(r, 1/lambda) * Ne;
   if(debug) cout << "Normal coalescence time  " << tn << endl;
   double t = tn;
   if(beta > 0){  // simulate exponential growth
       // t = (log(1 + beta * tn * exp(- beta * t_tot))) / beta;  // Formula in book "Gene Genealogies, Variation and Evolution", P99
       t = ((log(exp(beta*t_tot) + beta * tn)) / beta)- t_tot;    // Formula in CoalEvol7.3.5.c, equivalent to the one above
   }
   if(debug) cout << "exponential coalescence time  " << t << endl;
   t_tot += t;
   epoch_times.push_back(t_tot);

   // add in the time for the root and germline nodes
   times[nsample + 1] = t_tot;
   times[nsample] = t_tot;

   edges.push_back(nsample + 1);
   edges.push_back(node_count - 1);
   lengths.push_back(t);

   edges.push_back(nsample + 1);
   edges.push_back(nsample);
   lengths.push_back(0);

   // invert the times
   for(int i = 0; i < times.size(); ++i){
     times[i] = t_tot - times[i];
   }

   for(int l = 0; l < lengths.size(); ++l){
       lengths[l] = lengths[l] * gtime;
   }
   // cout << "total time of tree: " << t_tot << " : ";
   // for(int i = 0; i < epoch_times.size(); ++i) cout << "\t" << epoch_times[i];
   // cout << endl;
   //
   // // invert the times
   // cout << "times of nodes:" << endl;
   // for(int i = 0; i < times.size(); ++i){
   //   cout << i + 1 << "\t" << times[i] << endl;
   // }

   evo_tree ret(nsample + 1, edges, lengths);

   return ret;
 }



// randomly assign leaf edges to time points t0, t1, t2, t3, ...
// adjust terminal branch lengths accordingly
void assign_tip_times(double delta_t, int Ns, gsl_rng* r, vector<double>& tobs, const vector<int>& edges, vector<double>& lengths){
    tobs.clear(); // clear the sample time differences if there are some

    cout << "assigning temporal structure" << endl;

    bool assign0 = false;
    int bcount = 0;
    int num_diff = Ns; // maximum of tips with different times

    for(int l = 1; l < edges.size(); l = l + 2){
        if(edges[l] < Ns){
            int ind = 0;

            if(!assign0){
            ind = 0;
            assign0 = true;
            }else{
            ind = gsl_rng_uniform_int(r, num_diff);
            }

            double stime = ind * delta_t;
            cout << "\t sample / sampling time point :" << edges[l] + 1 << "\t" << stime << endl;
            lengths[bcount] = lengths[bcount] + stime;
            // Store the sampling time for time scaling by the age of a patient
            tobs.push_back(stime);
        }
        bcount++;
    }
}


evo_tree generate_random_tree(int Ns, gsl_rng* r, long unsigned (*fp_myrng)(long unsigned), int Ne, int age, double beta, double gtime, double delta_t, int cons, int debug){
    vector<int> edges;
    vector<double> lengths;
    vector<double> epoch_times;
    vector<double> node_times;

    cout << "Generating random coalescence trees" << endl;
    cout << " The effective population size is " << Ne << endl;
    cout << " The generation time is " << gtime << endl;
    if(beta > 0){
        cout << " The exponential growth rate is " << beta << endl;
    }
    generate_coal_tree(Ns, r, fp_myrng, edges, lengths, epoch_times, node_times, Ne, beta, gtime);

    if(debug){
      cout << "Initial coalescence tree: " << endl;
      evo_tree test_tree0(Ns + 1, edges, lengths);
      test_tree0.print();
    }

    vector<double> tobs(Ns, 0.0);
    if(delta_t > 0){
        assign_tip_times(delta_t, Ns, r, tobs, edges, lengths);
    }

    evo_tree test_tree(Ns + 1, edges, lengths);
    if(debug){
      cout << "Tree after assigning sampling times to tips: " << endl;
      test_tree.print();
    }
    // Ensure that the rescaled tree by age has a height smaller than the age of last sample
    // Need tobs to get the miminmal tree height
    double old_tree_height = get_tree_height(test_tree.get_node_times());
    if(cons && old_tree_height > age){
      double min_height = *max_element(tobs.begin(), tobs.end());
      double max_height = age + min_height;
      double tree_height = runiform(r, min_height, max_height);    
      double ratio = tree_height / old_tree_height;

      if(delta_t > 0){
        // Only scaling internal nodes to allow large differences at the tip
        test_tree.scale_time_internal(ratio);
      }else{
        test_tree.scale_time(ratio);
      }

      if(debug){
          cout << "Tree height before scaling " << old_tree_height << endl;
          cout << "Scaling ratio " << ratio << endl;
          test_tree.print();
      }
    }

    return test_tree;
}


// Create a new tree with the same topology as input tree but different branch lengths
evo_tree create_new_tree(gsl_vector* blens, evo_tree& rtree, const double& max_tobs, int cons){
    int nedge = 2 * rtree.nleaf - 2;

    vector<edge> enew;
    // cout << "copy original edges" << endl;
    enew.assign(rtree.edges.begin(), rtree.edges.end());

    if(cons){
        // cout << "branches constrained" << endl;
        int count = 0;
        for(int i = 0; i < nedge - 1; ++i){
            if(enew[i].end > rtree.nleaf - 1){
                // cout << "count " << count << endl;
                enew[i].length = gsl_vector_get(blens, count);
                count++;
            }else{
                enew[i].length = 0;
            }
        }

        // Time to first sample
        // cout << "time to 1st sample: " << total_time << endl;
        evo_tree new_tree(rtree.nleaf, enew, get_total_time(rtree.get_node_times(), max_tobs));
        new_tree.mu = rtree.mu;
        new_tree.dup_rate = rtree.dup_rate;
        new_tree.del_rate = rtree.del_rate;
        new_tree.chr_gain_rate = rtree.chr_gain_rate;
        new_tree.chr_loss_rate = rtree.chr_loss_rate;
        new_tree.wgd_rate = rtree.wgd_rate;

        return new_tree;
    }else{
        // cout << "branches unconstrained" << endl;
        for(int i = 0; i < nedge - 1; ++i){
            enew[i].length = gsl_vector_get(blens,i);
        }

        evo_tree new_tree(rtree.nleaf, enew);
        new_tree.mu = rtree.mu;
        new_tree.dup_rate = rtree.dup_rate;
        new_tree.del_rate = rtree.del_rate;
        new_tree.chr_gain_rate = rtree.chr_gain_rate;
        new_tree.chr_loss_rate = rtree.chr_loss_rate;
        new_tree.wgd_rate = rtree.wgd_rate;

        return new_tree;
    }
}


// Adjust the time of sample1 so that it is not after all other tips. Then the subsequent adjustment will be addition, not introducing negative values
void adjust_sample1(evo_tree& rtree, const vector<double>& tobs, int sample1){
    // Find the maximum time differences to the first sample
    double max_diff = 0.0;
    for(int i = 0; i < tobs.size(); i++){
        // current differences to first sample
        double delta = rtree.nodes[i].time - rtree.nodes[sample1].time;
        // If delta is negative, there is no need to adjust node time of first sample
        if(delta > max_diff){
            max_diff = delta;
        }
    }
    // cout << "Sample 1 time before " << rtree.nodes[sample1].time << endl;
    if(max_diff > 0){
        // double increament = max_diff - max_tobs;
        rtree.nodes[sample1].time =  rtree.nodes[sample1].time + max_diff;
        update_edge_len(rtree, sample1);
    }
    // cout << "Sample 1 time after " << rtree.nodes[sample1].time << endl;
}


// Change the time and branch length for one tip
void adjust_one_tip(evo_tree& rtree, const vector<double>& tobs, int i, int sample1){
    int debug = 0;
    int nedge = 2 * rtree.nleaf - 2;

    for(int j = 0; j < nedge; j++){
        if (rtree.edges[j].end == i){
            if(debug){
                cout << "Adjust sample " << i << endl;
                cout << "Sample time before " << rtree.nodes[i].time << endl;
            }
            rtree.nodes[i].time =  rtree.nodes[sample1].time + tobs[i];
            if(debug) cout << "Sample time after " << rtree.nodes[i].time << endl;
            rtree.edges[j].length = rtree.nodes[i].time - rtree.nodes[rtree.edges[j].start].time;

            break;
        }
    }
}


// Adjust all tips recursively so that tree height is smaller than age
void adjust_all_tips(evo_tree& rtree, const double& max_tobs, int age){
    assert(is_blen_valid(rtree));

    double ttime = get_total_time(rtree.get_node_times(), max_tobs);
    if(ttime < age) return;

    int debug = 0;
    if(debug){
        // cout << "Age of patient " << age << endl;
        cout << "Time until first sample before " << ttime << endl;
    }

    // Reduce to_reduce from all tips
    // Add BLEN_MIN keep some distance from age to ensure minimal branch length after reduction
    double to_reduce = ttime - age + BLEN_MIN;

    // Reduce some amount from terminal branches to keep all of them still positive
    vector<double> tip_blens;
    int nedge = 2 * rtree.nleaf - 2;

    for(int j = 0; j < nedge; j++){
       if (rtree.edges[j].end < rtree.nleaf - 1){
           tip_blens.push_back(rtree.edges[j].length);
       }
    }

    double min_tblen = *min_element(tip_blens.begin(), tip_blens.end());

    if(debug){
        cout << "Differences to reduce " << to_reduce << endl;
        cout << "Minimal terminal branch length " << min_tblen << endl;
    }

    if(min_tblen >= to_reduce + BLEN_MIN){
        if(debug) cout << "Only need to adjust tips " << endl;
        for(int j = 0; j < nedge; j++){
           if (rtree.edges[j].end < rtree.nleaf - 1){
               // Only need to update the time of tips
               rtree.nodes[rtree.edges[j].end].time =  rtree.nodes[rtree.edges[j].end].time - to_reduce;
               rtree.edges[j].length = rtree.nodes[rtree.edges[j].end].time - rtree.nodes[rtree.edges[j].start].time;
           }
        }
    }
    else{   // min_tblen < to_reduce + BLEN_MIN
        // Always need to adjust terminal branch lengths
        vector<double> internal_blens;
        for(int j = 0; j < nedge; j++){
            if (rtree.edges[j].end < rtree.nleaf - 1){
                // Reduce BLEN_MIN so that the shortest terminal branch has length BLEN_MIN after reduction
                rtree.edges[j].length = rtree.edges[j].length - (min_tblen - BLEN_MIN);
            }
            if (rtree.edges[j].end > rtree.nleaf - 1){
               internal_blens.push_back(rtree.edges[j].length);
            }
        }
        // If there are still residuals, find some internal edge with larger length and reduce it
        double to_reduce2 = to_reduce - (min_tblen - BLEN_MIN);
        assert(to_reduce2 > 0);
        if(debug) cout << "Differences to reduce further " << to_reduce2 << endl;
        double max_iblen = *max_element(internal_blens.begin(), internal_blens.end());
        if(max_iblen >= to_reduce2 + BLEN_MIN){
            if(debug){
                cout << "Adjust one internal branch length " << endl;
                cout << "Edge lengths before: ";
                for(int j = 0; j < nedge; j++){
                    cout << "\t" << rtree.edges[j].length;
                }
                cout << endl;
                cout << "Node times before: ";
                for(int i = 0; i < rtree.nodes.size(); ++i){
                    cout << "\t" << rtree.nodes[i].time;
                }
                cout << endl;
            }
            for(int j = 0; j < nedge; j++){
                if (rtree.edges[j].end > rtree.nleaf - 1 && rtree.edges[j].length >= to_reduce2 + BLEN_MIN ){
                    if(debug) cout << "Adjusting edge " << j << endl;
                    rtree.edges[j].length = rtree.edges[j].length - to_reduce2;
                    break;
                }
            }
            if(debug){
                cout << "Edge lengths after: ";
                for(int j = 0; j < nedge; j++){
                    cout << "\t" << rtree.edges[j].length;
                }
                cout << endl;
            }
        }else{  // max_iblen < to_reduce2 + BLEN_MIN, rarely needed
            if(debug) cout << "Adjust several internal branch lengths " << endl;
            // All internal branch lengths are smaller than to_reduce2 + BLEN_MIN, so reduction has to be done multiple times
            double delta = to_reduce2;
            for(int j = 0; j < nedge; j++){
                if (rtree.edges[j].end > rtree.nleaf - 1 && rtree.edges[j].length > 2 * BLEN_MIN){
                    if(delta >= (rtree.edges[j].length - BLEN_MIN)){   // need another reduction
                        rtree.edges[j].length = BLEN_MIN;
                        delta = delta - (rtree.edges[j].length - BLEN_MIN);
                    }
                    else{   // delta < (rtree.edges[j].length - BLEN_MIN
                        rtree.edges[j].length = rtree.edges[j].length - delta;
                        delta = -1;
                    }
                    if(delta <= 0) break;
                }
            }
        }

        rtree.calculate_node_times();
        if(debug){
            cout << "Node times after: ";
            for(int i = 0; i < rtree.nodes.size(); ++i){
                cout << "\t" << rtree.nodes[i].time;
            }
            cout << endl;
        }
    }
    assert(is_blen_valid(rtree));

    if(debug){
        // cout << "Age of patient " << age << endl;
        cout << "Time until first sample after " << get_total_time(rtree.get_node_times(), max_tobs) << endl;
    }
}


void adjust_blen(double& nx, double a, double b){
    // double a = BLEN_MIN;
    // double b = BLEN_MAX;
    // cout << "before adjusting " << nx << endl;
    double n, e;
    if(nx < a){
        e = a - nx;
        n = floor(e / (b - a));
        nx = a + (e - n * (b - a));
        // cout << "after adjusting lower bound " << nx << endl;
        assert(nx > a);
    }
    if(nx > b){
        e = nx - b;
        n = floor(e / (b - a));
        nx = b - (e - n * (b - a));
        // cout << "after adjusting upper bound " << nx << endl;
        assert(nx < b);
    }
    // cout << "after adjusting " << nx << endl;
    assert(nx >= a && nx <= b);
}


// Seem not work well
bool is_tip_valid(const evo_tree& rtree, const vector<double>& tobs, int sample1){
    int debug = 0;
    for(int i = 0; i < tobs.size(); i++){
        // current differences to first sample
        double delta = rtree.nodes[i].time - rtree.nodes[sample1].time;
        if(fabs(delta - tobs[i]) > SMALL_VAL){
            if(debug){
                cout << delta << endl;
                cout << tobs[i] << endl;
                cout << fabs(delta - tobs[i]) << endl;
                cout << "tip " << i << " is not valid" << endl;
            }
            return false;
        }
    }
    return true;
}

// Adjust the terminal branch lengths of a tree with fixed topolgy when the constraints are violated
// The initial or perturbed tree may not have tip nodes in accord with the sampling time information
// Ensure edge length keep positive
void adjust_tree_tips(evo_tree& rtree, const vector<double>& tobs, int age){
    assert(is_blen_valid(rtree));

    int debug = 0;
    // find the first sample
    int sample1 = 0;
    for(int i = 0; i < tobs.size(); i++){
        //
        if(tobs[i] == 0){
            sample1 = i;
            break;
        }
    }
    if(debug){
        cout << "A first sample is " << sample1 << endl;
    }

    // Check the node times for tip nodes
    adjust_sample1(rtree, tobs, sample1);
    for(int i = 0; i < tobs.size(); i++){
        if(i == sample1) continue;
        if(debug) cout << "Adjusting tip " << i << endl;
        adjust_one_tip(rtree, tobs, i, sample1);
    }
    double max_tobs = *max_element(tobs.begin(), tobs.end());
    adjust_all_tips(rtree, max_tobs, age);
}


// The initial branch lengths are always positive.
// But after optimization, some branch lengths may become negative.
// Increase all branch lengths with the same value to keep tip timing difference
// The result tree may break the height constraint
void adjust_tree_blens_all(evo_tree& rtree){
    int debug = 0;
    vector<double> blens;
    int nedge = 2 * rtree.nleaf - 2;

    for(int i = 0; i < nedge; ++i){
        blens.push_back(rtree.edges[i].length);
    }

    double min_blen = *min_element(blens.begin(), blens.end());
    if(min_blen < 0){
        double delta = min_blen - BLEN_MIN;
        if(debug){
            cout << "Estimated branch lengths have negative values!" << endl;
            cout << "Increase all branch lengths (except the branch to normal node) to eliminate negative values by " << -delta << endl;
        }
        for(int i = 0; i < blens.size(); i++){
            if(rtree.edges[i].start == rtree.nleaf && rtree.edges[i].end == rtree.nleaf - 1) continue;
            blens[i] = blens[i] - delta;
            rtree.nodes[rtree.edges[i].end].time = rtree.nodes[rtree.edges[i].end].time - delta;
        }
    }

    assert(blens.size() == nedge);

    for(int i = 0; i < nedge; ++i){
        rtree.edges[i].length = blens[i];
    }
}


// Only change negative branches. The tip differences will be adjusted later
void adjust_tree_blens(evo_tree& rtree){
    int debug = 0;
    int nedge = 2 * rtree.nleaf - 2;

    for(int i = 0; i < nedge; ++i){
        if(rtree.edges[i].length < 0){
            if(debug) cout << "Negative length for edge " << i << endl;
            rtree.edges[i].length = BLEN_MIN;
            rtree.nodes[rtree.edges[i].end].time = rtree.nodes[rtree.edges[i].start].time + BLEN_MIN;
        }
    }
}


// Scale the tree so that total tree height is in the lifetime of the patient
// This may cause tips nodes violating the given sampling time
void adjust_tree_height(evo_tree& rtree, gsl_rng* r, double min_height, double max_height, int scale){
    int debug = 0;
    double tree_height = get_tree_height(rtree.get_node_times());

    // cout << "ttime " << ttime << endl;
    if(tree_height > max_height){
        // Use smaller height to avoid new violations after tip adjustment
        double max_height_scaled = (max_height / scale > min_height) ? max_height / scale : min_height + scale;
        double new_tree_height = runiform(r, min_height, max_height_scaled);
        double ratio = new_tree_height / tree_height;
        if(debug){
            cout << "Estimated tree height larger than age at last sample!" << endl;
            // randomly change tree height until first sample to satifisty the constraints
            rtree.print();
            cout << "Adjusting tree height with ratio " << ratio << endl;
        }
        rtree.scale_time(ratio);
        if(debug){
            rtree.print();
        }
    }
    assert(is_blen_valid(rtree));
}


// Check whether the set of branch lengths is valid under the time constraints
bool is_tree_valid(evo_tree& rtree, const double& max_tobs, int age, int cons){
    for(int i = 0; i < rtree.edges.size(); ++i){
        if(rtree.edges[i].length < 0){
            cout << "Negative branch length!" << endl;
            return false;
        }
    }
    if(cons){
        double tot = get_total_time(rtree.get_node_times(), max_tobs);
        if(floor(tot) > age){
            cout << "Wrong total time (larger than age)!" << "\t" << tot << "\t" << age << endl;
            return false;
        }
    }
    return true;
}

// Save all branch lengths in a vector for optimization
// Assume lenvec is the size of all branches, using ID of end node to distinguish each edge
void save_branch_lengths(evo_tree& rtree, DoubleVector &lenvec, int startid, Node* node, Node* dad) {
    int debug = 0;

    if(!node){
        node = &(rtree.nodes[rtree.root_node_id]);   // root
        int branchNum = rtree.edges.size();
        if(lenvec.empty()) lenvec.resize(branchNum);

        if(debug){
            rtree.print();
            cout << "Root is " << node->id + 1 << " with " << branchNum << " edges "  << node->neighbors.size() << " neighbors " << endl;
        }
    }

    if(debug){
        cout << "saving neighboring branch lengths of node " << node->id + 1 << " to a vector" << endl;
        cout << "initial size of all branches:";
        for(int i = 0; i < lenvec.size(); i++){
            cout << "\t" << lenvec[i];
        }
        cout << endl;
        cout << " traversing neighbors of node " << node->id + 1 << endl;
    }

    FOR_NEIGHBOR_IT(node, dad, it){
        // cout << " neighbor " <<  (*it)->id + 1 << " with length " <<  (*it)->length << endl;
    	(*it)->getLength(lenvec, (*it)->id);
        // cout  << " length of neighbor " << (*it)->id + 1 << " is " << lenvec[(*it)->id] << endl;
    	save_branch_lengths(rtree, lenvec, startid, (Node* )(*it)->node, node);
    }

    if(debug){
        cout << "size of all branches:";
        for(int i = 0; i < lenvec.size(); i++){
            cout << "\t" << lenvec[i];
        }
        cout << endl;
    }
}


// Restoring branch lengths from a vector
void restore_branch_lengths(evo_tree& rtree, DoubleVector &lenvec, int startid, Node* node, Node* dad) {
    int debug = 0;
    if(debug) cout << "\nrestoring branch lengths from a vector" << endl;

    if(!node){
        node = &(rtree.nodes[rtree.root_node_id]);   // root
        assert(lenvec.size() == rtree.edges.size());
    }

    FOR_NEIGHBOR_IT(node, dad, it){
        int idx = (*it)->id + startid;   // edge id

    	(*it)->setLength(lenvec, idx);
        (*it)->node->findNeighbor(node)->setLength(lenvec, idx);

        rtree.edges[idx].length = lenvec[idx];

        if(debug) cout << "size of branch " << (*it)->node->id + 1 << ", " << (*it)->node->findNeighbor(node)->id + 1 << " after restoring is " << lenvec[(*it)->id + startid] << endl;

    	restore_branch_lengths(rtree, lenvec, startid, (Node* )(*it)->node, node);
    }
}


void save_mutation_rates(const evo_tree& rtree, DoubleVector& muvec){
    muvec.push_back(rtree.mu);
    muvec.push_back(rtree.dup_rate);
    muvec.push_back(rtree.del_rate);
    muvec.push_back(rtree.chr_gain_rate);
    muvec.push_back(rtree.chr_loss_rate);
    muvec.push_back(rtree.wgd_rate);
}


void restore_mutation_rates(evo_tree& rtree, const DoubleVector& muvec){
    rtree.mu = muvec[0];
    rtree.dup_rate = muvec[1];
    rtree.del_rate = muvec[2];
    rtree.chr_gain_rate = muvec[3];
    rtree.chr_loss_rate = muvec[4];
    rtree.wgd_rate = muvec[5];
}


// TODO: Build parsimony tree from copy number changes (breakpoints)
evo_tree build_parsimony_tree(int Ns, vector<vector<int>>& data){
    vector<int> nodes;
    vector<edge> edges;

    evo_tree ptree = evo_tree(Ns + 1, edges);

    return ptree;
}


evo_tree read_tree_info(const string& filename, const int& Ns, int debug){
  if(debug) cout << "\tread_tree_info" << endl;

  vector<edge> edges;
  int id = 0;

  ifstream infile(filename.c_str());
  if(infile.is_open()){
    std::string line;
    while(!getline(infile, line).eof()){
      if(line.empty()) continue;

      std::vector<std::string> split;
      std::string buf;
      stringstream ss(line);
      while(ss >> buf) split.push_back(buf);

      if(id > 0){
      	int start = atoi(split[0].c_str());
      	int end = atoi(split[1].c_str());
      	double length = atof(split[2].c_str());
        if(end == Ns + 1) length = 0;
        if( !(length > 0) && end != Ns + 1 ){
          cout << "Invalid branch length: " << id << "\t" << start << "\t" << end << "\t" << length << endl;
          exit(EXIT_FAILURE);
          // length = 1;        
        }

      	edges.push_back(edge(id - 1, start - 1, end - 1, length));
      }
      id++;
    }
  }else{
    std::cerr << "Error: open of tree data unsuccessful: " <<  filename << std::endl;
    exit(EXIT_FAILURE);
  }

  evo_tree new_tree(Ns + 1, edges);
  if(debug) new_tree.print();

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
//         exit(EXIT_FAILURE);
//     }
//
//     evo_tree new_tree(Ns + 1, edges);
//     //new_tree.print();
//
//     return new_tree;
// }


void adjust_tip_time(evo_tree& rtree, const vector<double>& tobs, int Ns, int same_tip, int debug){
    for(int i = 0; i < rtree.nodes.size(); i++){
        Node* node = &rtree.nodes[i];
        if(node->id < Ns){
            if(same_tip){
                node->time += node->age;
            }
            node->time += tobs[node->id];
        }
    }
    rtree.calculate_age_from_time();

    // update branch lengths based on node times
    for(int i = 0; i < rtree.edges.size(); i++){
        edge *e = &rtree.edges[i];
        if(e->end < Ns) e->length += tobs[e->end];
    }
}


// Read parsimony trees built by other tools as starting trees, assign timings to tip nodes and initialize mutation rates
evo_tree read_parsimony_tree(const string& tree_file, const int& Ns, const vector<double>& rates, const vector<double>& tobs, gsl_rng* r, int age, int cons){
    int debug = 0;
    if(debug)   cout << "reading from file " << tree_file << endl;

    evo_tree rtree = read_tree_info(tree_file, Ns);
    // The branch lengths in parsimony tree may be very large
    if(debug){
        cout << "original tree: " << rtree.make_newick() << endl;
        rtree.print();
    }

    double old_tree_height = get_tree_height(rtree.get_node_times());
    if(cons && old_tree_height > age){
        double min_height = *max_element(tobs.begin(), tobs.end());
        double max_height = age - min_height;  // allow adding extra time at the tips
        assert(min_height < max_height);
        double tree_height = runiform(r, min_height, max_height);       
        double ratio = tree_height / old_tree_height;
        rtree.scale_time(ratio);
    }

    bool non_zero = std::any_of(tobs.begin(), tobs.end(), [](double i) { return i > 0.0; });
    if(non_zero){
        adjust_tip_time(rtree, tobs, Ns, 1, debug);
    }

    if(debug){
        cout << "adjusted tree: " << rtree.make_newick() << endl;
        rtree.print();
    }

    restore_mutation_rates(rtree, rates);

    return rtree;
}


// vector<double> get_blens_from_intervals(evo_tree& rtree, double *x){
//     vector<double> intervals;
//     // The estimated value may be nan
//     for(int i = 0; i < rtree.nleaf - 1; i++){
//         double len = x[i + 1];
//         bool is_nan = std::isnan(len);
//         if ( is_nan || (!is_nan && len < BLEN_MIN)){
//            len = BLEN_MIN;
//         }
//         intervals.push_back(len);
//     }
//     // vector<double> intervals(x, x + sizeof x / sizeof x[0]);
//     if(debug){
//         rtree.print();
//         cout << "Current length of intervals: " << endl;
//         for(int i = 0; i < intervals.size(); i++){
//             cout << i + 1 << "\t" << "\t" << intervals[i] << endl;
//         }
//         cout << "Corresponding " << rtree.nleaf << " nodes: ";
//         for(int i = 0; i < rtree.top_tnodes.size(); i++){
//             cout << rtree.top_tnodes[i] + 1 << "\t" << "Prevous node time " << rtree.nodes[rtree.top_tnodes[i]].time << endl;
//         }
//     }
//
//     vector<double> blens = rtree.get_edges_from_interval(intervals, rtree.top_tnodes);
//
//     return blens;
// }


// Print branch lengths from time intervals
// void print_invl_blen(evo_tree& new_tree, string header){
//     cout << header << endl;
//     new_tree.print();
//
//     cout << "\tTop " << new_tree.nleaf - 1 << " time intervals: ";
//     for(int i = 0; i < new_tree.top_tinvls.size(); i++){
//         cout << "\t" << new_tree.top_tinvls[i];
//     }
//     cout<< endl;
//     cout << "\tCorresponding " << rtree.nleaf - 1+ 1 << " nodes: ";
//     for(int i = 0; i < new_tree.top_tnodes.size(); i++){
//         cout << new_tree.top_tnodes[i] + 1 << "\t" << "\t" << new_tree.node_times[new_tree.top_tnodes[i]] << endl;
//     }
//
//     vector<double> blens = new_tree.get_edges_from_interval(new_tree.top_tinvls, new_tree.top_tnodes);
//     cout << "\tNew branch lengths: ";
//     for(int i = 0; i < blens.size(); i++){
//         cout << i + 1 << "\t" << "\t" << blens[i] << endl;
//     }
// }



/*
vector<edge> create_edges_from_nodes(const vector<node>& nodes, const vector<double>& node_times ){
  vector<edge> enew;

  // Add leaf and internal nodes
  int root_id = 0;
  int id = 0;
  for(int i = 0; i < nodes.size(); ++i){
    if(nodes[i].parent == -1) root_id = nodes[i].id;
  }

  for(int i = 0; i < nodes.size(); ++i){
    if(nodes[i].id != root_id && nodes[i].id != root_id - 1 && nodes[i].parent != root_id){
      enew.push_back( edge(id, nodes[i].parent, nodes[i].id, node_times[nodes[i].id] - node_times[nodes[i].parent]) );
      id++;
    }
  }
  // Add root nodes
  for(int i = 0; i < nodes.size(); ++i){
  if( nodes[i].parent == root_id && nodes[i].id != root_id - 1){
    enew.push_back( edge(id, nodes[i].parent, nodes[i].id, node_times[nodes[i].id] - node_times[nodes[i].parent]) );
    id++;
  }
  }
  enew.push_back( edge(id, root_id, Ns, 0) );

  return enew;
}
*/
