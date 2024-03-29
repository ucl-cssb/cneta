#include "state.hpp"


void print_tree_state(const evo_tree& rtree, const vector<vector<int>>& S_sk_k, int nstate){
    cout << "\nStates so far:\n";

    int ntotn = 2 * rtree.nleaf - 1;
    for(int i = 0; i < ntotn; ++i){
        for(int j = 0; j < nstate; ++j){
          cout << "\t" << S_sk_k[i][j];
        }
        cout << endl;
    }
}


void extract_tree_ancestral_state(const evo_tree& rtree, const vector<int>& knodes, const set<vector<int>>& comps, const vector<vector<double>>& L_sk_k, const vector<vector<int>>& S_sk_k, int model, int cn_max, int is_total, int m_max, map<int, int> &asr_states){
    int debug = 0;

    if(debug){
        cout << "Get most likely joint estimation of ancestral nodes" << endl;
        // No computation on root
        // cout << "likelihood at root:";
        // for(int j = 0; j < L_sk_k[rtree.nleaf].size(); ++j){
        //   cout << "\t" << L_sk_k[rtree.nleaf][j];
        // }
        // cout << endl;
    }

    int parent_state;
    if(model == BOUNDA){
        parent_state = NORM_ALLElE_STATE;
    }else if(model == DECOMP){
        int k = 0;
        for (auto v : comps){
            // if(debug) cout << k << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << "\t" << v[3] << "\t" << v[4] << endl;
            bool zeros = all_of(v.begin(), v.end(), [](int i){ return i == 0; });
            if(zeros){
                break;
            }
            k++;
        }
        parent_state = k;
    }else{
        parent_state = NORM_PLOIDY;
    }
    asr_states[rtree.nleaf] = parent_state;   // for root

    // Traverse the tree from root to tips, need to know the parent of each node
    for(int i = knodes.size() - 2; i >= 0; i--){  // starting from node ID for MRCA
        int nid = knodes[i];
        // Find the parent node of current node
        int parent = rtree.nodes[nid].parent;
        if(asr_states.find(parent) == asr_states.end()){
            cout << "Cannot find state for the parent of node " << nid + 1 << ", " << parent + 1 << endl;
            exit(EXIT_FAILURE);
        }
        parent_state = asr_states[parent];

        int state = S_sk_k[nid][parent_state];
        asr_states[nid] = state;

        if(debug){
            cout << "\t\tnode " << nid + 1 << " with state " << state << " and parent " << parent + 1 << " whose state is " << parent_state << endl;
        }
    }

}


void initialize_asr_table(const vector<int>& obs, const evo_tree& rtree, const vector<double>& blens, const vector<double*>& pmat_per_blen, vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, int model, int nstate, int is_total){
    int debug = 0;
    if(debug){
        cout << "Initializing tables for reconstructing joint ancestral state" << endl;
        cout << "branch lengths so far:";
        for(auto b : blens){
            cout << "\t" << b;
        }
        cout << endl;
    }

    for(int i = 0; i < rtree.nleaf - 1; ++i){
        // cout << "node " << i + 1  << ", observed CN " << obs[i] << endl;
        // Find the parent node
        int parent = rtree.edges[rtree.nodes[i].e_in].start;
        double blen = rtree.edges[rtree.nodes[i].e_in].length;
        // cout << "parent " << parent + 1 << ", blen " << blen << endl;

        auto pi = equal_range(blens.begin(), blens.end(), blen);
        // assert(distance(pi.first, pi.second) == 1);
        int idx_blen = distance(blens.begin(), pi.first);

        if(debug){
            cout << "Pmatrix for branch length " << blen << " " << blens[idx_blen] << endl;
            r8mat_print(nstate, nstate, pmat_per_blen[idx_blen], "  P matrix:");
        }

        // Find the state(s) of current node
        vector<int> tip_states;
        if(model == BOUNDA){
          if(is_total){
              int si = (obs[i] * (obs[i] + 1)) / 2;
              int ei = si + obs[i];
              for(int k = si; k <= ei; k++){
                 tip_states.push_back(k);
              }
          }else{ // With haplotype-specific copy number, only the specific site needs to be filled
              tip_states.push_back(obs[i]);
          }
        }

        if(debug){
            cout << "There are " << tip_states.size() << " states for copy number " << obs[i] << endl;
        }

        for(int j = 0; j < nstate; ++j){  // For each possible parent state, find the most likely tip states
            vector<double> vec_li(nstate, SMALL_LNL);
            // another loop as there maybe multiple states for a specific total CN
            for(int m = 0; m < tip_states.size(); ++m){
                int k = tip_states[m];
                if(debug){
                    cout << "parent state " << j << ", child state " << k << endl;
                }
                double li = 0.0;
                if(model == MK){
                    li = get_transition_prob(rtree.mu, blen, j, k);
                }else{
                    li = pmat_per_blen[idx_blen][j  + k * nstate];  // assume parent has state j
                }
                // if(li > 0) li = log(li);
                // else li = SMALL_LNL;
                vec_li[k] = li;
            }

            int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
            double max_li = *max_element(vec_li.begin(), vec_li.end());
            assert(max_li == vec_li[max_i]);
            // cout << "for node: i " << i + 1 << ", parent state " << j << ", max state is " << max_i << " with probability " << exp(max_li) << endl;
            S_sk_k[i][j] = max_i;
            L_sk_k[i][j] = max_li;
        }
    }

    if(debug){
      print_lnl_at_tips(rtree, obs, L_sk_k, nstate);

      cout << "\nState vector for tips (should be observed data):\n";
      for(int i = 0; i < rtree.nleaf; ++i){
          for(int j = 0; j < nstate; ++j){
            cout << "\t" << S_sk_k[i][j];
          }
          cout << endl;
      }
    }

}


void initialize_asr_table_decomp(const vector<int>& obs, const evo_tree& rtree, const set<vector<int>>& comps, MAX_DECOMP& max_decomp, PMAT_DECOMP& pmat_decomp, vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, int nstate, int is_total){
    int debug = 0;
    if(debug) cout << "Initializing tables for reconstructing joint ancestral state" << endl;

    int s_wgd, s_chr, s_seg, s_chr2, s_seg2;
    int e_wgd, e_chr, e_seg, e_chr2, e_seg2;
    double prob_wgd, prob_chr, prob_seg, prob_chr2, prob_seg2, prob_chr_all, prob_seg_all;
    int dim_wgd = max_decomp.max_wgd + 1;
    int dim_chr = 2 * max_decomp.max_chr_change + 1;
    int dim_seg = 2 * max_decomp.max_site_change + 1;
    int m_max = max_decomp.m_max;

    int Ns = rtree.nleaf - 1;

    for(int i = 0; i < Ns; ++i){
        // cout << "node " << i + 1 << endl;
        // Find the parent node
        int parent = rtree.edges[rtree.nodes[i].e_in].start;
        double bli = rtree.edges[rtree.nodes[i].e_in].length;
        // cout << "parent " << parent + 1 << endl;
        // cout << "blen " << blen << endl;
        double *pbli_wgd;
        double *pbli_chr;
        double *pbli_seg;
        if(dim_wgd > 1){
            pbli_wgd = pmat_decomp.pmats_wgd[bli];
        }
        if(dim_chr > 1){
            pbli_chr = pmat_decomp.pmats_chr[bli];
        }
        if(dim_seg > 1){
            pbli_seg = pmat_decomp.pmats_seg[bli];
        }

        // Find the state(s) of current node
        vector<int> tip_states;
        int k = 0;
        for (auto c : comps){
            for(int m1 = 0; m1 <= m_max; m1++){
                for(int m2 = 0; m2 <= m_max; m2++){
                    int sum = pow(2, c[0] + 1) + m1 * c[1] + c[2] + 2 * m2 * c[3] + 2 * c[4];
                    if(sum==obs[i]){
                        tip_states.push_back(k);
                    }
                }
            }
            k++;
        }
        if(debug) cout << "There are " << tip_states.size() << " states for copy number " << obs[i] << endl;

        for(int j = 0; j < nstate; ++j){  // For each possible parent state, find the most likely tip states
            set<vector<int>>::iterator iter = comps.begin();
            // It will move forward the passed iterator by passed value
            std::advance(iter, j);
            vector<int> s = *iter;
            // get_decomposition(sk, s_wgd, s_chr, s_seg, decomp_table);
            s_wgd = s[0];
            s_chr = s[1];
            s_seg = s[2];
            s_chr2 = s[3];
            s_seg2 = s[4];
            // The indices for chromosome and segment matrix have to be ajusted
            int delta_chr = (dim_chr - 1)/2;
            int delta_seg = (dim_seg - 1)/2;
            if(debug){
                cout << "Starting state " << j << "\t" << s_wgd << "\t" << s_chr << "\t" << s_seg << "\t" << s_chr2 << "\t" << s_seg2 << "\n";
                cout << "   Offset for chr and segment matrices " << delta_chr << "\t" << delta_seg << "\n";
            }

            vector<double> vec_li(nstate, SMALL_LNL);
            for(int m = 0; m < tip_states.size(); ++m){
                int k = tip_states[m];
                if(debug) cout << "parent state " << j << ", child state " << k << endl;
                double li = 0;
                prob_wgd = 1;
                prob_chr = 1;
                prob_seg = 1;
                prob_chr2 = 1;
                prob_seg2 = 1;
                prob_chr_all = 1;
                prob_seg_all = 1;

                iter = comps.begin();
                std::advance(iter, k);
                vector<int> e = *iter;
                // get_decomposition(sk, s_wgd, s_chr, s_seg, decomp_table);
                e_wgd = e[0];
                e_chr = e[1];
                e_seg = e[2];
                e_chr2 = e[3];
                e_seg2 = e[4];
                if(dim_wgd > 1) prob_wgd = pbli_wgd[s_wgd + e_wgd * dim_wgd];
                if(dim_chr > 1){
                    prob_chr = pbli_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
                    prob_chr2 = pbli_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
                    prob_chr_all = prob_chr * prob_chr2;
                }
                if(dim_seg > 1){
                    prob_seg = pbli_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
                    prob_seg2 = pbli_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
                    prob_seg_all = prob_seg * prob_seg2;
                }
                li = prob_wgd * prob_chr_all * prob_seg_all;
                if(debug){
                    cout << "End state " << k << "\t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2  << "\n";
                    cout << "Prob for each event " << "\t" << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << li << "\n";
                }

                // if(li > 0) li = log(li);
                // else li = SMALL_LNL;
                vec_li[k] = li;
            }

            int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
            double max_li = *max_element(vec_li.begin(), vec_li.end());
            assert(max_li == vec_li[max_i]);
            // cout << "for node: i " << i + 1 << ", parent state " << j << ", max state is " << max_i << " with probability " << exp(max_li) << endl;
            S_sk_k[i][j] = max_i;
            L_sk_k[i][j] = max_li;
        }
    }

    if(debug){
      cout << "\nCNs at tips:\n";
      for(int i = 0; i < Ns; ++i){
          cout<< "\t" << obs[i];
      }
      cout << endl;
      cout << "\nLikelihood for tips:\n";
      for(int i = 0; i < rtree.nleaf; ++i){
          for(int j = 0; j < nstate; ++j){
            cout << "\t" << L_sk_k[i][j];
          }
          cout << endl;
      }
      cout << "\nState vector for tips:\n";
      for(int i = 0; i < rtree.nleaf; ++i){
          for(int j = 0; j < nstate; ++j){
            cout << "\t" << S_sk_k[i][j];
          }
          cout << endl;
      }
    }

}


double get_max_prob_children(const vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, double* pblen, int k, int nstate, int sp, int ni, int nj, int blen, int model){
    int debug = 0;

    vector<double> vec_li;
    double li = 0.0;
    // loop over possible si for a fixed state of parent (sp)
    for(int si = 0; si < nstate; ++si){
        if(model == MK){
          li = get_transition_prob(rtree.mu, blen, sp, si);
        }else{
          li = pblen[sp + si * nstate];
        }

        // if(debug){
        //     cout << "\tfor state " << si << endl;
        //     cout << "\t\tseparate likelihood "  << li << "\t" << L_sk_k[ni][si] << "\t" << L_sk_k[nj][si] << endl;
        // }

        li = li * L_sk_k[ni][si] * L_sk_k[nj][si];

        // if(debug) cout << "\t\tthe likelihood of the best reconstruction of subtree at " << sp + 1 << " is: " << li << endl;

        vec_li.push_back(li);
    }

    int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
    double max_li = *max_element(vec_li.begin(), vec_li.end());
    assert(max_li == vec_li[max_i]);
    S_sk_k[k][sp] = max_i;

    if(debug){
        cout << "all likelihoods:";
        for(int i = 0; i < vec_li.size(); i++){
            cout << "\t" << vec_li[i];
        }
        cout << endl;
        cout << "for node: k " << k + 1 << ", parent state " << sp << ", max state is " << max_i << " with probability " << max_li << endl;
    }

    return max_li;
}


double get_max_children_decomp2(const vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, int k, int nstate, double* pbli_wgd, double* pbli_chr, double* pbli_seg, DIM_DECOMP& dim_decomp, int sp, int ni, int nj, int blen){
    int debug = 0;
    int s_wgd, s_chr, s_seg, s_chr2, s_seg2;
    int e_wgd, e_chr, e_seg, e_chr2, e_seg2;
    double prob_wgd, prob_chr, prob_seg, prob_chr2, prob_seg2, prob_chr_all, prob_seg_all;
    // Each copy number is decomposed into a set of 3-tuples
    set<vector<int>>::iterator iter = comps.begin();
    std::advance(iter, sp);
    vector<int> s = *iter;
    s_wgd = s[0];
    s_chr = s[1];
    s_seg = s[2];
    s_chr2 = s[3];
    s_seg2 = s[4];
    int dim_wgd = dim_decomp.dim_wgd;
    int dim_chr = dim_decomp.dim_chr;
    int dim_seg = dim_decomp.dim_seg;
    // The indices for chromosome and segment matrix have to be ajusted
    int delta_chr = (dim_chr - 1)/2;
    int delta_seg = (dim_seg - 1)/2;

    if(debug){
        cout << "Starting state " << sp << "\t" << s_wgd << "\t" << s_chr << "\t" << s_seg << "\t" << s_chr2 << "\t" << s_seg2 << "\n";
        cout << "   Offset for chr and segment matrices " << delta_chr << "\t" << delta_seg << "\n";
    }

    vector<double> vec_li;
    double li = 0;
    int si = 0;
    for(auto e : comps){
        prob_wgd = 1;
        prob_chr = 1;
        prob_seg = 1;
        prob_chr2 = 1;
        prob_seg2 = 1;
        prob_chr_all = 1;
        prob_seg_all = 1;

        e_wgd = e[0];
        e_chr = e[1];
        e_seg = e[2];
        e_chr2 = e[3];
        e_seg2 = e[4];
        if(dim_wgd > 1){
          prob_wgd = pbli_wgd[s_wgd + e_wgd * dim_wgd];
        }
        if(dim_chr > 1){
            prob_chr = pbli_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
            prob_chr2 = pbli_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
            prob_chr_all = prob_chr * prob_chr2;
        }
        if(dim_seg > 1){
            prob_seg = pbli_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
            prob_seg2 = pbli_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
            prob_seg_all = prob_seg * prob_seg2;
        }
        li = prob_wgd * prob_chr_all * prob_seg_all;

        if(debug){
            cout << "End state " << si << "\t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2 << "\n";
            cout << "Prob for each event " << "\t" << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << li  << "\n";
        }

        // if(li > 0){
        //     li = log(li);
        // }
        // else li = SMALL_LNL;
        // li += (L_sk_k[ni][si]);
        // li += (L_sk_k[nj][si]);
        li = li * L_sk_k[ni][si] * L_sk_k[nj][si];

        // if(debug) cout << "\t\tscoring: Li\t" << li << endl;
        // if(std::isnan(li) || li < SMALL_LNL) li = SMALL_LNL;
        if(debug) cout << "\t\trevised scoring: Li\t" << li << endl;
        vec_li.push_back(li);

        si++;
     }
    int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
    double max_li = *max_element(vec_li.begin(), vec_li.end());
    assert(max_li == vec_li[max_i]);
    S_sk_k[k][sp] = max_i;

    if(debug){
        cout << "all likelihoods:";
        for(int i = 0; i < vec_li.size(); i++){
            cout << "\t" << vec_li[i];
        }
        cout << endl;
        cout << "for node: k " << k + 1 << ", parent state " << sp << ", max state is " << max_i << " with probability " << exp(max_li) << endl;
    }

    return max_li;
}


void get_ancestral_states_site(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, int nstate, int model){
  int debug = 0;
  if(debug){
      cout << "Getting ancestral state for one site" << endl;
  }

  int Ns = rtree.nleaf - 1;
  for(int kn = 0; kn < knodes.size() - 1; ++kn){   // not include root which is always normal
    int k = knodes[kn];
    int np = rtree.edges[rtree.nodes[k].e_in].start;
    double blen = rtree.edges[rtree.nodes[k].e_in].length;
    int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
    int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;

    auto pi = equal_range(blens.begin(), blens.end(), blen);
    // assert(distance(pi.first, pi.second) == 1);
    int idx_blen = distance(blens.begin(), pi.first);
    double* pblen = pmat_per_blen[idx_blen];

    if(debug) cout << "node:" << np + 1 << " -> " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , "  <<  nj + 1 << " , " << blen << endl;

    // loop over possible values of sk
    if(k == 2 * Ns){    // root node is always normal, for edge (Ns + 1, 2Ns)
        if(debug) cout << "Getting states for node MRCA " << k + 1 << endl;
        int sp = NORM_PLOIDY;
        if(model == BOUNDA) sp = NORM_ALLElE_STATE;
        if(debug) cout << "likelihood for state " << sp << endl;
        L_sk_k[k][sp] = get_max_prob_children(L_sk_k, S_sk_k, rtree, pblen, k, nstate, sp, ni, nj, blen, model);
    }else{
        for(int sp = 0; sp < nstate; ++sp){  // looping over all possible states of its parent
            if(debug) cout << "likelihood for state " << sp << endl;
            L_sk_k[k][sp] = get_max_prob_children(L_sk_k, S_sk_k, rtree, pblen, k, nstate, sp, ni, nj, blen, model);
        }
    }
  }

  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
    print_tree_state(rtree, S_sk_k, nstate);
  }
}


void get_ancestral_states_site_decomp(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const vector<int>& knodes, const set<vector<int>>& comps, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int nstate){
  int debug = 0;
  if(debug){
      cout << "Getting ancestral state for one site under independent Markov chain model" << endl;
      cout << dim_decomp.dim_wgd << "\t" << dim_decomp.dim_chr << "\t"  << dim_decomp.dim_seg << "\t"  << nstate << endl;
  }

  for(int kn = 0; kn < knodes.size() - 1; ++kn){
        int k = knodes[kn];
        int np = rtree.edges[rtree.nodes[k].e_in].start;
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blen = rtree.edges[rtree.nodes[k].e_in].length;

        double *pbli_wgd;
        double *pbli_chr;
        double *pbli_seg;

        if(dim_decomp.dim_wgd > 1){
            pbli_wgd = pmat_decomp.pmats_wgd[blen];
        }
        if(dim_decomp.dim_chr > 1){
            pbli_chr = pmat_decomp.pmats_chr[blen];
        }
        if(dim_decomp.dim_seg > 1){
            pbli_seg = pmat_decomp.pmats_seg[blen];
        }

        if(debug) cout << "node:" << np + 1 << " -> " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , "  <<  nj + 1 << " , " << blen << endl;

        // loop over possible observed states of start nodes
        if(k == rtree.nleaf){    // root node is always normal
            if(debug) cout << "Getting states for node MRCA " << k + 1 << endl;
            int sp = 0;
            for(auto v : comps){
                bool zeros = all_of(v.begin(), v.end(), [](int i) { return i == 0; });
                if(zeros){
                    L_sk_k[k][sp] = get_max_children_decomp2(L_sk_k, S_sk_k, rtree, comps, k, nstate, pbli_wgd, pbli_chr, pbli_seg, dim_decomp, sp, ni, nj, blen);
                    break;
                }
                sp++;
            }
        }
        else{
            for(int sp = 0; sp < nstate; ++sp){
                L_sk_k[k][sp] = get_max_children_decomp2(L_sk_k, S_sk_k, rtree, comps, k, nstate, pbli_wgd, pbli_chr, pbli_seg,  dim_decomp, sp, ni, nj, blen);
            }
        }
  }
  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
    print_tree_state(rtree, S_sk_k, nstate);
  }
}



void set_pmat(const evo_tree& rtree, int Ns, int nstate, int model, int cn_max, const vector<int>& knodes, vector<double>& blens, vector<double*>& pmat_per_blen, ofstream& fout){
  int debug = 0;

  int dim_mat = nstate * nstate;
  double *qmat = new double[dim_mat];
  memset(qmat, 0.0, dim_mat*sizeof(double));

  assert(model > 0);
  if(debug){
      cout << "Getting rate matrix" << endl;
  }
  if(model == BOUNDA){
      get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
  }else{
      get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
  }

  for(int kn = 0; kn < knodes.size(); ++kn){
    int k = knodes[kn];
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

    if(find(blens.begin(), blens.end(), bli) == blens.end()){
      double *pmati = new double[dim_mat];
      memset(pmati, 0.0, dim_mat*sizeof(double));
      get_transition_matrix_bounded(qmat, pmati, bli, nstate);
      pmat_per_blen.push_back(pmati);
      blens.push_back(bli);
    }
    if(find(blens.begin(), blens.end(), blj) == blens.end()){
      double *pmatj = new double[dim_mat];
      memset(pmatj, 0.0, dim_mat*sizeof(double));
      get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
      pmat_per_blen.push_back(pmatj);
      blens.push_back(blj);
    }
  }

  auto p = sort_permutation(blens, [&](const double& a, const double& b){ return a < b; });
  blens = apply_permutation(blens, p);
  pmat_per_blen = apply_permutation(pmat_per_blen, p);

  if(debug){
      cout << "branch lengths so far:";
      for(auto b : blens){
          cout << "\t" << b;
      }
      cout << endl;

      for(int i = 0; i < pmat_per_blen.size(); ++i){
          double blen = blens[i];
          cout << "Pmatrix for branch length " << blen << endl;
          r8mat_print(nstate, nstate, pmat_per_blen[i], "  P matrix:");
      }
  }

  delete [] qmat;
}


void set_pmat_decomp(const evo_tree& rtree, MAX_DECOMP& max_decomp, int nstate, const vector<int>& knodes, DIM_DECOMP& dim_decomp, PMAT_DECOMP& pmat_decomp, ofstream& fout){
    int debug = 0;

    string header = "node\tsite\tcn\tstate\tmax probability";
    for(int i = 0; i < nstate; i++){
        header += "\tprobablity " + to_string(i + 1);
    }
    fout << header << endl;

    // For WGD model
    int max_wgd = max_decomp.max_wgd;
    int max_chr_change = max_decomp.max_chr_change;
    int max_site_change = max_decomp.max_site_change;

    int dim_wgd = max_wgd + 1;
    int dim_chr = 2 * max_chr_change + 1;
    int dim_seg = 2 * max_site_change + 1;
    dim_decomp = {dim_wgd, dim_chr, dim_seg};

    double *qmat_wgd, *qmat_chr, *qmat_seg;
    double *pmat_wgd, *pmat_chr, *pmat_seg;
    // Find the transition probability matrix for each branch
    map<double, double*> pmats_wgd;
    map<double, double*> pmats_chr;
    map<double, double*> pmats_seg;
    double *pmati_wgd, *pmatj_wgd;
    double *pmati_chr, *pmatj_chr;
    double *pmati_seg, *pmatj_seg;

    if(max_wgd > 0){
        qmat_wgd = new double[dim_wgd * dim_wgd];  // WGD
        memset(qmat_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
        get_rate_matrix_wgd(qmat_wgd, rtree.wgd_rate, max_wgd);
    }

    if(max_chr_change > 0){
        qmat_chr = new double[(dim_chr)*(dim_chr)];   // chromosome gain/loss
        memset(qmat_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
        get_rate_matrix_chr_change(qmat_chr, rtree.chr_gain_rate, rtree.chr_loss_rate, max_chr_change);
    }

    if(max_site_change > 0){
        qmat_seg = new double[(dim_seg)*(dim_seg)];  // site duplication/deletion
        memset(qmat_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
        get_rate_matrix_site_change(qmat_seg, rtree.dup_rate, rtree.del_rate, max_site_change);
    }

    if(debug){
          cout << "Dimension of Qmat " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << "\n";
    }

  for(int kn = 0; kn < knodes.size(); ++kn){
        int k = knodes[kn];
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

        if(max_wgd > 0){
            if(pmats_wgd.find(bli) == pmats_wgd.end()){
                pmati_wgd = new double[(dim_wgd)*(dim_wgd)];
                memset(pmati_wgd, 0.0, (dim_wgd)*(dim_wgd)*sizeof(double));
                get_transition_matrix_bounded(qmat_wgd, pmati_wgd, bli, dim_wgd);
                pmats_wgd[bli] = pmati_wgd;
            }
            if(pmats_wgd.find(blj) == pmats_wgd.end()){
                pmatj_wgd = new double[(dim_wgd)*(dim_wgd)];
                memset(pmatj_wgd, 0.0, (dim_wgd)*(dim_wgd)*sizeof(double));
                get_transition_matrix_bounded(qmat_wgd, pmatj_wgd, blj, dim_wgd);
                pmats_wgd[blj] = pmatj_wgd;
            }
        }

        // For chr gain/loss
        if(max_chr_change > 0){
            if(pmats_chr.find(bli) == pmats_chr.end()){
                pmati_chr = new double[(dim_chr)*(dim_chr)];
                memset(pmati_chr, 0.0, (dim_chr)*(dim_chr)*sizeof(double));
                get_transition_matrix_bounded(qmat_chr, pmati_chr, bli, dim_chr);
                pmats_chr[bli] = pmati_chr;
            }
            if(pmats_chr.find(blj) == pmats_chr.end()){
                pmatj_chr = new double[(dim_chr)*(dim_chr)];
                memset(pmatj_chr, 0.0, (dim_chr)*(dim_chr)*sizeof(double));
                get_transition_matrix_bounded(qmat_chr, pmatj_chr, blj, dim_chr);
                pmats_chr[blj] = pmatj_chr;
            }
        }

        // For site duplication/deletion
        if(max_site_change > 0){
            if(pmats_seg.find(bli) == pmats_seg.end()){
                pmati_seg = new double[(dim_seg)*(dim_seg)];
                memset(pmati_seg, 0.0, (dim_seg)*(dim_seg)*sizeof(double));
                get_transition_matrix_bounded(qmat_seg, pmati_seg, bli, dim_seg);
                pmats_seg[bli] = pmati_seg;
            }
            if(pmats_seg.find(blj) == pmats_seg.end()){
                pmatj_seg = new double[(dim_seg)*(dim_seg)];
                memset(pmatj_seg, 0.0, (dim_seg)*(dim_seg)*sizeof(double));
                get_transition_matrix_bounded(qmat_seg, pmatj_seg, blj, dim_seg);
                pmats_seg[blj] = pmatj_seg;
            }
        }
    }
    if(debug){
        for(auto it = pmats_wgd.begin(); it != pmats_wgd.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print(dim_wgd, dim_wgd, it->second, "  P-WGD matrix:");
        }
        for(auto it = pmats_chr.begin(); it != pmats_chr.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print(dim_chr, dim_chr, it->second, "  P-CHR matrix:");
        }
        for(auto it = pmats_seg.begin(); it != pmats_seg.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print(dim_seg, dim_seg, it->second, "  P-SEG matrix:");
        }
    }

    pmat_decomp = {pmats_wgd, pmats_chr, pmats_seg};

    if(max_wgd > 0){
        delete [] qmat_wgd;
    }
    if(max_chr_change > 0){
        delete [] qmat_chr;
    }
    if(max_site_change > 0){
        delete [] qmat_seg;
    }

}


string get_prob_line(const vector<vector<double>>& L_sk_k, int nid, int nchr, int nc, int is_total, int cn_max){
    string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc + 1);
    //    + "\t" + to_string(cn) + "\t" + to_string(max_ln);
    // cout << line;

    // need to convert state probability to cn probability
    if(is_total){
        vector<double> Lsk_cn(cn_max + 1, 0.0);  // aggregate probabilities over the same total CN
        for(int i = 0; i < L_sk_k[nid].size(); i++){
            // cout << "\t" << to_string(L_sk_k[nid][i]);
            int cni = state_to_total_cn(i, cn_max);
            Lsk_cn[cni] += L_sk_k[nid][i];
        }
        // change to relative probability
        double sum_prob = accumulate(Lsk_cn.begin(), Lsk_cn.end(), 0.0);
        for(int i = 0; i < Lsk_cn.size(); i++){
            double rprob = Lsk_cn[i] / sum_prob;
            line += "\t" + to_string(rprob);
        }
    }else{
        double sum_prob = accumulate(L_sk_k[nid].begin(), L_sk_k[nid].end(), 0.0);
        for(int i = 0; i < L_sk_k[nid].size(); i++){
            double rprob = L_sk_k[nid][i] / sum_prob;
            line += "\t" + to_string(rprob);
        }
    }

    return line;
}


void get_site_cnp(const vector<vector<double>>& L_sk_k, int nid, int nchr, int nc, int is_total, int cn_max, copy_number& cnp){
    // need to convert state probability to cn probability
    if(is_total){
        vector<double> Lsk_cn(cn_max + 1, 0.0);
        for(int i = 0; i < L_sk_k[nid].size(); i++){
            int cni = state_to_total_cn(i, cn_max);
            Lsk_cn[cni] += L_sk_k[nid][i];
        }
        // get the most likely total CN
        int tcn = distance(Lsk_cn.begin(), max_element(Lsk_cn.begin(), Lsk_cn.end()));
        cnp[nchr][nc] = tcn;
    }else{
        // get the most likely haplotype-specific CN
        int state = distance(L_sk_k[nid].begin(), max_element(L_sk_k[nid].begin(), L_sk_k[nid].end()));
        cnp[nchr][nc] = state;
    }
}


void print_node_cnp(ofstream& fout, const copy_number& cnp, int nid, int cn_max, int is_total){
    for(auto site_cn : cnp){
        int nchr = site_cn.first;
        for(auto seg: site_cn.second){
            int nc = seg.first;
            // output CN for all sites: nid, chr, seg, cn
            string line_cn = to_string(nid + 1) + "\t" + to_string(nchr) + "\t" + to_string(nc + 1) + "\t";

            if(is_total){
                int tcn = seg.second;
                line_cn += to_string(tcn);
            }else{
                int cnA;
                int cnB;
                int state = seg.second;
                state_to_allele_cn(state, cn_max, cnA, cnB);
                line_cn += to_string(cnA) + "\t" + to_string(cnB);
            }
            fout << line_cn << endl;
        }
    }
}


// Infer the copy number of the MRCA given a tree at a site, assuming independent Markov chains
double reconstruct_marginal_ancestral_state_decomp(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, const vector<int>& knodes, const set<vector<int>>& comps, OBS_DECOMP& obs_decomp, int use_repeat, int infer_wgd, int infer_chr, int cn_max, string ofile, int is_total){
    int debug = 0;
    if(debug) cout << "\treconstruct marginal ancestral state with independent chain model" << endl;

    string ofile_mrca = ofile + ".mrca.state";
    ofstream fout(ofile_mrca);

    int Ns = rtree.nleaf - 1;
    vector<int> cn_mrca; // CNs for MRCA
    int nid = 2 * (Ns + 1) - 2;

    int nstate = comps.size();
    double logL = 0.0;    // for all chromosomes

    PMAT_DECOMP pmat_decomp;
    DIM_DECOMP dim_decomp;
    MAX_DECOMP max_decomp{obs_decomp.m_max, obs_decomp.max_wgd, obs_decomp.max_chr_change, obs_decomp.max_site_change};
    set_pmat_decomp(rtree, max_decomp, nstate, knodes, dim_decomp, pmat_decomp, fout);

    // Use a map to store computed log likelihood
    map<vector<int>, vector<vector<double>>> sites_lnl_map;

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug) cout << "Computing likelihood on Chr " << nchr << endl;
      double site_logL = 0.0;   // log likelihood for all sites on a chromosome
      for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
          // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k;
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  L_sk_k = initialize_lnl_table_decomp(obs, obs_decomp, nchr, rtree, comps, infer_wgd, infer_chr, cn_max, is_total);
                  get_likelihood_site_decomp(L_sk_k, rtree, comps, knodes, pmat_decomp, dim_decomp, cn_max, is_total);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  if(debug) cout << "\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              L_sk_k = initialize_lnl_table_decomp(obs, obs_decomp, nchr, rtree, comps, infer_wgd, infer_chr, cn_max, is_total);
              get_likelihood_site_decomp(L_sk_k, rtree, comps, knodes, pmat_decomp, dim_decomp, cn_max, is_total);
          }

          site_logL += extract_tree_lnl_decomp(L_sk_k, comps, Ns);

          // Get the likelihood table of MRCA node (with largest ID) in the tree from likelihood table
          int state = distance(L_sk_k[nid].begin(), max_element(L_sk_k[nid].begin(), L_sk_k[nid].end()));

          set<vector<int>>::iterator iter = comps.begin();
          // It will move forward the passed iterator by passed value
          std::advance(iter, state);
          vector<int> c = *iter;
          // There may be multiple possible copy numbers given different coeffients for chromosome-level changes
          int cn = pow(2, c[0] + 1) +  c[1] + c[2] + 2 * c[3] + 2 * c[4];

          cn_mrca.push_back(cn);
          // Print the state of MRCA at this site
          string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc + 1) + "\t" + to_string(cn);
          for(int i = 0; i < L_sk_k[nid].size(); i++){
              line += "\t" + to_string(L_sk_k[nid][i]);
          }
          fout << line << endl;

          if(debug){
              cout << "\t\tState with maximum likelihood " << state << " at chr " << nchr << " site " << nc << endl;
              cout << "\t\t" << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << c[3] << "\t" << c[4] << endl;
          }
      }
      logL += site_logL;
      if(debug){
          cout << "\nLikelihood for chromosome " << nchr << " is " << site_logL << endl;
      }
    } // for each chromosome
    if(debug){
        cout << "\nLikelihood without correcting acquisition bias: " << logL << endl;
        cout << "CNs at MRCA is: " << endl;
        for(int i = 0; i < cn_mrca.size(); i++){
            cout << cn_mrca[i] << endl;
        }
    }

    for(auto m : pmat_decomp.pmats_wgd){
        delete [] m.second;
    }
    for(auto m : pmat_decomp.pmats_chr){
        delete [] m.second;
    }
    for(auto m : pmat_decomp.pmats_seg){
        delete [] m.second;
    }

    return logL;
}


double reconstruct_marginal_ancestral_state(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, const vector<int>& knodes, int model, int cn_max, int use_repeat, int is_total, string ofile){
    int debug = 0;
    if(debug) cout << "\treconstruct marginal ancestral state" << endl;

    string ofile_mrca = ofile + ".mrca.state";
    ofstream fout(ofile_mrca);

    string ofile_mrca_cn = ofile + ".mrca.cn";
    ofstream fout_cn(ofile_mrca_cn);

    int Ns = rtree.nleaf - 1;
    int nid = 2 * Ns;    // node ID for MRCA
    // For copy number instantaneous changes
    int nstate = cn_max + 1;
    if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;

    string header="node\tsite";
    if(is_total){
        for(int i = 0; i < cn_max + 1; i++){
            header += "\tprobablity_" + to_string(i);
        }
    }else{
        for(int i = 0; i < nstate; i++){
            header += "\tprobablity_" + to_string(i);
        }
    }
    fout << header << endl;

    // Find the transition probability matrix for each branch
    vector<double> blens;
    vector<double*> pmat_per_blen;

    set_pmat(rtree, Ns, nstate, model, cn_max, knodes, blens, pmat_per_blen, fout);

    double logL = 0.0;    // for all chromosomes
    map<vector<int>, vector<vector<double>>> sites_lnl_map;
    copy_number cn_mrca;  // all CNs for MRCA

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug) cout << "Computing likelihood for chr " << nchr << " with  " << vobs[nchr].size() << " sites" << endl;
      double site_logL = 0.0;   // log likelihood for all sites on a chromosome

      for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
          bool is_repeated = false;

          if(use_repeat){
              if(debug) cout << " Use repeated site patterns on site " << nc << endl;
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  if(debug) cout << "sites first seen" << endl;
                  initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, 0, 0, model, nstate);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  if(debug) cout << "sites repeated on site " << nc << endl;
                  L_sk_k = sites_lnl_map[obs];
                  is_repeated = true;
              }
          }else{
              initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
              get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, 0, 0, model, nstate);
          }

          site_logL += extract_tree_lnl(L_sk_k, Ns, model);

          if(!is_repeated){
            string line = get_prob_line(L_sk_k, nid, nchr, nc, is_total, cn_max);
            // fout << setprecision(dbl::max_digits10) << line << endl;
            fout << line << endl;
          }
          get_site_cnp(L_sk_k, nid, nchr, nc, is_total, cn_max, cn_mrca);
      }

      logL += site_logL;

      if(debug){
          cout << "Site Likelihood for " << nchr << " is "  << site_logL << endl;
      }
    } // for each chromosome

    print_node_cnp(fout_cn, cn_mrca, nid, cn_max, is_total);

    fout.close();
    fout_cn.close();

    if(debug){
        cout << "\nLikelihood without correcting acquisition bias: " << logL << endl;
        // cout << "CNs at MRCA is: " << endl;
        // for(int i = 0; i < cn_mrca.size(); i++){
        //     cout << cn_mrca[i] << endl;
        // }
    }

    for_each(pmat_per_blen.begin(), pmat_per_blen.end(), DeleteObject());

    return logL;
}



// Infer the copy number of all internal nodes given a tree at a site, assuming independent Markov chains
void reconstruct_joint_ancestral_state_decomp(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs,  vector<int>& knodes, const set<vector<int>>& comps, MAX_DECOMP& max_decomp, int use_repeat, int cn_max, string ofile, int is_total){
    int debug = 0;
    if(debug) cout << "\treconstruct joint ancestral state with independent chain model" << endl;

    string ofile_joint = ofile + ".joint.state";
    ofstream fout(ofile_joint);

    int Ns = rtree.nleaf - 1;
    int max_id = 2 * Ns;   // node ID for MRCA

    int ntotn = 2 * rtree.nleaf - 1;
    int nstate = comps.size();

    PMAT_DECOMP pmat_decomp;
    DIM_DECOMP dim_decomp;

    set_pmat_decomp(rtree, max_decomp, nstate, knodes, dim_decomp, pmat_decomp, fout);

    map<vector<int>, vector<vector<double>>> sites_lnl_map;

    // for(int nchr = 1; nchr <= vobs.size(); nchr++){     // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with " << vobs[nchr].size() << " sites " << endl;

      for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
          if(debug) cout << "\tfor site " << nc << " on chromosome " << nchr << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k(ntotn, vector<double>(nstate, 0.0));
          vector<vector<int>> S_sk_k(ntotn, vector<int>(nstate, 0));
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  initialize_asr_table_decomp(obs, rtree, comps, max_decomp, pmat_decomp, L_sk_k, S_sk_k, nstate, is_total);
                  get_ancestral_states_site_decomp(L_sk_k, S_sk_k, rtree, knodes, comps, pmat_decomp, dim_decomp, nstate);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  if(debug) cout << "\t\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              initialize_asr_table_decomp(obs, rtree, comps, max_decomp, pmat_decomp, L_sk_k, S_sk_k, nstate, is_total);
              get_ancestral_states_site_decomp(L_sk_k, S_sk_k, rtree, knodes, comps, pmat_decomp, dim_decomp, nstate);
          }

          map<int, int> asr_states;
          extract_tree_ancestral_state(rtree, knodes, comps, L_sk_k, S_sk_k, DECOMP, cn_max, is_total, max_decomp.m_max, asr_states);
          for(int nid = max_id; nid > Ns + 1; nid--){
                int state = asr_states[nid];

                // Get original copy numbers from states
                set<vector<int>>::iterator iter = comps.begin();
                std::advance(iter, state);
                vector<int> c = *iter;
                set<int> cns;
                for(int m1 = 0; m1 <= max_decomp.m_max; m1++){
                    for(int m2 = 0; m2 <= max_decomp.m_max; m2++){
                        int cn = pow(2, c[0] + 1) + m1 * c[1] + c[2] + 2 * m2 * c[3] + 2 * c[4];
                        // only add unique copy number
                        cns.insert(cn);
                    }
                }
                string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc + 1);

                line += "\t" + to_string(state) + "\t" + to_string(pow(10, L_sk_k[nid][state]));
                for(int i = 0; i < L_sk_k[nid].size(); i++){
                    line += "\t" + to_string(pow(10, L_sk_k[nid][i]));
                }
                // fout << setprecision(dbl::max_digits10) << line << endl;
                fout << line << endl;
          }
      }
    } // for each chromosome

    for(auto m : pmat_decomp.pmats_wgd){
        delete [] m.second;
    }
    for(auto m : pmat_decomp.pmats_chr){
        delete [] m.second;
    }
    for(auto m : pmat_decomp.pmats_seg){
        delete [] m.second;
    }

}


// Infer the copy number of all internal nodes given a tree at a site, assuming only site duplication/deletion
void reconstruct_joint_ancestral_state(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, vector<int>& knodes, int model, int cn_max, int use_repeat, int is_total, int m_max, string ofile){
    int debug = 0;
    if(debug) cout << "\treconstruct joint ancestral state" << endl;

    string ofile_joint = ofile + ".joint.state";
    ofstream fout(ofile_joint);

    string header;
    if(is_total){
        header = "node\tsite\tcn";
    }else{
        header = "node\tsite\tcnA\tcnB";
    }
    fout << header << endl;

    string ofile_joint_cn = ofile + ".joint.cn";
    ofstream fout_cn(ofile_joint_cn);

    int Ns = rtree.nleaf - 1;
    // For copy number instantaneous changes
    int nstate = cn_max + 1;
    if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;
    int max_id = 2 * Ns;    // node ID for MRCA

    // Find the transition probability matrix for each branch
    vector<double> blens;
    vector<double*> pmat_per_blen;

    set_pmat(rtree, Ns, nstate, model, cn_max, knodes, blens, pmat_per_blen, fout);

    if(debug){
        cout << "branch lengths so far:";
        for(auto b : blens){
            cout << "\t" << b;
        }
        cout << endl;
    }

    int ntotn = 2 * rtree.nleaf - 1;
    map<vector<int>, vector<vector<double>>> sites_lnl_map;   // ignore chr ID as the likelihood is the same
    map<vector<int>, vector<vector<int>>> sites_state_map;
    // map<vector<int>, int> sites_duplicated;
    map<int, copy_number> cnps;  // all CNs for all internal nodes

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with " << vobs[nchr].size() << " sites" << endl;

      for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
          if(debug) cout << "\tfor site " << nc << " on chromosome " << nchr << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k(ntotn, vector<double>(nstate, 0.0));
          vector<vector<int>> S_sk_k(ntotn, vector<int>(nstate, 0));
          bool is_repeated = false;

          if(use_repeat){
              if(debug) cout << " Use repeated site patterns on site " << nc << endl;
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  if(debug) cout << "sites first seen" << endl;
                  initialize_asr_table(obs, rtree, blens, pmat_per_blen, L_sk_k, S_sk_k, model, nstate, is_total);
                  get_ancestral_states_site(L_sk_k, S_sk_k, rtree, knodes, blens, pmat_per_blen, nstate, model);
                  sites_lnl_map[obs] = L_sk_k;
                  sites_state_map[obs] = S_sk_k;
                //   sites_duplicated[obs] = 1;
              }else{
                  if(debug) cout << "sites repeated on site " << nc << endl;
                //   sites_duplicated[obs]++;
                  L_sk_k = sites_lnl_map[obs];
                  S_sk_k =  sites_state_map[obs];
                  is_repeated = true;
              }
          }else{
              initialize_asr_table(obs, rtree, blens, pmat_per_blen, L_sk_k, S_sk_k, model, nstate, is_total);
              get_ancestral_states_site(L_sk_k, S_sk_k, rtree, knodes, blens, pmat_per_blen, nstate, model);
          }

          if(debug){
            cout << " Get the likelihood table of all internal nodes" << endl;
            print_tree_lnl(rtree, L_sk_k, nstate);
            print_tree_state(rtree, S_sk_k, nstate);
          }

        map<int, int> asr_states;     // The state ID used in rate matrix
        set<vector<int>> comps;    // empty containtor for argument
        extract_tree_ancestral_state(rtree, knodes, comps, L_sk_k, S_sk_k, model, cn_max, is_total, m_max, asr_states);

        if(!is_repeated){
            // int best_state = asr_states[Ns + 1];
            // cout << "optimal state at node "  << Ns + 1 << " is " << best_state << endl;
            // double prob = L_sk_k[max_id][best_state];  // posterior probability of optimal reconstruction
            // cout << "The probability vector at node MRCA " << max_id + 1 << ":";
            // for(int i = 0; i < nstate; i++){
            //     cout << "\t" <<  L_sk_k[max_id][i];
            // }
            // cout << endl;
            // double prob_sum = accumulate(L_sk_k[max_id].begin(), L_sk_k[max_id].end(), 0.0);
            // double prob_rel = prob / prob_sum;
            // cout << prob << "\t" << prob_sum << "\t" << prob_rel << endl;

            for(int nid = max_id; nid > Ns + 1; nid--){
                int state = asr_states[nid];  // state assigned to nid when its parent is optimal
                string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc + 1);

                if(is_total){
                    int tcn = state_to_total_cn(state, cn_max);
                    line += "\t" + to_string(tcn);
                }else{
                    int cnA, cnB;
                    state_to_allele_cn(state, cn_max, cnA, cnB);
                    line += "\t" + to_string(cnA) + "\t" + to_string(cnB);
                }

                // line += "\t" + to_string(prob_rel);

                // fout << setprecision(dbl::max_digits10) << line << endl;
                fout << line << endl;
            }
        }
        // For all sites
        for(int nid = max_id; nid > Ns + 1; nid--){
            int state = asr_states[nid];

            if(is_total){
                int tcn = state_to_total_cn(state, cn_max);
                cnps[nid][nchr][nc] = tcn;
            }else{
                cnps[nid][nchr][nc] = state;
            }
        }
      }
    } // for each chromosome

     for(int nid = max_id; nid > Ns + 1; nid--){
        print_node_cnp(fout_cn, cnps[nid], nid, cn_max, is_total);
     }

    fout.close();
    fout_cn.close();

    if(debug){
        if(use_repeat){
            ofstream fout3("./pat_all");
            // for each chromosome
            for(auto vcn : vobs){
                int nchr = vcn.first;
                for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
                    // for each site of the chromosome (may be repeated)
                    vector<int> obs = vobs[nchr][nc];
                    for(auto c : obs){
                        fout3 << "\t" << c;
                    }
                    fout3 << endl;
                }
            }
            fout3.close();

            cout << "write unique patterns" << endl;
            ofstream fout1("./pat_uniq");
            for(auto sm : sites_lnl_map){
                vector<int> key = sm.first;
                for(auto c : key){
                    fout1 << "\t" << c;
                }
                fout1 << endl;
            }
            fout1.close();

            // cout << "write duplicated patterns" << endl;
            // ofstream fout2("./pat_dup");
            // for(auto sd : sites_duplicated){
            //     fout2 << sd.second;
            //     for(auto c : sd.first){
            //         fout2 << "\t" << c;
            //     }
            //     fout2 << endl;
            // }
            // fout2.close();
        }
    }

    for_each(pmat_per_blen.begin(), pmat_per_blen.end(), DeleteObject());
}
