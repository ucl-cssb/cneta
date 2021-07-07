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


map<int, vector<int>> extract_tree_ancestral_state(const evo_tree& rtree, const set<vector<int>>& comps, const vector<vector<double>>& L_sk_k, const vector<vector<int>>& S_sk_k, int model, int cn_max, int is_total, int m_max, map<int, int> &asr_states){
    int debug = 0;
    map<int, vector<int>> asr_cns;     // copy numbers for each internal node at one site, stored in the order of decreasing IDs
    int Ns = rtree.nleaf - 1;
    int max_id = 2 * Ns;
    // int model = lnl_type.model;
    // int cn_max = lnl_type.cn_max;

    if(debug){
        for(int j = 0; j < L_sk_k[Ns + 1].size(); ++j){
          cout << "\t" << L_sk_k[Ns + 1][j];
        }
        cout << endl;
    }

    int parent_state = 2;
    if(model == BOUNDA){
        parent_state = 4;
    }
    if(model == DECOMP){
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
    }
    asr_states[Ns + 1] = parent_state;
    if(is_total){
        vector<int> cns{2};
        asr_cns[Ns + 1] = cns;
    }
    else{
        vector<int> cns{1,1};
        asr_cns[Ns + 1] = cns;
    }

    // Traverse the tree from root to tips, need to know the parent of each node
    for(int nid = max_id; nid > Ns + 1; nid--){
        // Find the parent node of current node
        int parent = rtree.nodes[nid].parent;
        assert(asr_states.find(parent) != asr_states.end());
        parent_state = asr_states[parent];

        int state = S_sk_k[nid][parent_state];
        asr_states[nid] = state;
        if(debug){
            cout << "\t\tnode " << nid + 1 << " with state " << state  << " and parent " << parent + 1 << " whose state is " << parent_state << endl;
        }

        vector<int> cns;
        if(model == DECOMP){
            // Find the decomposition with state
            set<vector<int>>::iterator iter = comps.begin();
            advance(iter, state);
            vector<int> c = *iter;
            set<int> sums;
            for(int m1 = 0; m1 <= m_max; m1++){
                for(int m2 = 0; m2 <= m_max; m2++){
                    int cn = pow(2, c[0] + 1) + m1 * c[1] + c[2] + 2 * m2 * c[3] + 2 * c[4];
                    // only add unique copy number
                    sums.insert(cn);
                }
            }
            for(auto s : sums){
                cns.push_back(s);
            }
        }else{
            if(is_total){
                int cn = state_to_total_cn(state, cn_max);
                cns.push_back(cn);
            }else{
                int cnA, cnB;
                state_to_allele_cn(state, cn_max, cnA, cnB);
                cns.push_back(cnA);
                cns.push_back(cnB);
            }
        }
        asr_cns[nid] = cns;
    }

    return asr_cns;
}


// Create likelihood vectors and state vectors at the tip node for reconstructing joint ancestral state
// L_sk_k (S_sk_k) has one row for each tree node and one column for each possible state
void initialize_asr_table(const vector<int>& obs, const evo_tree& rtree, map<double, double*>& pmats, vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, int model, int nstate, int is_total){
    int debug = 0;
    if(debug) cout << "Initializing tables for reconstructing joint ancestral state" << endl;

    int Ns = rtree.nleaf - 1;
    for(int i = 0; i < Ns; ++i){
        // cout << "node " << i + 1 << endl;
        // Find the parent node
        int parent = rtree.edges[rtree.nodes[i].e_in].start;
        double blen = rtree.edges[rtree.nodes[i].e_in].length;
        // cout << "parent " << parent + 1 << endl;
        // cout << "blen " << blen << endl;
        double* pblen = pmats[blen];
        // Find the state(s) of current node
        vector<int> tip_states;
        if(model == BOUNDA){
          if(is_total){
              int si = (obs[i] * (obs[i] + 1))/2;
              int ei = si + obs[i];
              for(int k = si; k<=ei; k++){
                 tip_states.push_back(k);
              }
          }else{ // With allele-specific copy number, only the specific site needs to be filled
              tip_states.push_back(obs[i]);
          }
        }
        if(debug) cout << "There are " << tip_states.size() << " states for copy number " << obs[i] << endl;
        for(int j = 0; j < nstate; ++j){  // For each possible parent state, find the most likely tip states
            vector<double> vec_li(nstate, SMALL_LNL);
            for(int m = 0; m < tip_states.size(); ++m){
                int k = tip_states[m];
                if(debug) cout << "parent state " << j << ", child state " << k << endl;
                double li = 0;
                if (model == MK){
                    li = get_transition_prob(rtree.mu, blen, j, k);
                }
                else{
                    li =  pblen[j  + k * nstate];  // assume parent has state j
                }
                if(li > 0) li = log(li);
                else li = SMALL_LNL;
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


// Create likelihood vectors and state vectors at the tip node for reconstructing joint ancestral state, for independent chain model
// L_sk_k (S_sk_k) has one row for each tree node and one column for each possible state
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
            advance(iter, j);
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
                advance(iter, k);
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

                if(li > 0) li = log(li);
                else li = SMALL_LNL;
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


// Find the most likely state for a node under each possible state, assuming current node has state nsk and parent node (connected by branch of length blen) has state np
double get_max_prob_children(const vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, double* pblen, int k, int nstate, int sp, int ni, int nj, int blen, int model){
    int debug = 0;

    vector<double> vec_li;
    double li = 0;
    // loop over possible si for a fixed state of parent (sp)
    for(int si = 0; si < nstate; ++si){
        if (model == MK){
          li = get_transition_prob(rtree.mu, blen, sp, si);
        }
        else{
          li = pblen[sp + si * nstate];
        }

        if(debug){
            cout << "\tfor state " << si << endl;
            cout << "\t\tseparate likelihood "  << li << "\t" << L_sk_k[ni][si] << "\t" << L_sk_k[nj][si] << endl;
        }
        if(li > 0){
            li = log(li);
        }
        else li = SMALL_LNL;
        li += (L_sk_k[ni][si]);
        li += (L_sk_k[nj][si]);

        if(debug) cout << "\t\tscoring: Li\t" << li << endl;
        if(std::isnan(li) || li < SMALL_LNL) li = SMALL_LNL;
        if(debug) cout << "\t\trevised scoring: Li\t" << li << endl;
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
        cout << "for node: k " << k + 1 << ", parent state " << sp << ", max state is " << max_i << " with probability " << exp(max_li) << endl;
    }

    return max_li;
}



// Assume the likelihood table is for each combination of states
double get_max_children_decomp2(const vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, int k, int nstate, double* pbli_wgd, double* pbli_chr, double* pbli_seg, DIM_DECOMP& dim_decomp, int sp, int ni, int nj, int blen){
    int debug = 0;
    int s_wgd, s_chr, s_seg, s_chr2, s_seg2;
    int e_wgd, e_chr, e_seg, e_chr2, e_seg2;
    double prob_wgd, prob_chr, prob_seg, prob_chr2, prob_seg2, prob_chr_all, prob_seg_all;
    // Each copy number is decomposed into a set of 3-tuples
    set<vector<int>>::iterator iter = comps.begin();
    advance(iter, sp);
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

        if(li > 0){
            li = log(li);
        }
        else li = SMALL_LNL;
        li += (L_sk_k[ni][si]);
        li += (L_sk_k[nj][si]);

        if(debug) cout << "\t\tscoring: Li\t" << li << endl;
        if(std::isnan(li) || li < SMALL_LNL) li = SMALL_LNL;
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


// Get the most likely state on one site of a chromosome (assuming higher level events on nodes)
void get_ancestral_states_site(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const vector<int>& knodes, map<double, double*>& pmats, int nstate, int model){
  int debug = 0;
  if(debug){
      cout << "Getting ancestral state for one site" << endl;
  }

  int Ns = rtree.nleaf - 1;
  for(int kn = 0; kn < knodes.size(); ++kn){
        int k = knodes[kn];
        int np = rtree.edges[rtree.nodes[k].e_in].start;
        double blen = rtree.edges[rtree.nodes[k].e_in].length;
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double* pblen = pmats[blen];

        if(debug) cout << "node:" << np + 1 << " -> " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , "  <<  nj + 1 << " , " << blen << endl;

        //loop over possible values of sk
        if(k == 2 * Ns){    // root node is always normal, for edge (Ns + 1, 2Ns)
            if(debug) cout << "Getting states for node MRCA " << k + 1 << endl;
            int sp = 2;
            if(model == BOUNDA) sp = 4;
            if(debug) cout << "likelihood for state " << sp << endl;
            L_sk_k[k][sp] = get_max_prob_children(L_sk_k, S_sk_k, rtree, pblen, k, nstate, sp, ni, nj, blen, model);
        }
        else{
            for(int sp = 0; sp<nstate; ++sp){  // looping over all possible states of its parent
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



// Get the ancestral state on one site of a chromosome
// Assuming each observed copy number is composed of three type of events.
// Sum over all possible states for initial and final nodes
// Allow at most one WGD event along a branch
void get_ancestral_states_site_decomp(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const vector<int>& knodes, const set<vector<int>>& comps, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int nstate){
  int debug = 0;
  if(debug){
      cout << "Getting ancestral state for one site under independent Markov chain model" << endl;
      cout << dim_decomp.dim_wgd << "\t" << dim_decomp.dim_chr << "\t"  << dim_decomp.dim_seg << "\t"  << nstate << endl;
  }

  for(int kn = 0; kn < knodes.size(); ++kn){
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
            for(int sp = 0; sp<nstate; ++sp){
                L_sk_k[k][sp] = get_max_children_decomp2(L_sk_k, S_sk_k, rtree, comps, k, nstate, pbli_wgd, pbli_chr, pbli_seg,  dim_decomp, sp, ni, nj, blen);
            }
        }
  }
  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
    print_tree_state(rtree, S_sk_k, nstate);
  }
}


// Infer the copy number of the MRCA given a tree at a site, assuming independent Markov chains
double reconstruct_marginal_ancestral_state_decomp(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, const set<vector<int>>& comps, OBS_DECOMP& obs_decomp, int use_repeat, int infer_wgd, int infer_chr, int cn_max, ofstream& fout, double min_asr, int is_total){
    int debug = 0;
    if(debug) cout << "\treconstruct marginal ancestral state with independent chain model" << endl;

    // if(!is_tree_valid(rtree, tobs, age, cons)){
    //     return SMALL_LNL;
    // }
    int Ns = rtree.nleaf - 1;
    vector<int> cn_mrca; // CNs for MRCA
    int nid = 2 * (Ns + 1) - 2;

    string header="node\tsite\tcn";
    for(int i = 0; i < comps.size(); i++){
        header += "\tprobablity " + to_string(i + 1);
    }
    fout << header << endl;

    // For WGD model
    int max_wgd = obs_decomp.max_wgd;
    int max_chr_change = obs_decomp.max_chr_change;
    int max_site_change = obs_decomp.max_site_change;

    int dim_wgd = max_wgd + 1;
    int dim_chr = 2 * max_chr_change + 1;
    int dim_seg = 2 * max_site_change + 1;
    DIM_DECOMP dim_decomp = {dim_wgd, dim_chr, dim_seg};

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
        qmat_seg = new double[(dim_seg)*(dim_seg)];  // segment duplication/deletion
        memset(qmat_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
        get_rate_matrix_site_change(qmat_seg, rtree.dup_rate, rtree.del_rate, max_site_change);
    }

    if(debug){
          cout << "Dimension of Qmat " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << "\n";
    }

    //create a list of nodes to loop over (only internal nodes), making sure the root is last
    vector<int> knodes;
    int ntotn = 2 * rtree.nleaf - 1;
    for(int k = Ns+2; k < ntotn; ++k) knodes.push_back(k);
    knodes.push_back(Ns + 1);

    for(int kn = 0; kn < knodes.size(); ++kn){
           int k = knodes[kn];
           double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
           double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

           // For WGD
           if(max_wgd > 0){
               pmati_wgd = new double[(dim_wgd)*(dim_wgd)];
               pmatj_wgd = new double[(dim_wgd)*(dim_wgd)];
               memset(pmati_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
               memset(pmatj_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
               if(pmats_wgd.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_wgd, pmati_wgd, bli, dim_wgd);
                   pmats_wgd[bli] = pmati_wgd;
               }
               if(pmats_wgd.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_wgd, pmatj_wgd, blj, dim_wgd);
                   pmats_wgd[blj] = pmatj_wgd;
               }
           }

           // For chr gain/loss
           if(max_chr_change > 0){
               pmati_chr = new double[(dim_chr)*(dim_chr)];
               pmatj_chr = new double[(dim_chr)*(dim_chr)];
               memset(pmati_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
               memset(pmatj_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
               if(pmats_chr.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_chr, pmati_chr, bli, dim_chr);
                   pmats_chr[bli] = pmati_chr;
               }
               if(pmats_chr.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_chr, pmatj_chr, blj, dim_chr);
                   pmats_chr[blj] = pmatj_chr;
               }
           }

           // For segment duplication/deletion
           if(max_site_change > 0){
               pmati_seg = new double[(dim_seg)*(dim_seg)];
               pmatj_seg = new double[(dim_seg)*(dim_seg)];
               memset(pmati_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
               memset(pmatj_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
               if(pmats_seg.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_seg, pmati_seg, bli, dim_seg);
                   pmats_seg[bli] = pmati_seg;
               }
               if(pmats_seg.count(blj) == 0){
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

    int nstate = comps.size();
    double logL = 0;    // for all chromosmes

    PMAT_DECOMP pmat_decomp = {pmats_wgd, pmats_chr, pmats_seg};
    for(int nchr = 1; nchr <= vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "Computing likelihood on Chr " << nchr << endl;
      double site_logL = 0;   // log likelihood for all sites on a chromosome
      // Use a map to store computed log likelihood
      map<vector<int>, vector<vector<double>>> sites_lnl_map;
      // cout << " chromosome number change is " << 0 << endl;
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
          // site_logL += extract_tree_lnl(L_sk_k, Ns, model);
          double lnl = extract_tree_lnl_decomp(L_sk_k, comps, Ns);
          site_logL += lnl;

          // Get the likelihood table of MRCA node (with largest ID) in the tree from likelihood table
          int state = distance(L_sk_k[nid].begin(), max_element(L_sk_k[nid].begin(), L_sk_k[nid].end()));

          set<vector<int>>::iterator iter = comps.begin();
          // It will move forward the passed iterator by passed value
          advance(iter, state);
          vector<int> c = *iter;
          // There may be multiple possible copy numbers given different coeffients for chromosome-level changes
          int cn = pow(2, c[0] + 1) +  c[1] + c[2] + 2 * c[3] + 2 * c[4];

          cn_mrca.push_back(cn);
          // Print the state of MRCA at this site
          string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc) + "\t" + to_string(cn);
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
    if(max_wgd > 0){
        free(qmat_wgd);
    }
    if(max_chr_change > 0){
        free(qmat_chr);
    }
    if(max_site_change > 0){
        free(qmat_seg);
    }
    for(auto m : pmats_wgd){
        free(m.second);
    }
    for(auto m : pmats_chr){
        free(m.second);
    }
    for(auto m : pmats_seg){
        free(m.second);
    }

    return logL;
}

// Infer the copy number of the MRCA given a tree at a site, assuming only segment duplication/deletion
double reconstruct_marginal_ancestral_state(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, int model, int cn_max, int use_repeat, int is_total, ofstream& fout, double min_asr){
    int debug = 0;
    // if(!is_tree_valid(rtree, tobs, age, cons)){
    //     return SMALL_LNL;
    // }

    int Ns = rtree.nleaf - 1;

    double logL = 0;    // for all chromosmes
    // For copy number instantaneous changes
    int nstate = cn_max + 1;
    if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;

    vector<int> cn_mrca; // CNs for MRCA
    int nid = 2 * (Ns + 1) - 2;
    if(debug) cout << "\treconstruct marginal ancestral state" << endl;

    string header="node\tsite\tcn\tmax probability";
    for(int i = 0; i < nstate; i++){
        header += "\tprobablity " + to_string(i + 1);
    }
    fout << header << endl;

    double *qmat = new double[(nstate)*(nstate)];
    memset(qmat, 0, (nstate)*(nstate)*sizeof(double));
    if(model > 0){
      if(debug){
          cout << "Getting rate matrix" << endl;
      }
      if(model == BOUNDA){
          get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }else{
          get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }
    }

    //create a list of nodes to loop over (only internal nodes), making sure the root is last
    vector<int> knodes;
    int ntotn = 2 * rtree.nleaf - 1;
    for(int k= Ns+2; k < ntotn; ++k) knodes.push_back(k);
    knodes.push_back(Ns + 1);

    // Find the transition probability matrix for each branch
    map<double, double*> pmats;
    for(int kn = 0; kn < knodes.size(); ++kn){
    int k = knodes[kn];
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
    if(model == 1 || model == BOUNDA){
        double *pmati = new double[(nstate)*(nstate)];
        double *pmatj = new double[(nstate)*(nstate)];
        memset(pmati, 0, (nstate)*(nstate)*sizeof(double));
        memset(pmatj, 0, (nstate)*(nstate)*sizeof(double));

        if(pmats.count(bli) == 0){
            get_transition_matrix_bounded(qmat, pmati, bli, nstate);
            pmats[bli] = pmati;
        }
        if(pmats.count(blj) == 0){
            get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
            pmats[blj] = pmatj;
        }
     }
    }

    if(debug){
      for(auto it = pmats.begin(); it != pmats.end(); ++it )
      {
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print(nstate, nstate, it->second, "  P matrix:");
      }
    }

    for(int nchr = 1; nchr <= vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "\tComputing likelihood on Chr " << nchr << endl;
      double site_logL = 0;   // log likelihood for all sites on a chromosome
      map<vector<int>, vector<vector<double>>> sites_lnl_map;

      for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
          if(debug) cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k;
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  L_sk_k = initialize_lnl_table(obs, rtree, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, pmats, 0, 0, model, nstate);
              }else{
                  // cout << "sites repeated" << end1;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              L_sk_k = initialize_lnl_table(obs, rtree, model, nstate, is_total);
              get_likelihood_site(L_sk_k, rtree, knodes, pmats, 0, 0, model, nstate);
          }
          double lnl = extract_tree_lnl(L_sk_k, Ns, model);
          site_logL += lnl;

          // Get the likelihood table of MRCA node (with largest ID) in the tree from likelihood table
          double max_ln = *max_element(L_sk_k[nid].begin(), L_sk_k[nid].end());
          int state = distance(L_sk_k[nid].begin(), max_element(L_sk_k[nid].begin(), L_sk_k[nid].end()));
          if(debug) cout << "\t\tState with maximum likelihood " << state << " at chr " << nchr << " site " << nc << endl;
          int cn = state;
          if(model == BOUNDA){
              cn = state_to_total_cn(state, cn_max);
          }
          cn_mrca.push_back(cn);

          string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc) + "\t" + to_string(cn) + "\t" + to_string(max_ln);
          for(int i = 0; i < L_sk_k[nid].size(); i++){
              line += "\t" + to_string(L_sk_k[nid][i]);
          }
          // fout << setprecision(dbl::max_digits10) << line << endl;
          fout << line << endl;
      }
      logL += site_logL;
      if(debug){
          cout << "Site Likelihood for " << nchr << " is "  << site_logL << endl;
      }
    } // for each chromosome
    if(debug){
        cout << "\nLikelihood without correcting acquisition bias: " << logL << endl;
        cout << "CNs at MRCA is: " << endl;
        for(int i = 0; i < cn_mrca.size(); i++){
            cout << cn_mrca[i] << endl;
        }
    }

    free(qmat);
    for(auto m : pmats){
        free(m.second);
    }

    return logL;
}



// Infer the copy number of all internal nodes given a tree at a site, assuming independent Markov chains
void reconstruct_joint_ancestral_state_decomp(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, const set<vector<int>>& comps, MAX_DECOMP& max_decomp, int use_repeat, int cn_max, ofstream& fout, double min_asr, int is_total){
    int debug = 0;
    if(debug) cout << "\treconstruct joint ancestral state with independent chain model" << endl;

    // if(!is_tree_valid(rtree, tobs, age, cons)){
    //     return;
    // }
    int Ns = rtree.nleaf - 1;
    int max_id = 2 * Ns;

    // For WGD model
    int m_max = max_decomp.m_max;
    int max_wgd = max_decomp.max_wgd;
    int max_chr_change = max_decomp.max_chr_change;
    int max_site_change = max_decomp.max_site_change;

    int dim_wgd = max_wgd + 1;
    int dim_chr = 2 * max_chr_change + 1;
    int dim_seg = 2 * max_site_change + 1;
    DIM_DECOMP dim_decomp = {dim_wgd, dim_chr, dim_seg};

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
        qmat_seg = new double[(dim_seg)*(dim_seg)];  // segment duplication/deletion
        memset(qmat_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
        get_rate_matrix_site_change(qmat_seg, rtree.dup_rate, rtree.del_rate, max_site_change);
    }

    if(debug){
          cout << "Dimension of Qmat " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << "\n";
    }

    //create a list of nodes to loop over (only internal nodes), making sure the root is last
    vector<int> knodes;
    int ntotn = 2 * rtree.nleaf - 1;
    for(int k = Ns+2; k < ntotn; ++k) knodes.push_back(k);
    knodes.push_back(Ns + 1);

    for(int kn = 0; kn < knodes.size(); ++kn){
           int k = knodes[kn];
           double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
           double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

           // For WGD
           if(max_wgd > 0){
               pmati_wgd = new double[(dim_wgd)*(dim_wgd)];
               pmatj_wgd = new double[(dim_wgd)*(dim_wgd)];
               memset(pmati_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
               memset(pmatj_wgd, 0, (dim_wgd)*(dim_wgd)*sizeof(double));
               if(pmats_wgd.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_wgd, pmati_wgd, bli, dim_wgd);
                   pmats_wgd[bli] = pmati_wgd;
               }
               if(pmats_wgd.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_wgd, pmatj_wgd, blj, dim_wgd);
                   pmats_wgd[blj] = pmatj_wgd;
               }
           }

           // For chr gain/loss
           if(max_chr_change > 0){
               pmati_chr = new double[(dim_chr)*(dim_chr)];
               pmatj_chr = new double[(dim_chr)*(dim_chr)];
               memset(pmati_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
               memset(pmatj_chr, 0, (dim_chr)*(dim_chr)*sizeof(double));
               if(pmats_chr.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_chr, pmati_chr, bli, dim_chr);
                   pmats_chr[bli] = pmati_chr;
               }
               if(pmats_chr.count(blj) == 0){
                   get_transition_matrix_bounded(qmat_chr, pmatj_chr, blj, dim_chr);
                   pmats_chr[blj] = pmatj_chr;
               }
           }

           // For segment duplication/deletion
           if(max_site_change > 0){
               pmati_seg = new double[(dim_seg)*(dim_seg)];
               pmatj_seg = new double[(dim_seg)*(dim_seg)];
               memset(pmati_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
               memset(pmatj_seg, 0, (dim_seg)*(dim_seg)*sizeof(double));
               if(pmats_seg.count(bli) == 0){
                   get_transition_matrix_bounded(qmat_seg, pmati_seg, bli, dim_seg);
                   pmats_seg[bli] = pmati_seg;
               }
               if(pmats_seg.count(blj) == 0){
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

    int nstate = comps.size();
    double logL = 0;    // for all chromosmes

    string header = "node\tsite\tcn\tstate\tmax probability";
    for(int i = 0; i < nstate; i++){
        header += "\tprobablity " + to_string(i + 1);
    }
    fout << header << endl;

    knodes.pop_back();  // no need to reconstruct root which is always normal

    PMAT_DECOMP pmat_decomp = {pmats_wgd, pmats_chr, pmats_seg};

    for(int nchr = 1; nchr <= vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with " << vobs[nchr].size() << " sites " << endl;
      map<vector<int>, vector<vector<double>>> sites_lnl_map;
      for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
          if(debug) cout << "\tfor site " << nc << " on chromosome " << nchr << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k(ntotn, vector<double>(nstate,0));
          vector<vector<int>> S_sk_k(ntotn, vector<int>(nstate,0));
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  initialize_asr_table_decomp(obs, rtree, comps, max_decomp, pmat_decomp, L_sk_k, S_sk_k, nstate, is_total);
                  get_ancestral_states_site_decomp(L_sk_k, S_sk_k, rtree, knodes, comps, pmat_decomp, dim_decomp, nstate);
              }else{
                  if(debug) cout << "\t\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              initialize_asr_table_decomp(obs, rtree, comps, max_decomp, pmat_decomp, L_sk_k, S_sk_k, nstate, is_total);
              get_ancestral_states_site_decomp(L_sk_k, S_sk_k, rtree, knodes, comps, pmat_decomp, dim_decomp, nstate);
          }

          map<int, int> asr_states;
          map<int, vector<int>> asr_cns = extract_tree_ancestral_state(rtree, comps, L_sk_k, S_sk_k, DECOMP, cn_max, is_total, m_max, asr_states);
          for(int nid = max_id; nid > Ns + 1; nid--){
              vector<int> cns = asr_cns[nid];
              int state = asr_states[nid];
              string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc);
              if(cns.size()>1){
                  line += "\t" + to_string(cns[0]) + "," + to_string(cns[1]);
              }else{
                line += "\t" + to_string(cns[0]);
              }
              line += "\t" + to_string(state) + "\t" + to_string(pow(10, L_sk_k[nid][state]));
              for(int i = 0; i < L_sk_k[nid].size(); i++){
                  line += "\t" + to_string(pow(10, L_sk_k[nid][i]));
              }
              // fout << setprecision(dbl::max_digits10) << line << endl;
              fout << line << endl;
          }
      }
    } // for each chromosome

    if(max_wgd > 0){
        free(qmat_wgd);
    }
    if(max_chr_change > 0){
        free(qmat_chr);
    }
    if(max_site_change > 0){
        free(qmat_seg);
    }
    for(auto m : pmats_wgd){
        free(m.second);
    }
    for(auto m : pmats_chr){
        free(m.second);
    }
    for(auto m : pmats_seg){
        free(m.second);
    }

}


// Infer the copy number of all internal nodes given a tree at a site, assuming only segment duplication/deletion
void reconstruct_joint_ancestral_state(const evo_tree& rtree, map<int, vector<vector<int>>>& vobs, int model, int cn_max, int use_repeat, int is_total, int m_max, ofstream &fout, double min_asr){
    int debug = 0;
    if(debug) cout << "\treconstruct joint ancestral state" << endl;

    // if(!is_tree_valid(rtree, tobs, age, cons)){
    //     return;
    // }

    int Ns = rtree.nleaf - 1;

    // For copy number instantaneous changes
    int nstate = cn_max + 1;
    if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;
    int max_id = 2 * Ns;

    string header="node\tsite\tcn\tstate\tmax probability";
    for(int i = 0; i < nstate; i++){
        header += "\tprobablity " + to_string(i + 1);
    }
    fout << header << endl;

    double *qmat = new double[(nstate)*(nstate)];
    memset(qmat, 0, (nstate)*(nstate)*sizeof(double));
    if(model > 0){
        if(debug){
            cout << "Getting rate matrix" << endl;
        }
        if(model == BOUNDA){
            get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
        }else{
            get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
        }
    }

    //create a list of nodes to loop over (only internal nodes)
    vector<int> knodes;
    int ntotn = 2 * rtree.nleaf - 1;
    for(int k = Ns+2; k < ntotn; ++k) knodes.push_back(k);
    knodes.push_back(Ns + 1);     // needed for computing P-matrix

    // Find the transition probability matrix for each branch
    map<double, double*> pmats;
    for(int kn = 0; kn < knodes.size(); ++kn){
      int k = knodes[kn];
      double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
      double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
      if(model == 1 || model == BOUNDA){
          double *pmati = new double[(nstate)*(nstate)];
          double *pmatj = new double[(nstate)*(nstate)];
          memset(pmati, 0, (nstate)*(nstate)*sizeof(double));
          memset(pmatj, 0, (nstate)*(nstate)*sizeof(double));

          if(pmats.count(bli) == 0){
              get_transition_matrix_bounded(qmat, pmati, bli, nstate);
              pmats[bli] = pmati;
          }
          if(pmats.count(blj) == 0){
              get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
              pmats[blj] = pmatj;
          }
       }
    }

    if(debug){
        for(auto it = pmats.begin(); it != pmats.end(); ++it )
        {
            double key = it->first;
            cout << "Get Pmatrix for branch length " << key << endl;
            r8mat_print(nstate, nstate, it->second, "  P matrix:");
        }
    }

    knodes.pop_back();  // no need to reconstruct root which is always normal

    for(int nchr = 1; nchr <= vobs.size(); nchr++){     // for each chromosome
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with " << vobs[nchr].size() << " sites " << endl;
      map<vector<int>, vector<vector<double>>> sites_lnl_map;
      for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
          if(debug) cout << "\tfor site " << nc << " on chromosome " << nchr << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k(ntotn, vector<double>(nstate,0));
          vector<vector<int>> S_sk_k(ntotn, vector<int>(nstate,0));
          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  initialize_asr_table(obs, rtree, pmats, L_sk_k, S_sk_k, model, nstate, is_total);
                  get_ancestral_states_site(L_sk_k, S_sk_k, rtree, knodes, pmats, nstate, model);
              }else{
                  if(debug) cout << "\t\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              initialize_asr_table(obs, rtree, pmats, L_sk_k, S_sk_k, model, nstate, is_total);
              get_ancestral_states_site(L_sk_k, S_sk_k, rtree, knodes, pmats, nstate, model);
          }

          map<int, int> asr_states;     // The state ID used in rate matrix
          set<vector<int>> comps;  // empty containtor for argument
          map<int, vector<int>> asr_cns = extract_tree_ancestral_state(rtree, comps, L_sk_k, S_sk_k, model, cn_max, is_total, m_max, asr_states);
          for(int nid = max_id; nid > Ns + 1; nid--){
              vector<int> cns = asr_cns[nid];
              int state = asr_states[nid];
              string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc);
              if(cns.size()>1){
                  line += "\t" + to_string(cns[0]) + "," + to_string(cns[1]);
              }else{
                line += "\t" + to_string(cns[0]);
              }
              line += "\t" + to_string(state) + "\t" + to_string(pow(10, L_sk_k[nid][state]));
              for(int i = 0; i < L_sk_k[nid].size(); i++){
                  line += "\t" + to_string(pow(10, L_sk_k[nid][i]));
              }
              fout << line << endl;
          }
      }
    } // for each chromosome

    free(qmat);
    for(auto m : pmats){
        free(m.second);
    }

}
