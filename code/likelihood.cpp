#include "likelihood.hpp"



void print_lnl_at_tips(const evo_tree& rtree, const vector<int>& obs, const vector<vector<double>>& L_sk_k, int nstate){
    cout << "\nCNs at tips:\n";
    for(int i = 0; i < rtree.nleaf - 1; ++i){
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
}


void initialize_lnl_table(vector<vector<double>>& L_sk_k, const vector<int>& obs, const evo_tree& rtree, int model, int nstate, int is_total){
    int debug = 0;

    int Ns = rtree.nleaf - 1;
    if(model > 1){
        for(int i = 0; i < Ns; ++i){
            // For total copy number, all the possible combinations have to be considered.
            // Set related allele specific cases to be 1, with index from obs[i] * (obs[i] + 1)/2 to obs[i] * (obs[i] + 1)/2 + obs[i]. The index is computed based on pre-specified order.
            if(is_total){
                int si = (obs[i] * (obs[i] + 1)) / 2;
                int ei = si + obs[i];
                for(int k = si; k <= ei; k++){
                    L_sk_k[i][k] = 1.0;
                }
            }else{ // With haplotype-specific copy number, only the specific site needs to be filled
                L_sk_k[i][obs[i]] = 1.0;
            }
        }
        // set unaltered 1/1
        L_sk_k[Ns][NORM_ALLElE_STATE] = 1.0;
    }else{
        for(int i = 0; i < Ns; ++i){
          for(int j = 0; j < nstate; ++j){
    	         if(j == obs[i]) L_sk_k[i][j] = 1.0;
          }
        }
        // set unaltered
        L_sk_k[Ns][NORM_PLOIDY] = 1.0;
    }

    if(debug){
      print_lnl_at_tips(rtree, obs, L_sk_k, nstate);
    }
}


// L_sk_k has one row for each tree node and one column for each possible state; chr starting from 1
// This function is critical in obtaining correct likelihood. If one tip is not initialized, the final likelihood will be 0.
vector<vector<double>> initialize_lnl_table_decomp(vector<int>& obs, OBS_DECOMP& obs_decomp, int chr, const evo_tree& rtree, const set<vector<int>>& comps, int infer_wgd, int infer_chr, int cn_max, int is_total){
    int debug = 0;
    // construct a table for each state of each node
    int ntotn = 2 * rtree.nleaf - 1;
    int nstate = comps.size();
    vector<vector<double>> L_sk_k(ntotn, vector<double>(nstate, 0.0));

    for(int i = 0; i < rtree.nleaf - 1; ++i){
        // For total copy number, all the possible combinations have to be considered.
        // Set related allele specific cases to be 1, with index from obs[i] * (obs[i] + 1)/2 to obs[i] * (obs[i] + 1)/2 + obs[i]. The index is computed based on pre-specified order.
        int cn = obs[i];
        int num_change = 0;
        if(infer_chr && chr > 0){
            num_change = obs_decomp.obs_change_chr[i][chr-1];
        }
        if(!is_total){    // changing input allel-specific copy number to total
            cn = state_to_total_cn(obs[i], cn_max);
        }
        if(debug) cout << "\nStart filling likelihood table for sample "  << i + 1 << " chromosome " << chr  << " copy number " << cn << endl;
        // Fill all the possible state combinations
        int k = 0;
        for (auto c : comps){
            if(debug){
                cout << "\tcn vector:";
                for(int k = 0; k < c.size(); k++){
                    cout << "\t" << c[k];
                }
                cout << endl;
            }
            int alpha = c[0];  // wgd component
            // If a sample has one WGD event, the correponding component must be the same
            if(infer_wgd && alpha != obs_decomp.obs_num_wgd[i]) continue;

            // WGD may occur before, at least one chromosome gain before or after WGD
            if(num_change >= 1 && (c[1] <= 0 && c[3] <= 0)){
                if(debug){
                    cout << "\t\tpotential chromosome gain in sample "  << i + 1 << " chromosome " << chr << " is " << num_change << endl;
                }
                continue;
            }
            if(num_change <= -1 && (c[1] >= 0 && c[3] >= 0)){
                if(debug){
                    cout << "\t\tpotential chromosome loss in sample "  << i + 1 << " chromosome " << chr << " is " << num_change << endl;
                }
                continue;
            }
            // int sum = pow(2, alpha + 1) + c[5] * c[1] + c[2] + 2 * c[6] * c[3] + 2 * c[4];
            // if(sum == cn){
            //     if(debug) cout << "\t\tfilling 1 here" << endl;
            //     L_sk_k[i][k] = 1.0;
            // }
            // assuming m_max >= 1. It is likely that all copies of a segment is lost before chromosome gain/loss
            for(int m1 = 0; m1 <= obs_decomp.m_max; m1++){
                for(int m2 = 0; m2 <= obs_decomp.m_max; m2++){
                    int sum = pow(2, alpha + 1) + m1 * c[1] + c[2] + 2 * m2 * c[3] + 2 * c[4];
                    if(sum == cn){
                        if(debug) cout << "\t\tfilling 1 here" << endl;
                        L_sk_k[i][k] = 1.0;
                    }
                }
            }
            k++;
        }
        // each row should have one entry being 1
        int sum = 0;
        for(int j = 0; j < L_sk_k[i].size(); j++){
            sum += L_sk_k[i][j];
        }
        if(sum < 1){
            cout << "Error in filling table for copy number " << cn << " in sample " << i + 1 << " chromosome " << chr << endl << endl;
        }
    }
    // set likelihood for normal sample
    // Fill all the possible state combinations
    for(int j = 0; j < comps.size(); j++){
        int k = 0;
        for(auto v : comps){
            bool zeros = all_of(v.begin(), v.end(), [](int i) { return i == 0; });
            if(zeros){
                L_sk_k[rtree.nleaf - 1][k] = 1.0;
                break;
            }
            k++;
        }
    }

    if(debug){
        print_lnl_at_tips(rtree, obs, L_sk_k, nstate);
    }

    return L_sk_k;
}


// Assume the likelihood table is for each haplotype-specific copy number
double get_prob_children_decomp(vector<vector<double>>& L_sk_k, const evo_tree& rtree, map<int, set<vector<int>>>& decomp_table, int sk, int cn_max, int nstate, PROB_DECOMP& prob_decomp, DIM_DECOMP& dim_decomp, int ni, int nj, int bli, int blj, int is_total){
    int debug = 0;
    int s_wgd, s_chr, s_seg;
    int e_wgd, e_chr, e_seg;
    // Each copy number is decomposed into a set of 3-tuples
    set<vector<int>> comp_start = decomp_table[sk];


    double Li = 0.0;
    for(auto s : comp_start){
        s_wgd = s[0];
        s_chr = s[1] + s[3];
        s_seg = s[2] + s[4];
        // loop over possible si
        for(int si = 0; si < nstate; ++si){
            // get the start and end state for each type
            int cn = state_to_total_cn(si, cn_max);
            set<vector<int>> comp_end = decomp_table[cn];

            for(auto e : comp_end){
                e_wgd = e[0];
                e_chr = e[1] + e[3];
                e_seg = e[2] + e[4];
                double prob_wgd = prob_decomp.pbli_wgd[s_wgd + e_wgd * dim_decomp.dim_wgd];
                double prob_chr = prob_decomp.pbli_chr[s_chr + e_chr * dim_decomp.dim_chr];
                double prob_seg = prob_decomp.pbli_seg[s_seg + e_seg * dim_decomp.dim_seg];
                double prob = prob_wgd * prob_chr * prob_seg * L_sk_k[ni][si];
                Li += prob;
                if(debug) cout << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob << "\n";
            }
          if(debug) cout << "\tscoring: Li\t" << Li << endl;
        }
    }

    double Lj = 0.0;
    for(auto s : comp_start){
        s_wgd = s[0];
        s_chr = s[1] + s[3];
        s_seg = s[2] + s[4];
        // loop over possible sj
        for(int sj = 0; sj < nstate; ++sj){
            int cn = state_to_total_cn(sj, cn_max);
            set<vector<int>> comp_end = decomp_table[cn];

            for(auto e : comp_end){
                e_wgd = e[0];
                e_chr = e[1] + e[3];
                e_seg = e[2] + e[4];
                double prob_wgd = prob_decomp.pblj_wgd[s_wgd + e_wgd * dim_decomp.dim_wgd];
                double prob_chr = prob_decomp.pblj_chr[s_chr + e_chr * dim_decomp.dim_chr];
                double prob_seg = prob_decomp.pblj_seg[s_seg + e_seg * dim_decomp.dim_seg];
                double prob = prob_wgd * prob_chr * prob_seg * L_sk_k[nj][sj];
                Lj += prob;
                if(debug) cout << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob << "\n";
            }
        }
        if(debug) cout << "\tscoring: Lj\t" << Lj << endl;
    }

    return Li * Lj;
}


// Assume the likelihood table is for each combination of states
double get_prob_children_decomp2(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, int sk, int cn_max, int nstate, PROB_DECOMP& prob_decomp, DIM_DECOMP& dim_decomp, int ni, int nj, int bli, int blj, int is_total){
    int debug = 0;
    int s_wgd, s_chr, s_seg, s_chr2, s_seg2;
    int e_wgd, e_chr, e_seg, e_chr2, e_seg2;
    double prob_wgd, prob_chr, prob_seg, prob_chr2, prob_seg2, prob_chr_all, prob_seg_all;
    // Each copy number is decomposed into a set of 3-tuples
    set<vector<int>>::iterator iter = comps.begin();
    // It will move forward the passed iterator by passed value
    advance(iter, sk);
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
    int delta_chr = (dim_chr - 1) / 2;
    int delta_seg = (dim_seg - 1) / 2;

    if(debug){
        cout << "Starting state " << sk << "\t" << s_wgd << "\t" << s_chr << "\t" << s_seg << "\t" << s_chr2 << "\t" << s_seg2 << "\n";
        cout << "   Offset for chr and segment matrices " << delta_chr << "\t" << delta_seg << "\n";
    }


    double Li = 0.0;
    int si = 0;
    for(auto e : comps){
        double prob = 0;
        prob_wgd = 1;
        prob_chr = 1;
        prob_seg = 1;
        prob_chr2 = 1;
        prob_seg2 = 1;
        prob_chr_all = 1;
        prob_seg_all = 1;
        if(L_sk_k[ni][si] > 0){
            e_wgd = e[0];
            e_chr = e[1];
            e_seg = e[2];
            e_chr2 = e[3];
            e_seg2 = e[4];
            if(dim_wgd > 1) prob_wgd = prob_decomp.pbli_wgd[s_wgd + e_wgd * dim_wgd];
            if(dim_chr > 1){
                prob_chr = prob_decomp.pbli_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
                prob_chr2 = prob_decomp.pbli_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
                prob_chr_all = prob_chr * prob_chr2;
            }
            if(dim_seg > 1){
                prob_seg = prob_decomp.pbli_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
                prob_seg2 = prob_decomp.pbli_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
                prob_seg_all = prob_seg * prob_seg2;
            }
            prob = prob_wgd * prob_chr_all * prob_seg_all * L_sk_k[ni][si];
            if(debug){
                cout << "End state " << si << "\t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2 << "\n";
                cout << "Prob for each event " << "\t" << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << prob << "\n";
            }
        }
        Li += prob;
        si++;
    }
    if(debug) cout << "\tscoring: Li\t" << Li << endl;

    double Lj = 0.0;
    int sj = 0;
    for(auto e : comps){
        double prob = 0;
        prob_wgd = 1;
        prob_chr = 1;
        prob_seg = 1;
        prob_chr2 = 1;
        prob_seg2 = 1;
        prob_chr_all = 1;
        prob_seg_all = 1;
        if(L_sk_k[nj][sj] > 0){
            e_wgd = e[0];
            e_chr = e[1];
            e_seg = e[2];
            e_chr2 = e[3];
            e_seg2 = e[4];
            if(dim_wgd > 1) prob_wgd = prob_decomp.pblj_wgd[s_wgd + e_wgd * dim_wgd];
            if(dim_chr > 1){
                prob_chr = prob_decomp.pblj_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
                prob_chr2 = prob_decomp.pblj_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
                prob_chr_all = prob_chr * prob_chr2;
            }
            if(dim_seg > 1){
                prob_seg = prob_decomp.pblj_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
                prob_seg2 = prob_decomp.pblj_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
                prob_seg_all = prob_seg * prob_seg2;
            }
            prob = prob_wgd * prob_chr_all * prob_seg_all * L_sk_k[nj][sj];
            if(debug){
                cout << "End state " << sj << "\t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2 << "\n";
                cout << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << prob << "\n";
            }
        }
        Lj += prob;
        sj++;
    }
    if(debug) cout << "\tscoring: Lj\t" << Lj << endl;

    return Li * Lj;
}



// Get the likelihood on one site of a chromosome (assuming higher level events on nodes)
// z: possible changes in copy number caused by chromosome gain/loss
void get_likelihood_site(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, const int& has_wgd, const int& z, const int& model, const int& nstate){
  int debug = 0;
  if(debug){
      cout << "Computing likelihood for one site" << endl;
  }

  for(int kn = 0; kn < knodes.size(); ++kn){
    int k = knodes[kn];
    int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

    // find index of bli and blj in blens
    auto pi = std::equal_range(blens.begin(), blens.end(), bli);
    // assert(distance(pi.first, pi.second) == 1);
    int idx_bli = std::distance(blens.begin(), pi.first);
    auto pj = std::equal_range(blens.begin(), blens.end(), blj);
    // assert(distance(pj.first, pj.second) == 1);
    int idx_blj = std::distance(blens.begin(), pj.first);

    double* pbli = pmat_per_blen[idx_bli];
    double* pblj = pmat_per_blen[idx_blj];

    if(debug){
      cout << "branch lengths so far:";
      for(auto b : blens){
          cout << "\t" << b;
      }
      cout << endl;
      cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;
      cout << "Get Pmatrix for branch length " << bli << "\t" << blens[idx_bli] << endl;
      r8mat_print(nstate, nstate, pmat_per_blen[idx_bli], "  P matrix:");
      cout << "Get Pmatrix for branch length " << blj << "\t" << blens[idx_blj] << endl;
      r8mat_print(nstate, nstate, pmat_per_blen[idx_blj], "  P matrix:");
    }

    //loop over possible values of sk
    if(k == rtree.nleaf){    // root node is always normal
        if(debug) cout << "Getting likelihood for root node " << k << endl;
        int nsk = NORM_PLOIDY;
        if(model == BOUNDA) nsk = NORM_ALLElE_STATE;
        L_sk_k[k][nsk] = get_prob_children(L_sk_k, rtree, pbli, pblj, nsk, ni, nj, bli, blj, model, nstate);
    }else{
        for(int sk = 0; sk < nstate; ++sk){
            int nsk = sk;  // state after changes by other large scale events
            if(has_wgd) nsk = 2 * sk;
            nsk += z;
            if(debug) cout << "likelihood for state " << nsk << endl;
            if(nsk < 0 || nsk >= nstate) continue;
            // cout << " getting likelihood of children nodes " << endl;
            L_sk_k[k][nsk] = get_prob_children(L_sk_k, rtree, pbli, pblj, nsk, ni, nj, bli, blj, model, nstate);
        }
    }
  }
  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
  }
}



// Get the likelihood on one site of a chromosome
// Assuming each observed copy number is composed of three type of events.
// Sum over all possible states for initial and final nodes
// Allow at most one WGD event along a branch
void get_likelihood_site_decomp(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, const vector<int>& knodes, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int cn_max, int is_total){
  int debug = 0;

  int dim_wgd = dim_decomp.dim_wgd;
  int dim_chr = dim_decomp.dim_chr;
  int dim_seg = dim_decomp.dim_seg;
  int nstate = comps.size();

  if(debug){
      cout << "Computing likelihood for one site" << endl;
      cout << dim_wgd << "\t" << dim_chr << "\t"  << dim_seg << "\t"  << nstate << endl;
  }

  for(int kn = 0; kn < knodes.size(); ++kn){
        int k = knodes[kn];
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

        double *pbli_wgd, *pblj_wgd;
        double *pbli_chr, *pblj_chr;
        double *pbli_seg, *pblj_seg;

        if(dim_wgd > 1){
            pbli_wgd = pmat_decomp.pmats_wgd[bli];
            pblj_wgd = pmat_decomp.pmats_wgd[blj];
        }
        if(dim_chr > 1){
            pbli_chr = pmat_decomp.pmats_chr[bli];
            pblj_chr = pmat_decomp.pmats_chr[blj];
        }
        if(dim_seg > 1){
            pbli_seg = pmat_decomp.pmats_seg[bli];
            pblj_seg = pmat_decomp.pmats_seg[blj];
        }

        PROB_DECOMP prob_decomp;
        prob_decomp.pbli_wgd = pbli_wgd;
        prob_decomp.pblj_wgd = pblj_wgd;
        prob_decomp.pbli_chr = pbli_chr;
        prob_decomp.pblj_chr = pblj_chr;
        prob_decomp.pbli_seg = pbli_seg;
        prob_decomp.pblj_seg = pblj_seg;


        if(debug) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;

        // loop over possible observed states of start nodes
        if(k == rtree.nleaf){    // root node is always normal
            // int sk = 4;
            // L_sk_k[k][sk] = get_prob_children_decomp(L_sk_k, rtree, decomp_table, sk, cn_max, nstate, prob_decomp, dim_decomp, ni, nj, bli, blj, is_total);
            if(debug) cout << "Getting likelihood for root node " << k << endl;
            int sk = 0;
            for(auto v : comps){
                bool zeros = all_of(v.begin(), v.end(), [](int i) { return i == 0; });
                if(zeros){
                    L_sk_k[k][sk] = get_prob_children_decomp2(L_sk_k, rtree, comps, sk, cn_max, nstate, prob_decomp, dim_decomp, ni, nj, bli, blj, is_total);
                    break;
                }
                sk++;
            }
        }else{
            for(int sk = 0; sk < nstate; ++sk){
                // cout << " getting likelihood of children nodes " << endl;
                // L_sk_k[k][sk] = get_prob_children_decomp(L_sk_k, rtree, decomp_table, sk, cn_max, nstate, prob_decomp, dim_decomp, ni, nj, bli, blj, is_total);
                L_sk_k[k][sk] = get_prob_children_decomp2(L_sk_k, rtree, comps, sk, cn_max, nstate, prob_decomp, dim_decomp, ni, nj, bli, blj, is_total);
            }
        }
  }
  if(debug > 0){
    print_tree_lnl(rtree, L_sk_k, nstate);
  }
}


double get_likelihood_chr(map<int, vector<vector<int>>>& vobs, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, const int& has_wgd, const int& only_seg, const int& use_repeat, const int& model, const int& nstate, const int& is_total){
    int debug = 0;
    double logL = 0.0;    // for all chromosmes
    double chr_gain = 0.0;
    double chr_loss = 0.0;
    // Use a map to store computed log likelihood to save computation on duplicated site patterns
    map<vector<int>, vector<vector<double>>> sites_lnl_map;

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with  " << vobs[nchr].size() << " sites " << endl;
      double chr_logL = 0.0;  // for one chromosome
      double chr_logL_normal = 0.0, chr_logL_gain = 0.0, chr_logL_loss = 0.0;
      double site_logL = 0.0;   // log likelihood for all sites on a chromosome
      int z = 0;    // no chr gain/loss
      // cout << " chromosome number change is " << 0 << endl;

      for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs[nchr][nc];
          vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));

          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, has_wgd, z, model, nstate);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  // cout << "sites repeated" << end1;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
              get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, has_wgd, z, model, nstate);
          }

          site_logL += extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);

          if(debug){
              // cout << "\nLikelihood for site " << nc << " is " << lnl << endl;
              print_tree_lnl(rtree, L_sk_k, nstate);
          }
      }

      double chr_normal = 1.0;
      if(!only_seg){
          chr_gain = rtree.chr_gain_rate;
          chr_loss = rtree.chr_loss_rate;

          if(debug){
              cout << "chromosome gain rate " << chr_gain << endl;
              cout << "chromosome loss rate " << chr_loss << endl;
              cout << "Number of chr so far " << vobs.size() << endl;
          }

          if(fabs(chr_loss) > SMALL_VAL){
             chr_normal -= chr_loss;
          }
          if(fabs(chr_gain) > SMALL_VAL){
             chr_normal -= chr_gain;
          }
      }

      chr_logL_normal = log(chr_normal) + site_logL;
      chr_logL += chr_logL_normal;

      if(debug){
          cout << "Likelihood without chr gain/loss for " << nchr << " is "  << chr_normal << endl;
          cout << "Site Likelihood for " << nchr << " is "  << site_logL << endl;
          cout << "Likelihood without chr gain/loss: " << chr_logL_normal << endl;
      }

      if(!only_seg){
          if(fabs(chr_loss) > SMALL_VAL){
              z = -1;
              double site_logL = 0.0;   // log likelihood for all sites on a chromosome
              // cout << " chromosome number change is " << z << endl;
              for(int nc = 0; nc < vobs[nchr].size(); nc++){
                  // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
                  // for each site of the chromosome
                  vector<int> obs = vobs[nchr][nc];
                  vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
                  initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);

                  get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, has_wgd, z, model, nstate);
                  site_logL += extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);

                  if(debug){
                      print_tree_lnl(rtree, L_sk_k, nstate);
                  }
              } // for all sites on a chromosome

              chr_logL_loss = log(chr_loss) + site_logL;
              chr_logL += log(1 + exp(chr_logL_loss - chr_logL_normal));
              if(debug){
                  cout << "\nLikelihood before chr loss for " << nchr << " is " << site_logL << endl;
                  cout << "\nLikelihood after chr loss: " << chr_logL_loss << endl;
              }
          } // for all chromosome loss

          if(fabs(chr_gain) > SMALL_VAL){
              z = 1;
              double site_logL = 0.0;   // log likelihood for all sites on a chromosome
              // cout << " chromosome number change is " << z << endl;
              for(int nc = 0; nc < vobs[nchr].size(); nc++){
                  // cout << "Number of sites for this chr " << vobs[nchr].size() << endl;
                  // for each site of the chromosome
                  vector<int> obs = vobs[nchr][nc];
                  vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
                  initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, has_wgd, z, model, nstate);
                  site_logL += extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);

                  if(debug){
                      print_tree_lnl(rtree, L_sk_k, nstate);
                  }
              } // for all sites on a chromosome

              chr_logL_gain = log(chr_gain) + site_logL;
              if(chr_logL_loss > 0){
                  chr_logL += log(1 + 1 / (exp(chr_logL_normal - chr_logL_gain) + exp(chr_logL_loss - chr_logL_gain)));
              }else{
                  chr_logL += log(1 + exp(chr_logL_gain - chr_logL_normal));
              }

              if(debug){
                  cout << "\nLikelihood before chr gain for " << nchr << " is " << site_logL << endl;
                  cout << "\nLikelihood after chr gain: " << chr_logL_gain << endl;
              }
          } // for all chromosome loss
          // chr_logL = chr_logL_normal + log(1 + exp(chr_logL_loss-chr_logL_normal)) + log(1 + 1 / (exp(chr_logL_normal-chr_logL_gain) + exp(chr_logL_loss-chr_logL_gain)));
      }

      logL += chr_logL;

      if(debug){
          cout << "\nLikelihood after considering chr gain/loss for  " << nchr << " is " << logL << endl;
      }
    } // for each chromosome

    if(debug){
        cout << "\nLikelihood with chr gain/loss for all chromosmes: " << logL << endl;
    }

    return logL;
}


double get_likelihood_chr_decomp(map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const evo_tree& rtree, const set<vector<int>>& comps, const vector<int>& knodes, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int infer_wgd, int infer_chr, int use_repeat, int cn_max, int is_total){
    int debug = 0;
    double logL = 0.0;    // for all chromosmes
    int nstate = comps.size();
    // Use a map to store computed log likelihood
    map<vector<int>, vector<vector<double>>> sites_lnl_map;

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug){
        cout << "Computing likelihood on Chr " << nchr <<  " with " << vobs[nchr].size() << "sites" << endl;
      }
      double site_logL = 0.0;   // log likelihood for all sites on a chromosome

      // cout << " chromosome number change is " << 0 << endl;
      for(int nc = 0; nc < vobs[nchr].size(); nc++){    // for each segment on the chromosome
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

          site_logL += extract_tree_lnl_decomp(L_sk_k, comps, rtree.nleaf - 1);

          if(debug){
              // cout << "Crtree.nleaf - 1 at this site: ";
              for(int i = 0; i < obs.size(); i++){
                  cout << "\t" << obs[i];
              }
              cout << endl;
              print_tree_lnl(rtree, L_sk_k, nstate);
          }
      }

      logL += site_logL;

      if(debug){
          cout << "\nLikelihood for chromosome " << nchr << " is " << site_logL << endl;
      }
    } // for each chromosome

    return logL;
}


double get_likelihood_revised(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, LNL_TYPE& lnl_type){
  // int debug = 0;
  // if(debug) cout << "\tget_likelihood by matrix exponential" << endl;

  if(!is_tree_valid(rtree, lnl_type.max_tobs, lnl_type.patient_age, lnl_type.cons)){
       return SMALL_LNL;
   }

  int model = lnl_type.model;
  int cn_max = lnl_type.cn_max;
  int is_total = lnl_type.is_total;

  // For copy number instantaneous changes
  int nstate = cn_max + 1;
  if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;

  // may be different for each tree due to change of mutation rate
  int dim_mat = nstate * nstate;
  double *qmat = new double[dim_mat];
  memset(qmat, 0.0, dim_mat * sizeof(double));

  // model == BOUNDT || model == BOUNDA
  assert(model > 0);
  if(model == BOUNDA){
      get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
  }else{
      get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
  }

  // if(debug){
  //   r8mat_print(nstate, nstate, qmat, "  Q matrix:" );
  //   check_matrix_row_sum(qmat, nstate);
  // }

  // Find the transition probability matrix for each branch
  vector<int> knodes = lnl_type.knodes;
  vector<double> blens;
  vector<double*> pmat_per_blen;

  for(int kn = 0; kn < knodes.size(); ++kn){
    int k = knodes[kn];
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

    if(find(blens.begin(), blens.end(), bli) == blens.end()){
        double *pmati = new double[dim_mat];
        memset(pmati, 0.0, dim_mat * sizeof(double));
        get_transition_matrix_bounded(qmat, pmati, bli, nstate);
        pmat_per_blen.push_back(pmati);
        blens.push_back(bli);
        // if(debug){
        //     cout << "Get Pmatrix for branch length " << bli << endl;
        //     r8mat_print(nstate, nstate, pmati, "  P matrix:" );
        // }
    }
    if(find(blens.begin(), blens.end(), blj) == blens.end()){
        double *pmatj = new double[dim_mat];
        memset(pmatj, 0.0, dim_mat * sizeof(double));
        get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
        pmat_per_blen.push_back(pmatj);
        blens.push_back(blj);
        // if(debug){
        //     cout << "Get Pmatrix for branch length " << blj << endl;
        //     r8mat_print(nstate, nstate, pmatj, "  P matrix:" );
        // }
    }
  }

  // sort pmats according to branch lengths
  auto p = sort_permutation(blens, [&](const double& a, const double& b){ return a < b; });
  blens = apply_permutation(blens, p);
  pmat_per_blen = apply_permutation(pmat_per_blen, p);

  // if(debug){
  //     for(int i = 0; i < pmat_per_blen.size(); ++i){
  //         double blen = blens[i];
  //         cout << "Get Pmatrix for branch length " << blen << endl;
  //         r8mat_print(nstate, nstate, pmat_per_blen[i], "  P matrix:");
  //     }
  // }

  double logL = 0.0;

  if(lnl_type.only_seg){
      // if(debug) cout << "Computing the likelihood without consideration of WGD" << endl;
      logL += get_likelihood_chr(vobs, rtree, knodes, blens, pmat_per_blen, 0, lnl_type.only_seg, lnl_type.use_repeat, model, nstate, is_total);
  }else{
      // if(debug) cout << "Computing the likelihood with consideration of WGD" << endl;
      logL += (1 - rtree.wgd_rate) * get_likelihood_chr(vobs, rtree, knodes, blens, pmat_per_blen, 0, lnl_type.only_seg, lnl_type.use_repeat, model, nstate, is_total);
      logL += rtree.wgd_rate * get_likelihood_chr(vobs, rtree, knodes, blens, pmat_per_blen, 1, lnl_type.only_seg, lnl_type.use_repeat, model, nstate, is_total);
  }

  // if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;

  if(lnl_type.correct_bias){
      // if(debug) cout << "Correcting for the skip of invariant sites" << endl;

      // Compute the likelihood of dummy sites consisting entirely of 2s for the tree
      double lnl_invar = 0.0;
      // Suppose the value is 2 for all samples
      int normal_cn = 2;
      if(!is_total){
          normal_cn = 4;    // state ID for haplotype-specific CN 1/1
      }

      vector<int> obs(rtree.nleaf - 1, normal_cn);
      vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
      initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);

      for(int kn = 0; kn < knodes.size(); ++kn){
          int k = knodes[kn];
          int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
          double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
          int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
          double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

          DBIterPair pi = equal_range(blens.begin(), blens.end(), bli);
          // assert(distance(pi.first, pi.second) == 1);
          int idx_bli = distance(blens.begin(), pi.first);
          DBIterPair pj = equal_range(blens.begin(), blens.end(), blj);
          // assert(distance(pj.first, pj.second) == 1);
          int idx_blj = distance(blens.begin(), pj.first);

          //loop over possible values of sk
          for(int sk = 0; sk < nstate; ++sk){
            double Li = 0.0;
            for(int si = 0; si < nstate; ++si){
                if(model == MK){
                  Li += get_transition_prob(rtree.mu, bli, sk, si) * L_sk_k[ni][si];
                }else{
                  Li += pmat_per_blen[idx_bli][sk + si * nstate] * L_sk_k[ni][si];
                }
            }
            double Lj = 0.0;
            for(int sj = 0; sj < nstate; ++sj){
                if(model == MK){
                     Lj += get_transition_prob(rtree.mu, blj, sk, sj) * L_sk_k[nj][sj];
                }else{
                     Lj += pmat_per_blen[idx_blj][sk + sj * nstate] * L_sk_k[nj][sj];
                }
            }

            L_sk_k[k][sk] = Li * Lj;
         }
      }

      lnl_invar = extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);

      double bias = lnl_type.num_invar_bins * lnl_invar;
      logL = logL + bias;

      // if(debug){
      //     cout << "Likelihood of an invariant bin: " << lnl_invar << endl;
      //     cout << "Number of invariant bins " << lnl_type.num_invar_bins << endl;
      //     cout << "Bias to correct " << bias << endl;
      //     cout << "Final likelihood after correcting acquisition bias: " << logL << endl;
      // }
  }

  if(std::isnan(logL) || logL < SMALL_LNL) logL = SMALL_LNL;
  // if(debug){
  //     cout << "Final likelihood: " << logL << endl;
  // }

  delete [] qmat;
  for_each(pmat_per_blen.begin(), pmat_per_blen.end(), DeleteObject());

  return logL;
}

// Computing likelihood when WGD and chr gain/loss are incorporated
// Assume likelihood is for haplotype-specific information
double get_likelihood_decomp(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type){
  int debug = 0;
  if(debug) cout << "\tget_likelihood from multiple chains" << endl;

  if(!is_tree_valid(rtree, lnl_type.max_tobs, lnl_type.patient_age, lnl_type.cons)){
       return SMALL_LNL;
   }

  int model = lnl_type.model;
  int cn_max = lnl_type.cn_max;
  int is_total = lnl_type.is_total;

  int max_wgd = obs_decomp.max_wgd;
  int max_chr_change = obs_decomp.max_chr_change;
  int max_site_change = obs_decomp.max_site_change;

  // For WGD model
  int dim_wgd = max_wgd + 1;
  int dim_chr = 2 * max_chr_change + 1;
  int dim_seg = 2 * max_site_change + 1;
  int dim_mat_wgd = dim_wgd * dim_wgd;
  int dim_mat_chr = dim_chr * dim_chr;
  int dim_mat_seg = dim_seg * dim_seg;

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
      qmat_wgd = new double[dim_mat_wgd];  // WGD
      memset(qmat_wgd, 0.0, dim_mat_wgd * sizeof(double));
      get_rate_matrix_wgd(qmat_wgd, rtree.wgd_rate, max_wgd);
  }

  if(max_chr_change > 0){
      qmat_chr = new double[dim_mat_chr];   // chromosome gain/loss
      memset(qmat_chr, 0.0, dim_mat_chr * sizeof(double));
      get_rate_matrix_chr_change(qmat_chr, rtree.chr_gain_rate, rtree.chr_loss_rate, max_chr_change);
  }

  if(max_site_change > 0){
      qmat_seg = new double[dim_mat_seg];  // segment duplication/deletion
      memset(qmat_seg, 0.0, dim_mat_seg * sizeof(double));
      get_rate_matrix_site_change(qmat_seg, rtree.dup_rate, rtree.del_rate, max_site_change);
  }

  if(debug){
    cout << "Dimension of Qmat " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << "\n";
  }

  vector<int> knodes = lnl_type.knodes;
  for(int kn = 0; kn < knodes.size(); ++kn){
         int k = knodes[kn];
         double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
         double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

         // For WGD
         if(max_wgd > 0){
             if(pmats_wgd.find(bli) == pmats_wgd.end()){
                pmati_wgd = new double[dim_mat_wgd];
                memset(pmati_wgd, 0.0, dim_mat_wgd*sizeof(double));
                get_transition_matrix_bounded(qmat_wgd, pmati_wgd, bli, dim_wgd);
                pmats_wgd[bli] = pmati_wgd;
             }
             if(pmats_wgd.find(blj) == pmats_wgd.end()){
                pmatj_wgd = new double[dim_mat_wgd];
                memset(pmatj_wgd, 0.0, dim_mat_wgd*sizeof(double));
                get_transition_matrix_bounded(qmat_wgd, pmatj_wgd, blj, dim_wgd);
                pmats_wgd[blj] = pmatj_wgd;
             }
         }

         // For chr gain/loss
         if(max_chr_change > 0){
             if(pmats_chr.find(bli) == pmats_chr.end()){
                 pmati_chr = new double[dim_mat_chr];
                 memset(pmati_chr, 0.0, dim_mat_chr*sizeof(double));
                 get_transition_matrix_bounded(qmat_chr, pmati_chr, bli, dim_chr);
                 pmats_chr[bli] = pmati_chr;
             }
             if(pmats_chr.find(blj) == pmats_chr.end()){
                 pmatj_chr = new double[dim_mat_chr];
                 memset(pmatj_chr, 0.0, dim_mat_chr*sizeof(double));
                 get_transition_matrix_bounded(qmat_chr, pmatj_chr, blj, dim_chr);
                 pmats_chr[blj] = pmatj_chr;
             }
         }

         // For segment duplication/deletion
         if(max_site_change > 0){
             if(pmats_seg.find(bli) == pmats_seg.end()){
                 pmati_seg = new double[dim_mat_seg];
                 memset(pmati_seg, 0.0, dim_mat_seg*sizeof(double));
                 get_transition_matrix_bounded(qmat_seg, pmati_seg, bli, dim_seg);
                 pmats_seg[bli] = pmati_seg;
             }
             if(pmats_seg.find(blj) == pmats_seg.end()){
                 pmatj_seg = new double[dim_mat_seg];
                 memset(pmatj_seg, 0.0, dim_mat_seg*sizeof(double));
                 get_transition_matrix_bounded(qmat_seg, pmatj_seg, blj, dim_seg);
                 pmats_seg[blj] = pmatj_seg;
            }
        }
  }

  if(debug){
      for(auto it = pmats_wgd.begin(); it != pmats_wgd.end(); ++it){
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print(dim_wgd, dim_wgd, it->second, "  P-WGD matrix:");
      }
      for(auto it = pmats_chr.begin(); it != pmats_chr.end(); ++it){
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print(dim_chr, dim_chr, it->second, "  P-CHR matrix:");
      }
      for(auto it = pmats_seg.begin(); it != pmats_seg.end(); ++it){
          double key = it->first;
          cout << "Get Pmatrix for branch length " << key << endl;
          r8mat_print(dim_seg, dim_seg, it->second, "  P-SEG matrix:");
      }
  }

  double logL = 0.0;

  PMAT_DECOMP pmat_decomp;
  pmat_decomp.pmats_wgd = pmats_wgd;
  pmat_decomp.pmats_chr = pmats_chr;
  pmat_decomp.pmats_seg = pmats_seg;

  DIM_DECOMP dim_decomp;
  dim_decomp.dim_wgd = dim_wgd;
  dim_decomp.dim_chr = dim_chr;
  dim_decomp.dim_seg = dim_seg;

  // cout << "Number of states is " << nstate << endl;
  logL = get_likelihood_chr_decomp(vobs, obs_decomp, rtree, comps, knodes, pmat_decomp, dim_decomp, lnl_type.infer_wgd, lnl_type.infer_chr, lnl_type.use_repeat, cn_max, is_total);

  if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;
  if(lnl_type.correct_bias){
    if(debug) cout << "Correcting for the skip of invariant sites" << endl;
    vector<int> obs(rtree.nleaf - 1, NORM_PLOIDY);
    vector<vector<double>> L_sk_k = initialize_lnl_table_decomp(obs, obs_decomp, 0, rtree, comps, lnl_type.infer_wgd, lnl_type.infer_chr, cn_max, is_total);
    get_likelihood_site_decomp(L_sk_k, rtree, comps, knodes, pmat_decomp, dim_decomp, cn_max, is_total);
    double lnl_invar = extract_tree_lnl_decomp(L_sk_k, comps, rtree.nleaf - 1);

    double bias = lnl_type.num_invar_bins * lnl_invar;
    logL = logL + bias;
    if(debug){
        cout << "Likelihood of an invariant bin " << lnl_invar << endl;
        cout << "Number of invariant bins " << lnl_type.num_invar_bins << endl;
        cout << "Bias to correct " << bias << endl;
        cout << "Final likelihood after correcting acquisition bias: " << logL << endl;
    }
  }

  if(std::isnan(logL) || logL < SMALL_LNL) logL = SMALL_LNL;
  if(debug){
      cout << "Final likelihood: " << logL << endl;
      cout << "Free memory" << endl;
  }

  if(max_wgd > 0){
      delete [] qmat_wgd;
      for(auto m : pmats_wgd){
          delete [] m.second;
      }
  }
  if(max_chr_change > 0){
      delete [] qmat_chr;
      for(auto m : pmats_chr){
          delete [] m.second;
      }
  }
  if(max_site_change > 0){
      delete [] qmat_seg;
      for(auto m : pmats_seg){
          delete [] m.second;
      }
  }

  return logL;
}


// Compute the likelihood without grouping sites by chromosome, only considering segment duplication/deletion
/*
vobs: the observed data matrix
rtree: the given tree
model: model of evolution
*/
double get_likelihood(const vector<vector<int>>& vobs, evo_tree& rtree, const vector<int>& knodes, int model, int cn_max, int is_total){
  int debug = 0;
  if(debug) cout << "\tget_likelihood" << endl;

  assert(vobs.size() > 1);
  int Nchar = vobs[0].size();  // Nchar: number of characters for each sample

  int nstate = cn_max + 1;
  if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;

  int dim_mat = nstate * nstate;
  double *qmat = new double[dim_mat];
  memset(qmat, 0.0, dim_mat * sizeof(double));

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

  double *pmati = new double[dim_mat];
  double *pmatj = new double[dim_mat];
  memset(pmati, 0.0, dim_mat * sizeof(double));
  memset(pmatj, 0.0, dim_mat * sizeof(double));

  double logL = 0;
  for(int nc = 0; nc < Nchar; ++nc){
    vector<int> obs = vobs[nc];
    if(debug){
      cout << "char: " << nc << endl;
      for(int i = 0; i < rtree.nleaf - 1; ++i) cout << "\t" << obs[i];
      cout << endl;
    }

    vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
    initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
    if(debug){
        print_tree_lnl(rtree, L_sk_k, nstate);
    }

    for(int kn = 0; kn < knodes.size(); ++kn){
      int k = knodes[kn];
      int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
      double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
      int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
      double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
      if(model > 0){
           get_transition_matrix_bounded(qmat, pmati, bli, nstate);
           get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
           if(debug){
               cout << "Get Pmatrix for branch length " << bli << endl;
               r8mat_print(nstate, nstate, pmati, "  P matrix after change:" );
               cout << "Get Pmatrix for branch length " << blj << endl;
               r8mat_print(nstate, nstate, pmatj, "  P matrix after change:" );
           }
      }

      if(debug) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;

      //loop over possible values of sk
      for(int sk = 0; sk < nstate; ++sk){
    	  double Li = 0;
    	  // loop over possible si
    	  for(int si = 0; si < nstate; ++si){
            if (model == MK){
                Li += get_transition_prob(rtree.mu, bli, sk, si) * L_sk_k[ni][si];
            }
            else{
                Li += pmati[sk + si * nstate] * L_sk_k[ni][si];
            }
    	       //cout << "\tscoring: Li\t" << li << "\t" << get_transition_prob(mu, bli, sk, si ) << "\t" << L_sk_k[ni][si] << endl;
        }

    	  double Lj = 0;
    	  // loop over possible sj
    	  for(int sj = 0; sj < nstate; ++sj){
            if (model == MK){
    	         Lj += get_transition_prob(rtree.mu, blj, sk, sj) * L_sk_k[nj][sj];
            }
            else{
               Lj += pmatj[sk + sj * nstate] * L_sk_k[nj][sj];
            }
    	  }

	      //cout << "scoring: sk" << sk << "\t" <<  Li << "\t" << Lj << endl;
	      L_sk_k[k][sk] = Li * Lj;
      }

      if(debug){
    	  print_tree_lnl(rtree, L_sk_k, nstate);
      }
    }

    logL += extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);
  }

  if(debug) cout << "Final likelihood: " << logL << endl;

  delete [] qmat;
  delete [] pmati;
  delete [] pmatj;

  return logL;
}


double extract_tree_lnl(vector<vector<double>>& L_sk_k, int Ns, int model){
    int debug = 0;

    if(debug){
        for(int j = 0; j < L_sk_k[Ns + 1].size(); ++j){
          cout << "\t" << L_sk_k[Ns + 1][j];
        }
        cout << endl;
    }

    if(model == BOUNDA){
        // The index is changed from 2 to 4 (1/1)
        if(debug) cout << "Likelihood for root is " << L_sk_k[Ns + 1][NORM_ALLElE_STATE] << endl;

        if(L_sk_k[Ns + 1][NORM_ALLElE_STATE] > 0) return log(L_sk_k[Ns + 1][NORM_ALLElE_STATE]);
        else return LARGE_LNL;
    }else{
        if(debug) cout << "Likelihood for root is " << L_sk_k[Ns + 1][NORM_PLOIDY] << endl;

        if(L_sk_k[Ns + 1][NORM_PLOIDY] > 0) return log(L_sk_k[Ns + 1][NORM_PLOIDY]);
        else return LARGE_LNL;
    }
}




// Get the likelihood of the tree from likelihood table of state combinations
double extract_tree_lnl_decomp(vector<vector<double>>& L_sk_k, const set<vector<int>>& comps, int Ns){
    int debug = 0;
    if(debug) cout << "Extracting likelihood for the root" << endl;

    double likelihood = 0;
    int k = 0;

    for (auto v : comps){
        // if(debug) cout << k << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << "\t" << v[3] << "\t" << v[4] << endl;
        bool zeros = all_of(v.begin(), v.end(), [](int i){ return i == 0; });
        if(zeros){
            likelihood = L_sk_k[Ns + 1][k];
            break;
        }
        k++;
    }

    if(debug){
        for(int j = 0; j < comps.size(); ++j){
          cout << "\t" << L_sk_k[Ns + 1][j];
        }
        cout << endl;
    }

    if(likelihood > 0) return log(likelihood);
    else return LARGE_LNL;
}


void print_tree_lnl(const evo_tree& rtree, vector<vector<double>>& L_sk_k, int nstate){
    cout << "\nLikelihood so far:\n";

    int ntotn = 2 * rtree.nleaf - 1;
    for(int i = 0; i < ntotn; ++i){
        for(int j = 0; j < nstate; ++j){
          cout << "\t" << L_sk_k[i][j];
        }
        cout << endl;
    }
}
