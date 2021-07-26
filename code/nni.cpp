#include "nni.hpp"

// Two NNIs are considered conflicting if they operate on the same inner branch or adjacent branches.
void get_compatible_NNIs(vector<NNIMove>& nniMoves, vector<NNIMove>& compatibleNNIs){
    int debug = 0;
    compatibleNNIs.clear();

  	for(vector<NNIMove>::iterator it1 = nniMoves.begin(); it1 != nniMoves.end(); it1++){
  		bool select = true;
  		for(vector<NNIMove>::iterator it2 = compatibleNNIs.begin(); it2 != compatibleNNIs.end(); it2++){
  			if((*it1).node1 == (*(it2)).node1
  					|| (*it1).node2 == (*(it2)).node1
  					|| (*it1).node1 == (*(it2)).node2
  					|| (*it1).node2 == (*(it2)).node2){
  		        select = false;
              break;
        }
      }
  		if(select){
        compatibleNNIs.push_back(*it1);
      }
    }

    if(debug){
        cout << "There are " << compatibleNNIs.size() << " compatible NNIs" << endl;
        for(int i = 0; i < compatibleNNIs.size(); i++){
            cout << compatibleNNIs[i].node1->id + 1 << "\t" << compatibleNNIs[i].node2->id + 1 << endl;
        }
    }
}


// check age of the sibling of node
// find the top node u with one child v and another child c, only feasible when t(v) > t(c), according to P293, Yang, 2014
bool is_valid_NNI(const evo_tree& rtree, const Branch& curBranch){
    int debug = 0;

    int id_v = curBranch.second->id;  // chosen branch for NNI (u,v)

    // find another child of node
    int id_c = 0;
    for(auto child : curBranch.first->daughters){
      if(child != id_c){
        id_c = child;
        break;
      }
    }

    double tv = rtree.nodes[id_v].age;
    double tc = rtree.nodes[id_c].age;

    if(debug){
      cout << "checking neighbors of " << curBranch.first->id << endl;
      cout << id_c << "\t" << id_v << endl;
      cout << tc << "\t" << tv << endl;
    }

    if(tv - tc > BLEN_MIN){
      return true;
    }else{
      return false;
    }
}


// Get all branches where NNI is feasible, namely all internal branches
// NNI is only feasible when it does not violate parent-child age constraint
void get_NNI_branches(const evo_tree& rtree, Branches& nniBranches, Node* node, Node* dad){
    int debug = 0;
    assert(node != NULL);

    // find neighbor of "node" which is not "dad"
    FOR_NEIGHBOR_IT(node, dad, it){
        if(is_inner_branch((*it)->node, node)){
            Branch curBranch(node, (*it)->node);
            assert((*it)->direction == AWAYFROM_ROOT);  // node is parent

            bool is_valid = is_valid_NNI(rtree, curBranch);

            if(is_valid){
              if(debug) cout << "Only add NNI branch when tc > tv" << endl;
              int branchID = pairInteger(curBranch.first->id, curBranch.second->id);
              nniBranches.insert(pair<int, Branch>(branchID, curBranch));
            }
        }

        get_NNI_branches(rtree, nniBranches, (*it)->node, node);
    }
}


// Find branches at most "depth" branches away from the tagged branch
void get_neighbor_inner_branches(const evo_tree& rtree, Node* node, Node* dad, int depth, Branches& surrBranches){
    int debug = 0;

    if(debug) cout << "Find branches at most " << depth << " branches away from the tagged branch" << endl;

    if(depth == 0)
      return;

      FOR_NEIGHBOR_IT(node, dad, it){
          if(!(*it)->node->is_leaf() && is_inner_branch((*it)->node, node)){
              Branch curBranch(node, (*it)->node);

              if((*it)->direction == TOWARD_ROOT){
                Node* tmp = curBranch.first;
                curBranch.first = curBranch.second;
                curBranch.second = tmp;
              }

              bool is_valid = is_valid_NNI(rtree, curBranch);

              int branchID = pairInteger(node->id, (*it)->node->id);
              if(is_valid && surrBranches.find(branchID) == surrBranches.end())
                  surrBranches.insert(pair<int, Branch>(branchID, curBranch));

              get_neighbor_inner_branches(rtree, (*it)->node, node, depth - 1, surrBranches);
          }
      }
}


// Find branches at most two branches away from the tagged branch, used when doing reduced NNI
void filter_NNI_branches(const evo_tree& rtree, vector<NNIMove>& appliedNNIs, Branches& nniBranches){
    int debug = 0;

    if(debug){
        cout << "Find branches at most two branches away from the tagged branch" << endl;
    }

    for(vector<NNIMove>::iterator it = appliedNNIs.begin(); it != appliedNNIs.end(); it++){
        if(debug) cout << "applied NNI " << it->node1->id + 1 << ", " << it->node2->id + 1 << endl;

        Branch curBranch(it->node1, it->node2);

        if(is_valid_NNI(rtree, curBranch)){
          int branchID = pairInteger(it->node1->id, it->node2->id);
          if(nniBranches.find(branchID) == nniBranches.end())
              nniBranches.insert(pair<int, Branch>(branchID, curBranch));
        }

        // find neighbors of node1 except node2
        get_neighbor_inner_branches(rtree, it->node1, it->node2, 2, nniBranches);
        get_neighbor_inner_branches(rtree, it->node2, it->node1, 2, nniBranches);
    }
}


// Apply one NNI move
// Only affect topology, excepting to adjust branch length of related branches to satifisty time constaints when cons is true
// Assume neighbors have been generated for the tree
void do_one_NNI(evo_tree& rtree, NNIMove& move, int cons){
    int debug = 0;

    // (node1, node2) is the chosen NNI branch, (u, v)
    Node* node1 = move.node1;
    Node* node2 = move.node2;
    assert(node1->degree() == 3 && node2->degree() == 3);

    if(debug){
        cout << "tree address for node " << node1->id + 1 << " is " << &rtree.nodes[node1->id] << ", move address is " << node1 << endl;
        cout << "tree address for node " << node2->id + 1 << " is " << &rtree.nodes[node2->id] << ", move address is " << node2 << endl;
        cout << " Apply one NNI move " << node1->id + 1 << "," << node2->id + 1 << endl;
    }

    NeighborVec::iterator node1Nei_it = move.node1Nei_it;
    NeighborVec::iterator node2Nei_it = move.node2Nei_it;
    Neighbor* node1Nei = *(node1Nei_it);
    Neighbor* node2Nei = *(node2Nei_it);

    NeighborVec::iterator it;
    if(debug){
        cout << " old neighbors of node " << node1->id + 1 << endl;
        FOR_NEIGHBOR(node1, node2, it)
            cout << "\t" << (*it)->node->id + 1 << "\t" << (*it)->length;
        cout << endl;
        cout << " old neighbors of node " << node2->id + 1 << endl;
        FOR_NEIGHBOR(node2, node1, it)
            cout << "\t" << (*it)->node->id + 1 << "\t" << (*it)->length;
        cout << endl;
        cout << "tree before NNI" << endl;
        // rtree.print();
        cout << rtree.make_newick() << endl;
        rtree.print_neighbors();
    }

    // do the NNI swap
    node1->updateNeighbor(node1Nei_it, node2Nei);
    node2Nei->node->updateNeighbor(node2, node1);
    node2->updateNeighbor(node2Nei_it, node1Nei);
    node1Nei->node->updateNeighbor(node1, node2);

    // update edges (u, c) to (u, a)
    int eid_uc = rtree.get_edge_id(node1->id, node1Nei->node->id);
    double blen_uc = rtree.edges[eid_uc].length;
    rtree.edges[eid_uc].end = node2Nei->node->id;
    rtree.nodes[node2Nei->node->id].e_in = eid_uc;
    // update node daughters and parent
    rtree.nodes[node1->id].daughters.erase(std::remove(rtree.nodes[node1->id].daughters.begin(), rtree.nodes[node1->id].daughters.end(),
        node1Nei->node->id), rtree.nodes[node1->id].daughters.end());
    rtree.nodes[node1->id].daughters.push_back(node2Nei->node->id);
    rtree.nodes[node2Nei->node->id].parent = node1->id;

    // update edges (v, a) to (v, c)
    int eid_va = rtree.get_edge_id(node2->id, node2Nei->node->id);
    double blen_va = rtree.edges[eid_va].length;
    rtree.edges[eid_va].end = node1Nei->node->id;
    rtree.nodes[node1Nei->node->id].e_in = eid_va;
    // update node daughters and parent
    rtree.nodes[node2->id].daughters.erase(std::remove(rtree.nodes[node2->id].daughters.begin(), rtree.nodes[node2->id].daughters.end(), 
        node2Nei->node->id), rtree.nodes[node2->id].daughters.end());
    rtree.nodes[node2->id].daughters.push_back(node1Nei->node->id);
    rtree.nodes[node1Nei->node->id].parent = node2->id;

    // node times and ages do not change, change length accordingly
    // branch length should be positive if NNI is valid
    if(cons){
        int eid_uv = rtree.get_edge_id(node1->id, node2->id);
        double blen_uv = rtree.edges[eid_uv].length;
        double blen_vc = blen_uc - blen_uv;
        double blen_ua = blen_uv + blen_va;
        rtree.edges[eid_va].length = blen_vc;
        rtree.edges[eid_uc].length = blen_ua;

        // update neighbor lengths for (v, c)
        Neighbor* node1_node2_nei = (Neighbor*) node2->findNeighbor(node1Nei->node);
        Neighbor* node2_node1_nei = (Neighbor*) node1Nei->node->findNeighbor(node2);
        node1_node2_nei->setLength(blen_vc);
        node2_node1_nei->setLength(blen_vc);

         // update neighbor lengths for (u, a)
        node1_node2_nei = (Neighbor*) node1->findNeighbor(node2Nei->node);
        node2_node1_nei = (Neighbor*) node2Nei->node->findNeighbor(node1);
        node1_node2_nei->setLength(blen_ua);
        node2_node1_nei->setLength(blen_ua);       
    }

    if(debug){
        cout << " current neighbors of node " << node1->id + 1 << endl;
        FOR_NEIGHBOR(node1, node2, it)
            cout << "\t" << (*it)->node->id + 1 << "\t" << (*it)->length;
        cout << endl;
        cout << " current neighbors of node " << node2->id + 1 << endl;
        FOR_NEIGHBOR(node2, node1, it)
            cout << "\t" << (*it)->node->id + 1 << "\t" << (*it)->length;
        cout << endl;
        cout << " tree after topology NNI move " << endl;
        // rtree.print();
        cout << rtree.make_newick() << endl;
        // rtree.print_neighbors();
    }

    assert(rtree.is_blen_valid());
}


// Update branch lengths related to the NNI move
// include four adjacent branches when nni5 is true
// used in do_all_NNIs
void change_NNI_Brans(evo_tree& rtree, NNIMove& nnimove, bool nni5){
  int debug = 0;
  if(debug) cout << "  Update branch lengths related to the NNI move" << endl;

  Node* node1 = nnimove.node1;
  Node* node2 = nnimove.node2;

  NeighborVec::iterator it;
  Neighbor* node1_node2_nei = (Neighbor*) node1->findNeighbor(node2);
  Neighbor* node2_node1_nei = (Neighbor*) node2->findNeighbor(node1);
  node1_node2_nei->setLength(nnimove.newLen[0]);
  node2_node1_nei->setLength(nnimove.newLen[0]);
  // cout << node1_node2_nei->node->id <<  "\t" << node1_node2_nei->length << endl;

  // update node times and ages later as a whole
  int eid = rtree.get_edge_id(node1->id, node2->id);
  double delta = nnimove.newLen[0][0] - rtree.edges[eid].length;
//   rtree.update_node_age(node1->id, delta);
//   rtree.update_node_time(node2->id, delta);
  rtree.edges[eid].length = nnimove.newLen[0][0];

  if(debug)
  {
      cout << " old neighbors of node " << node1->id + 1 << endl;
      FOR_NEIGHBOR(node1, node2, it)
          cout << "\t" << (*it)-> node->id + 1 << "\t" << (*it)->id + 1 << "\t" << (*it)->length;
      cout << endl;
      cout << " old neighbors of node " << node2->id + 1 << endl;
      FOR_NEIGHBOR(node2, node1, it)
          cout << "\t" << (*it)-> node->id + 1 << "\t" << (*it)->id + 1 << "\t" << (*it)->length;
      cout << endl;
      // cout << "branch length difference " << delta << endl;
  }

  if(nni5){
    int i = 1;
    Neighbor* nei;
    Neighbor* nei_back;
    NeighborVec::iterator it;

    FOR_NEIGHBOR(node1, node2, it)
    {
    	nei = (*it)->node->findNeighbor(node1);
    	nei_back = (node1)->findNeighbor((*it)->node);
    	nei->setLength(nnimove.newLen[i]);
    	nei_back->setLength(nnimove.newLen[i]);

        eid = rtree.get_edge_id(node1->id, (*it)->node->id);
        rtree.edges[eid].length = nnimove.newLen[i][0];
        if(debug){
            cout << "update edge " << eid << " with node " << node1->id+1 << "\t" << (*it)->node->id+1 << " by length " << nnimove.newLen[i][0] << endl;
        }
        // // update node times and ages
        // double delta = nnimove.newLen[i][0] - rtree.edges[eid].length;
        // if(fabs(delta) > SMALL_VAL){
        //     rtree.update_node_age(node1->id, delta);
        //     rtree.update_node_time((*it)->node->id, delta);
        // }

    	i++;
    }

    FOR_NEIGHBOR(node2, node1, it)
    {
    	nei = (*it)->node->findNeighbor(node2);
    	nei_back = (node2)->findNeighbor((*it)->node);
    	nei->setLength(nnimove.newLen[i]);
    	nei_back->setLength(nnimove.newLen[i]);

        eid = rtree.get_edge_id(node2->id, (*it)->node->id);
        rtree.edges[eid].length = nnimove.newLen[i][0];
        if(debug) cout << "update edge " << eid << " with node " << node2->id+1 << "\t" << (*it)->node->id+1 << " by length " << nnimove.newLen[i][0] << endl;
        // // update node times and ages
        // double delta = nnimove.newLen[i][0] - rtree.edges[eid].length;
        // if(fabs(delta) > SMALL_VAL){
        //     // update node times and ages
        //     double delta = nnimove.newLen[0][0] - rtree.edges[eid].length;
        //     rtree.update_node_age(node1->id, delta);
        //     rtree.update_node_time((*it)->node->id, delta);
        // }

    	i++;
    }
  }

  if(debug){
      cout << " current neighbors of node " << node1->id + 1 << endl;
      FOR_NEIGHBOR(node1, node2, it)
          cout << "\t" << (*it)->node->id + 1 << "\t" << (*it)->id + 1 << "\t" << (*it)->length;
      cout << endl;
      cout << " current neighbors of node " << node2->id + 1 << endl;
      FOR_NEIGHBOR(node2, node1, it)
          cout << "\t" << (*it)->node->id + 1 << "\t" << (*it)->id + 1 << "\t" << (*it)->length;
      cout << endl;
      cout << " tree after branch length change " << rtree.make_newick() <<  endl;
      // rtree.print();
  }
}


// Simultaneously apply all NNIs, assigning new branch lengths to related branches
// Mutation rates are estimated at the same time if maxj = 1
void do_all_NNIs(evo_tree& rtree, vector<NNIMove>& compatibleNNIs, bool changeBran, bool nni5, int cons){
    int debug = 0;
    if(debug) cout << "\nSimultaneously apply all NNIs (" << compatibleNNIs.size() << ")" << endl;

    for(vector<NNIMove>::iterator it = compatibleNNIs.begin(); it != compatibleNNIs.end(); it++){
      if(debug){
        cout << " applying NNI Move " << it->node1->id + 1 << ", " << it->node2->id + 1 << endl;
      }

	  do_one_NNI(rtree, *it, cons);

      if(changeBran){
            // apply new branch lengths obtained from optimization of individual branches
            if(debug){
                cout << " apply new branch lengths obtained from optimization of individual branches" << endl;
            }
		    change_NNI_Brans(rtree, *it, nni5);
            // update age and time to ensure all variables are consistent
            rtree.calculate_node_times();
            rtree.calculate_age_from_time();
      }
    }
}


NNIMove get_random_NNI(Branch& branch, gsl_rng* r){
    assert(is_inner_branch(branch.first, branch.second));
    // for rooted tree
    if(((Neighbor*)branch.first->findNeighbor(branch.second))->direction == TOWARD_ROOT){
        // swap node1 and node2 if the direction is not right, only for nonreversible models
        Node* tmp = branch.first;
        branch.first = branch.second;
        branch.second = tmp;
    }
    NNIMove nni;
    nni.node1 = (Node*) branch.first;
    nni.node2 = (Node*) branch.second;

    FOR_NEIGHBOR_IT(branch.first, branch.second, node1NeiIt){
      if(((Neighbor*)*node1NeiIt)->direction != TOWARD_ROOT){
          nni.node1Nei_it = node1NeiIt;
          break;
      }
    }

    int randInt = gsl_rng_uniform_int(r, branch.second->neighbors.size()-1);
    int cnt = 0;
    FOR_NEIGHBOR_IT(branch.second, branch.first, node2NeiIt){
        // if this loop, is it sure that direction is away from root because node1->node2 is away from root
        if(cnt == randInt){
            nni.node2Nei_it = node2NeiIt;
            break;
        }else{
            cnt++;
        }
    }
	  assert(*nni.node1Nei_it != NULL && *nni.node2Nei_it != NULL);
    assert(((Neighbor*)*nni.node1Nei_it)->direction != TOWARD_ROOT && ((Neighbor*)*nni.node2Nei_it)->direction != TOWARD_ROOT);
    nni.newloglh = 0.0;

    return nni;
}


// do stochastic NNI
void do_random_NNIs(evo_tree& rtree, gsl_rng* r, int cons){
    int debug = 0;
    if(debug) cout << "doing stochastic NNIs" << endl;

    int cntNNI = 0;
    int numRandomNNI;
    Branches nniBranches;
    // half the number of internal branches in a rooted tree
    numRandomNNI = floor((rtree.nleaf - 2) * 0.5);

    while (cntNNI < numRandomNNI){
        nniBranches.clear();
        get_NNI_branches(rtree, nniBranches, &rtree.nodes[rtree.nleaf], NULL);

        if(nniBranches.size() == 0) break;

        // Convert the map data structure Branches to vector of Branch
        vector<Branch> vectorNNIBranches;
        int i = 0;
        for(Branches::iterator it = nniBranches.begin(); it != nniBranches.end(); ++it){
            vectorNNIBranches.push_back(it->second);
            i++;
            // cout << "NNI branch " << i << " is " << it->second->first->id << ", " << it->second->second->id << endl;
        }

        int randInt = gsl_rng_uniform_int(r, vectorNNIBranches.size());
        NNIMove randNNI = get_random_NNI(vectorNNIBranches[randInt], r);

        do_one_NNI(rtree, randNNI, cons);

        cntNNI++;
    }
    if(debug)
	    cout << "Tree perturbation: number of random NNI performed = " << cntNNI << endl;

}


/**
   Search for the best NNI move corresponding to the chosen inner branch
   @return NNIMove the best NNI, this NNI could be worse than the current tree
   according to the evaluation scheme in use
   @param node1 1 of the 2 nodes on the branch
   @param node2 1 of the 2 nodes on the branch
   @param nniMoves (IN/OUT) detailed information of the 2 NNIs
   adapted from IQ-TREE package, phylotree.cpp
 */
NNIMove get_best_NNI_for_bran(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, 
    LNL_TYPE& lnl_type, OPT_TYPE& opt_type, NNIMove* nniMoves, bool nni5){
    int debug = 0;

    if(debug){
        cout << "\nComputing the approximate likelihood of tree " << rtree.make_newick() << endl;
        // rtree.print_neighbors();
    }

    NeighborVec::iterator it, node1_it, node2_it;

    bool newNNIMoves = false;
    NNIMove* nniMoves2 = new NNIMove[2];    // moves on replicated trees
    if(!nniMoves){
        // Initialize the 2 NNI moves
        newNNIMoves = true;
        nniMoves = new NNIMove[2];
        nniMoves[0].node1 = nniMoves[1].node1 = NULL;
        nniMoves2[0].node1 = nniMoves2[1].node1 = NULL;
    }

    // create NNI moves from original tree, so that the pointers are for the original tree
    Node* node1 = &rtree.nodes[rtree.node1->id];
    Node* node2 = &rtree.nodes[rtree.node2->id];
    // assert(!node1->is_leaf() && !node2->is_leaf());
    assert(node1->degree() == 3 && node2->degree() == 3);
    if(((Neighbor*)node1->findNeighbor(node2))->direction == TOWARD_ROOT){
        // swap node1 and node2 if the direction is not right, only for nonreversible models
        // cout << "swapping node " << node1->id + 1 << ", " << node2->id + 1 << endl;
        Node* tmp = node1;
        node1 = node2;
        node2 = tmp;
    }
    // cout << "(node1) has id " << node1->id + 1 << " address " << node1 << " neighbor " << node1->neighbors.size() << endl;
    // cout << "(node2) has id " << node2->id + 1 << " address " << node2 << " neighbor " << node2->neighbors.size() << endl;

    int cnt_move = 0;
    if(nniMoves[0].node1 && nniMoves[1].node1){
        // assuming that node1Nei_it and node2Nei_it are defined in nniMoves structure
        for(cnt_move = 0; cnt_move < 2; cnt_move++){
            // sanity check
            if(!node1->findNeighbor((*nniMoves[cnt_move].node1Nei_it)->node) ||
            !node2->findNeighbor((*nniMoves[cnt_move].node2Nei_it)->node)){
              cout << "Error in neighbors" << endl;
            }
            // copy original nniMoves
            nniMoves2[cnt_move] = nniMoves[cnt_move];
        }
    }else{
        if(debug) cout << "\ncreating NNI move" << endl;
        FOR_NEIGHBOR_IT(node1, node2, node1_it){
            if(((Neighbor*)*node1_it)->direction != TOWARD_ROOT){
                cnt_move = 0;
                FOR_NEIGHBOR_IT(node2, node1, node2_it){
                    //  Initialize the 2 NNI moves
                    nniMoves[cnt_move].node1Nei_it = node1_it;
                    nniMoves[cnt_move].node2Nei_it = node2_it;
                    cnt_move++;
                }
                break;
            }
        }
        assert(cnt_move == 2);
    }

    // Initialize node1 and node2 in nniMoves
    nniMoves[0].node1 = nniMoves[1].node1 = node1;
    nniMoves[0].node2 = nniMoves[1].node2 = node2;
    nniMoves[0].newloglh = nniMoves[1].newloglh = -DBL_MAX;

    if(debug){
        cout << " NNI move 0 " << nniMoves[0].node1->id + 1  << "\t" << nniMoves[0].node2->id + 1 << "\t" << ((Neighbor*)(*nniMoves[0].node1Nei_it))->node->id + 1 << "\t" << ((Neighbor*)(*nniMoves[0].node2Nei_it))->node->id + 1 << endl;
        cout << " NNI move 1 " << nniMoves[1].node1->id + 1  << "\t" << nniMoves[1].node2->id + 1 << "\t" << ((Neighbor*)(*nniMoves[1].node1Nei_it))->node->id + 1 << "\t" << ((Neighbor*)(*nniMoves[1].node2Nei_it))->node->id + 1 << endl;
    }

    // make two copies to avoid changing original tree. Not updating node times
    evo_tree rtree1(rtree);     // for 1st move
    rtree1.generate_neighbors();
    evo_tree rtree2(rtree);     // for 2nd move
    rtree2.generate_neighbors();
    vector<evo_tree*> trees{&rtree1, &rtree2};

    int cons = lnl_type.cons;
    int maxj = opt_type.maxj;
    double tolerance = opt_type.tolerance;
   
    // 2 moves associated with the same branch (corresponding to rtree1 and rtree2 respectively)
    for(int cnt = 0; cnt < 2; cnt++){
        if(debug){
            cout << "\nstart a move " << cnt << endl;
            // trees[cnt]->print();
            cout << "tree before move " << trees[cnt]->make_newick() << endl;
        }

        node1 = &trees[cnt]->nodes[rtree.node1->id];
        node2 = &trees[cnt]->nodes[rtree.node2->id];

        if(((Neighbor*)node1->findNeighbor(node2))->direction == TOWARD_ROOT){
            // swap node1 and node2 if the direction is not right, only for nonreversible models
            // cout << "swapping node " << node1->id + 1 << ", " << node2->id + 1 << endl;
            Node* tmp = node1;
            node1 = node2;
            node2 = tmp;
        }

        // create nniMoves for this tree
        if(nniMoves2[cnt].node1 == NULL){
            if(debug) cout << "\ncreating copied NNI move " << cnt << endl;
            FOR_NEIGHBOR_IT(node1, node2, node1_it){
                // cout << "checking node 1" << endl;
                // Node* curr = ((Neighbor*)*node1_it)->node;
                // cout << curr << " has id " << curr->id + 1 << endl;
                if(((Neighbor*)*node1_it)->direction != TOWARD_ROOT)
                {
                    cnt_move = 0;
                    FOR_NEIGHBOR_IT(node2, node1, node2_it){
                        // cout << "checking node 2" << endl;
                        // Initialize the 2 NNI moves
                        // Node* curr2 = ((Neighbor*)*node2_it)->node;
                        // cout << curr2 << " has id " << curr2->id + 1 << endl;
                        if(cnt_move != cnt){
                            cnt_move++;
                            continue;
                        }
                        nniMoves2[cnt_move].node1Nei_it = node1_it;
                        nniMoves2[cnt_move].node2Nei_it = node2_it;
                        cnt_move++;
                    }
                    break;
                }
            }
            assert(cnt_move == 2);
        }

        // Initialize node1 and node2 in nniMoves
      	nniMoves2[cnt].node1 = node1;
      	nniMoves2[cnt].node2 = node2;
        nniMoves2[cnt].newloglh = -DBL_MAX;

        bool is_valid = true;
        if(cons){
          Branch nni_bran(node1, node2);
          is_valid = is_valid_NNI(*trees[cnt], nni_bran);
        }

        if(debug){
            cout << " NNI move " << cnt << "\t" << nniMoves2[cnt].node1->id + 1  << "\t" << nniMoves2[cnt].node2->id + 1 << "\t" << ((Neighbor*)(*nniMoves2[cnt].node1Nei_it))->node->id + 1 << "\t" << ((Neighbor*)(*nniMoves2[cnt].node2Nei_it))->node->id + 1 << endl;
            cout << " \ndoing NNI move " << cnt << " between node " << node1->id + 1 << " and " << node2->id + 1 << endl;
            cout << "tree before NNI " << trees[cnt]->make_newick() << endl;

            // node1_it = nniMoves2[cnt].node1Nei_it;
            // node2_it = nniMoves2[cnt].node2Nei_it;
            // Neighbor* node1_nei = *node1_it;
            // Neighbor* node2_nei = *node2_it;
            // cout << " old neighbors of node " << node1->id + 1 << endl;
            // FOR_NEIGHBOR(node1, node2, it)
            //     cout << "\t" << (*it)->node->id + 1;
            // cout << endl;
            // cout << " old neighbors of node " << node2->id + 1 << endl;
            // FOR_NEIGHBOR(node2, node1, it)
            //     cout << "\t" << (*it)->node->id + 1;
            // cout << endl;
            // cout << " old length of node1, node2 is " << node1->findNeighbor(node2)->length << endl;
            // cout << node1_nei->node->id + 1 << "\t" << node2_nei->node->id + 1 << endl;
        }

        // do the NNI swap if the move is valid
        // neighbors may change after NNI
        if(is_valid){
            do_one_NNI(*trees[cnt], nniMoves2[cnt], cons);
            vector<int> inodes;
            Node* root = &(trees[cnt]->nodes[trees[cnt]->root_node_id]);
            trees[cnt]->get_inodes_postorder(root, inodes);
            lnl_type.knodes = inodes;
        }else{
            cout << "NNI not valid, checking next one ......" << endl;
        }

        if(debug){
            cout << "after NNI move " << cnt << endl;
            // cout << " current neighbors of node " << node1->id + 1<< endl;
            // for(it = node1->neighbors.begin(); it != node1->neighbors.end(); it++)
            //     cout << "\t" << (*it)->node->id + 1;
            // cout << endl;
            // cout << " current neighbors of node " << node2->id + 1 << endl;
            // for(it = node2->neighbors.begin(); it != node2->neighbors.end(); it++)
            //     cout << "\t" << (*it)->node->id + 1;
            // cout << endl;
            // trees[cnt]->print();
            cout << "tree after NNI " << trees[cnt]->make_newick() << endl;
        }

        int nni5_num_eval = 1;   // Number of steps for the loop evaluating 5 branches around NNI, fix for now for convenience
        double score = 0.0;
        for(int step = 0; step < nni5_num_eval; step++){
      	    int i = 1;

            if(nni5){
                if(debug) cout << "\noptimizing neighbors of node 1 (excluding node2): " << i << endl;
                FOR_NEIGHBOR(node1, node2, it){
                    if(debug){
                        // Node* n1_nei = (*it)->node;
                        // cout << node1->id + 1 << "\t" << n1_nei->id + 1 << endl;
                        // cout << "tree address before " << &trees[cnt] << endl;
                        // trees[cnt]->print();
                        cout << "tree before " << trees[cnt]->make_newick() << endl;
                    }
                    if(cons){
                        if(is_inner_branch(node1, (Node* )(*it)->node))
                            score = optimize_one_branch_BFGS(*trees[cnt], vobs, obs_decomp, comps, lnl_type, opt_type, node1, (Node* )(*it)->node);
                    }else{
                        if(!(node1->id == trees[cnt]->nleaf - 1 || ((Node* )(*it)->node)->id == trees[cnt]->nleaf - 1))
                            score = optimize_one_branch(*trees[cnt], vobs, obs_decomp, comps, lnl_type, tolerance, node1, (Node* )(*it)->node);
                    }

                    // if(cons && !is_tip_age_valid(trees[cnt]->get_node_ages(), opt_type.tobs)){
                    //   trees[cnt]->print();
                    //   cout << trees[cnt]->make_newick() << endl;
                    //   cout << "Tip timings incorrect after optimizing neighbor of " << node1->id + 1 << endl;
                    //   exit(1);
                    // }
                    if(debug){
                        cout << "tree after " << trees[cnt]->make_newick() << endl;
                    }
                    node1->findNeighbor((*it)->node)->getLength(nniMoves2[cnt].newLen[i]);
                    nniMoves[cnt].newLen[i] = nniMoves2[cnt].newLen[i];
                    i++;
                }
            }

            if(debug) cout << "\noptimizing (node 1, node2): " << node1->id + 1 << ", " << node2->id + 1 << endl;

            if(cons){
                assert(is_inner_branch(node1, node2));
                score = optimize_one_branch_BFGS(*trees[cnt], vobs, obs_decomp, comps, lnl_type, opt_type, node1, node2);
            }else{
                score = optimize_one_branch(*trees[cnt], vobs, obs_decomp, comps, lnl_type, tolerance, node1, node2);
            }

            // if(cons && !is_tip_age_valid(trees[cnt]->get_node_ages(), opt_type.tobs)){
            //   trees[cnt]->print();
            //   cout << trees[cnt]->make_newick() << endl;
            //   cout << "Tip timings incorrect after optimizing neighbor of " << node1->id + 1 << ", " << node2->id + 1 << endl;
            //   exit(1);
            // }

            node1->findNeighbor(node2)->getLength(nniMoves2[cnt].newLen[0]);
            nniMoves[cnt].newLen[0] = nniMoves2[cnt].newLen[0];

            if(nni5){
                if(debug) cout << "\noptimizing neighbors of node 2 (excluding node1): " << i << endl;
                FOR_NEIGHBOR(node2, node1, it){
                    if(cons){
                        if(is_inner_branch(node2, (Node* )(*it)->node))
                            score = optimize_one_branch_BFGS(*trees[cnt], vobs, obs_decomp, comps, lnl_type, opt_type, node2, (Node* )(*it)->node);
                    }else{
                        if(!(node2->id == trees[cnt]->nleaf - 1 || ((Node* )(*it)->node)->id == trees[cnt]->nleaf - 1))
                            score = optimize_one_branch(*trees[cnt], vobs, obs_decomp, comps, lnl_type, tolerance, node2, (Node* )(*it)->node);
                    }

                    // if(cons && !is_tip_age_valid(trees[cnt]->get_node_ages(), opt_type.tobs)){
                    //   trees[cnt]->print();
                    //   cout << trees[cnt]->make_newick() << endl;
                    //   cout << "Tip timings incorrect after optimizing neighbor of " << node2->id + 1 << endl;
                    //   exit(1);
                    // }

                    node2->findNeighbor((*it)->node)->getLength(nniMoves2[cnt].newLen[i]);
                    nniMoves[cnt].newLen[i] = nniMoves2[cnt].newLen[i];

                    i++;
                }
            }
        }

        // double score = -DBL_MAX;
        // if(lnl_type.model == DECOMP){
        //     score = get_likelihood_decomp(*trees[cnt], vobs, obs_decomp, comps, lnl_type);
        // }else{
        //     score = get_likelihood_revised(*trees[cnt], vobs, lnl_type);
        // }
        nniMoves[cnt].newloglh = score;

        if(debug){
            cout << "\nfinish NNI " << node1->id + 1 << " - " << node2->id + 1 << ", approximate likelihood " << nniMoves[cnt].newloglh << ", branch lengths: ";
            for(int j = 0; j < 5; j++){
                cout << "\t" << nniMoves[cnt].newLen[j][0];
            }
            cout << endl;
            // trees[cnt]->print();
            cout << "tree after NNI " << trees[cnt]->make_newick() << endl;
        }
     }

    NNIMove res;
    if(nniMoves[0].newloglh > nniMoves[1].newloglh){
      res = nniMoves[0];
    }else{
      res = nniMoves[1];
    }

  	if(newNNIMoves){
  		delete [] nniMoves;
  	}
    delete [] nniMoves2;

    rtree1.delete_neighbors();
    rtree2.delete_neighbors();
    // cout << "tree address for node " << node1->id+1 << " is " << &rtree.nodes[node1->id] << ", move address is " << res.node1 << endl;
    // cout << "tree address for node " << node2->id+1 << " is " << &rtree.nodes[node2->id] << ", move address is " << res.node2 << endl;

	return res;
}


// Find NNI increasing likelihood of current tree
void evaluate_NNIs(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, Branches& nniBranches, vector<NNIMove>& positiveNNIs, double curScore){
    int debug = 0;

    for(Branches::iterator it = nniBranches.begin(); it != nniBranches.end(); it++){
        rtree.node1 = (Node* )it->second.first;
        rtree.node2 = (Node* )it->second.second;
        NNIMove nni = get_best_NNI_for_bran(rtree, vobs, obs_decomp, comps, lnl_type, opt_type);

        if(debug) cout << "\n   evaluating NNI: " << nni.node1->id + 1 << ", " << nni.node2->id + 1 <<" with score " << nni.newloglh << endl;

        if(nni.newloglh > curScore){
            positiveNNIs.push_back(nni);
        }
    }

    if(debug){
        cout << "\nThere are " << positiveNNIs.size() << " positive NNIs" << endl;
        for(int i = 0; i < positiveNNIs.size(); i++){
            cout << "\t" << positiveNNIs[i].node1->id + 1 << ", " << positiveNNIs[i].node2->id + 1 << ", " << positiveNNIs[i].newloglh << endl;
        }
        cout << endl;
    }
}


// Apply hill climbing perturbation to obtain a locally optimal tree (by NNI)
// score used in this function is log likelihood, the larger the better
void do_hill_climbing_NNI(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double loglh_epsilon, int speed_nni, bool nni5){
    int debug = 0;

    int totalNNIApplied = 0;
    int numSteps = 0;
    int max_steps = rtree.nleaf - 1;
    int cons = lnl_type.cons;
    int maxj = opt_type.maxj;
    int miter = opt_type.miter;
    double tolerance = opt_type.tolerance;

    Branches nniBranches;
    vector<NNIMove> positiveNNIs;
    vector<NNIMove> appliedNNIs;
    bool changeBran = true;
    if(cons)    changeBran = false;  // not change branch lengths in NNIs when doing constrained optimization for convenience

    double curScore = rtree.score;
    double curBestScore = curScore;
    double originalScore = curScore;
    if(debug){
        cout << "Newick String for start tree is " << rtree.make_newick() << endl;
        cout << "score at the beginning of hill climbing NNIs " << curScore << endl;
    }

    for(numSteps = 1; numSteps <= max_steps; numSteps++){
        if(debug){
            cout << "\nStep " << numSteps << endl;
            // rtree.print();
            cout << "Newick String for current tree is " << rtree.make_newick() << endl;
        }

        double oldScore = curScore;

        bool startSpeedNNI = false;
        if(speed_nni && !appliedNNIs.empty()){
            startSpeedNNI = true;
        }

        // Do reduced NNI search.
        if(startSpeedNNI){
            if(debug) cout << "Doing reduced NNI search." << endl;
            // speedNNI option: only evaluate NNIs that are 2 branches away from the previously applied NNI
            Branches filteredNNIBranches;
            filter_NNI_branches(rtree, appliedNNIs, filteredNNIBranches);

            for(Branches::iterator it = filteredNNIBranches.begin(); it != filteredNNIBranches.end(); it++){
                Branch curBranch = it->second;

                if(debug) cout << "filtered NNI " << curBranch.first->id + 1 << ", " << curBranch.second->id + 1 << endl;

                Neighbor* nei = (Neighbor*) curBranch.first->findNeighbor(curBranch.second);
                int branchID = pairInteger(curBranch.first->id, curBranch.second->id);
                nniBranches.insert(pair<int, Branch>(branchID, curBranch));
            }
        }else{
            if(debug) cout << "Doing complete NNI search." << endl;
            get_NNI_branches(rtree, nniBranches, &rtree.nodes[rtree.nleaf], NULL);
        }

        if(debug){
            cout << "\nThere are " <<  nniBranches.size() << " NNI branches" << endl;
            for(auto p: nniBranches){
                int branchID = p.first;
                Branch curBranch = p.second;
                cout << branchID << "\t" << curBranch.first->id + 1 << "," << curBranch.second->id + 1 << endl;
            }
        }

        // Only consider NNIs that increase the likelihood of current tree
        positiveNNIs.clear();
        evaluate_NNIs(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, nniBranches, positiveNNIs, curScore);

        if(debug){
            cout << "\nSearching NNIs increasing likelihood of current tree with score " << curScore << endl;
            cout << "Newick String for current tree after evaluating NNI (same as original tree) is " << rtree.make_newick() << endl;
        }

        if(positiveNNIs.empty()){
            break;
        }

        // if(debug){
            // cout << "tree address for node  " << positiveNNIs[0].node1->id+1 << " is " << &rtree.nodes[positiveNNIs[0].node1->id] << endl;
            // cout << "tree address for node  " << positiveNNIs[0].node2->id+1 << " is " << &rtree.nodes[positiveNNIs[0].node2->id] << endl;
        // }

        /* sort all positive NNI moves (ASCENDING) */
        sort(positiveNNIs.begin(), positiveNNIs.end());
        /* remove conflicting NNIs */
        appliedNNIs.clear();
        if(debug) cout << "\nGetting compatible NNIs" << endl;
        get_compatible_NNIs(positiveNNIs, appliedNNIs);

        if(debug){
            cout << "\nThere are " << appliedNNIs.size() << " applied NNIs" << endl;
            for(int i = 0; i < appliedNNIs.size(); i++){
                cout << "\t" << appliedNNIs[i].node1->id + 1 << ", " << appliedNNIs[i].node2->id + 1 << ", " << appliedNNIs[i].newloglh << ", branch lengths:";
                for(int j = 0; j < 5; j++){
                    cout << "\t" << appliedNNIs[i].newLen[j][0];
                }
                cout << endl;
            }
            cout << endl;
        }
       
        nniBranches.clear();

        // save all current branch lengths (and mutation rates)
        DoubleVector lenvec;
        save_branch_lengths(rtree, lenvec);    // both edges and neighbors      
        DoubleVector muvec;
        if(maxj){
          save_mutation_rates(rtree, muvec);
        }

        // do non-conflicting positive NNIs
        do_all_NNIs(rtree, appliedNNIs, changeBran, nni5, cons);
        // get the correct node order for likelihood computation 
        vector<int> inodes;
        Node* root = &(rtree.nodes[rtree.root_node_id]);
        rtree.get_inodes_postorder(root, inodes);
        lnl_type.knodes = inodes;

        if(debug){
            double tscore = 0.0;
            if(lnl_type.model == DECOMP){
                tscore = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
            }else{
                tscore = get_likelihood_revised(rtree, vobs, lnl_type);
            }           
            cout << "score after applying all NNIs " << tscore << endl;
            // rtree.print();
            cout << "Newick String for current tree after applying all NNIs is " << rtree.make_newick() << endl;
        }

        if(cons){  // optimization of all branches under time constraints
            if(debug) cout << "\n1st global optimization" << endl;
            double min_nlnl = MAX_NLNL;
            opt_type.maxj = 0;
            while(-min_nlnl < curScore){
                min_nlnl = MAX_NLNL;
                max_likelihood_BFGS(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, min_nlnl);
            } 
            opt_type.maxj = 1;           
            curScore = -min_nlnl;

            // need to keep neighbours updated in NNI after updating lengths of all edges, which will be used in NNIs
            // if(debug){
            //     cout << "current tree after 1st global optimization " << rtree.make_newick() << endl;
            //     cout << "neighbors before updating lengths (may not be consistent with newick string)" << endl;
            //     rtree.print_neighbors();       
            // }           
            rtree.update_neighbor_lengths();
            // if(debug){
            //     cout << "neighbors after updating lengths (should be consistent with newick string)" << endl;
            //     rtree.print_neighbors();       
            // }             
        }else{
            // achieved by optimizing each branch in a specific order, 2 iterations by default
            curScore = optimize_all_branches(rtree, vobs, obs_decomp, comps, lnl_type, 2, loglh_epsilon);
        }

        if(maxj){
            curScore = optimize_mutation_rates(rtree, vobs, obs_decomp, comps, lnl_type, tolerance);
        }

        if(debug){
            cout << "score after global optimization " << curScore << endl;
            cout << "Newick String for current tree after global optimization is " << rtree.make_newick() << endl;
            cout << "Current estimation of mutation rates are: " << endl;
            rtree.print_mutation_rates(lnl_type.model, lnl_type.only_seg);
        }

        if(curScore < appliedNNIs.at(0).newloglh - loglh_epsilon){
            if(debug) cout << "Tree getting worse: curScore = " << curScore << " / best score = " <<  appliedNNIs.at(0).newloglh << endl;
            // tree cannot be worse if only 1 NNI is applied
            if(appliedNNIs.size() > 1){
                // revert all applied NNIs (topolgy will be restored but not branch lengths or mutation rates)
                do_all_NNIs(rtree, appliedNNIs, false, nni5, cons);

                if(debug){
                    cout << "tree after reverting all applied NNIs " << rtree.make_newick() << endl;
                    // rtree.print();
                    cout << "original branch lengths:";
                    for(auto len : lenvec){
                        cout << "\t" << len;
                    }
                    cout << endl;
                }

                restore_branch_lengths(rtree, lenvec);   // both edges and neigbors are updated 

                rtree.calculate_node_times();
                rtree.calculate_age_from_time();             

                if(debug){
                    cout << "tree after reverting recomputing branch lengths and node times " << rtree.make_newick() << endl;
                    // rtree.print();
                }

                if(cons && !is_age_time_consistent(rtree.get_node_times(), rtree.get_node_ages())){
                    cout << "Wrong node age or time after reverting NNIs" << endl;
                    rtree.print();
                    cout << rtree.make_newick() << endl;
                    exit(1);
                }
         
                if(maxj) restore_mutation_rates(rtree, muvec);

                // rtree.generate_neighbors();
                // only do the best NNI
                if(debug){
                    cout << "only do the best NNI for the original tree " << rtree.make_newick() << endl;
                    // rtree.print();
                }

                appliedNNIs.resize(1);
                do_all_NNIs(rtree, appliedNNIs, changeBran, nni5, cons);
                vector<int> inodes;
                Node* root = &(rtree.nodes[rtree.root_node_id]);
                rtree.get_inodes_postorder(root, inodes);
                lnl_type.knodes = inodes;

                if(debug){
                    double tscore = 0.0;
                    if(lnl_type.model == DECOMP){
                        tscore = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
                    }else{
                        tscore = get_likelihood_revised(rtree, vobs, lnl_type);
                    }                   
                    cout << "score after applying best NNI " << tscore << endl;
                    // rtree.print();
                    cout << "Newick String for current tree after applying best NNI is " << rtree.make_newick() << endl;
                }


                if(cons){
                    if(debug) cout << "2nd global optimization" << endl;
                    double min_nlnl = MAX_NLNL;
                    opt_type.maxj = 0;
                    while(-min_nlnl < appliedNNIs.at(0).newloglh){
                        min_nlnl = MAX_NLNL;
                        max_likelihood_BFGS(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, min_nlnl);
                    }
                    opt_type.maxj = 1;
                    curScore = -min_nlnl;
                    // if(debug){
                    //     cout << "current tree " << rtree.make_newick() << endl;
                    //     cout << "neighbors before updating lengths" << endl;
                    //     rtree.print_neighbors();       
                    // }           
                    rtree.update_neighbor_lengths();
                    // if(debug){
                    //     cout << "neighbors after updating lengths" << endl;
                    //     rtree.print_neighbors();       
                    // }  
                }else{
                    curScore = optimize_all_branches(rtree, vobs, obs_decomp, comps, lnl_type, 2, loglh_epsilon);
                }

                if(maxj){
                    curScore = optimize_mutation_rates(rtree, vobs, obs_decomp, comps, lnl_type, tolerance);
                }                

                if(debug){
                    cout << "current score " << curScore << ", best NNI score " << appliedNNIs.at(0).newloglh << endl;
                    cout << "Current estimation of mutation rates are: " << endl;
                    rtree.print_mutation_rates(lnl_type.model, lnl_type.only_seg);                    
                }
                assert(curScore > appliedNNIs.at(0).newloglh - 0.1);
            }else{
                if(debug) cout << "Applied all NNIs successfully " << endl;
                // failed sometimes
                if(debug) cout << "current score " << curScore << ", best NNI score " << appliedNNIs.at(0).newloglh << endl;
                assert(curScore > appliedNNIs.at(0).newloglh - 0.1 && "Using one NNI reduces LogL");
            }
            totalNNIApplied++;

            // update applied NNIs to point to the new tree
            appliedNNIs[0].node1 = &rtree.nodes[appliedNNIs[0].node1->id];
            appliedNNIs[1].node2 = &rtree.nodes[appliedNNIs[1].node2->id];
        }else{
            totalNNIApplied += appliedNNIs.size();
        }

        if(curScore < oldScore - loglh_epsilon){    // error
            cout << "$$$$$$$$: " << curScore << "\t" << oldScore << "\t" << curScore - oldScore << endl;

        }
        if(curScore - oldScore <  loglh_epsilon)   // no improvement
            break;
        if(curScore > curBestScore + 0.1){
            curBestScore = curScore;
        }
    }

    if(totalNNIApplied == 0 && debug){
        cout << "NOTE: Input tree is already NNI-optimal" << endl;
    }
    if(numSteps == max_steps){
        cout << "WARNING: NNI search needs unusual large number of steps (" << numSteps << ") to converge!" << endl;
    }
    if(curScore < originalScore - loglh_epsilon){ // error
        cout << "AAAAAAAAAAAAAAAAAAA: " << curScore << "\t" << originalScore << "\t" << curScore - originalScore << endl;
    }

    // Get the new tree (already changed when applying NNIs)
    rtree.score = curScore;

    if(debug){
        cout << "new tree after hill climbing NNI" << endl;
        rtree.print();
        string newick = rtree.make_newick();
        cout << "Newick String for current tree is " << newick << endl;
    }

}



/**************************** NNI operations used in MCMC approach ****************************/

// Use NNI for rooted (clock) trees
evo_tree do_random_NNI(evo_tree& rtree, gsl_rng* r,int debug){
    if(debug){
        printf ("Before:\n");
        rtree.print();
        // Output to a file to facilitate ploting
        // ofstream out_tree("./tree-before-NNI.txt");
        // rtree.write(out_tree);
        // out_tree.close();
    }

    evo_tree ntree(rtree);

    /* pick an interior branch (excluding the top edge connecting to root rtree.nleaf), around which it is possible to make an NNI */
    double tv = -1;
    double tc = -2;
    int i, u, v, c, a, b;
    edge *e_uv, *e_uc, *e_va;
    vector<int> children;
    double rnum;
    int count = 0;
    double blen_uc;
    double blen_va;
    set<double> chosen; // Stop if all the internal edges have been tested
    vector<edge*> intedges = rtree.get_internal_edges();
    int nintedge = rtree.nleaf - 2;

    while (tc <= tv && chosen.size() < nintedge - 1){
        count++;
        // Exclude the edges connecting to unaltered genome
        assert(intedges[nintedge - 1]->start == rtree.nleaf);
        i = gsl_rng_uniform_int(r, nintedge - 1);
        chosen.insert(i);
        if(debug){
            cout << "Picking " << i + 1 << "th internal edge: " << intedges[i]->start + 1 << ", " << intedges[i]->end + 1 << endl;
        }
        e_uv = intedges[i];
        // Find the children of the start node
        u = e_uv->start;
        children = rtree.nodes[u].daughters;
        // Select the child with subtree as v
        v = *max_element(children.begin(), children.end());
        c = *min_element(children.begin(), children.end());

        assert(v > rtree.nleaf);

        tv = rtree.nodes[v].time;
        tc = rtree.nodes[c].time;

        // Find the children of v
        children = rtree.nodes[v].daughters;
        assert(children.size() == 2);
        rnum = runiform(r, 0, 1);
        if(rnum > 0.5){
            a = children[0];
        }
        else{
            a = children[1];
        }
    }

    if(tc <= tv){
        if(debug) cout << "Cannot update topology!" << endl;
        return ntree;
    }
    if(debug){
        cout << "Picked an edge after " << count << " iteration(s)" << endl;
        cout << "u, v, a, c " << u + 1 << "\t" << v + 1 << "\t" << a + 1 << "\t" << c + 1 << endl;
    }

    /* record branch lengths */
    // oldALength = e_va->length;
    // oldCLength = e_uc->length;

    /* make topology change */
    if(debug){
        cout << "Making topology change" << endl;
    }
    // (u,c) --> (u,a); (v,a) --> (v,c)
    // move the subtree below
    // Find the length of edge uc and va
    for(int j = 0; j < ntree.edges.size(); j++){
        edge *e = &ntree.edges[j];
        if(e->start == u && e->end == c){
            blen_uc = e->length;
            e->end = a;
            e_uc = e;
        }
        if(e->start == v && e->end == a){
            blen_va = e->length;
            e->end = c;
            e_va = e;
        }
    }

    // update branch length to keep node ages unchanged
    e_va->length = blen_uc - e_uv->length;
    e_uc->length = e_uv->length + blen_va;
    // cout << " edge v->a length " << e_va->length << endl;
    // cout << " edge u->c length " << e_uc->length << endl;

    adjust_blen(e_va->length, BLEN_MIN, BLEN_MAX);
    adjust_blen(e_uc->length, BLEN_MIN, BLEN_MAX);


    if(debug){
        printf("After:\n");
        ntree.print();
        // Output to a file to facilitate ploting
        // ofstream out_tree("./test/tree-after-NNI.txt");
        // rtree.write(out_tree);
        // out_tree.close();
    }
    return ntree;
}


// Do NNI on an internal branch i, according to P293 (Yang, 2014)
vector<evo_tree> do_NNI_on_branch(evo_tree& rtree, int i, vector<NNIMove>& moves, int debug){
    if(debug){
        printf("Before:\n");
        rtree.print();
        // Output to a file to facilitate ploting
        // ofstream out_tree("./tree-before-NNI.txt");
        // rtree.write(out_tree);
        // out_tree.close();
    }
    vector<evo_tree> nni_trees;
    evo_tree tree1(rtree);
    evo_tree tree2(rtree);
    /* pick an interior branch (excluding the top edge connecting to root rtree.nleaf), around which it is possible to make an NNI */
    double tv, tc;
    int u, v, c, a1, a2, b;
    edge *e_uv, *e_uc, *e_va;
    vector<int> children;
    double blen_uc, blen_va;
    vector<edge*> intedges = rtree.get_internal_edges();
    tc = -2;
    tv = -1;

    if(debug){
        cout << "Picking " << i + 1 << "th internal edge: " << intedges[i]->start + 1 << ", " << intedges[i]->end + 1 << endl;
    }
    e_uv = intedges[i];
    // Find the children of the start node
    u = e_uv->start;
    children = rtree.nodes[u].daughters;
    // Select the child with subtree as v
    v = *max_element(children.begin(), children.end());
    c = *min_element(children.begin(), children.end());
    assert(v > rtree.nleaf);
    tv = rtree.nodes[v].time;
    tc = rtree.nodes[c].time;
    if(tc <= tv){
        cout << "Cannot update topology!" << endl;
        return nni_trees;
    }
    // Find the children of v
    children = rtree.nodes[v].daughters;
    assert(children.size() == 2);
    a1 = children[0];
    a2 = children[1];

    if(debug){
        cout << "u, v, a1, a2, c " << u + 1 << "\t" << v + 1 << "\t" << a1 + 1 << "\t" << a2 + 1 << "\t" << c + 1 << endl;
    }

    /* make topology change */
    // (u,c) --> (u,a); (v,a) --> (v,c)
    // move the subtree below
    // Find the length of edge uc and va
    for(int j = 0; j < tree1.edges.size(); j++){
        edge *e = &tree1.edges[j];
        if(e->start == u && e->end == c){
            blen_uc = e->length;
            e->end = a1;
            e_uc = e;
        }
        if(e->start == v && e->end == a1){
            blen_va = e->length;
            e->end = c;
            e_va = e;
        }
    }
    // update branch length to keep node ages unchanged
    e_va->length = blen_uc - e_uv->length;
    e_uc->length = e_uv->length + blen_va;
    adjust_blen(e_va->length, BLEN_MIN, BLEN_MAX);
    adjust_blen(e_uc->length, BLEN_MIN, BLEN_MAX);

    // NNIMove m1;
    // m1.swap_id = i;
    // m1.node1 = u;
    // m1.node2 = v;

    for(int j = 0; j < tree2.edges.size(); j++){
        edge *e = &tree2.edges[j];
        if(e->start == u && e->end == c){
            blen_uc = e->length;
            e->end = a2;
            e_uc = e;
        }
        if(e->start == v && e->end == a2){
            blen_va = e->length;
            e->end = c;
            e_va = e;
        }
    }
    // update branch length to keep node ages unchanged
    e_va->length = blen_uc - e_uv->length;
    e_uc->length = e_uv->length + blen_va;
    adjust_blen(e_va->length, BLEN_MIN, BLEN_MAX);
    adjust_blen(e_uc->length, BLEN_MIN, BLEN_MAX);

    // NNIMove m2;
    // m2.swap_id = i;
    // m2.node1 = u;
    // m2.node2 = v;

    if(debug){
        printf ("After:\n");
        tree1.print();
        tree2.print();
        // Output to a file to facilitate ploting
        // ofstream out_tree("./test/tree-after-NNI.txt");
        // rtree.write(out_tree);
        // out_tree.close();
    }

    nni_trees.push_back(tree1);
    nni_trees.push_back(tree2);

    // moves.push_back(m1);
    // moves.push_back(m2);

    return nni_trees;
}



// Find all the NNI trees for a given tree
map<int, vector<evo_tree>> get_all_NNI_trees(evo_tree& rtree, vector<NNIMove>& moves, int debug){
    map<int, vector<evo_tree>> nni_trees_all;
    // vector<evo_tree> nni_trees_all;
    vector<edge*> intedges = rtree.get_internal_edges();

    for(int i = 0; i < intedges.size(); i++){
        vector<evo_tree> nni_trees = do_NNI_on_branch(rtree, i, moves, debug);
        nni_trees_all[i] = nni_trees;
        // nni_trees_all.insert(nni_trees_all.end(), nni_trees.begin(), nni_trees.end());
    }

    return nni_trees_all;
}
