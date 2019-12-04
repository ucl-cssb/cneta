//evo_tree.hpp

#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <stack>
#include <numeric>      // std::iota
#include <string>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;

const double BLEN_MIN = 1e-10;
const double BLEN_MAX = 100;
const double RATE_MIN = 1e-10;
const double RATE_MAX = 1;
const double RATE_MIN_LOG = -10;
const double RATE_MAX_LOG = 0;


typedef vector<double> DoubleVector;

/**
 * direction of a Neighbor from the root, for rooted tree only
 */
enum RootDirection {UNDEFINED_DIRECTION, TOWARD_ROOT, AWAYFROM_ROOT};

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {
  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

class Node;

/**
    Neighbor list of a node in the tree, from node.h (IQTREE)
 */
class Neighbor {
public:

    /**
        the other end of the branch
     */
    Node* node;

    /**
        branch length
     */
    double length;

    /**
        branch ID
     */
    int id;

    /**
     * direction of the Neighbor in a rooted tree
     */
    RootDirection direction;

    /** size of subtree below this neighbor in terms of number of taxa */
    int size;

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
     */
    Neighbor(Node *anode, double alength) {
        node = anode;
        length = alength;
        id = -1;
        // split = NULL;
        direction = UNDEFINED_DIRECTION;
        size = 0;
    }

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
        @param id branch ID
     */
    Neighbor(Node *anode, double alength, int aid) {
        node = anode;
        length = alength;
        id = aid;
        // split = NULL;
        direction = UNDEFINED_DIRECTION;
        size = 0;
    }

    /**
        construct class with another Neighbor
        @param nei another Neighbor
     */
    Neighbor(Neighbor *nei) {
        node = nei->node;
        length = nei->length;
        id = nei->id;
    }

    /**
     * true if this Neighbor is directed towards the root
     */
    bool isTowardsRoot() {
        assert(direction != UNDEFINED_DIRECTION);
        return (direction == TOWARD_ROOT);
    }

    int getSize() {
        return size;
    }

    void getLength(DoubleVector &vec) {
        vec.resize(1);
        vec[0] = length;
    }

    void getLength(DoubleVector &vec, int start_pos) {
        assert(start_pos < vec.size());
        vec[start_pos] = length;
    }

    void setLength(DoubleVector &vec) {
        assert(vec.size() == 1);
        length = vec[0];
    }

    void setLength(DoubleVector &vec, int start_pos) {
        assert(start_pos < vec.size());
        length = vec[start_pos];
    }

    void setLength(Neighbor *nei) {
        length = nei->length;
    }
};

/**
    Neighbor vector
 */
typedef vector<Neighbor*> NeighborVec;

/*
    some macros to transverse neighbors of a node
 */
#define FOR_NEIGHBOR(mynode, mydad, it) \
	for (it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))

#define FOR_NEIGHBOR_IT(mynode, mydad, it) \
	for (NeighborVec::iterator it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))

#define FOR_NEIGHBOR_DECLARE(mynode, mydad, it) \
	NeighborVec::iterator it; \
	for (it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))



struct NNIMove {
    // Two nodes representing the central branch
    Node *node1, *node2;

    // Roots of the two subtree that are swapped
    NeighborVec::iterator node1Nei_it, node2Nei_it;

    // log-likelihood of the tree after applying the NNI
    double newloglh;

    int swap_id;

    // new branch lengths of 5 branches corresponding to the NNI
    DoubleVector newLen[5];
    // pattern likelihoods
    // double *ptnlh;

    bool operator<(const NNIMove & rhs) const {
        return newloglh > rhs.newloglh;
    }
};


class Node {
public:
  int id;
  int isRoot;
  int isLeaf;
  int parent;

  int e_in;
  vector<int> e_ot;
  vector<int> daughters;
  NeighborVec neighbors;

  double height;

  Node(const int& _id){
    id = _id;
    isRoot = 0;
    isLeaf = 0;
    parent = -1;
    e_in = -1;
  }

  Node(const int& _id, const int& _isRoot, const int& _isLeaf){
    id = _id;
    isRoot = _isRoot;
    isLeaf = _isLeaf;
    parent = -1;
    e_in = -1;
  }

  Node(const Node& _n2){
    id     = _n2.id;
    isRoot = _n2.isRoot;
    isLeaf = _n2.isLeaf;
    parent = _n2.parent;
    e_in   = _n2.e_in;

    e_ot.clear();
    e_ot.insert(e_ot.end(), _n2.e_ot.begin(), _n2.e_ot.end() );
    daughters.clear();
    daughters.insert(daughters.end(), _n2.daughters.begin(), _n2.daughters.end() );
    neighbors.clear();
    neighbors.insert(neighbors.end(), _n2.neighbors.begin(), _n2.neighbors.end() );
  }

  bool is_leaf() {
      return neighbors.size() <= 1;
  }

  int degree() {
      return neighbors.size();
  }


  Neighbor *findNeighbor(Node *node) {
  	int size = neighbors.size();
      for (int i = 0; i < size; i++)
          if (neighbors[i]->node == node) return neighbors[i];
      /*
      for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it ++)
              if ((*it)->node == node)
                      return (*it);*/
      cout << "ERROR : Could not find neighbors of node " << node->id + 1 << endl;
      assert(0);
      return NULL;
  }

  bool isNeighbor(Node* node) {
      int size = neighbors.size();
      for (int i = 0; i < size; i++)
          if (neighbors[i]->node == node) return true;
      return false;
  }

  /**
      @param node the target node
      @return the iterator to the neighbor that has the node. If not found, return neighbors.end()
   */
  NeighborVec::iterator findNeighborIt(Node *node) {
      for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++)
          if ((*it)->node == node)
              return it;
      assert(0);
      return neighbors.end();
  }

  void addNeighbor(Node *node, double length, int id) {
      neighbors.push_back(new Neighbor(node, length, id));
  }

  void addNeighbor(Node *node, DoubleVector &length, int id) {
  //	assert(!length.empty());
      if (length.empty())
          addNeighbor(node, -1.0, id);
      else
          addNeighbor(node, length[0], id);
  }

  void updateNeighbor(NeighborVec::iterator nei_it, Neighbor *newnei) {
      assert(nei_it != neighbors.end());
      *nei_it = newnei;
  }

  void updateNeighbor(NeighborVec::iterator nei_it, Neighbor *newnei, double newlen) {
      assert(nei_it != neighbors.end());
      *nei_it = newnei;
      newnei->length = newlen;
  }

  void updateNeighbor(Node *node, Neighbor *newnei) {
      NeighborVec::iterator nei_it = findNeighborIt(node);
      assert(nei_it != neighbors.end());
      *nei_it = newnei;
  }

  void updateNeighbor(Node *node, Neighbor *newnei, double newlen) {
      NeighborVec::iterator nei_it = findNeighborIt(node);
      assert(nei_it != neighbors.end());
      *nei_it = newnei;
      newnei->length = newlen;
  }

  void updateNeighbor(Node* node, Node *newnode, double newlen) {
      for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++)
          if ((*it)->node == node) {
              (*it)->node = newnode;
              (*it)->length = newlen;
              break;
          }
  }

  double updateNeighbor(Node* node, Node *newnode) {
      for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++)
          if ((*it)->node == node) {
              (*it)->node = newnode;
              return (*it)->length;
          }
      return -1;
  }
};


typedef vector<Node*> NodeVector;
typedef pair<Node*, Node*> Branch;
typedef map<int, Branch> Branches;


class edge {
public:
  int id;
  int start;
  int end;
  double length;
  int parent;
  int nmuts;

  edge(const int& _id, const int& _start, const int& _end, const double& _length){
    id = _id;
    start = _start;
    end = _end;
    length = _length;
    nmuts = 0;
  }

  edge( const edge& _e2 ){
    id     = _e2.id;
    start  = _e2.start;
    end    = _e2.end;
    length = _e2.length;
    parent = _e2.parent;
    nmuts  = _e2.nmuts;
  }
};



class evo_tree {
public:
  int nleaf;
  int nnode;
  int nedge;
  int ntotn;
  int nintedge;
  vector<edge>   edges;
  vector<double> lengths;
  // vector<int> muts;
  vector<Node>   nodes;
  vector<double> node_times;
  // pointers to internal edges
  vector<edge*> intedges;
  int root_node_id;

  // vector<vector<int>> chars;
  vector<double> tobs;

  double score = 0;
  double mu = 0;
  double dup_rate = 0;
  double del_rate = 0;
  double chr_gain_rate = 0;
  double chr_loss_rate = 0;
  double wgd_rate = 0;
  double tree_height;
  double total_time;
  vector<double> top_tinvls;    // Top Ns+1 intervals for optimization with constaints on tip times
  vector<int> top_tnodes;

  // The chosen branch
  int branch_start, branch_end;
  Neighbor* current_it;
  Neighbor* current_it_back;

  evo_tree(){}
  evo_tree(const int& _nleaf, const vector<int>& _edges, const vector<double>& _lengths, int gen_node = 1);
  evo_tree(const int& _nleaf, const vector<edge>& _edges, int gen_node = 1);
  // reparameterized tree: expects terminal edges to all have length 0
  evo_tree(const int& _nleaf, const vector<edge>& _edges, const double& total_time, const vector<double>& _tobs);
  evo_tree(const evo_tree& _t2);

  void   get_nodes_preorder(Node* root, vector<Node*>& nodes_preorder);
  string make_newick(int precision);
  string make_newick_nmut(int precision, vector<int> nmuts);
  void   write_nexus(string newick, ofstream& fout);
  void   scale_time(double ratio);
  void   scale_time_internal(double ratio);
  void   generate_nodes();
  void   calculate_node_times();
  void   generate_int_edges();
  vector<int> get_internal_lengths();

  double get_height2sample(const int& node_id);
  void get_ntime_interval(int k, int& sample1);
  vector<int> get_ntime_ratio();
  pair<int, double> is_edge(int _start, int _end);
  pair<int, double> is_interval(int i, const map<pair<int, int>, double>& slens);
  double find_interval_len(int& _start, int _end, map<pair<int, int>, double>& slens);
  vector<double> get_edges_from_interval(const vector<double>& intervals, vector<int>& tnodes);
  vector<double> get_edges_from_ratio(const vector<double>& ratios);
  // The time from beginning to the time of first sample
  double get_total_time(){ return *max_element(node_times.begin(), node_times.end()) - *max_element(tobs.begin(), tobs.end()); }
  double get_tree_height(){ return *max_element(node_times.begin(), node_times.end()); }


  vector<int> get_ancestral_nodes(const int& node_id) const {
    vector<int> ret;
    //cout << "\tnode_id " << node_id+1 << " :";

    int id = node_id;
    while(nodes[id].isRoot == 0){
      id = nodes[id].parent;
      //cout << "\t" << id+1;
      ret.push_back(id);
    }
    //cout << endl;

    return ret;
  }

  vector<int> get_ancestral_edges(const int& node_id) const {
    vector<int> ret;
    //cout << "\tnode_id " << node_id+1 << " :";

    int eid = nodes[node_id].e_in;
    ret.push_back(eid);
    //cout << "\t" << eid+1 << " (" << edges[eid].start+1 << " -> " << edges[eid].end+1 << ") ";
    while(nodes[ edges[eid].start ].isRoot == 0){
      eid = nodes[ edges[eid].start].e_in;
      //cout << "\t" << eid+1 << " (" << edges[eid].start+1 << " -> " << edges[eid].end+1 << ") ";
      ret.push_back(eid);
    }
    //cout << endl;

    return ret;
  }

  // Find the tip nodes below one node
  vector<int> get_tips_below(int node_id){
      assert(node_id > nleaf-1);
      vector<int> tips;
      for(int i=0; i< nleaf - 1; i++){
          vector<int> anodes = get_ancestral_nodes(i);
          if(find(anodes.begin(), anodes.end(), node_id) != anodes.end()){
              tips.push_back(i);
          }
      }
      return tips;
  }

  double get_descendants_min_time(int node_id){
      if(node_id < nleaf)
        return node_times[node_id];

      vector<int> descendants = get_tips_below(node_id);
      assert(descendants.size() >= 1);
      double min_time = node_times[descendants[0]];
      for(int j = 1; j < descendants.size(); j++){
          double tj = node_times[descendants[j]];
          if(tj < min_time) min_time = tj;
      }
      return min_time;
  }


  // Find the number of mutations on each branch
  vector<int> get_nmuts(vector<double> mu){
      vector<int> nmuts;
      for(int i=0; i< lengths.size(); i++){
          int mut = mu[i] * lengths[i];
          nmuts.push_back(mut);
      }
      return nmuts;
  }

  void print(){
    cout << "EDGES:" << endl;
    cout << "\tid\tstart\tend\tlength" << endl;
    for(int i=0; i<nedge; ++i){
      cout << "\t" << edges[i].id+1 << "\t" << edges[i].start+1 << "\t" << edges[i].end+1 << "\t" << edges[i].length << endl;
    }
    cout << "NODES:" << endl;
    cout << "\tid\tparent\td1\td2\tisLeaf\tisRoot\te_in\te_o1\te_o2\ttime" << endl;
    for(int i=0; i<nodes.size(); ++i){
      cout << "\t" << nodes[i].id+1 << "\t" <<  nodes[i].parent+1;
      if(nodes[i].daughters.size() == 2){
        cout << "\t" << nodes[i].daughters[0]+1 << "\t" << nodes[i].daughters[1]+1;
      }
      else{
        cout << "\t0\t0";
      }

      cout << "\t" << nodes[i].isLeaf  << "\t" << nodes[i].isRoot << "\t" << nodes[i].e_in+1;
      if(nodes[i].e_ot.size() == 2){
        cout<< "\t" << nodes[i].e_ot[0]+1 << "\t" << nodes[i].e_ot[1]+1;
      }
      else{
        cout << "\t0\t0";
      }
      cout << "\t" << node_times[i];
      cout << endl;
    }

  }

  void write(ofstream& of){
    of << "start\tend\tlength" << endl;
    for(int i=0; i<nedge; ++i){
      of << edges[i].start+1 << "\t" << edges[i].end+1 << "\t" << edges[i].length << endl;
    }
  }

  void write_with_mut(ofstream& of, vector<int> nmuts){
    of << "start\tend\tlength\tnmut" << endl;
    for(int i=0; i<nedge; ++i){
      of << edges[i].start+1 << "\t" << edges[i].end+1 << "\t" << edges[i].length << "\t" << nmuts[i] << endl;
    }
  }

  void print_ancestral_edges(const int& node_id) const{
    vector<int> aedges = get_ancestral_edges(node_id);
    reverse(aedges.begin(),aedges.end());
    for(int j=0; j<aedges.size(); ++j) cout << "\t" << edges[aedges[j]].id+1
					    << " (" << edges[aedges[j]].start+1 << " -> " << edges[aedges[j]].end+1 << ") ";
    cout << endl;
  }


  void print_ancestral_edges(const int& node_id, ostream& stream) const{
    stream << "ANCESTRY    (" << node_id+1 << ")";
    vector<int> aedges = get_ancestral_edges(node_id);
    reverse(aedges.begin(),aedges.end());
    for(int j=0; j<aedges.size(); ++j) stream << "\t" << edges[aedges[j]].id+1
					    << " (" << edges[aedges[j]].start+1 << " -> " << edges[aedges[j]].end+1 << ") ";
    stream << endl;
  }

  void compute_branch_direction(Node *node = NULL, Node *dad = NULL) {
    if (!node) {
       node = &(nodes[root_node_id]);
    }
    if (dad)
       ((Neighbor*)node->findNeighbor(dad))->direction = TOWARD_ROOT;
    // cout << "Generating branch directions for edge ending at " << node->id + 1 << endl;
    NeighborVec::iterator it;
    FOR_NEIGHBOR_IT(node, dad, it) {
       // do not update if direction was already computed
       assert(((Neighbor*)*it)->direction != TOWARD_ROOT);
       if (((Neighbor*)*it)->direction != UNDEFINED_DIRECTION)
           continue;
       // otherwise undefined.
       ((Neighbor*)*it)->direction = AWAYFROM_ROOT;
       compute_branch_direction((Node*)(*it)->node, node);
    }
  }

  void generate_neighbors(){
      int debug = 0;
      for(int i = 0; i < nodes.size(); i++){
          Node *n = &nodes[i];
          if(debug) cout << "Generating neighbors for node " << n->id + 1 << endl;
          n->neighbors.clear();

          Node *dad = NULL;
          if(n->e_in >= 0){  // root has no incoming edge
              edge e = edges[n->e_in];
              dad = &nodes[e.start];
              n->addNeighbor(dad, e.length, e.id);
              if(debug) cout << "\tadding node " << e.start + 1 << " with edge " << e.id << endl;
          }

          for(int j = 0; j < n->e_ot.size(); j++){
              edge e = edges[n->e_ot[j]];
              n->addNeighbor(&nodes[e.end], e.length, e.id);
              if(debug) cout << "\tadding node " << e.end + 1 << " with edge " << e.id << endl;
          }
      }

      compute_branch_direction();

      if(debug){
          for(int i = 0; i < nodes.size(); i++){
              Node *n = &nodes[i];
              cout << "Neighbor of node " << n->id + 1 << ":";
              for(int i = 0; i < n->neighbors.size(); i++){
                  cout << "\t" << (n->neighbors[i])->node->id + 1 << ", " << (n->neighbors[i])-> length << ", " << (n->neighbors[i])-> id << ", " << (n->neighbors[i])-> direction;
              }
              cout << endl;
          }
      }
  }

  Node* find_farthest_leaf(Node *node = NULL, Node *dad = NULL) {
      if (!node)
          node = &(nodes[root_node_id]);

      if (dad && node->is_leaf()) {
          node->height = 0.0;
          return node;
      }
      Node *res = NULL;
      node->height = 0.0;
      FOR_NEIGHBOR_IT(node, dad, it) {
          Node *leaf = find_farthest_leaf((*it)->node, node);
          if (node->height < (*it)->node->height+1) {
              node->height = (*it)->node->height+1;
              res = leaf;
          }
      }
      return res;
  }

  void get_preorder_branches(NodeVector &nodes, NodeVector &nodes2, Node *node, Node *dad = NULL) {
      if (dad) {
          nodes.push_back(node);
          nodes2.push_back(dad);
      }

      NeighborVec neivec = node->neighbors;
      NeighborVec::iterator i1, i2;
      for (i1 = neivec.begin(); i1 != neivec.end(); i1++)
          for (i2 = i1+1; i2 != neivec.end(); i2++)
              if ((*i1)->node->height > (*i2)->node->height) {
                  Neighbor *nei = *i1;
                  *i1 = *i2;
                  *i2 = nei;
              }
      for (i1 = neivec.begin(); i1 != neivec.end(); i1++)
          if ((*i1)->node != dad)
              get_preorder_branches(nodes, nodes2, (*i1)->node, node);
  //    FOR_NEIGHBOR_IT(node, dad, it)
  //        get_preorder_branches(nodes, nodes2, (*it)->node, node);
  }

};


evo_tree::evo_tree(const evo_tree& _t2) {
  nleaf = _t2.nleaf;
  nnode = _t2.nnode;
  nedge = _t2.nedge;
  ntotn = _t2.ntotn;
  nintedge = _t2.nintedge;

  root_node_id = _t2.root_node_id;
  score = _t2.score;
  tobs = _t2.tobs;

  mu = _t2.mu;
  dup_rate = _t2.dup_rate;
  del_rate = _t2.del_rate;
  chr_gain_rate = _t2.chr_gain_rate;
  chr_loss_rate = _t2.chr_loss_rate;
  wgd_rate = _t2.wgd_rate;

  top_tinvls = _t2.top_tinvls;
  top_tnodes = _t2.top_tnodes;

  edges.clear();
  edges.insert(edges.end(), _t2.edges.begin(), _t2.edges.end() );
  lengths.clear();
  lengths.insert(lengths.end(), _t2.lengths.begin(), _t2.lengths.end() );
  nodes.clear();
  nodes.insert(nodes.end(), _t2.nodes.begin(), _t2.nodes.end() );
  node_times.clear();
  node_times.insert(node_times.end(), _t2.node_times.begin(), _t2.node_times.end() );

  // chars.clear();
  // for(int i=0; i< _t2.chars.size(); ++i) chars.push_back( _t2.chars[i] );

  intedges.clear();
  generate_int_edges();

  generate_neighbors();
}


evo_tree::evo_tree(const int& _nleaf, const vector<int>& _edges, const vector<double>& _lengths, int gen_node){
  nleaf = _nleaf;
  nnode = nleaf - 1; // internal nodes
  nedge = 2*nleaf - 2;
  ntotn = 2*nleaf - 1; // all nodes
  nintedge = nedge - nleaf;

  mu = 0;
  dup_rate = 0;
  del_rate = 0;
  chr_gain_rate = 0;
  chr_loss_rate = 0;
  wgd_rate = 0;

  // create list of edges
  int count = 0;
  lengths.clear();
  for(int i=0; i<nedge; ++i){
    edges.push_back( edge(i,_edges[count], _edges[count+1], _lengths[i]) );
    lengths.push_back(_lengths[i]);
    count = count + 2;
  }

  if(gen_node){
    generate_nodes();
    calculate_node_times();
    generate_int_edges();
    generate_neighbors();
  }
  score = 0;
}

evo_tree::evo_tree(const int& _nleaf, const vector<edge>& _edges, int gen_node){
  nleaf = _nleaf;
  nnode = nleaf - 1; // internal nodes
  nedge = 2*nleaf - 2;
  ntotn = 2*nleaf - 1; // all nodes
  nintedge = nedge - nleaf;

  mu = 0;
  dup_rate = 0;
  del_rate = 0;
  chr_gain_rate = 0;
  chr_loss_rate = 0;
  wgd_rate = 0;

  edges.clear();
  edges.insert(edges.end(), _edges.begin(), _edges.end());
  lengths.clear();
  for(int i=0; i<nedge; ++i){
      lengths.push_back(_edges[i].length);
  }

  if(gen_node){
    generate_nodes();
    calculate_node_times();
    generate_int_edges();
    generate_neighbors();
  }
  score = 0;
}

evo_tree::evo_tree(const int& _nleaf, const vector<edge>& _edges, const double& total_time, const vector<double>& _tobs){
  int debug = 0;
  //cout << "creating constrained tree" << endl;
  nleaf = _nleaf;
  nnode = nleaf - 1; // internal nodes
  nedge = 2*nleaf - 2;
  ntotn = 2*nleaf - 1; // all nodes
  nintedge = nedge - nleaf;

  mu = 0;
  dup_rate = 0;
  del_rate = 0;
  chr_gain_rate = 0;
  chr_loss_rate = 0;
  wgd_rate = 0;

  edges.clear();
  edges.insert(edges.end(), _edges.begin(), _edges.end() );

  nodes.clear();
  generate_nodes();

  tobs = _tobs;

  // Fill external edge lengths by looping over nodes
  for(int i=0; i<nleaf-1; ++i){
    vector<int> es = get_ancestral_edges( nodes[i].id );
    reverse(es.begin(),es.end());

    if(debug){
        cout << "node id / edges: \t" << nodes[i].id+1 << " : ";
        for(int j=0; j<es.size(); ++j){
          cout << "\t" << edges[es[j]].id+1;
        }
        cout << endl;
        cout << "edge " << es.back() << "; tip " <<  nodes[i].id+1 <<  "; total time " << total_time << "; offset " << tobs[ nodes[i].id ] << endl;
    }
    edges[ es.back() ].length = total_time + tobs[ nodes[i].id ];
    for(int j=0; j<es.size()-1; ++j){
      edges[ es.back() ].length -= edges[es[j]].length;
    }
  }

  calculate_node_times();
  generate_int_edges();
  generate_neighbors();

  lengths.clear();
  for(int i=0; i<nedge; ++i){
      lengths.push_back(edges[i].length);
  }
}

// Scale all the branches (easy to implement)
void evo_tree::scale_time(double ratio) {
    // Scale the tree height and times so that it is less than the age of the patient
    // cout << "Current age of patient " << curr_age << endl;
    for(int i=0; i < edges.size(); i++)
    {
        // cout << "Previous length: " << lengths[i] << endl;
        edges[i].length = edges[i].length * ratio;
        lengths[i] = lengths[i] * ratio;
        // cout << "Scaled length: " << lengths[i] << endl;
        // change the time of tip node to reflect the sampling point?
    }
    calculate_node_times();
}


// Scale non-tip nodes only
// Node times starts with 0 (root)
void evo_tree::scale_time_internal(double ratio) {
    // Scale the tree height and times so that it is less than the age of the patient
    // cout << "Current age of patient " << curr_age << endl;
    for(int i=0; i<nedge; ++i){
      if( edges[i].end >= nleaf ){
          // cout << "Previous length: " << lengths[i] << endl;
          edges[i].length = edges[i].length * ratio;
          lengths[i] = lengths[i] * ratio;
          // cout << "Scaled length: " << lengths[i] << endl;
          // change the time of tip node to reflect the sampling point?
      }
    }
    calculate_node_times();
}

void evo_tree::generate_int_edges(){
  //cout << "generate_int_edges" << endl;
  for(int i=0; i<nedge; ++i){
    if( edges[i].end >= nleaf ){
      intedges.push_back( &edges[i] );
      //cout << "internal edge: " << (intedges.back()->id)+1 << endl;
    }
  }
}

// Find the sum of internal lengths along the path to each tip to get more constrained optimisation
vector<int> evo_tree::get_internal_lengths(){
    int debug = 0;
    vector<int> internal_lens;
    // Fill external edge lengths by looping over nodes
    for(int i=0; i<nleaf-1; ++i){
      double len = 0;
      vector<int> es = get_ancestral_edges( nodes[i].id );
      // reverse(es.begin(),es.end());

      // Exclude the terminal edge (index 0)
      for(int j=1; j<es.size(); ++j){
        len += edges[es[j]].length;
      }

      if(debug){
          cout << "node id / edges: \t" << nodes[i].id+1 << " : " << len << endl;
          for(int j=0; j<es.size(); ++j){
            cout << "\t" << edges[es[j]].id+1;
          }
          cout << endl;
          cout << "edge " << es.back() << "; tip " <<  nodes[i].id+1 <<  "; offset " << tobs[ nodes[i].id ] << endl;
      }
      internal_lens.push_back(len);
    }
    return internal_lens;
}


// Compute the tree height until a certain tip node
double evo_tree::get_height2sample(const int& node_id){
      // Find the ancestral edges of the node
      double height = 0;
      vector<int> aedges = get_ancestral_edges(node_id);
      for(int i = 0; i < aedges.size(); i++){
          height += edges[aedges[i]].length;
      }

      return height;
}


// Find the top k time intervals with only one interval ending at a tip node
// If the tree is not of the disired topolgy, sample1 may not be the first sample
void evo_tree::get_ntime_interval(int k, int& sample1){
    int debug = 0;

    assert(k < node_times.size());
    top_tinvls.clear();
    top_tnodes.clear();
    bool has_tip = false;
    // Record the index of the sorted elements
    for (auto i: sort_indexes(node_times)) {
        if(top_tnodes.size() >= k+1){
            break;
        }
        if(i == nleaf-1 || (i < nleaf - 1 && tobs[i] > 0)) continue;   // Skip the normal leaf node
        if(has_tip && i < nleaf - 1) continue;  // only allow one tip node with tobs 0
        if(!has_tip && i < nleaf - 1) {
            has_tip = true;
            sample1 = i;
        }
        // cout << i+1 << "\t" << node_times[i] << endl;
        top_tnodes.push_back(i);
    }

    for(int i = 0; i < k; i++){
        double invl = node_times[top_tnodes[i+1]] - node_times[top_tnodes[i]];
        top_tinvls.push_back(invl);
    }

    if(debug){
        cout << "Top " << top_tinvls.size() << " intervals " << endl;
        for(int i = 0; i < top_tinvls.size(); i++){
            cout << i+1 << "\t" << top_tinvls[i] << endl;
        }
    }
}


// Convert the time constraints among nodes into a set of ratios for bounded estimation
vector<int> evo_tree::get_ntime_ratio(){
    vector<int> ratios;
    ratios.push_back(node_times[root_node_id]);

    for (int i = root_node_id; i < ntotn; i++){
        // double max_ti = get_descendants_min_time(i);
        // Find its children
        Node ni = nodes[i];
        for(int j = 0; j < ni.daughters.size(); j++){
            double max_tj = get_descendants_min_time(j);
            double ratio = (node_times[i] - max_tj) / (node_times[j] - max_tj);
            ratios.push_back(ratio);
        }
    }

    return ratios;
}

// Check whether if an interval is in the tree. If yes, return the ID.
pair<int, double> evo_tree::is_edge(int _start, int _end){
  for(int i = 0; i<edges.size(); i++){
      int start = edges[i].start;
      int end = edges[i].end;
      if (start == _start && end == _end){
          return pair<int, double>(edges[i].id, edges[i].length);
      }
  }
  return pair<int, double>(-1,0);
}

// Check whether if an edge is in the intervals. If yes, return the ID.
pair<int, double> evo_tree::is_interval(int i, const map<pair<int, int>, double>& slens){
    int _start = edges[i].start;
    int _end = edges[i].end;

    for(auto it : slens){
      pair<int, int> seg = it.first;
      int start = seg.first;
      int end = seg.second;

      // Sometimes, when two interval nodes have the same node times, the direction can be either way
      if ((start == _start && end == _end) || (start == _end && end == _start)){
          return pair<int, double>(edges[i].id, it.second);
      }
    }
    return pair<int, double>(-1,0);
}

// Given the start point (and/or end point) of an interval, find the interval starting with it
double evo_tree::find_interval_len(int& _start, int _end, map<pair<int, int>, double>& slens){
    for(auto it : slens){
        pair<int, int> seg = it.first;
        int start = seg.first;
        if(start != _start) continue;

        // multiple intervals starts from a tip node
        if(start < nleaf - 1){
            for (auto it2 : slens){ // Find the exact tip node
                pair<int, int> seg2 = it2.first;
                if(seg2.first == start && (seg2.second == _end || seg2.second > nleaf)){
                    _start = seg2.second;
                    return it2.second;
                }
                // Sometimes the interval end points may be reversed in a tie
                if(seg2.second == start && (seg2.first == _end)){
                    _start = seg2.first;
                    return it2.second;
                }
            }
        }
        else{
            _start = seg.second;
            return it.second;
        }
        // cout << "Changing start to " << _start + 1 << endl;
    }

    return 0;
}


// Find branch lengths from ratios of node times
vector<double> evo_tree::get_edges_from_ratio(const vector<double>& ratios){
    int debug = 1;
    vector<double> blens;
    // for root
    node_times[root_node_id] = ratios[0];
    for (int i = 1; i < ratios.size(); i++){
        double ratio = ratios[i];
        int nid = i + root_node_id;  // corresponding node ID
        // double max_ti = get_descendants_min_time(nid);
        // Find its children
        Node ni = nodes[nid];
        for(int j = 0; j < ni.daughters.size(); j++){
            double max_tj = get_descendants_min_time(j);
            double ntime = ((node_times[nid] - max_tj) / ratio) + max_tj;
            node_times[j] = ntime;
        }
    }

    // Compute branch lengths from time intervals

    if(debug){
        cout << "branch lengths from ratios:";
        for(int i = 0; i < blens.size(); i++){
            cout << "\t" << blens[i];
        }
        cout << endl;
    }

    return  blens;
}


// Find branch lengths from top k time intervals and tobs
vector<double> evo_tree::get_edges_from_interval(const vector<double>& intervals, vector<int>& tnodes){
    int debug = 0;
    vector<double> blens;
    map<int, double> blens_map; // Ensure branches match with a map
    map<pair<int, int>, double> slens;
    // Find the time for the top k + 1 nodes first
    int sample1 = 0;
    bool is_below = true; // Whether or not the first sample is below all interval nodes. If yes, sample1 should be the last element in tnodes.

    // Find the lengths of all intervals
    for(int i=0; i<intervals.size(); i++){
        if(tnodes[i+1] < nleaf - 1){
            sample1 = tnodes[i+1];
        }
        if(tnodes[i] < nleaf - 1){
            is_below = false;
        }
        pair<int, int> invl(tnodes[i], tnodes[i+1]);
        pair<pair<int, int>, double> s(invl, intervals[i]);
        slens.insert(s);
    }

    if(debug)
    {
        cout << "The 1st sample is " << sample1 + 1 << endl;
        if(!is_below){
            cout << "The first sample is above some interval nodes" << endl;
        }
    }

    for(int i=0; i < tobs.size(); i++){
        // cout << "\t" << tobs[i];
        if(i == sample1) { continue; }
        pair<int, int> invl(sample1, i);
        pair<pair<int, int>, double> s(invl, tobs[i]);
        slens.insert(s);
    }
    // cout << endl;

    if(debug){
        for(auto it: slens){
            pair<int, int> seg(it.first);
            cout << "segment " << seg.first + 1 << "\t" << seg.second + 1 << "\t" << it.second << endl;
        }
    }

    // Find the branch length of each edge
    for (int i=0; i< edges.size(); i++){
        // Skip the edge from LUCA to normal genome
        if(edges[i].start==nleaf && edges[i].end==nleaf-1){
            pair<int, double> id_len(i, edges[i].length);
            blens_map.insert(id_len);
            continue;
        }
        // check if the edge exists in the graph
        pair<int, double> id_len = is_interval(i, slens);
        if(id_len.first >= 0){
            blens_map.insert(id_len);
        }
        else{
            int _start = edges[i].start;
            int _end = edges[i].end;
            double _len = 0;    // length of the edge
            // cout << "Computing length for edge " << _start + 1 << ", " << _end + 1 << endl;
            // Find the index of sample1 and _start in tnodes to decide the direction of the branch
            vector<int>::iterator idx_sample1 = find(tnodes.begin(), tnodes.end(), sample1);
            vector<int>::iterator idx_estart = find(tnodes.begin(), tnodes.end(), _start);
            int dist = distance(idx_estart, idx_sample1);
            if( dist > 0 || (_end > nleaf)){  // When the first sample is below all internal nodes, the branch is the concontation of consecutive line segments
                // cout << "Length as sum of segments" << endl;
                int start = _start;
                while(start != _end){
                    // cout << "\tstart node " << start << endl;
                    _len += find_interval_len(start, _end, slens);
                    // cout << "\tadd segment ending at " << start + 1  << "\t" << _len << endl;
                }
            }
            else{   // When the first sample is on top of some internal nodes, the branch is the differences of delta_t - (segments from first sample until the common ancestor of the tip node)
                // cout << "Length as difference of segments" << endl;
                assert(_end < nleaf - 1);
                double len1 = slens[pair<int, int>(sample1, _end)];
                // cout << "Length for interval "  << sample1 + 1 << ", " << _end + 1 << " is " << len1 << endl;
                // Find the length from the first sample until the common ancestor of the tip node
                double len2 = 0;
                int start = sample1;
                while(start != _start){
                    // cout << "\tstart node " << start << endl;
                    len2 += find_interval_len(start, _start, slens);
                    // cout << "\tadd segment ending at " << start + 1  << " to get length \t" << len2 << endl;
                }
                // cout << "Length for interval "  << sample1 + 1 << ", " << _start + 1 << " is " << len2 << endl;
                _len = len1 - len2;
            }
            // cout << "Length for edge "  << _start + 1 << ", " << _end + 1 << " is " << _len << endl;
            pair<int, double> id_len(i, _len);
            blens_map.insert(id_len);
        }
    }

    // cout << "branch lengths from intervals:" << "\n";
    for(auto it : blens_map){
        // cout << it.first << "\t" << it.second << endl;
        blens.push_back(it.second);
    }

    if(debug){
        cout << "branch lengths from intervals:" << "\n";
        for(auto it : blens_map){
            cout << it.first << "\t" << it.second << endl;
        }
    }

    return  blens;
}

void evo_tree::generate_nodes(){
  // create list of nodes
  for(int i=0; i<ntotn; ++i){
    if(i < nleaf){
      nodes.push_back( Node(i,0,1));
    }
    else{
      nodes.push_back( Node(i,0,0));
    }
  }

  // fill in node details
  for(int i=0; i<nedge; ++i){
    nodes[ edges[i].end ].parent = edges[i].start;
    nodes[ edges[i].start ].daughters.push_back(edges[i].end);

    // edges[i].start are the e_ot of the nodes
    // edges[i].end are the e_in of the nodes
    nodes[ edges[i].end ].e_in = edges[i].id;
    nodes[ edges[i].start ].e_ot.push_back(edges[i].id);
  }

  // locate root node
  for(int i=0; i<nodes.size(); ++i){
    if(nodes[i].parent == -1){
      nodes[i].isRoot = 1;
      root_node_id = nodes[i].id;
    }
  }
  //cout << "finished generating nodes " << endl;
  //for(int i=0; i<ntotn; ++i){
  //cout << "\t" << nodes[i].id << endl;
  //}
}

// assuming branch lengths are times then calculate the times of each node with MCRA at t=0
void evo_tree::calculate_node_times(){
    node_times.clear();
    for(int i=0; i<nodes.size(); ++i){
      node_times.push_back(0);

      if( nodes[i].id != root_node_id){
        vector<int> aedges = get_ancestral_edges(i);
        reverse(aedges.begin(),aedges.end());
        double time = 0;
        for(int j=0; j<aedges.size(); ++j){
           time += edges[aedges[j]].length;
        }
        node_times[ nodes[i].id ] = time;
      }
      else{
	    node_times[ nodes[i].id ] = 0;
      }
      // cout << "node times:" << nodes[i].id+1  << "\t" << node_times[i] << endl;
    }
  }


// Find the preorder of nodes in the tree
void evo_tree::get_nodes_preorder(Node* root, vector<Node*>& nodes_preorder){
    nodes_preorder.push_back(root);
    for(int j=0; j<root->daughters.size();j++){
            get_nodes_preorder(&nodes[root->daughters[j]], nodes_preorder);
    }
}


string evo_tree::make_newick(int precision){
    string newick;
    const boost::format tip_node_format(boost::str(boost::format("%%d:%%.%df") % precision));
    const boost::format internal_node_format(boost::str(boost::format(")%%d:%%.%df") % precision));
    const boost::format root_node_format(boost::str(boost::format(")%%d")));
    stack<Node*> node_stack;
    vector<Node*> nodes_preorder;
    Node* root;

    for(int i=0; i<nodes.size(); ++i){
      if(nodes[i].isRoot){
          root = &nodes[i];
          break;
      }
    }
    // cout << "root " << root->id + 1 << endl;
    get_nodes_preorder(root, nodes_preorder);

    // Traverse nodes in preorder
    for (int i = 0; i<nodes_preorder.size(); i++)
    {
        Node* nd = nodes_preorder[i];
        // cout << nd->id + 1 << endl;
        if (nd->daughters.size()>0) // internal nodes
        {
            newick += "(";
            node_stack.push(nd);
        }
        else
        {
            newick += boost::str(boost::format(tip_node_format) % (nd->id + 1) % lengths[nd->e_in]);
            if (nd->id == nodes[nd->parent].daughters[0])   //left child
                newick += ",";
            else
            {
                Node* popped = (node_stack.empty() ? 0 : node_stack.top());
                while (popped && popped->parent>0 && popped->id == nodes[popped->parent].daughters[1]) // right sibling of the previous node
                {
                    node_stack.pop();
                    newick += boost::str(boost::format(internal_node_format) % (popped->id + 1) % lengths[popped->e_in]);
                    popped = node_stack.top();
                }
                if (popped && popped->parent>0 && popped->id == nodes[popped->parent].daughters[0]) // left child, with another sibling
                {
                    node_stack.pop();
                    newick += boost::str(boost::format(internal_node_format) % (popped->id + 1) % lengths[popped->e_in]);
                    newick += ",";
                }
                if (node_stack.empty())
                {
                    newick += ")";
                }
            }
        }
    }
    newick +=  boost::str(boost::format(root_node_format) % (root->id + 1));
    return newick;
}

string evo_tree::make_newick_nmut(int precision, vector<int> nmuts){
    string newick;
    const boost::format tip_node_format(boost::str(boost::format("%%d:%%.%df") % precision));
    const boost::format internal_node_format(boost::str(boost::format(")%%d:%%.%df") % precision));
    const boost::format root_node_format(boost::str(boost::format(")%%d")));
    stack<Node*> node_stack;
    vector<Node*> nodes_preorder;
    Node* root;

    for(int i=0; i<nodes.size(); ++i){
      if(nodes[i].isRoot){
          root = &nodes[i];
          break;
      }
    }
    // cout << "root " << root->id + 1 << endl;
    get_nodes_preorder(root, nodes_preorder);

    // Traverse nodes in preorder
    for (int i = 0; i<nodes_preorder.size(); i++)
    {
        Node* nd = nodes_preorder[i];
        // cout << nd->id + 1 << endl;
        if (nd->daughters.size()>0) // internal nodes
        {
            newick += "(";
            node_stack.push(nd);
        }
        else
        {
            newick += boost::str(boost::format(tip_node_format) % (nd->id + 1) % nmuts[nd->e_in]);
            if (nd->id == nodes[nd->parent].daughters[0])   //left child
                newick += ",";
            else
            {
                Node* popped = (node_stack.empty() ? 0 : node_stack.top());
                while (popped && popped->parent>0 && popped->id == nodes[popped->parent].daughters[1]) // right sibling of the previous node
                {
                    node_stack.pop();
                    newick += boost::str(boost::format(internal_node_format) % (popped->id + 1) % nmuts[popped->e_in]);
                    popped = node_stack.top();
                }
                if (popped && popped->parent>0 && popped->id == nodes[popped->parent].daughters[0]) // left child, with another sibling
                {
                    node_stack.pop();
                    newick += boost::str(boost::format(internal_node_format) % (popped->id + 1) % nmuts[popped->e_in]);
                    newick += ",";
                }
                if (node_stack.empty())
                {
                    newick += ")";
                }
            }
        }
    }
    newick +=  boost::str(boost::format(root_node_format) % (root->id + 1));
    return newick;
}

// Print the tree in nexus format to be used in other tools for downstream analysis
void evo_tree::write_nexus(string newick, ofstream& fout){
    fout << "#nexus" << endl;
    fout << "begin trees;" << endl;
    // string newick = make_newick(precision);
    fout << "tree 1 = " << newick << ";" << endl;
    fout << "end;" << endl;
}


// The string will be different for different labeled histories
string create_tree_string( evo_tree tree ){
  stringstream sstm;
  for(int i=0; i<tree.ntotn; ++i){
    sstm << tree.nodes[i].id+1;
    if(tree.nodes[i].daughters.size() == 2){
      sstm << ";" << tree.nodes[i].daughters[0]+1 << ";" << tree.nodes[i].daughters[1]+1;
    }
    sstm << ":";
  }
  return sstm.str();
}


// The string will be unique for different topolies
// The key is to name the internal node based on tips
string create_tree_string_uniq(evo_tree tree){
  int debug = 0;
  stringstream sstm;
  // Use a dictory to store the names of internal nodes;
  map<int, string> names;
  for(int i=0; i<tree.ntotn; ++i){
    int nid = tree.nodes[i].id + 1;
    if(nid <= tree.nleaf){
        sstm << nid;
    }
    else if(nid == tree.nleaf + 1){  // root
        sstm << nid << ";" << tree.nodes[i].daughters[0]+1 << ";" << tree.nodes[i].daughters[1]+1;
    }
    else{    // internal nodes
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
            }else if(d2 > tree.nleaf && d1 <= tree.nleaf){  // d2 > Ns
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
      tree.print();
      for(auto it: names){
          cout << it.first << "\t" << it.second << endl;
      }
  }
  return sstm.str();
}

string order_tree_string( string tree ){
  stringstream sstm;
  vector<string> split1;
  boost::split(split1, tree, [](char c){return c == ':';});

  for(int i=0; i<split1.size()-1; ++ i){     // split creates an empty string at the end
    //sstm << split1[i];
    //cout << "\t" << split1[i] << endl;
    vector<string> split2;
    boost::split(split2, split1[i], [](char c){return c == ';';});

    if( split2.size() == 1){
      sstm << split1[i];
    }
    else{
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

// Sort the parent nodes at first
string order_tree_string_uniq( string tree ){
  stringstream sstm;
  vector<string> split1;
  boost::split(split1, tree, [](char c){return c == ':';});
  // Sort split1
  sort(split1.begin(), split1.end());
  for(int i=0; i<split1.size()-1; ++ i){     // split creates an empty string at the end
    //sstm << split1[i];
    //cout << "\t" << split1[i] << endl;
    vector<string> split2;
    boost::split(split2, split1[i], [](char c){return c == ';';});

    if( split2.size() == 1){
      sstm << split1[i];
    }
    else{
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


void test_evo_tree(const evo_tree& tree){
  // generate the list of ancestral nodes belonging to leaf nodes
  cout << "leaf ancestral nodes:" << endl;
  for(int i=0; i<tree.nleaf-1; ++i){
    vector<int> anodes = tree.get_ancestral_nodes( tree.nodes[i].id );
    cout << "\tnode " << tree.nodes[i].id+1;
    for(int j=0; j<anodes.size(); ++j) cout << "\t" << tree.nodes[anodes[j]].id+1;
    cout << endl;
  }

  // generate the list of ancestral edges belonging to leaf nodes
  cout << "leaf ancestral edges:" << endl;
  for(int i=0; i<tree.nleaf-1; ++i){
    vector<int> aedges = tree.get_ancestral_edges( tree.nodes[i].id );
    reverse(aedges.begin(),aedges.end());
    cout << "\tnode " << tree.nodes[i].id+1;
    for(int j=0; j<aedges.size(); ++j) cout << "\t" << tree.edges[aedges[j]].id+1
					    << " (" << tree.edges[aedges[j]].start+1 << " -> " << tree.edges[aedges[j]].end+1 << ") ";
    cout << endl;
  }
}
