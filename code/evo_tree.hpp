#ifndef EVO_TREE_HPP    // To make sure you don't declare the function more than once by including the header multiple times.
#define EVO_TREE_HPP

#include "common.hpp"


const double BLEN_MIN = 1e-3;   // The shortest branch length allowed (in year), about 9h (<1 cell cycle)
const double BLEN_MAX = 100;
const double RATE_MIN = 1e-5;
const double RATE_MAX = 1;
const double RATE_MIN_LOG = -10;
const double RATE_MAX_LOG = 0;
const double SMALL_DIFF_BRANCH = 1.0e-10;   // used to compare branch length
const double SMALL_AGE = 1.0e-5;  // used when converting branch lengths to ratios
const double MIN_RATIO = 1e-2;    // The minimum age ratio allowed
const double MAX_RATIO = 0.99;
const int NRATE = 6;  // type of mutation rates considered

typedef vector<double> DoubleVector;

/**
 * direction of a Neighbor from the root, for rooted tree only
 */
enum RootDirection {UNDEFINED_DIRECTION, TOWARD_ROOT, AWAYFROM_ROOT};


/*
For convenience, the IDs of nodes follow a specified order
leaf nodes: from 1 to Ns
normal node: Ns + 1
root node: Ns + 2
internal node: Ns + 3 to 2 * Ns + 1, small to larger from bottom to up
*/


class Node;


/**
    Neighbor list of a node in the tree, from node.h (IQTREE), used for NNI
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
    Neighbor(Node* anode, double alength);

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
        @param id branch ID
     */
    Neighbor(Node* anode, double alength, int aid);

    /**
        construct class with another Neighbor
        @param nei another Neighbor
     */
    Neighbor(Neighbor* nei);

    // ~Neighbor();

    /**
     * true if this Neighbor is directed towards the root
     */
    bool isTowardsRoot();

    int getSize();

    void getLength(DoubleVector &vec);

    void getLength(DoubleVector &vec, int start_pos);

    void setLength(DoubleVector &vec);

    void setLength(const double& blen);

    void setLength(DoubleVector &vec, int start_pos);

    void setLength(Neighbor* nei);

};


/**
    Neighbor vector
 */
// typedef vector<shared_ptr<Neighbor>> NeighborVec;
typedef vector<Neighbor*> NeighborVec;

/*
    some macros to transverse neighbors of a node
 */
// find neighbors of node which is not dad
#define FOR_NEIGHBOR(mynode, mydad, it) \
	for (it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))

#define FOR_NEIGHBOR_IT(mynode, mydad, it) \
	for (NeighborVec::iterator it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))

// #define FOR_NEIGHBOR_DECLARE(mynode, mydad, it) \
// 	NeighborVec::iterator it; \
// 	for (it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); it++) \
// 		if ((*it)->node != (mydad))


struct NNIMove{
    // Two nodes representing the central branch
    Node* node1, *node2;

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

  int parent;   // node id of its parent
  int e_in;
  vector<int> e_ot;
  vector<int> daughters;  // node ids of its children

  int height;  // used in preorder traverse
  double time;
  double age;

  NeighborVec neighbors;   // vector of pointers, used in NNI


  Node(const int& _id);
  Node(const int& _id, const int& _isRoot, const int& _isLeaf);
  Node(const Node& _n2);

  void deleteNeighbors();
  ~Node();

  bool is_leaf();

  int degree();

  Neighbor* findNeighbor(Node* node);

  bool isNeighbor(Node* node);

  /**
      @param node the target node
      @return the iterator to the neighbor that has the node. If not found, return neighbors.end()
   */
  NeighborVec::iterator findNeighborIt(Node* node);

  void addNeighbor(Node* node, double length, int id);
  void addNeighbor(Node* node, DoubleVector &length, int id);

  void updateNeighbor(NeighborVec::iterator nei_it, Neighbor* newnei);
  void updateNeighbor(NeighborVec::iterator nei_it, Neighbor* newnei, double newlen);
  void updateNeighbor(Node* node, Neighbor* newnei);
  void updateNeighbor(Node* node, Neighbor* newnei, double newlen);
  void updateNeighbor(Node* node, Node* newnode, double newlen);
  double updateNeighbor(Node* node, Node* newnode);
};


typedef vector<Node*> NodeVector;
typedef pair<Node*, Node*> Branch;
typedef map<int, Branch> Branches;


class edge {
public:
  int id;
  int start;
  int end;
  double length;   // time
  int parent;
  int nmuts;

  edge(const int& _id, const int& _start, const int& _end, const double& _length);
  edge(const edge& _e2);
};



class evo_tree {
private:
  // used in converting tree to newick format
  void get_nodes_preorder(Node* root, vector<Node*>& nodes_preorder);

  // used in generating node neighbors
  void compute_branch_direction(Node* node = NULL, Node* dad = NULL);

  // used in ratio computing
  int get_node_height(int node_id); // get node depth (start from root)
  vector<int> get_tips_below(int node_id);  // Find the tip nodes below one node
  int get_node_dist(int ni, int nj);
  double get_tips_max_age(int node_id, int incl_depth = 1);  // Find the maximum age of tips below a node
  void set_edge_length(int start, int end, double blen);

public:
  int nleaf;

  int root_node_id;   // initialized when generating all nodes
  double score;     // likelihood of the tree

  // int nnode = nleaf - 1; // internal nodes
  // int nedge = 2 * nleaf - 2;
  // int ntotn = 2 * nleaf - 1; // all nodes
  // int nintedge = nleaf - 2;

  vector<edge> edges;
  vector<Node> nodes;

  // vector<double> node_times;    // time start from 0 at root
  // vector<double> node_ages;    // count ages from present (recent sample has age 0), used to transform time constraints into ratios
  // vector<edge*> intedges;     // pointers to internal edges
  // double* ratios;  // branch ratios used for constrained optimization

  double mu;   // overall mutation rate
  double dup_rate;
  double del_rate;
  double chr_gain_rate;
  double chr_loss_rate;
  double wgd_rate;

  // The chosen branch for NNI (used in Brent optimization)
  int current_eid;
  Neighbor* current_it;
  Neighbor* current_it_back;
  Node* node1, *node2;  // current branch involved in NNI

  evo_tree();
  evo_tree(const int& _nleaf, const vector<int>& _edges, const vector<double>& _lengths, int gen_node = 1);
  evo_tree(const int& _nleaf, const vector<edge>& _edges, int gen_node = 1);
  evo_tree(const int& _nleaf, const vector<edge>& _edges, const vector<double>& tobs, const double& total_time); // reparameterized tree: expects terminal edges to all have length 0, not used
  evo_tree(const evo_tree& _t2);
  // ~evo_tree();
  evo_tree& operator=(const evo_tree& _t2);


  // used in initialization
  void generate_nodes();
  void generate_neighbors();  // neighbor for one internal node: 1 incoming edges and 2 outgoing edges
  void update_neighbor_lengths(Node* node = NULL, Node* dad = NULL);   // Update neighbor lengths based on updated edge lengths
  void delete_neighbors();

  void calculate_node_times();
  void calculate_age_from_time(bool keep_tip = false);

  vector<edge*> get_internal_edges();
  vector<double> get_node_times();
  vector<double> get_node_ages();

  void scale_time(double ratio);
  void scale_time_internal(double ratio);   // used in sveta

  // update node ages of a node and its ancestor (after the relevant edge is updated)
  void update_node_age(int node_id, double delta);
  // update node times of a node and its descendants (after the relevant edge is updated)
  void update_node_time(int node_id, double delta);

  // used in optimization branch by branch
  void get_preorder_branches(NodeVector &nodes, NodeVector &nodes2, Node* node, Node* dad = NULL);
  // Get postorder of internal nodes for likelihood computation
  void get_inodes_postorder(Node* node, vector<int>& inodes);
  int get_edge_id(int start, int end);   // used in NNI and optimization
  edge* get_edge(int start, int end);

  // Functions related to convertion of branch length ratios and node age
  vector<double> get_ratio_from_age();
  void update_edges_from_ratios(const vector<double>& ratios, const vector<int>& knodes);
  void update_edge_from_ratio(double ratio, int eid);    // not used for now, issues remain on updating dependent edges

  vector<int> get_ancestral_nodes(const int& node_id) const;   // used in testing and internal calling
  vector<int> get_ancestral_edges(const int& node_id) const;   // used in testing and internal calling

  inline bool is_blen_valid() const{
    for(auto e : this->edges){
      if(e.length < 0){
        return false;
      }
    }
    return true;
  }

  void print();
  void print_neighbors() const;
  void print_ancestral_edges(const int& node_id) const;
  void print_ancestral_edges(const int& node_id, ostream& stream) const;
  void print_mutation_rates(int model, int only_seg) const;

  void write(ofstream& of) const;
  void write_with_mut(ofstream& of, const vector<int>& nmuts) const;

  string make_newick(int precision = PRINT_PRECISION);
  string make_newick_nmut(int precision, const vector<int>& nmuts);
  void write_nexus(const string& newick, ofstream& fout) const;

  vector<int> get_nmuts(const vector<double>& mu);    // Find the number of mutations on each branch, used in ML tree building
  Node* find_farthest_leaf(Node* node = NULL, Node* dad = NULL);   // used in optimization

  // functions used in constrained optimization (deprecated due to insufficiency)
  // pair<int, double> is_edge(int _start, int _end);     // Check whether if an interval is in the tree. If yes, return the ID. not used.
  // pair<int, double> is_interval(int i, const map<pair<int, int>, double>& slens);  // Check whether if an edge is in the intervals. If yes, return the ID. not used.
  // vector<int> get_internal_lengths(); // Find the sum of internal lengths along the path to each tip to get more constrained optimisation, not use

};


#endif
