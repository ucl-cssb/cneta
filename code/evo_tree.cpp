#include "evo_tree.hpp"


/**
    Neighbor list of a node in the tree, from node.h (IQTREE), used for NNI
 */

/**
    construct class with a node and length
    @param anode the other end of the branch
    @param alength length of branch
 */
Neighbor::Neighbor(Node* anode, double alength):
node(anode), length(alength), id(-1), direction(UNDEFINED_DIRECTION), size(0){
}

/**
    construct class with a node and length
    @param anode the other end of the branch
    @param alength length of branch
    @param id branch ID
 */
Neighbor::Neighbor(Node* anode, double alength, int aid):
node(anode), length(alength), id(aid), direction(UNDEFINED_DIRECTION), size(0){
}

/**
    construct class with another Neighbor
    @param nei another Neighbor
 */
Neighbor::Neighbor(Neighbor* nei){
  node = (nei->node);
  length = (nei->length);
  id = (nei->id);
  direction = (nei->direction);
  size = (nei->size);
}


/**
 * true if this Neighbor is directed towards the root
 */
bool Neighbor::isTowardsRoot(){
    assert(direction != UNDEFINED_DIRECTION);
    return (direction == TOWARD_ROOT);
}

int Neighbor::getSize(){
    return size;
}

void Neighbor::getLength(DoubleVector& vec){
    vec.resize(1);
    vec[0] = length;
}

void Neighbor::getLength(DoubleVector& vec, int start_pos){
    assert(start_pos < vec.size());
    vec[start_pos] = length;
}

void Neighbor::setLength(const double& blen){
    length = blen;
}

void Neighbor::setLength(DoubleVector& vec){
    assert(vec.size() == 1);
    length = vec[0];
}

void Neighbor::setLength(DoubleVector& vec, int start_pos){
    assert(start_pos < vec.size());
    length = vec[start_pos];
}

void Neighbor::setLength(Neighbor* nei){
    length = nei->length;
}




Node::Node(const int& _id):
id(_id), isRoot(0), isLeaf(0), parent(-1), e_in(-1), height(0), time(0.0), age(0.0){
}


Node::Node(const int& _id, const int& _isRoot, const int& _isLeaf):
id(_id), isRoot(_isRoot), isLeaf(_isLeaf), parent(-1), e_in(-1), height(0), time(0.0), age(0.0){
}


Node::Node(const Node& _n2){
  id = (_n2.id);

  isRoot = (_n2.isRoot);
  isLeaf = (_n2.isLeaf);

  parent = (_n2.parent);
  e_in = (_n2.e_in);
  e_ot.assign(_n2.e_ot.begin(), _n2.e_ot.end());
  daughters.assign(_n2.daughters.begin(), _n2.daughters.end());

  height = (_n2.height);
  time = (_n2.time);
  age = (_n2.age);

  // not copy neighbors as the pointers are for Neighbor objects in _n2
  // regenerate all neighbors in evo_tree
}


bool Node::is_leaf(){
    return neighbors.size() <= 1;
}


int Node::degree(){
    return neighbors.size();
}


Neighbor* Node::findNeighbor(Node* node){
	  int size = neighbors.size();
    for(int i = 0; i < size; i++)
        if(neighbors[i]->node == node) return neighbors[i];
    /*
    for(NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it ++)
            if((*it)->node == node)
                    return (*it);*/
    cout << "ERROR : Could not find neighbor node " << node->id + 1 << " at " << node << " for node " << id + 1 << " on tree " << this << endl;
    assert(0);
    return NULL;
}

bool Node::isNeighbor(Node* node){
    int size = neighbors.size();
    for(int i = 0; i < size; i++)
        if(neighbors[i]->node == node) return true;
    return false;
}

/**
    @param node the target node
    @return the iterator to the neighbor that has the node. If not found, return neighbors.end()
 */
NeighborVec::iterator Node::findNeighborIt(Node* node){
    for(NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        if((*it)->node == node)
            return it;
    assert(0);
    return neighbors.end();
}

void Node::addNeighbor(Node* node, double length, int id){
    neighbors.push_back(new Neighbor(node, length, id));
}

void Node::addNeighbor(Node* node, DoubleVector& length, int id){
//	assert(!length.empty());
    if(length.empty())
        addNeighbor(node, -1.0, id);
    else
        addNeighbor(node, length[0], id);
}

void Node::updateNeighbor(NeighborVec::iterator nei_it, Neighbor* newnei){
    assert(nei_it != neighbors.end());
    // cout << "update neighbor to " << newnei->node->id + 1 << endl;
    *nei_it = newnei;
}

void Node::updateNeighbor(NeighborVec::iterator nei_it, Neighbor* newnei, double newlen){
    assert(nei_it != neighbors.end());
    // cout << "update neighbor to " << newnei->node->id + 1 << " with length " << newlen << endl;
    *nei_it = newnei;
    newnei->length = newlen;
}

void Node::updateNeighbor(Node* node, Neighbor* newnei){
    NeighborVec::iterator nei_it = findNeighborIt(node);
    assert(nei_it != neighbors.end());
    // cout << "update neighbor of " << node->id + 1 << " to " << newnei->node->id + 1 << endl;
    *nei_it = newnei;
}

void Node::updateNeighbor(Node* node, Neighbor* newnei, double newlen){
    NeighborVec::iterator nei_it = findNeighborIt(node);
    assert(nei_it != neighbors.end());
    // cout << "update neighbor of " << node->id + 1 << " to " << newnei->node->id + 1 << " with length " << newlen << endl;
    *nei_it = newnei;
    newnei->length = newlen;
}

void Node::updateNeighbor(Node* node, Node* newnode, double newlen){
    for(NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        if((*it)->node == node){
            (*it)->node = newnode;
            (*it)->length = newlen;
            break;
        }
}

double Node::updateNeighbor(Node* node, Node* newnode){
    for(NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        if((*it)->node == node){
            (*it)->node = newnode;
            return (*it)->length;
        }
    return -1;
}

void Node::deleteNeighbors(){
  // NeighborVec::reverse_iterator it;
  // for(it = neighbors.rbegin(); it != neighbors.rend(); it++){
  //   delete (*it);
  //   (*it) = NULL;
  // }
  for(auto nei: neighbors){
    delete nei;
    nei = NULL;
  }
  neighbors.clear();
}

Node::~Node(){
  deleteNeighbors();
}


edge::edge(const int& _id, const int& _start, const int& _end, const double& _length):
id(_id), start(_start), end(_end), length(_length), parent(-1), nmuts(0){
}

edge::edge(const edge& _e2){
  id = (_e2.id);
  start = (_e2.start);
  end = (_e2.end);
  length = (_e2.length);
  parent = (_e2.parent);
  nmuts = (_e2.nmuts);
}




evo_tree::evo_tree():nleaf(0), current_eid(-1), mu(0.0), dup_rate(0.0), del_rate(0.0), chr_gain_rate(0.0), chr_loss_rate(0.0), wgd_rate(0.0), score(0.0){}


// used in generate_coal_tree
evo_tree::evo_tree(const int& _nleaf, const vector<int>& _edges, const vector<double>& _lengths, int gen_node):
nleaf(_nleaf), current_eid(-1), mu(0.0), dup_rate(0.0), del_rate(0.0), chr_gain_rate(0.0), chr_loss_rate(0.0), wgd_rate(0.0), score(0.0){
  // create list of edges
  int count = 0;
  int nedge = 2 * _nleaf - 2;
  for(int i = 0; i < nedge; ++i){
    edges.push_back(edge(i, _edges[count], _edges[count + 1], _lengths[i]));
    count = count + 2;
  }

  if(gen_node){
    generate_nodes();
    calculate_node_times();
    calculate_age_from_time();
    // generate_neighbors();
  }
}


evo_tree::evo_tree(const int& _nleaf, const vector<edge>& _edges, int gen_node):
nleaf(_nleaf), current_eid(-1), mu(0.0), dup_rate(0.0), del_rate(0.0), chr_gain_rate(0.0), chr_loss_rate(0.0), wgd_rate(0.0), score(0.0){
  edges.assign(_edges.begin(), _edges.end());

  if(gen_node){
    generate_nodes();
    calculate_node_times();
    calculate_age_from_time();
    // generate_neighbors();
  }
}


// update edge lengths by total_time and tobs
// total_time: time to 1st sample
evo_tree::evo_tree(const int& _nleaf, const vector<edge>& _edges, const vector<double>& tobs, const double& total_time):
nleaf(_nleaf), current_eid(-1), mu(0.0), dup_rate(0.0), del_rate(0.0), chr_gain_rate(0.0), chr_loss_rate(0.0), wgd_rate(0.0), score(0.0){
  int debug = 0;

  edges.assign(_edges.begin(), _edges.end());

  generate_nodes();

  // Fill external edge lengths by looping over nodes
  for(int i = 0; i < nleaf - 1; ++i){
    vector<int> es = get_ancestral_edges(nodes[i].id);
    reverse(es.begin(), es.end());

    if(debug){
        cout << "node id / edges: \t" << nodes[i].id + 1 << " : ";
        for(int j = 0; j < es.size(); ++j){
          cout << "\t" << edges[es[j]].id + 1;
        }
        cout << endl;
        cout << "edge " << es.back() << "; tip " <<  nodes[i].id + 1 <<  "; total time " << total_time << "; offset " << tobs[ nodes[i].id ] << endl;
    }

    edges[es.back()].length = total_time + tobs[ nodes[i].id ];
    for(int j = 0; j < es.size()-1; ++j){
      edges[es.back()].length -= edges[es[j]].length;
    }
  }

  calculate_node_times();
  calculate_age_from_time();
}


evo_tree::evo_tree(const evo_tree& _t2){
  nleaf = _t2.nleaf;

  root_node_id = _t2.root_node_id;
  score = _t2.score;

  mu = _t2.mu;
  dup_rate = _t2.dup_rate;
  del_rate = _t2.del_rate;
  chr_gain_rate = _t2.chr_gain_rate;
  chr_loss_rate = _t2.chr_loss_rate;
  wgd_rate = _t2.wgd_rate;

  current_eid = _t2.current_eid;

  edges.assign(_t2.edges.begin(), _t2.edges.end());
  nodes.assign(_t2.nodes.begin(), _t2.nodes.end());
}


// evo_tree::~evo_tree(){
//   delete_neighbors();
// }


evo_tree& evo_tree::operator=(const evo_tree& _t2){
    nleaf = _t2.nleaf;

    root_node_id = _t2.root_node_id;
    score = _t2.score;

    mu = _t2.mu;
    dup_rate = _t2.dup_rate;
    del_rate = _t2.del_rate;
    chr_gain_rate = _t2.chr_gain_rate;
    chr_loss_rate = _t2.chr_loss_rate;
    wgd_rate = _t2.wgd_rate;

    current_eid = _t2.current_eid;

    edges.assign(_t2.edges.begin(), _t2.edges.end());
    nodes.assign(_t2.nodes.begin(), _t2.nodes.end());

    return *this;
}


// create nodes from edges
void evo_tree::generate_nodes(){
  this->nodes.clear();

  // create list of nodes
  int ntotn = 2 * this->nleaf - 1;
  for(int i = 0; i < ntotn; ++i){
    if(i < this->nleaf){
      this->nodes.push_back(Node(i,0,1));
    }else{
      this->nodes.push_back(Node(i,0,0));
    }
  }

  // fill in node details
  int nedge = 2 * this->nleaf - 2;
  for(int i = 0; i < nedge; ++i){
    int ne = this->edges[i].end;
    int ns = this->edges[i].start;

    this->nodes[ne].parent = ns;
    this->nodes[ns].daughters.push_back(ne);

    // edges[i].start are the e_ot of the nodes
    // edges[i].end are the e_in of the nodes
    this->nodes[ne].e_in = edges[i].id;
    this->nodes[ns].e_ot.push_back(edges[i].id);
  }

  // locate root node after all edges are parsed
  for(int i = 0; i < nodes.size(); ++i){
    if(nodes[i].parent == -1){
      this->nodes[i].isRoot = 1;
      this->root_node_id = nodes[i].id;
    }
  }
}


void evo_tree::print_neighbors() const{
  for(int i = 0; i < this->nodes.size(); i++){
    const Node* n = &this->nodes[i];
    cout << "Neighbor of node " << n->id + 1 << ":";
    for(int i = 0; i < n->neighbors.size(); i++){
        cout << "\t" << (n->neighbors[i])->node->id + 1 << ", " << (n->neighbors[i])->length << ", " << (n->neighbors[i])->id << ", " << (n->neighbors[i])->direction;
    }
    cout << endl;

    if(n->e_in >= 0){  // root has no incoming edge
        const edge *e = &this->edges[n->e_in];
        cout << "\tnode " << n->id + 1 << " has neighbor " << e->start + 1 << " with incoming edge " << e->id+1 << endl;
    }

    for(int j = 0; j < n->e_ot.size(); j++){
        const edge *e = &this->edges[n->e_ot[j]];
        cout << "\tnode " << n->id + 1 << " has neighbor " << e->end + 1 << " with outgoing edge " << e->id+1 << endl;
    }
  }
}


// neighbor for one internal node: 1 incoming edges and 2 outgoing edges
void evo_tree::generate_neighbors(){
    int debug = 0;

    for(int i = 0; i < this->nodes.size(); i++){
        Node* n = &this->nodes[i];
        if(debug) cout << "Generating neighbors for node " << n->id + 1 << " at " << n << " with original size " << n->neighbors.size() << endl;
        n->deleteNeighbors();

        Node* dad = NULL;
        if(n->e_in >= 0){  // root has no incoming edge
            edge *e = &this->edges[n->e_in];
            dad = &nodes[e->start];
            n->addNeighbor(dad, e->length, e->id);
            if(debug) cout << "\tadding node " << e->start + 1 << " with incoming edge " << e->id << endl;
        }

        for(int j = 0; j < n->e_ot.size(); j++){
            edge *e = &edges[n->e_ot[j]];
            n->addNeighbor(&nodes[e->end], e->length, e->id);
            if(debug) cout << "\tadding node " << e->end + 1 << " with outgoing edge " << e->id << endl;
        }
    }

    if(debug){
        print_neighbors();
    }

    if(debug) cout << "begin computing direction" << endl;
    compute_branch_direction();

    if(debug){
        print_neighbors();
        cout << "finish generating neighbors" << endl;
    }
}


// Update neighbor lengths based on updated edge lengths
void evo_tree::update_neighbor_lengths(Node* node, Node* dad){
    if(!node){
        node = &(this->nodes[this->root_node_id]);   // root
    }

    FOR_NEIGHBOR_IT(node, dad, it){
        Node* nnei = (Node* )(*it)->node; 
        // edge id of a neighbor may change after NNI
        int idx = get_edge_id(node->id, nnei->id);
        (*it)->id = idx;
        double blen = this->edges[idx].length;

    	  (*it)->setLength(blen);
        (*it)->node->findNeighbor(node)->setLength(blen);

        update_neighbor_lengths((Node* )(*it)->node, node);
    }   
}


void evo_tree::delete_neighbors(){
  for(auto n : this->nodes){
    n.deleteNeighbors();
  }
}


// only use branch lengths (may not obay time constaints at tips)
// assuming branch lengths are times then calculate the times of each node with MCRA at t=0
void evo_tree::calculate_node_times(){
  for(int i = 0; i < this->nodes.size(); ++i){
    if(nodes[i].id != this->root_node_id){
      vector<int> aedges = get_ancestral_edges(i);
      reverse(aedges.begin(), aedges.end());

      double time = 0.0;
      for(int j = 0; j < aedges.size(); ++j){
         time += this->edges[aedges[j]].length;
      }
      nodes[i].time = time;
    }else{
      nodes[i].time = 0.0;
    }
    // cout << "node times:" << nodes[i].id + 1  << "\t" << node_times[i] << endl;
  }
}


// Using only node times
// Node age start from 0 at the lowest node
void evo_tree::calculate_age_from_time(bool keep_tip){
    int debug = 0;
    if(debug) cout << "generate ages from node times" << endl;

    vector<double> node_times = get_node_times();
    double max_time = *max_element(node_times.begin(), node_times.end());

    for(int i = 0; i < this->nodes.size(); i++){
      Node* n = &this->nodes[i];
      if(keep_tip && n->id < this->nleaf - 1){  // fix tip ages
        continue;
      }
      n->age = max_time - node_times[n->id];
    }

    if(debug){
        cout << "ages of nodes: ";
        vector<double> node_ages = get_node_ages();
        for(int i = 0; i < node_ages.size(); i++){
            cout << i + 1 << "\t" << node_ages[i] << endl;
        }
    }
}


vector<int> evo_tree::get_ancestral_nodes(const int& node_id) const{
  vector<int> ret;
  //cout << "\tnode_id " << node_id + 1 << " :";

  int id = node_id;
  while(!this->nodes[id].isRoot){
    id = this->nodes[id].parent;
    //cout << "\t" << id + 1;
    ret.push_back(id);
  }
  //cout << endl;

  return ret;
}


vector<int> evo_tree::get_ancestral_edges(const int& node_id) const{
  vector<int> ret;
  //cout << "\tnode_id " << node_id + 1 << " :";

  int eid = this->nodes[node_id].e_in;
  ret.push_back(eid);
  //cout << "\t" << eid + 1 << " (" << edges[eid].start + 1 << " -> " << edges[eid].end + 1 << ") ";
  while(!this->nodes[this->edges[eid].start].isRoot){
    eid = this->nodes[this->edges[eid].start].e_in;
    //cout << "\t" << eid + 1 << " (" << edges[eid].start + 1 << " -> " << edges[eid].end + 1 << ") ";
    ret.push_back(eid);
  }
  //cout << endl;

  return ret;
}


// Find the tip nodes below one node
vector<int> evo_tree::get_tips_below(int node_id){
    assert(node_id > this->nleaf - 1);
    vector<int> tips;

    for(int i = 0; i < this->nleaf - 1; i++){
        vector<int> anodes = get_ancestral_nodes(i);
        if(find(anodes.begin(), anodes.end(), node_id) != anodes.end()){
            tips.push_back(i);
        }
    }

    return tips;
}


// same edge ID may refer to different edges after topolgy change
int evo_tree::get_edge_id(int start, int end){
    for(int i = 0; i < this->edges.size(); i++){
        edge* e = &this->edges[i];
        if((e->start == start && e->end == end) || (e->start == end && e->end == start)){
            return e->id;
        }
    }
    return -1;
}


edge* evo_tree::get_edge(int start, int end){
    for(int i = 0; i < this->edges.size(); i++){
        edge* e = &this->edges[i];
        if((e->start == start && e->end == end) || (e->start == end && e->end == start)){
            return e;
        }
    }
    return NULL;
}

void evo_tree::set_edge_length(int start, int end, double blen){
    for(int i = 0; i < this->edges.size(); i++){
        edge *e = &this->edges[i];
        if(e->start == start && e->end == end){
            e->length = blen;
        }
    }
}


// get node depth (leaf node has depth 0), used in computing ratios
int evo_tree::get_node_height(int node_id){
  if(node_id < this->nleaf){
    return 0;
  }else{
    int max_d = 0;
    for(int i = 0; i < this->nodes[node_id].daughters.size(); i++){
      int nd = this->nodes[node_id].daughters[i];
      int di = get_node_height(nd);
      if(di > max_d){
        max_d = di;
      }
    }
    max_d = max_d + 1;
    return max_d;
  }
}


// assume ni is ancestor of nj, 
// dist: number of branches between ni and nj
int evo_tree::get_node_dist(int ni, int nj){ 
  int dist = 1;
  int njp = nodes[nj].parent;
  while(njp != ni){
    dist = dist + 1;
    njp =  nodes[njp].parent;
  }
  return dist;
}



// Find the maximum age of tips below a node
double evo_tree::get_tips_max_age(int node_id, int incl_depth){
    if(node_id < this->nleaf){
    }

    vector<int> tips = get_tips_below(node_id);
    assert(tips.size() >= 1);

    double max_age = nodes[tips[0]].age;
    if(incl_depth){     
      int dist = get_node_dist(node_id, tips[0]);
      max_age = max_age + dist * BLEN_MIN;
      // cout << "tip " << tips[0] + 1 << ", distance from  " << node_id + 1 << " is " << dist << endl;
    }
    for(int j = 1; j < tips.size(); j++){      
        double tj = nodes[tips[j]].age;      
        if(incl_depth){
          int dist = get_node_dist(node_id, tips[j]);           
          tj = tj + dist * BLEN_MIN;
          // cout << "tip " << tips[j] + 1 << ", distance from  " << node_id + 1 << " is " << dist << endl;
        }
        if(tj > max_age) max_age = tj;
    }

    return max_age;
}


// update node ages of a node and its ancestor (after the relevant edge is updated)
void evo_tree::update_node_age(int node_id, double delta){
  int debug = 0;

  nodes[node_id].age = nodes[node_id].age + delta;

  vector<int> aedges = get_ancestral_edges(node_id);

  for(int j = 0; j < aedges.size(); ++j){
    int nid = edges[aedges[j]].start;
    nodes[nid].age = nodes[nid].age + delta;
  }

  nodes[nleaf - 1].age = nodes[nleaf].age;

  if(debug){
    cout << "update nodes above " << node_id  << " by " << delta << endl;
    for(int j = 0; j < aedges.size(); ++j){
      int nid = edges[aedges[j]].start;
      cout << " node " << nid << " with new age " << nodes[nid].age << endl;
    }
    cout << endl;
  }

}


// update node times of a node and its descendants (after the relevant edge is updated)
void evo_tree::update_node_time(int node_id, double delta){
  int debug = 0;
  nodes[node_id].time = nodes[node_id].time + delta;
  if(debug){
    cout << "update time of node " << node_id << " to " << nodes[node_id].time << endl;
  }

  for(int i = 0; i < nodes[node_id].daughters.size(); i++){
    update_node_time(nodes[node_id].daughters[i], delta);
  }
}


// mu: overall mutation rates for all branches
vector<int> evo_tree::get_nmuts(const vector<double>& mu){
    vector<int> nmuts;
    for(int i = 0; i < this->edges.size(); i++){
        int mut = mu[i] * this->edges[i].length;
        nmuts.push_back(mut);
    }
    return nmuts;
}


void evo_tree::compute_branch_direction(Node* node, Node* dad){
  if(!node){
     node = &(nodes[root_node_id]);
  }
  // cout << " get direction for node " << node->id + 1 << " at " << node << endl;
  if(dad){
      // cout << " dad is " << dad->id + 1 << " at " << dad << endl;
      ((Neighbor*)node->findNeighbor(dad))->direction = TOWARD_ROOT;
  }
  // cout << "Generating branch directions for edge ending at " << node->id + 1 << endl;
  NeighborVec::iterator it;
  FOR_NEIGHBOR_IT(node, dad, it){
     // do not update if direction was already computed
     // cout << " neighbor of node " << (*it)->node->id + 1 << " at " << (*it)->node << endl;
     assert(((Neighbor*)*it)->direction != TOWARD_ROOT);
     if(((Neighbor*)*it)->direction != UNDEFINED_DIRECTION)
         continue;
     // otherwise undefined.
     ((Neighbor*)*it)->direction = AWAYFROM_ROOT;
     compute_branch_direction((Node*)(*it)->node, node);
  }
}


Node* evo_tree::find_farthest_leaf(Node* node, Node* dad){
    if(!node)
        node = &(nodes[root_node_id]);

    if(dad && node->is_leaf()){
        node->height = 0;
        return node;
    }
    Node* res = NULL;
    node->height = 0;
    FOR_NEIGHBOR_IT(node, dad, it){
        Node* leaf = find_farthest_leaf((*it)->node, node);
        if(node->height < (*it)->node->height + 1){
            node->height = (*it)->node->height + 1;
            res = leaf;
        }
    }
    return res;
}

void evo_tree::get_preorder_branches(NodeVector &nodes, NodeVector &nodes2, Node* node, Node* dad){
    if(dad){
        nodes.push_back(node);
        nodes2.push_back(dad);
    }

    NeighborVec neivec = node->neighbors;
    NeighborVec::iterator i1, i2;
    for(i1 = neivec.begin(); i1 != neivec.end(); i1++)
        for(i2 = i1 + 1; i2 != neivec.end(); i2++)
            if((*i1)->node->height > (*i2)->node->height){
                Neighbor* nei = *i1;
                *i1 = *i2;
                *i2 = nei;
            }
    for(i1 = neivec.begin(); i1 != neivec.end(); i1++)
        if((*i1)->node != dad)
            get_preorder_branches(nodes, nodes2, (*i1)->node, node);
}


// Scale all the branches (easy to implement)
void evo_tree::scale_time(double ratio){
    // Scale the tree height and times so that it is less than the age of the patient
    // cout << "Current age of patient " << curr_age << endl;
    for(int i = 0; i < edges.size(); i++)
    {
        // cout << "Previous length: " << lengths[i] << endl;
        edges[i].length = edges[i].length * ratio;
        // lengths[i] = lengths[i] * ratio;
        // cout << "Scaled length: " << lengths[i] << endl;
        // change the time of tip node to reflect the sampling point?
    }
    calculate_node_times();
    calculate_age_from_time();
}


// Scale non-tip nodes only
// Node times starts with 0 (root)
void evo_tree::scale_time_internal(double ratio){
    // Scale the tree height and times so that it is less than the age of the patient
    // cout << "Current age of patient " << curr_age << endl;
    int nedge = 2 * nleaf - 2;
    for(int i = 0; i < nedge; ++i){
      if(edges[i].end >= nleaf){
          // cout << "Previous length: " << lengths[i] << endl;
          edges[i].length = edges[i].length * ratio;
          // lengths[i] = lengths[i] * ratio;
          // cout << "Scaled length: " << lengths[i] << endl;
          // change the time of tip node to reflect the sampling point?
      }
    }
    calculate_node_times();
    calculate_age_from_time();
}


vector<edge*> evo_tree::get_internal_edges(){
  vector<edge*> intedges;
  //cout << "get_internal_edges" << endl;
  int nedge = 2 * nleaf - 2;
  for(int i = 0; i < nedge; ++i){
    if(edges[i].end >= nleaf){
      intedges.push_back(&edges[i]);
      //cout << "internal edge: " << (intedges.back()->id) + 1 << endl;
    }
  }
  return intedges;
}


vector<double> evo_tree::get_node_times(){
  vector<double> node_times(this->nodes.size(), 0.0);

  for(int i = 0; i < this->nodes.size(); i++){
    Node* n = &this->nodes[i];
    node_times[n->id] = n->time;
  }

  return node_times;
}


vector<double> evo_tree::get_node_ages(){
  vector<double> node_ages(this->nodes.size(), 0.0);

  for(int i = 0; i < this->nodes.size(); i++){
    Node* n = &this->nodes[i];
    node_ages[n->id] = n->age;
  }

  return node_ages;
}


// // Check whether if an interval is in the tree. If yes, return the ID.
// pair<int, double> evo_tree::is_edge(int _start, int _end){
//   for(int i = 0; i < this->edges.size(); i++){
//     edge* e = &edges[i];
//     if(e->start == _start && e->end == _end){
//         return pair<int, double>(e->id, e->length);
//     }
//   }
//   return pair<int, double>(-1,0);
// }


// // Check whether if an edge is in the intervals. If yes, return the ID.
// pair<int, double> evo_tree::is_interval(int i, const map<pair<int, int>, double>& slens){
//     int _start = edges[i].start;
//     int _end = edges[i].end;
//
//     for(auto it : slens){
//       pair<int, int> seg = it.first;
//       int start = seg.first;
//       int end = seg.second;
//
//       // Sometimes, when two interval nodes have the same node times, the direction can be either way
//       if((start == _start && end == _end) || (start == _end && end == _start)){
//           return pair<int, double>(edges[i].id, it.second);
//       }
//     }
//     return pair<int, double>(-1,0);
// }

// Given the start point (and/or end point) of an interval, find the interval starting with it
// double evo_tree::find_interval_len(int& _start, int _end, map<pair<int, int>, double>& slens){
//     for(auto it : slens){
//         pair<int, int> seg = it.first;
//         int start = seg.first;
//         if(start != _start) continue;
//
//         // multiple intervals starts from a tip node
//         if(start < nleaf - 1){
//             for(auto it2 : slens){ // Find the exact tip node
//                 pair<int, int> seg2 = it2.first;
//                 if(seg2.first == start && (seg2.second == _end || seg2.second > nleaf)){
//                     _start = seg2.second;
//                     return it2.second;
//                 }
//                 // Sometimes the interval end points may be reversed in a tie
//                 if(seg2.second == start && (seg2.first == _end)){
//                     _start = seg2.first;
//                     return it2.second;
//                 }
//             }
//         }
//         else{
//             _start = seg.second;
//             return it.second;
//         }
//         // cout << "Changing start to " << _start + 1 << endl;
//     }
//
//     return 0;
// }


// Convert the time constraints among nodes into a set of ratios for bounded estimation
// assume node ages are available, which were computed based on current branch lengths
// number of ratios equal to #internal nodes
vector<double> evo_tree::get_ratio_from_age(){
    int debug = 0;

    if(debug){
        cout << "getting ratio from age by edges " << endl;
        this->print();
        cout << "Newick String for current tree is " << make_newick(8) << endl;
    }

    vector<double> ratios(this->nleaf - 1, 0.0);

    int ntotn = 2 * this->nleaf - 1;
    ratios[0] = this->nodes[root_node_id].age;
    for(int ni = root_node_id; ni < ntotn; ni++){
        // Find its children
        Node* nodei = &this->nodes[ni];
        if(debug) cout << "node " << ni + 1 << ", node age " << nodei->age << endl;

        for(int j = 0; j < nodei->daughters.size(); j++){
            int nj = nodei->daughters[j];
            if(nj < root_node_id) continue;

            double max_tj = get_tips_max_age(nj);    // node depths are considered
            double t1 = nodei->age - max_tj - 2 * BLEN_MIN;   // for parent node, age > max_tj + 2 * BLEN_MIN
            if(fabs(t1) < SMALL_AGE)  t1 = SMALL_AGE;     // t1 and t2 may become very small values, or even negative
            double t2 = this->nodes[nj].age - max_tj - BLEN_MIN;   // for child node, age > max_tj + BLEN_MIN
            if(fabs(t2) < SMALL_AGE)  t2 = SMALL_AGE / 10;      // more likely to have both t1 and t2 smaller than SMALL_AGE
            
            double ratio =  t2 / t1;
            // transform values to keep range of ratio in [0.01, 0.99]
            if(ratio < MIN_RATIO){
              ratio = MIN_RATIO;
            }else if(ratio > MAX_RATIO){
              ratio = MAX_RATIO;
            }else{

            }     
           
            if(debug){
                cout << " child node " << nj + 1 << ", node age " << this->nodes[nj].age << ", max age " << max_tj << ", ratio " << ratio  << endl;
            }
            // parent node always has larger node age
            // assert(ratio < 1 && ratio > 0);
            if(ratio >= 1 || ratio <= 0){
                cout << "Invalid ratio " << ratio << " for node " << ni + 1 << "(" << t1 << ") and " << nj + 1 << "(" << t2 << ") " << "with max age of descendants being " << max_tj << endl;
                cout << "all the ratios so far: ";
                for(auto r : ratios){
                  cout << "\t" << r;
                }
                cout << endl;
                cout << "all the node ages so far: ";
                for(auto n : nodes){
                  cout << "\t" << n.age;
                }
                cout << endl;
                print();
                cout << make_newick() << endl;
                exit(EXIT_FAILURE);
            }
            ratios[nj - root_node_id] = ratio;
        }
    }

    return ratios;
}


// Find branch lengths from ratios of node times. Called each time a new tree is returned during tree search
// Assume tip times are given by tobs
// Have to update branch lengths from top to bottom
void evo_tree::update_edges_from_ratios(const vector<double>& ratios, const vector<int>& knodes){
    int debug = 0;
    int nnode = nleaf - 1;

    if(debug){
        cout << "\ngetting edges back from ratios below: " << endl;
        for(int i = 0; i < nnode; i++){
            cout << i + 1 << "\t" << "\t" << ratios[i] << endl;
        }
    }

    // for root and normal node
    // cout << "root time before " << node_ages[root_node_id] << endl;
    nodes[root_node_id].age = ratios[0];
    nodes[root_node_id - 1].age = ratios[0];     // normal node

    // another daughter of root (MRCA of tumor samples), fixed node ID
    int nj = nnode + root_node_id - 1;
    if(debug) cout << "root daughter node is " << nj + 1 << endl;
    double max_tj = get_tips_max_age(nj);
    // If node age for root is the same as max_tj, all other nodes will be equal to max_tj
    double ntime = (nodes[root_node_id].age - max_tj - 2 * BLEN_MIN) * ratios[nnode - 1] + max_tj + BLEN_MIN;
    nodes[nj].age = ntime;
    double blen = nodes[root_node_id].age - nodes[nj].age;
    if(debug){
        cout << "maximum time below root daughter node " << nj + 1 << " is " << max_tj << endl;
        cout << "new age for node " << nj + 1 << " is " << ntime << ", converted from ratio " << ratios[nnode - 1] << endl;
        cout << " parent node " << root_node_id + 1 << " child node " << nj + 1 << ", parent time after " << nodes[root_node_id].age << ", child time after " << nodes[nj].age << ", blen after " << blen << endl;
    }

    // blen < BLEN_MIN, due to precision, it may appear 0.001 smaller than 0.001
    if(blen < BLEN_MIN && fabs(blen - BLEN_MIN) > SMALL_VAL){ 
        cout << "\nBranch length (" << root_node_id + 1 << ", " << nj + 1 << ") " << blen << " smaller than " << BLEN_MIN << endl;
        this->print();
        cout << make_newick() << endl;
        exit(EXIT_FAILURE);  
    }

    set_edge_length(root_node_id, nj, blen);

    // get preorder of internal nodes
    vector<int> inodes(knodes.rbegin(), knodes.rend());

    for(int i = 1; i < inodes.size(); i++){
        int ni = inodes[i];  // corresponding node ID
        // cout << "updating node " << ni + 1 << endl;
        // Find its children
        Node* nodei = &nodes[ni];
        for(int j = 0; j < nodei->daughters.size(); j++){
            int nj = nodei->daughters[j];
            if(debug){
                cout << " child node " << nj + 1 << " with age " << nodes[nj].age << endl;
            }
            if(nj > root_node_id){  // when child is not a leaf, update child age
                double max_tj = get_tips_max_age(nj);
                double ntime = (nodei->age - max_tj - 2 * BLEN_MIN) * ratios[nj - root_node_id] + max_tj + BLEN_MIN;
                nodes[nj].age = ntime;

                if(debug){
                  cout << " maximum time below node " << nj + 1 << " is " << max_tj << endl;
                  cout << " new age for node " << nj + 1 << " is " << ntime << ", converted from ratio " << ratios[nj - root_node_id] << endl;
                }
            }
            double blen = nodei->age - nodes[nj].age;  // The age of parent node should be no less than the age of its child
            if(debug){
                cout << " parent node " << ni + 1 << " child node " << nj + 1 << ", parent age after " << nodei->age << ", child age after " << nodes[nj].age << ", blen after " << blen << endl;
            }

            if(blen < BLEN_MIN && fabs(blen - BLEN_MIN) > SMALL_VAL){ 
              cout << "\nBranch length (" << ni + 1 << ", " << nj + 1 << ") " << blen << " smaller than " << BLEN_MIN << endl;
              this->print();
              cout << make_newick() << endl;
              cout << "current set of ratios: ";
              for(auto r : ratios){
                  cout << "\t" << r;
              }
              cout << endl;              
              exit(EXIT_FAILURE);
            }
            set_edge_length(ni, nj, blen);
        }
    }

    if(debug){
        cout << "node ages from ratios:";
        for(int i = 0; i < nodes.size(); i++){
            cout << "\t" << nodes[i].id << ":" << nodes[i].age;
        }
        cout << endl;
    }

    // update node time and age to be consistent with edge length
    calculate_node_times();
    calculate_age_from_time(true);

    // check node times and node ages are consistent

    // Compute branch lengths from time intervals
    if(debug){
        cout << "branch lengths from ratios:";
        for(int i = 0; i < edges.size(); i++){
            cout << "\t" << edges[i].length;
        }
        cout << endl;

        cout << "node ages after recomputing from branch lengths and node times:";
        for(int i = 0; i < nodes.size(); i++){
            cout << "\t" << nodes[i].id << ":" << nodes[i].age;
        }
        cout << endl;

        cout << "Current tree after updating from new ratios is " << endl;
        print();
        string newick = make_newick(5);
        cout << newick << endl;
    }
}


// Find a single branch length from the ratio of end node times
void evo_tree::update_edge_from_ratio(double ratio, int eid){
    int debug = 0;
    if(eid < 0 || ratio <= 0 || ratio >= 1) return;

    edge *e = &edges[eid];
    double max_tj = get_tips_max_age(e->end);
    double ntime = ((nodes[e->start].age - max_tj - 2 * BLEN_MIN) * ratio) + max_tj + BLEN_MIN;
    nodes[e->end].age = ntime;

    // Compute branch lengths from time intervals
    if(debug){
        cout << "updating edge " << eid + 1 << " from ratio " << ratio << " with starting node at age " << nodes[e->start].age << " max tip time at " << max_tj << endl;
        cout << "node age from ratio:" << "\t" << nodes[e->end].age << endl;
    }

    double blen = nodes[e->start].age - nodes[e->end].age;
    assert(blen >= BLEN_MIN);
    double delta = blen - e->length;
    e->length = blen;

    // node age is bottom up
    update_node_age(e->start, delta);
    // node time is top down
    update_node_time(e->end, delta);

    if(debug){
        cout << "branch length from ratio: " << edges[eid].length << endl;
    }
}


// Find branch lengths from top k time intervals and tobs
// vector<double> evo_tree::get_edges_from_interval(const vector<double>& intervals, vector<int>& tnodes){
//     int debug = 0;
//     vector<double> blens;
//     map<int, double> blens_map; // Ensure branches match with a map
//     map<pair<int, int>, double> slens;
//     // Find the time for the top k + 1 nodes first
//     int sample1 = 0;
//     bool is_below = true; // Whether or not the first sample is below all interval nodes. If yes, sample1 should be the last element in tnodes.
//
//     // Find the lengths of all intervals
//     for(int i = 0; i < intervals.size(); i++){
//         if(tnodes[i + 1] < nleaf - 1){
//             sample1 = tnodes[i + 1];
//         }
//         if(tnodes[i] < nleaf - 1){
//             is_below = false;
//         }
//         pair<int, int> invl(tnodes[i], tnodes[i + 1]);
//         pair<pair<int, int>, double> s(invl, intervals[i]);
//         slens.insert(s);
//     }
//
//     if(debug)
//    {
//         cout << "The 1st sample is " << sample1 + 1 << endl;
//         if(!is_below){
//             cout << "The first sample is above some interval nodes" << endl;
//         }
//     }
//
//     for(int i = 0; i < tobs.size(); i++){
//         // cout << "\t" << tobs[i];
//         if(i == sample1){ continue; }
//         pair<int, int> invl(sample1, i);
//         pair<pair<int, int>, double> s(invl, tobs[i]);
//         slens.insert(s);
//     }
//     // cout << endl;
//
//     if(debug){
//         for(auto it: slens){
//             pair<int, int> seg(it.first);
//             cout << "segment " << seg.first + 1 << "\t" << seg.second + 1 << "\t" << it.second << endl;
//         }
//     }
//
//     // Find the branch length of each edge
//     for(int i = 0; i <  edges.size(); i++){
//         // Skip the edge from LUCA to normal genome
//         if(edges[i].start==nleaf && edges[i].end==nleaf-1){
//             pair<int, double> id_len(i, edges[i].length);
//             blens_map.insert(id_len);
//             continue;
//         }
//         // check if the edge exists in the graph
//         pair<int, double> id_len = is_interval(i, slens);
//         if(id_len.first >= 0){
//             blens_map.insert(id_len);
//         }
//         else{
//             int _start = edges[i].start;
//             int _end = edges[i].end;
//             double _len = 0;    // length of the edge
//             // cout << "Computing length for edge " << _start + 1 << ", " << _end + 1 << endl;
//             // Find the index of sample1 and _start in tnodes to decide the direction of the branch
//             vector<int>::iterator idx_sample1 = find(tnodes.begin(), tnodes.end(), sample1);
//             vector<int>::iterator idx_estart = find(tnodes.begin(), tnodes.end(), _start);
//             int dist = distance(idx_estart, idx_sample1);
//             if( dist > 0 || (_end > nleaf)){  // When the first sample is below all internal nodes, the branch is the concontation of consecutive line segments
//                 // cout << "Length as sum of segments" << endl;
//                 int start = _start;
//                 while(start != _end){
//                     // cout << "\tstart node " << start << endl;
//                     _len += find_interval_len(start, _end, slens);
//                     // cout << "\tadd segment ending at " << start + 1  << "\t" << _len << endl;
//                 }
//             }
//             else{   // When the first sample is on top of some internal nodes, the branch is the differences of delta_t - (segments from first sample until the common ancestor of the tip node)
//                 // cout << "Length as difference of segments" << endl;
//                 assert(_end < nleaf - 1);
//                 double len1 = slens[pair<int, int>(sample1, _end)];
//                 // cout << "Length for interval "  << sample1 + 1 << ", " << _end + 1 << " is " << len1 << endl;
//                 // Find the length from the first sample until the common ancestor of the tip node
//                 double len2 = 0;
//                 int start = sample1;
//                 while(start != _start){
//                     // cout << "\tstart node " << start << endl;
//                     len2 += find_interval_len(start, _start, slens);
//                     // cout << "\tadd segment ending at " << start + 1  << " to get length \t" << len2 << endl;
//                 }
//                 // cout << "Length for interval "  << sample1 + 1 << ", " << _start + 1 << " is " << len2 << endl;
//                 _len = len1 - len2;
//             }
//             // cout << "Length for edge "  << _start + 1 << ", " << _end + 1 << " is " << _len << endl;
//             pair<int, double> id_len(i, _len);
//             blens_map.insert(id_len);
//         }
//     }
//
//     // cout << "branch lengths from intervals:" << "\n";
//     for(auto it : blens_map){
//         // cout << it.first << "\t" << it.second << endl;
//         blens.push_back(it.second);
//     }
//
//     if(debug){
//         cout << "branch lengths from intervals:" << "\n";
//         for(auto it : blens_map){
//             cout << it.first << "\t" << it.second << endl;
//         }
//     }
//
//     return  blens;
// }



// Find the preorder of nodes in the tree
void evo_tree::get_nodes_preorder(Node* root, vector<Node*>& nodes_preorder){
    nodes_preorder.push_back(root);
    for(int j = 0; j < root->daughters.size();j++){
        get_nodes_preorder(&nodes[root->daughters[j]], nodes_preorder);
    }
}



// Get postorder of internal nodes for likelihood computation
void evo_tree::get_inodes_postorder(Node* node, vector<int>& inodes){
  // cout << "visiting node " << node->id << endl;
  for(auto child: node->daughters){
    get_inodes_postorder(&this->nodes[child], inodes);
  }

  if(!node->isLeaf){
    // cout << "non leaf for " << node->id << endl;
    inodes.push_back(node->id);
  }
}



/**************** general printing functions ****************************/

// Ensure branch lengths, node times and node ages are always consistent
void evo_tree::print(){
  // update node times in case branch lengths are changed somewhere
  // calculate_node_times();
  // calculate_age_from_time();

  cout << "EDGES:" << endl;
  cout << "\tid\tstart\tend\tlength" << endl;
  for(int i = 0; i < edges.size(); ++i){
    cout << "\t" << edges[i].id + 1 << "\t" << edges[i].start + 1 << "\t" << edges[i].end + 1 << "\t" << edges[i].length << endl;
  }

  cout << "NODES:" << endl;
  cout << "\tid\tparent\td1\td2\tisLeaf\tisRoot\te_in\te_o1\te_o2\ttime\tage" << endl;
  for(int i = 0; i < nodes.size(); ++i){
    cout << "\t" << nodes[i].id + 1 << "\t" <<  nodes[i].parent + 1;
    if(nodes[i].daughters.size() == 2){
      cout << "\t" << nodes[i].daughters[0] + 1 << "\t" << nodes[i].daughters[1] + 1;
    }
    else{
      cout << "\t0\t0";
    }
    cout << "\t" << nodes[i].isLeaf << "\t" << nodes[i].isRoot << "\t" << nodes[i].e_in + 1;
    if(nodes[i].e_ot.size() == 2){
      cout<< "\t" << nodes[i].e_ot[0] + 1 << "\t" << nodes[i].e_ot[1] + 1;
    }
    else{
      cout << "\t0\t0";
    }
    cout << "\t" << nodes[i].time << "\t" << nodes[i].age;
    cout << endl;
  }
  // cout << "root node is " << this->root_node_id << endl;
}


void evo_tree::write(ofstream& of) const{
  of << "start\tend\tlength\n";
  int nedge = 2 * nleaf - 2;
  for(int i = 0; i < nedge; ++i){
    of << edges[i].start + 1 << "\t" << edges[i].end + 1 << "\t" << edges[i].length << endl;
  }
}


void evo_tree::write_with_mut(ofstream& of, const vector<int>& nmuts) const{
  of << "start\tend\tlength\tnmut" << endl;
  int nedge = 2 * nleaf - 2;
  for(int i = 0; i < nedge; ++i){
    of << edges[i].start + 1 << "\t" << edges[i].end + 1 << "\t" << edges[i].length << "\t" << nmuts[i] << endl;
  }
}


void evo_tree::print_ancestral_edges(const int& node_id) const{
  vector<int> aedges = get_ancestral_edges(node_id);
  reverse(aedges.begin(),aedges.end());
  for(int j = 0; j < aedges.size(); ++j) cout << "\t" << edges[aedges[j]].id + 1
				    << " (" << edges[aedges[j]].start + 1 << " -> " << edges[aedges[j]].end + 1 << ") ";
  cout << endl;
}


void evo_tree::print_ancestral_edges(const int& node_id, ostream& stream) const{
  stream << "ANCESTRY    (" << node_id + 1 << ")";
  vector<int> aedges = get_ancestral_edges(node_id);
  reverse(aedges.begin(),aedges.end());
  for(int j = 0; j < aedges.size(); ++j) stream << "\t" << edges[aedges[j]].id + 1
				    << " (" << edges[aedges[j]].start + 1 << " -> " << edges[aedges[j]].end + 1 << ") ";
  stream << endl;
}


void evo_tree::print_mutation_rates(int model, int only_seg) const{
    if(model == MK){
      cout << "mutation rate (per allele per site per year):  " << mu << endl;
    }else{
      cout << "segment duplication rate (per allele per site per year):  " << dup_rate << endl;
      cout << "segment deletion rate (per allele per site per year):  " << del_rate << endl;
      if(!only_seg){
          cout << "chromosome gain rate (per chromosome per year):  " << chr_gain_rate << endl;
          cout << "chromosome loss rate (per chromosome per year):  " << chr_loss_rate << endl;
          cout << "whole genome doubling rate (per year):  " << wgd_rate << endl;
      }
    }
}

/****************** convert tree to newick format ************************/

string evo_tree::make_newick(int precision){
    string newick;
    const boost::format tip_node_format(boost::str(boost::format("%%d:%%.%df") % precision));
    const boost::format internal_node_format(boost::str(boost::format(")%%d:%%.%df") % precision));
    const boost::format root_node_format(boost::str(boost::format(")%%d")));
    stack<Node*> node_stack;
    vector<Node*> nodes_preorder;
    Node* root;

    for(int i = 0; i < nodes.size(); ++i){
      if(nodes[i].isRoot){
          root = &nodes[i];
          break;
      }
    }
    // cout << "root " << root->id + 1 << endl;
    get_nodes_preorder(root, nodes_preorder);

    // Traverse nodes in preorder
    for(int i = 0; i < nodes_preorder.size(); i++)
    {
        Node* nd = nodes_preorder[i];
        // cout << nd->id + 1 << endl;
        if(nd->daughters.size() > 0) // internal nodes
        {
            newick += "(";
            node_stack.push(nd);
        }
        else
        {
            newick += boost::str(boost::format(tip_node_format) % (nd->id + 1) % edges[nd->e_in].length);
            if(nd->id == nodes[nd->parent].daughters[0])   //left child
                newick += ",";
            else
           {
                Node* popped = (node_stack.empty() ? 0 : node_stack.top());
                while (popped && popped->parent > 0 && popped->id == nodes[popped->parent].daughters[1]) // right sibling of the previous node
               {
                    node_stack.pop();
                    newick += boost::str(boost::format(internal_node_format) % (popped->id + 1) % edges[popped->e_in].length);
                    popped = node_stack.top();
                }
                if(popped && popped->parent > 0 && popped->id == nodes[popped->parent].daughters[0]) // left child, with another sibling
               {
                    node_stack.pop();
                    newick += boost::str(boost::format(internal_node_format) % (popped->id + 1) % edges[popped->e_in].length);
                    newick += ",";
                }
                if(node_stack.empty())
               {
                    newick += ")";
                }
            }
        }
    }
    newick +=  boost::str(boost::format(root_node_format) % (root->id + 1));
    return newick;
}


// branch length as number of mutations
string evo_tree::make_newick_nmut(int precision, const vector<int>& nmuts){
    string newick;
    const boost::format tip_node_format(boost::str(boost::format("%%d:%%.%df") % precision));
    const boost::format internal_node_format(boost::str(boost::format(")%%d:%%.%df") % precision));
    const boost::format root_node_format(boost::str(boost::format(")%%d")));
    stack<Node*> node_stack;
    vector<Node*> nodes_preorder;
    Node* root;

    for(int i = 0; i < nodes.size(); ++i){
      if(nodes[i].isRoot){
          root = &nodes[i];
          break;
      }
    }
    // cout << "root " << root->id + 1 << endl;
    get_nodes_preorder(root, nodes_preorder);

    // Traverse nodes in preorder
    for(int i = 0; i < nodes_preorder.size(); i++)
    {
        Node* nd = nodes_preorder[i];
        // cout << nd->id + 1 << endl;
        if(nd->daughters.size() > 0) // internal nodes
        {
            newick += "(";
            node_stack.push(nd);
        }
        else
        {
            newick += boost::str(boost::format(tip_node_format) % (nd->id + 1) % nmuts[nd->e_in]);
            if(nd->id == nodes[nd->parent].daughters[0])   //left child
                newick += ",";
            else
            {
                Node* popped = (node_stack.empty() ? 0 : node_stack.top());
                while (popped && popped->parent > 0 && popped->id == nodes[popped->parent].daughters[1]) // right sibling of the previous node
               {
                    node_stack.pop();
                    newick += boost::str(boost::format(internal_node_format) % (popped->id + 1) % nmuts[popped->e_in]);
                    popped = node_stack.top();
                }
                if(popped && popped->parent > 0 && popped->id == nodes[popped->parent].daughters[0]) // left child, with another sibling
                {
                    node_stack.pop();
                    newick += boost::str(boost::format(internal_node_format) % (popped->id + 1) % nmuts[popped->e_in]);
                    newick += ",";
                }
                if(node_stack.empty())
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
void evo_tree::write_nexus(const string& newick, ofstream& fout) const{
    fout << "#nexus" << endl;
    fout << "begin trees;" << endl;
    fout << "tree 1 = " << newick << ";" << endl;
    fout << "end;" << endl;
}


/********************** unused funtions (commented out for now ) ***********************/

// vector<int> evo_tree::get_internal_lengths(){
//     int debug = 0;
//     vector<int> internal_lens;
//     // Fill external edge lengths by looping over nodes
//     for(int i = 0; i < nleaf-1; ++i){
//       double len = 0;
//       vector<int> es = get_ancestral_edges( nodes[i].id );
//       // reverse(es.begin(),es.end());
//
//       // Exclude the terminal edge (index 0)
//       for(int j=1; j < es.size(); ++j){
//         len += edges[es[j]].length;
//       }
//
//       if(debug){
//           cout << "node id / edges: \t" << nodes[i].id + 1 << " : " << len << endl;
//           for(int j = 0; j < es.size(); ++j){
//             cout << "\t" << edges[es[j]].id + 1;
//           }
//           cout << endl;
//           cout << "edge " << es.back() << "; tip " <<  nodes[i].id + 1 <<  "; offset " << tobs[ nodes[i].id ] << endl;
//       }
//       internal_lens.push_back(len);
//     }
//     return internal_lens;
// }
//

// Compute the tree height until a certain tip node
// double evo_tree::get_height2sample(const int& node_id){
//       // Find the ancestral edges of the node
//       double height = 0;
//       vector<int> aedges = get_ancestral_edges(node_id);
//       for(int i = 0; i < aedges.size(); i++){
//           height += edges[aedges[i]].length;
//       }
//
//       return height;
// }


// template <typename T>
// vector<size_t> sort_indexes(const vector<T> &v){
//   // initialize original index locations
//   vector<size_t> idx(v.size());
//   iota(idx.begin(), idx.end(), 0);
//
//   // sort indexes based on comparing values in v
//   sort(idx.begin(), idx.end(),
//        [&v](size_t i1, size_t i2){return v[i1] < v[i2];});
//
//   return idx;
// }


// Find the top k time intervals with only one interval ending at a tip node
// If the tree is not of the desired topolgy, sample1 may not be the first sample
// deprecated, due to complexity s
// void evo_tree::get_ntime_interval(int k, int& sample1){
//     int debug = 0;
//
//     assert(k < node_times.size());
//     top_tinvls.clear();
//     top_tnodes.clear();
//     bool has_tip = false;
//     // Record the index of the sorted elements
//     for(auto i: sort_indexes(node_times)){
//         if(top_tnodes.size() >= k + 1){
//             break;
//         }
//         if(i == nleaf-1 || (i < nleaf - 1 && tobs[i] > 0)) continue;   // Skip the normal leaf node
//         if(has_tip && i < nleaf - 1) continue;  // only allow one tip node with tobs 0
//         if(!has_tip && i < nleaf - 1){
//             has_tip = true;
//             sample1 = i;
//         }
//         // cout << i + 1 << "\t" << node_times[i] << endl;
//         top_tnodes.push_back(i);
//     }
//
//     for(int i = 0; i < k; i++){
//         double invl = node_times[top_tnodes[i + 1]] - node_times[top_tnodes[i]];
//         top_tinvls.push_back(invl);
//     }
//
//     if(debug){
//         cout << "Top " << top_tinvls.size() << " intervals " << endl;
//         for(int i = 0; i < top_tinvls.size(); i++){
//             cout << i + 1 << "\t" << top_tinvls[i] << endl;
//         }
//     }
// }
