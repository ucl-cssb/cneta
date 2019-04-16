//evo_tree.hpp

#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

using namespace std;

class node {
public:
  int id;
  int isRoot;
  int isLeaf;
  int parent;

  int e_in;
  vector<int> e_ot;
  vector<int> daughters;

  node(const int& _id, const int& _isRoot, const int& _isLeaf){
    id = _id;
    isRoot = _isRoot;
    isLeaf = _isLeaf;
    parent = -1;
    e_in = -1;
  }

  node( const node& _n2 ){
    id     = _n2.id;
    isRoot = _n2.isRoot;
    isLeaf = _n2.isLeaf;
    parent = _n2.parent;
    e_in   = _n2.e_in;

    e_ot.clear();
    e_ot.insert(e_ot.end(), _n2.e_ot.begin(), _n2.e_ot.end() );
    daughters.clear();
    daughters.insert(daughters.end(), _n2.daughters.begin(), _n2.daughters.end() );
  }
};

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
  vector<node>   nodes;
  vector<double> node_times;
  // pointers to internal edges
  vector<edge*> intedges;
  int root_node_id;
  vector<vector<int> > chars;
  vector<double> tobs;
  double score;
  double mu;

  evo_tree(){}
  evo_tree(const int& _nleaf, const vector<int>& _edges, const vector<double>& _lengths, int gen_node = 1);
  evo_tree(const int& _nleaf, const vector<edge>& _edges, int gen_node = 1);

  // reparameterized tree: expects terminal edges to all have length 0
  evo_tree(const int& _nleaf, const vector<edge>& _edges, const double& total_time, const vector<double>& _tobs);

  evo_tree(const evo_tree& _t2);

  void   generate_nodes();
  void   calculate_node_times();
  void   generate_int_edges();
  // The time from beginning to the time of first sample
  double get_total_time(){ return *max_element(node_times.begin(), node_times.end()) - *max_element(tobs.begin(), tobs.end()); }

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
};

evo_tree::evo_tree(const int& _nleaf, const vector<int>& _edges, const vector<double>& _lengths, int gen_node){
  nleaf = _nleaf;
  nnode = nleaf - 1; // internal nodes
  nedge = 2*nleaf - 2;
  ntotn = 2*nleaf - 1; // all nodes
  nintedge = nedge - nleaf;

  // create list of edges
  int count = 0;
  for(int i=0; i<nedge; ++i){
    edges.push_back( edge(i,_edges[count], _edges[count+1], _lengths[i]) );
    count = count + 2;
  }

  if(gen_node){
    generate_nodes();
    calculate_node_times();
    generate_int_edges();
  }
  score = 0;
}

evo_tree::evo_tree(const int& _nleaf, const vector<edge>& _edges, int gen_node){
  nleaf = _nleaf;
  nnode = nleaf - 1; // internal nodes
  nedge = 2*nleaf - 2;
  ntotn = 2*nleaf - 1; // all nodes
  nintedge = nedge - nleaf;

  edges.clear();
  edges.insert(edges.end(), _edges.begin(), _edges.end() );

  if(gen_node){
    generate_nodes();
    calculate_node_times();
    generate_int_edges();
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


void evo_tree::generate_nodes(){
  // create list of nodes
  for(int i=0; i<ntotn; ++i){
    if(i < nleaf){
      nodes.push_back( node(i,0,1));
    }
    else{
      nodes.push_back( node(i,0,0));
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

  edges.clear();
  edges.insert(edges.end(), _t2.edges.begin(), _t2.edges.end() );
  lengths.clear();
  lengths.insert(lengths.end(), _t2.lengths.begin(), _t2.lengths.end() );
  nodes.clear();
  nodes.insert(nodes.end(), _t2.nodes.begin(), _t2.nodes.end() );
  node_times.clear();
  node_times.insert(node_times.end(), _t2.node_times.begin(), _t2.node_times.end() );

  chars.clear();
  for(int i=0; i< _t2.chars.size(); ++i) chars.push_back( _t2.chars[i] );

  intedges.clear();
  generate_int_edges();
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
