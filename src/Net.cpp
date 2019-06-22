//
//  Net.cpp
//
//
//  Created by Guagnlei on 03/06/2017.
//

#include <gravity/Net.h>
#include <algorithm>
#include <map>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <queue>
#include <sstream>
#include <string>
//#define USEDEBUG
#ifdef USEDEBUG
#define Debug(x) cout << x
#else
#define Debug(x)
#endif
#define DebugOn(x) cout << x
#define DebugOff(x)

using namespace std;

static int max_line_len;
static char* line = nullptr;

Net::Net() {
  // horton_net = nullptr;
  _clone = nullptr;
}

/* returns true if an arc is already present between the given nodes */
bool Net::duplicate(int n1, int n2, int id1) {
  int id2 = get_arc(n1, n2)->id;
  if (id2 < id1)
    return true;
  else
    return false;
}

Net* Net::clone() {
  Net* copy_net = new Net();
  Node* node = NULL;

  for (int i = 0; i < nodes.size(); i++) {
    node = this->nodes[i];
    copy_net->add_node(node->clone());
  }

  Arc* arc = NULL;
  for (int i = 0; i < arcs.size(); i++) {
    /* ignores if the arc is a paralel line to an already existing arc */
    // if (duplicate(arcs[i]->_src->_name, arcs[i]->_dest->_name, arcs[i]->_id))
    // {
    //            continue;
    //        }
    arc = arcs[i]->clone();

    /* Update the source and destination to the new nodes in copy_net */
    arc->src = copy_net->get_node(arc->src->_name);
    arc->dest = copy_net->get_node(arc->dest->_name);

    /* Add the new arc to the list of arcs */
    copy_net->add_arc(arc);

    /* Connects it to its source and destination */
    arc->connect();
  }
  return copy_net;
}

const bool node_compare(const Node* n1, const Node* n2) {
  return n1->fill_in > n2->fill_in;
}

void Net::add_node(Node* node) {
  node->ID = (int)nodes.size();  // depends on the size.

  //  nodeID.insert(pair<string,Node*>(_name,node)).

  if (!nodeID.insert(pair<string, Node*>(node->_name, node)).second) {
    cerr << "ERROR: adding the same node twice!";
  }
  nodes.push_back(node);
}

void Net::print_nodes(){
  for (auto iter = nodes.begin(); iter != nodes.end()-1; iter++){
    cout << to_string((*iter)->ID) << ",";
  }
  cout << to_string(nodes.back()->ID) << endl;
  for (auto iter = nodes.begin(); iter != nodes.end()-1; iter++){
    cout << to_string((*iter)->fill_in) << ",";
  }
  cout << to_string(nodes.back()->fill_in) << endl;
}

Node* Net::get_node(string name) {
  if (nodeID.find(name) != nodeID.end())
    return nodeID.find(name)->second;
  else
    return nullptr;
}

/* returns the arc formed by node ids n1 and n2 */
Arc* Net::get_arc(Node* n1, Node* n2) {
  string src, dest, key, inv_key;
  src = n1->_name;
  dest = n2->_name;
  key.clear();
  inv_key.clear();
  key.append(src);
  inv_key.append(dest);
  key.append(",");
  inv_key.append(",");
  key.append(dest);
  inv_key.append(src);
  map<string, set<Arc*>*>::iterator it = lineID.find(key);
  if (it != lineID.end()) {
    for (auto a : *it->second) {
      //   if (!a->parallel) {
      return a;
      // }
    }
  }
  it = lineID.find(inv_key);
  if (it != lineID.end()) {
    for (auto a : *it->second) {
      //   if (!a->parallel) {
      return a;
      // }
    }
  }
  return nullptr;
}

/* returns the Id of the arc formed by node ids n1 and n2 */
Arc* Net::get_arc(int n1, int n2) {
  string src, dest, key, inv_key;
  src = to_string(n1);
  dest = to_string(n2);
  key.clear();
  inv_key.clear();
  key.append(src);
  inv_key.append(dest);
  key.append(",");
  inv_key.append(",");
  key.append(dest);
  inv_key.append(src);
  map<string, set<Arc*>*>::iterator it = lineID.find(key);
  if (it != lineID.end()) {
    for (auto a : *it->second) {
      //  if (!a->parallel) {
      return a;
      //}
    }
  }
  it = lineID.find(inv_key);
  if (it != lineID.end()) {
    for (auto a : *it->second) {
      // if (!a->parallel) {
      return a;
      //}
    }
  }
  return nullptr;
}

bool Net::add_arc(Arc* a) {
  bool parallel = false;
  set<Arc*>* s = NULL;
  string src, dest, key, inv_key;
  src = a->src->_name;
  dest = a->dest->_name;
  key.clear();
  inv_key.clear();
  key.append(src);
  inv_key.append(dest);
  key.append(",");
  inv_key.append(",");
  key.append(dest);
  inv_key.append(src);
  if (lineID.find(key) == lineID.end() &&
      lineID.find(inv_key) == lineID.end()) {
    s = new set<Arc*>;
    s->insert(a);
    lineID.insert(pair<string, set<Arc*>*>(key, s));
  } else {
    if (lineID.find(key) != lineID.end()) s = lineID[key];
    if (lineID.find(inv_key) != lineID.end()) s = lineID[inv_key];
    s->insert(a);
    cerr << "\nWARNING: adding another line between same nodes! \n Node name: "
         << src << " and Node name: " << dest << endl;
    exit(1);
  }
  arcs.push_back(a);
  return parallel;
}

void Net::remove_arc(Arc* a) {
  // arcs.erase(arcs.at(a->id));
  arcs[a->id] = nullptr;
  lineID.erase(a->src->_name + "," + a->dest->_name);
  // remove arc from their branches.
  auto u = a->src;
  auto v = a->dest;
  u->removeArc(a);
  v->removeArc(a);
  //    auto it = std::find(u->branches.begin(), u->branches.end(), a);
  //    if (it != u->branches.end()){
  //        u->branches.erase(it);
  //    }
  //    it = std::find(v->branches.begin(), v->branches.end(), a);
  //    if (it != v->branches.end()){
  //        v->branches.erase(it);
  //    }
}

// Reading files
char* Net::readline(FILE* input) {
  size_t len;
  // line, max_line_len have been declared
  if (std::fgets(line, max_line_len, input) == NULL) return NULL;

  while (strrchr(line, '\n') == NULL) {
    max_line_len *= 2;
    line = (char*)realloc(line, max_line_len);
    len = strlen(line);
    if (fgets(line + len, max_line_len - len, input) == NULL) break;
  }
  return line;
}

void Net::exit_input_error(int line_num) {
  fprintf(stderr, "Wrong input format at line %d\n", line_num);
  exit(1);
}

// readFile: just read a matrix, nothing new!
void Net::readFile(string fn) {
  auto fname = fn.c_str();
  FILE* fp = fopen(fname, "r");
  if (fp == NULL) {
    fprintf(stderr, "can’t open input file %s\n", fname);
    exit(1);
  }

  size_t max_line_len = 1024;
  char* line = new char[max_line_len];

  vector<vector<int>> matrix;
  int temp;

  stringstream linestream(line);

  while (readline(fp) != NULL) {
    vector<int> row;
    stringstream linestream(line);
    while (linestream >> temp) row.push_back(temp);
    matrix.push_back(row);
    // cout <<matrix.size()<< endl;
  }
  rewind(fp);
  delete[] line;
  fclose(fp);
}

// read rudy
void Net::readrudy(string fn) {
  auto fname = fn.c_str();
  int Num_nodes = 0;
  int Num_edges = 0;
  ifstream infile(fname);
  string sLine;

  if (infile.good()) {
    getline(infile, sLine);
    istringstream iss(sLine);
    iss >> Num_nodes;
    iss >> Num_edges;
  } else {
    fprintf(stderr, "can’t open input file %s\n", fname);
    exit(1);
  }
  string name;
  _clone = new Net();
  _chordalextension = new Net();
  Node* node = nullptr;
  Node* node_clone = nullptr;
  Node* node_chordal = nullptr;

  for (int i = 0; i < Num_nodes; i++) {
    name = to_string(i + 1);
    node = new Node(name, i);
    node_clone = new Node(name, i);
    node_chordal = new Node(name, i);
    add_node(node);
    _clone->add_node(node_clone);
    _chordalextension->add_node(node_chordal);
    // cout << "size of chordal extension is: " <<
    // _chordalextension->nodes.size() << endl;
  }

  // get arcs
  Arc* arc = NULL;
  Arc* arc_clone = NULL;
  Arc* arc_chordal = NULL;

  // note that src, dest are names of nodes.
  string src, dest;
  double weight;
  while (getline(infile, sLine, '\n')) {
    istringstream iss(sLine);
    iss >> src >> dest >> weight;
    // cout << src  << ", " << dest << ", " << weight << endl;

    name = src + "," + dest;  //

    arc = new Arc(name);
    arc_clone = new Arc(name);
    arc_chordal = new Arc(name);

    arc->id = (int)arcs.size();
    arc_clone->id = (int)_clone->arcs.size();
    arc_chordal->id = (int)_chordalextension->arcs.size();

    arc->src = get_node(src);
    arc->dest = get_node(dest);
    arc->weight = weight;
    add_arc(arc);
    arc->connect();

    arc_clone->src = _clone->get_node(src);
    arc_clone->dest = _clone->get_node(dest);
    arc_clone->weight = weight;
    _clone->add_arc(arc_clone);
    arc_clone->connect();

    arc_chordal->src = _chordalextension->get_node(src);
    arc_chordal->dest = _chordalextension->get_node(dest);
    arc_chordal->weight = weight;
    _chordalextension->add_arc(arc_chordal);
    arc_chordal->connect();
  }
  infile.close();
}

// populates the graph
void Net::topology(string fn, bool complement) {
  auto fname = fn.c_str();
  FILE* fp = fopen(fname, "r");
  if (fp == NULL) {
    fprintf(stderr, "can’t open input file %s\n", fname);
    exit(1);
  }

  max_line_len = 1024;
  line = new char[max_line_len];

  vector<vector<int>> matrix;
  int temp;
  while (readline(fp) != NULL) {
    vector<int> row;
    stringstream linestream(line);
    while (linestream >> temp) row.push_back(temp);
    matrix.push_back(row);
  }
  int n = 0;
  n = matrix.size();

  string name;
  int id = 0;
  _clone = new Net();

  Node* node = NULL;
  Node* node_clone = NULL;
  for (int i = 0; i < n; i++) {
    name = to_string(i);
    node = new Node(name, i);
    node_clone = new Node(name, i);
    add_node(node);
    _clone->add_node(node_clone);
  }

  Arc* arc = NULL;
  Arc* arc_clone = NULL;
  string src, dest;

  if (complement) {
    for (int i = 0; i < (n - 1); i++)
      for (int j = i + 1; j < n; j++) {
        if (matrix[i][j] == 0) {
          src = to_string(i);
          dest = to_string(j);

          id = (int)arcs.size();
          arc = new Arc(to_string(id));
          arc_clone = new Arc(to_string(id));
          arc->id = id;
          arc_clone->id = id;
          arc->src = get_node(src);
          arc->dest = get_node(dest);
          arc_clone->src = _clone->get_node(src);
          arc_clone->dest = _clone->get_node(dest);
          add_arc(arc);
          arc->connect();
          _clone->add_arc(arc_clone);
          arc_clone->connect();
        }
      }
  } else {
    for (int i = 0; i < (n - 1); i++)
      for (int j = i + 1; j < n; j++) {
        if (matrix[i][j] > 0) {
          src = to_string(i);
          dest = to_string(j);

          id = (int)arcs.size();
          arc = new Arc(to_string(id));
          arc_clone = new Arc(to_string(id));
          arc->id = id;
          arc_clone->id = id;
          arc->src = get_node(src);
          arc->dest = get_node(dest);
          arc_clone->src = _clone->get_node(src);
          arc_clone->dest = _clone->get_node(dest);
          add_arc(arc);
          arc->connect();
          _clone->add_arc(arc_clone);
          arc_clone->connect();
        }
      }
  }
  delete[] line;
  fclose(fp);
  cout << "Edges: " << arcs.size() << endl;
}

Net Net::get_complement(string fn) {
  Net complement;
  auto fname = fn.c_str();
  FILE* fp = fopen(fname, "r");
  if (fp == NULL) {
    fprintf(stderr, "can’t open input file %s\n", fname);
    exit(1);
  }

  max_line_len = 1024;
  line = new char[max_line_len];

  vector<vector<int>> matrix;
  int temp;
  while (readline(fp) != NULL) {
    vector<int> row;
    stringstream linestream(line);
    while (linestream >> temp) row.push_back(temp);
    matrix.push_back(row);
  }
  int n = 0;
  n = matrix.size();

  string name;
  int id = 0;
  _clone = new Net();

  Node* node = NULL;
  Node* node_clone = NULL;
  for (int i = 0; i < n; i++) {
    name = to_string(i);
    node = new Node(name, i);
    node_clone = new Node(name, i);
    complement.add_node(node);
    complement._clone->add_node(node_clone);
  }

  Arc* arc = NULL;
  Arc* arc_clone = NULL;
  string src, dest;

  for (int i = 0; i < (n - 1); i++)
    for (int j = i + 1; j < n; j++) {
      if (matrix[i][j] == 0) {
        src = to_string(i);
        dest = to_string(j);

        id = (int)arcs.size();
        arc = new Arc(to_string(id));
        arc_clone = new Arc(to_string(id));
        arc->id = id;
        arc_clone->id = id;
        arc->src = get_node(src);
        arc->dest = get_node(dest);
        arc_clone->src = _clone->get_node(src);
        arc_clone->dest = _clone->get_node(dest);
        complement.add_arc(arc);
        arc->connect();
        complement._clone->add_arc(arc_clone);
        arc_clone->connect();
      }
    }
  delete[] line;
  fclose(fp);
  return complement;
}

/*  @brief Remove node and all incident arcs from the network
 @note Does not remove the incident arcs from the list of arcs in the network!
 @return the id of the node removed
 */
string Net::remove_end_node() {
  Node* n = nodes.back();
  Node* nn = nullptr;
  string n_id = n->_name;
  for (auto a : n->branches) {
    nn = a->neighbour(n);
    nn->removeArc(a);
    for (auto aa : nn->branches) {
      if (!aa->neighbour(nn)->is_connected(n)) {
        nn->fill_in--;
        assert(nn->fill_in >= 0);
      }
    }
  }
  //    delete nodes.back();
  nodes.pop_back();
  return n_id;
}

void Net::remove_node(Node* n) {
  nodes[n->ID] = nullptr;
  nodeID.erase(n->_name);
  Node* nn = nullptr;
  for (auto a : n->branches) {
    nn = a->neighbour(n);
    remove_arc(a);
    for (auto aa : nn->branches) {
      if (!aa->neighbour(nn)->is_connected(n)) {
        nn->fill_in--;
        assert(nn->fill_in >= 0);
      }
    }
  }
}

const bool bag_compare(const vector<Node*>& a, const vector<Node*>& b) {
  return a.size() > b.size();
}

// find cliques of this graph.
void Net::get_tree_decomp_bags(bool print_bags) {
  Node* n = nullptr;
  Node* u = nullptr;
  Node* nn = nullptr;
  Arc* arc = nullptr;
  Arc* arc_chordal = nullptr;

  Node* u_chordal = nullptr;
  Node* nn_chordal = nullptr;

  string name = "";
  string name_chordal = "";
  _chordalextension = clone();

  Net* graph_clone = clone();

  /** cliques with 1 nodes are useless for us.*/
  //sort(graph_clone->nodes.begin(), graph_clone->nodes.end(), node_compare);
  //graph_clone->print_nodes();
  while (graph_clone->nodes.size() > 1) {
    sort(graph_clone->nodes.begin(), graph_clone->nodes.end(), node_compare);
    // last element has the minimum fill-in.
    n = graph_clone->nodes.back();
    vector<Node*> bag_copy;
    vector<Node*> bag; // store original nodes

    vector<Node*> N = n->get_neighbours();
    for (auto iter = N.begin(); iter != N.end(); ++iter) {
      nn = *iter;
      bag_copy.push_back(nn);
      bag.push_back(get_node(nn->_name));
    }
    graph_clone->remove_end_node();
    bag_copy.push_back(n);
    bag.push_back(get_node(n->_name));
    sort(bag_copy.begin(), bag_copy.end(),
         [](const Node* a, const Node* b) -> bool { return a->ID < b->ID; });
    sort(bag.begin(), bag.end(),
         [](const Node* a, const Node* b) -> bool { return a->ID < b->ID; });

    // update graph_graph and construct chordal extension.
    for (int i = 0; i < bag_copy.size() - 1; i++) {
      u = bag_copy.at(i);
      u_chordal = _chordalextension->get_node(u->_name);
      for (int j = i + 1; j < bag_copy.size(); j++) {
        nn = bag_copy.at(j);
        nn_chordal = _chordalextension->get_node(nn->_name);
        if (u->is_connected(nn)) {
          continue;
        }
        name = to_string((int)graph_clone->arcs.size() + 1);
        name_chordal = to_string((int)_chordalextension->arcs.size() + 1);

        arc = new Arc(name);
        arc_chordal = new Arc(name_chordal);

        arc->id = arcs.size();
        arc->src = u;
        arc->dest = nn;
        arc->connect();
        graph_clone->add_arc(arc);

        arc_chordal->id = _chordalextension->arcs.size();
        arc_chordal->src = u_chordal;
        arc_chordal->dest = nn_chordal;
        arc_chordal->connect();
        _chordalextension->add_arc(arc_chordal);
      }
    }
    _bags.push_back(bag);
    delete n;
  }
  // sort the bags by its size (descending order)
  sort(_bags.begin(), _bags.end(), bag_compare);
  printf(
      "With the greedy fill-in algirithm, the chordal graph added  %lu edges "
      "\n",
      (_chordalextension->arcs.size() - arcs.size()));
  printf("the number of bags is %lu \n", _bags.size());

  delete graph_clone;
}

// get cliques from the tree decomposition
// Two methods
// first one: check the inclusion relationship
// second one: use the RIP property of the tree decomposition, thus just need to
// check every leaf
void Net::get_cliquebags(bool print) {
  get_tree_decomp_bags(print);
  for (unsigned i = 0; i < _bags.size() - 1; i++) {
    for (unsigned j = i + 1; j < _bags.size();) {
      if (std::includes(_bags[i].begin(), _bags[i].end(), _bags[j].begin(),
                        _bags[j].end())) {
        _bags.erase(_bags.begin() + j);
      } else
        j++;
    }
  }
  cout << "Number of maximal cliques of the chordal extension = "
       << _bags.size() << endl;
}
/* Find the root of a node */
Node* find_root(Node* u) {
  /* Make the parent of the nodes in the path
   from u--> parent[u] point to parent[u] */
  if (u->ID != u->predecessor->ID) u->predecessor = find_root(u->predecessor);
  return u->predecessor;
}

/* check if two nodes are reachable by checking if their roots conincides*/
bool is_reachable(Node* n1, Node* n2) {
  return (find_root(n1)->ID == find_root(n2)->ID);
}
std::vector<Arc*> Net::minimal_spanning_tree_kruskal() {
  std::vector<Arc*> MST;
  // initilize predecessors
  int* rank = new int[nodes.size() + 1];
  for (auto n : nodes) {
    n->predecessor =
        n;  // set predecessor to itself (thus there are n disjoint components.)
    rank[n->ID] = 0;
  }

  // sort edges in an increasing order on basis of cost
  sort(arcs.begin(), arcs.end(), [](const Arc* a, const Arc* b) -> bool {
    return a->weight < b->weight;
  });
  double MST_weight = 0.0;
  /* Check if adding the current edge will form a cycle (using disjoint set) */
  for (auto a : arcs) {
    if (!is_reachable(a->src, a->dest)) {
      MST_weight += fabs(a->weight);
      MST.push_back(a);
      // cout << a->src->_name << "-" << a->dest->_name <<endl;
      // make them reachable by letting their predecessors reachable
      auto n1 = a->src->predecessor;
      auto n2 = a->dest->predecessor;
      if (rank[n2->ID] < rank[n1->ID])
        n1->predecessor = n2;
      else
        n2->predecessor = n1;
      if (rank[n2->ID] == rank[n1->ID]) rank[n2->ID]++;
    }
  }
  MST.resize(MST.size());
  cout << "MST weight = " << MST_weight << endl;
  return MST;
}

Net* Net::get_clique_tree_kruskal() {
  Net* graph_intersection = new Net();
  Net* cliquetree = new Net();
  string name;

  for (int i = 0; i < _bags.size(); i++) {
    auto node = new Node(to_string(i), i);
    auto node_ct = node->clone();
    graph_intersection->add_node(node);
    cliquetree->add_node(node_ct);
    sort(_bags[i].begin(), this->_bags[i].end(),
         [](const Node* a, const Node* b) -> bool { return a->ID < b->ID; });
  }
  unsigned weights_interesction = 0;
  int nb_cliques = _bags.size();

  vector<Node*> v3;
  for (int i = 0; i < nb_cliques; i++) {
    for (int j = i + 1; j < nb_cliques; j++) {
      v3.clear();
      set_intersection(
          _bags.at(i).begin(), _bags.at(i).end(), _bags.at(j).begin(),
          _bags.at(j).end(), back_inserter(v3),
          [](const Node* a, const Node* b) -> bool { return a->ID < b->ID; });
      if ((int)v3.size() > 0) {
        auto arc = new Arc(graph_intersection->nodes[i],
                           graph_intersection->nodes[j], -(int)(v3.size()));
        arc->connect();
        for (auto nn : v3) {
          arc->_intersection.push_back(nn);
        }
        graph_intersection->add_arc(arc);
        weights_interesction += (int)v3.size();
      }
    }
  }
  DebugOn("weights of the intersection graph: " << weights_interesction
                                                << endl);
  std::vector<Arc*> MST = graph_intersection->minimal_spanning_tree_kruskal();
  DebugOff("Print the total " << spanning_tree.size()
                              << " edges in the clique tree:" << endl);
  //////////CLIQUE TREE /////////////////////////////
  double total_weight = 0;

  for (auto e : MST) {
    int n1 = e->src->ID;
    int n2 = e->dest->ID;
    total_weight -= e->weight;
    Debug(e->src->_name << " <--> " << e->dest->_name
                        << " weight: " << -e->weight << endl);
    auto arc_ct =
        new Arc(cliquetree->nodes[n1], cliquetree->nodes[n2], -e->weight);
    arc_ct->id = cliquetree->arcs.size();
    arc_ct->_intersection.insert(arc_ct->_intersection.begin(),
                                 e->_intersection.begin(),
                                 e->_intersection.end());
    arc_ct->connect();
    cliquetree->add_arc(arc_ct);
  }
  DebugOn("total weight of the clique tree is: " << total_weight << endl);
  return cliquetree;
}

void Net::depth_first_vist(Node* root) {
  for (auto n : root->get_neighbours()) {
    if (root->predecessor != n) {
      n->predecessor = root;
      depth_first_vist(n);
    }
  }
}

/* Destructors */
Net::~Net() {
  if (!nodes.empty()) {
    for (vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++) {
      delete (*it);
    }
    nodes.clear();
  }
  if (!arcs.empty()) {
    for (vector<Arc*>::iterator it = arcs.begin(); it != arcs.end(); it++) {
      if (*it) delete (*it);
    }
    arcs.clear();
  }

  if (!cycle_basis.empty()) {
    for (vector<Path*>::iterator it = cycle_basis.begin();
         it != cycle_basis.end(); it++) {
      delete (*it);
    }
    arcs.clear();
  }
  for (pair<string, set<Arc*>*> it : lineID) {
    delete it.second;
  }
  delete _clone;
}

int Net::test() {
  string name;
  int id = 0;
  _clone = new Net();

  Node* node = NULL;
  Node* node_clone = NULL;

  for (int i = 0; i < 3; i++) {
    name = to_string(i);
    id = i;
    //      cout << "name " << name << " ID: " << i << endl;
    node = new Node(name, id);
    node_clone = new Node(name, i);
    add_node(node);
    _clone->add_node(node_clone);
  }

  Arc* arc = NULL;
  Arc* arc_clone = NULL;
  string src, dest;
  for (int i = 0; i < 2; i++) {
    src = to_string(i);
    dest = to_string(i + 1);
    id = (int)arcs.size();
    // arc = new Arc(src+dest);
    // arc_clone = new Arc(src+dest);
    arc = new Arc(to_string(id));
    arc_clone = new Arc(to_string(id));
    arc->id = id;
    arc_clone->id = id;
    arc->src = get_node(src);
    arc->dest = get_node(dest);
    arc_clone->src = _clone->get_node(src);
    arc_clone->dest = _clone->get_node(dest);
    add_arc(arc);
    arc->connect();
    _clone->add_arc(arc_clone);
    arc_clone->connect();
  }
  id = (int)arcs.size();
  // arc = new Arc(src+dest);
  // arc_clone = new Arc(src+dest);
  arc = new Arc(to_string(id));
  arc_clone = new Arc(to_string(id));
  arc->id = id;
  arc_clone->id = id;
  arc->src = get_node("2");
  arc->dest = get_node("0");
  arc_clone->src = _clone->get_node("2");
  arc_clone->dest = _clone->get_node("0");
  add_arc(arc);
  arc->connect();
  _clone->add_arc(arc_clone);
  arc_clone->connect();

  //    get_tree_decomp_bags();
  return 0;
}
