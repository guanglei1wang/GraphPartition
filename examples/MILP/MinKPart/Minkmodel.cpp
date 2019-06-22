//
//  MinkMinkmodel.cpp
//  Gravity
//
//  Created by Guanglei Wang on 19/6/17.
//
//

#include "Minkmodel.hpp"
//#define USEDEBUG
#ifdef USEDEBUG
#define Debug(x) cout << x
#else
#define Debug(x)
#endif
#define DebugOn(x) cout << x
#define DebugOff(x)
#define SIGMA 0.1

using namespace std;

Minkmodel::Minkmodel(){};

Minkmodel::~Minkmodel(){};

Minkmodel::Minkmodel(ModelType type, Net* graph, double K)
    : _type(type), _solver(cplex), _K(K), _graph(graph) {
  _cliqueid = make_shared<map<string, vector<unsigned>>>();
};

Minkmodel::Minkmodel(ModelType type, Net* graph, double K, SolverType solver)
    : _type(type), _solver(solver), _K(K), _graph(graph) {
  _cliqueid = make_shared<map<string, vector<unsigned>>>();
};

void Minkmodel::build() {
  switch (_type) {
    case MIP:
      add_vars_origin();
      add_triangle();
      add_clique();
      break;
    case SDP:
      add_vars_lifted();
      break;
    case MIP_tree:
      add_vars_origin_tree();
      cliquetree_decompose();
      add_triangle_tree();
      add_clique_tree();
      break;
    case SDP_tree:
      // add_vars_lifted_tree();
      add_vars_lifted_tree_compact();
      break;
    case SDP_tree_global:
      add_vars_lifted_tree();
      break;
    case Node_edge:
      node_edge_formulation();
      break;
    default:
      break;
  }
}
void Minkmodel::reset(){};
void Minkmodel::add_vars_origin() {
  var<bool> zij("zij");
  _model.add_var(zij ^ (_graph->nodes.size() * (_graph->nodes.size() - 1) / 2));

  func_ obj_MIP;
  int i = 0, j = 0;
  for (auto a : _graph->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    if (i <= j)
      obj_MIP += (a->weight) * zij(i, j);
    else
      obj_MIP += (a->weight) * zij(j, i);
  }
  _model.set_objective(min(obj_MIP));
}

void Minkmodel::add_vars_origin_tree() {
  var<bool> zij("zij");
  // the number of arcs in the chordal extension
  _model.add_var(zij ^ ((_graph->_chordalextension)->arcs.size()));
  func_ obj_MIP;
  int i = 0, j = 0;
  for (auto a : _graph->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    if (i <= j)
      obj_MIP += (a->weight) * zij(i, j);
    else
      obj_MIP += (a->weight) * zij(j, i);
  }
  _model.set_objective(min(obj_MIP));
}

void Minkmodel::add_vars_lifted() {
  sdpvar<double> X("X");
  // ^ specifies that X._symdim = _graph->nodes.size();
  _model.add_var(X ^ (_graph->nodes.size()));

  int i = 0, j = 0;
  func_ obj;
  for (auto a : _graph->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    if (i <= j)
      obj += a->weight * ((_K - 1) * X(i, j) + 1) / _K;
    else
      obj += a->weight * ((_K - 1) * X(j, i) + 1) / _K;
  }
  _model.set_objective(min(obj));

  for (i = 0; i < _graph->nodes.size(); i++) {
    Constraint diag("(" + to_string(i) + ")");
    diag = X(i, i) - 1;
    _model.add_constraint(diag = 0);
  }

  for (i = 0; i < _graph->nodes.size(); i++)
    for (j = i + 1; j < _graph->nodes.size(); j++) {
      Constraint bound("(" + to_string(i) + "," + to_string(j) + ")");
      bound = X(i, j) + 1 / (_K - 1);
      _model.add_constraint(bound >= 0);
    }
}

bool node_comparator(Node* a, Node* b) { return a->ID < b->ID; }

template <typename _cont, typename _func>
inline _cont my_set_intersection(_cont a, _cont b, _func c) {
  _cont res;
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        back_inserter(res), c);
  return res;
}

template <typename _cont, typename _func>
inline _cont my_set_union(_cont a, _cont b, _func c) {
  _cont res;
  std::set_union(a.begin(), a.end(), b.begin(), b.end(), back_inserter(res), c);
  return res;
}

/* This is an implmentation of merging tricks proposed by Nakata 2003 in
 * MathProg **/
void brother_merge(Net* graph, Net* T, Node* b1, Node* b2) {
  // merge corresponding bags and let bag b2 empty
  auto merged = my_set_union(graph->_bags.at(b1->ID), graph->_bags.at(b2->ID),
                             node_comparator);
  graph->_bags[b1->ID].clear();
  auto it = graph->_bags[b1->ID].begin();
  graph->_bags[b1->ID].insert(it, merged.begin(), merged.end());
  graph->_bags.at(b2->ID).clear();

  // add new arcs to b2' neighbours
  auto children = b2->get_neighbours();
  for (auto nn : children) {
    if (nn != b1->predecessor) {
      // remove old arc
      auto A = T->get_arc(nn, b2);
      T->remove_arc(A);
      // according CIP, the weight does not change
      Arc* a = new Arc(b1, nn, A->weight);
      a->id = T->arcs.size();
      // a->_intersection.insert(a->_intersection.begin(),A->_intersection.begin(),
      // A->_intersection.end());
      nn->predecessor = b1;
      T->add_arc(a);
    }
  }
  // remove node b2 in the clique graph
  T->remove_node(b2);
  // T->update_children(b1->predecessor); // update the list of unvisted
  // nodes.
}

void merge_parent_brother(Net* graph, Net* ct, Node* root, Node* n) {
  // n is child.
  // merge
  vector<Node*> merged;
  std::set_union(graph->_bags.at(n->ID).begin(), graph->_bags.at(n->ID).end(),
                 graph->_bags.at(root->ID).begin(),
                 graph->_bags.at(root->ID).end(), back_inserter(merged),
                 [](const Node* a, Node* b) -> bool { return a->ID < b->ID; });
  // merge
  graph->_bags[root->ID].clear();
  auto it = graph->_bags[root->ID].begin();
  graph->_bags[root->ID].insert(it, merged.begin(), merged.end());
  graph->_bags.at(n->ID).clear();
  // remove arc betten root and n
  ct->remove_arc(ct->get_arc(n, root));
  // making the corresponding bag empty.
  graph->_bags[n->ID].clear();
  // add new arcs
  for (auto nn : n->get_neighbours()) {
    if (nn != root) {
      // remove old arc
      auto A = ct->get_arc(nn, n);
      ct->remove_arc(A);
      // add new arc with new weight
      // according CIP, the weight does not change
      Arc* a = new Arc(root, nn, A->weight);
      a->id = ct->arcs.size();
      a->connect();
      a->_intersection.insert(a->_intersection.begin(),
                              A->_intersection.begin(), A->_intersection.end());
      nn->predecessor = root;
      ct->add_arc(a);
    }
  }
  // remove node n in the clique graph.
  ct->remove_node(n);
}
void depth_first_merge(Net* graph, Net* ct, Node* root) {
  std::vector<Node*> merged;
  auto ns = root->get_neighbours();
  while (ns.size() > 0) {
    auto n = ns.back();
    ns.pop_back();
    if (n != root->predecessor) {
      double a = ct->get_arc(root, n)->weight / graph->_bags[n->ID].size();
      double b = ct->get_arc(root, n)->weight / graph->_bags[root->ID].size();
      if (min(a, b) > SIGMA) {
        merge_parent_brother(graph, ct, root, n);
        ns = root->get_neighbours();
      } else {
        // next level;
        n->predecessor = root;
        depth_first_merge(graph, ct, n);
      }
    }
  }
}

// merging algorithm due to Nakata 2003.
void merge_cliques(Net* graph, Net* cliquetree, Node* root) {
  root->predecessor = root;
  cliquetree->depth_first_vist(root);
  depth_first_merge(graph, cliquetree, root);
  // clean clique tree by removing empty nodes and arcs
  Net* ct = new Net();
  vector<vector<Node*>> new_bags;
  for (auto n : cliquetree->nodes) {
    if (n != nullptr) {
      new_bags.push_back(graph->_bags.at(n->ID));
      // note that n id is updated
      ct->add_node(n);
    }
  }
  for (auto a : cliquetree->arcs) {
    if (a != nullptr) {
      a->id = ct->arcs.size();
      ct->add_arc(a);
    }
  }
  cliquetree->nodes.clear();
  cliquetree->nodes.insert(cliquetree->nodes.begin(), ct->nodes.begin(),
                           ct->nodes.end());
  cliquetree->arcs.clear();
  cliquetree->arcs.insert(cliquetree->arcs.begin(), ct->arcs.begin(),
                          ct->arcs.end());
  graph->_bags.clear();
  graph->_bags.insert(graph->_bags.begin(), new_bags.begin(), new_bags.end());
}

void Minkmodel::add_vars_lifted_tree() {
  if (_graph->_bags.size() == 0) {
    cerr << "ERROR: bags of cliques are not generated." << endl;
    exit(1);
  }
  auto cliquetree = _graph->get_clique_tree_kruskal();
  // merging algorithm due to Nakata 2003.
  merge_cliques(_graph, cliquetree, cliquetree->nodes.front());

  // Adding variables X.
  var<double> X("X");
  _model.add_var(X ^ (_graph->nodes.size() * (_graph->nodes.size() + 1) * 0.5));

  int clique_id = 0;
  for (auto clique : _graph->_bags) {
    // adding SDP constraint for each clique
    sdpvar<double> X_sub("X_sub_" + to_string(clique_id));
    _model.add_var(X_sub ^ (clique.size()));
    // linking X and X_sub
    for (int i = 0; i < clique.size(); i++) {
      for (int j = i; j < clique.size(); j++) {
        int index1 = clique[i]->ID;
        int index2 = clique[j]->ID;
        Constraint Linking_constraint("Linking_" + to_string(clique_id) + "(" +
                                      to_string(i) + "," + to_string(j) + ")");
        if (index1 < index2)
          Linking_constraint = X(index1, index2) - X_sub(i, j);
        else
          Linking_constraint = X(index2, index1) - X_sub(i, j);
        _model.add_constraint(Linking_constraint = 0);
      }
    }
    clique_id++;
  }

  int i = 0, j = 0;
  func_ obj;
  for (auto a : _graph->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    if (i < j)
      obj += a->weight * ((_K - 1) * X(i, j) + 1) / _K;
    else
      obj += a->weight * ((_K - 1) * X(j, i) + 1) / _K;
  }
  _model.set_objective(min(obj));

  for (i = 0; i < _graph->nodes.size(); i++) {
    Constraint diag("(" + to_string(i) + ")");
    diag = X(i, i) - 1;
    _model.add_constraint(diag = 0);
  }

  for (auto a : _graph->_chordalextension->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    Constraint bound("(" + to_string(i) + "," + to_string(j) + ")");
    if (i < j) {
      bound = X(i, j) + 1 / (_K - 1);
    } else
      bound = X(j, i) + 1 / (_K - 1);
    _model.add_constraint(bound >= 0);
  }
}

void Minkmodel::add_vars_lifted_tree_compact() {
  if (_graph->_bags.size() == 0) {
    cerr << "ERROR: bags of cliques are not generated." << endl;
    exit(1);
  }
  // find the clique tree structure
  auto cliquetree = _graph->get_clique_tree_kruskal();
  // merging algorithm due to Nakata 2003.
  merge_cliques(_graph, cliquetree, cliquetree->nodes.front());

  // For each (i, j), we find all maximal cliques containing i,j.
  std::map<pair<int, int>, vector<int>> index_maps2_bag;
  vector<int> bags;
  for (int i = 0; i < _graph->nodes.size(); i++) {
    auto key = make_pair(i, i);
    bags.clear();
    for (auto j = 0; j < _graph->_bags.size(); j++) {
      auto bag = _graph->_bags.at(j);
      auto node = _graph->nodes.at(i);
      if (std::find(bag.begin(), bag.end(), node) != bag.end()) {
        bags.push_back(j);
      }
    }
    index_maps2_bag.insert(make_pair<>(key, bags));
  }

  vector<int> inter;
  for (auto a : _graph->_chordalextension->arcs) {
    auto i = a->src->ID, j = a->dest->ID;
    auto key = pair<int, int>(i, j);
    auto keyi = pair<int, int>(i, i);
    auto keyj = pair<int, int>(j, j);
    auto i_bags = index_maps2_bag[keyi];
    auto j_bags = index_maps2_bag[keyj];
    inter.clear();
    set_intersection(i_bags.begin(), i_bags.end(), j_bags.begin(), j_bags.end(),
                     back_inserter(inter));
    index_maps2_bag.insert(make_pair<>(key, inter));
  }

  vector<map<int, int>> node_mapsto_index;
  for (auto bag : _graph->_bags) {
    map<int, int> node_index;
    for (int i = 0; i < bag.size(); i++) {
      node_index.insert(make_pair<>(bag.at(i)->ID, i));
    }
    node_mapsto_index.push_back(node_index);
  }

  int clique_id = 0;
  vector<sdpvar<double>> X_subs;
  for (auto clique : _graph->_bags) {
    // adding SDP constraint for each clique
    sdpvar<double> X_sub("X_sub_" + to_string(clique_id));
    _model.add_var(X_sub ^ (clique.size()));
    X_subs.push_back(X_sub);
    clique_id++;
  }
  // Let common elements of two neighbouring blocks equal
  for (auto edge : cliquetree->arcs) {
    int n1 = edge->src->ID;
    int n2 = edge->dest->ID;
    for (int i = 0; i < edge->_intersection.size(); i++) {
      int n1_i = node_mapsto_index.at(n1)[edge->_intersection.at(i)->ID];
      int n2_i = node_mapsto_index.at(n2)[edge->_intersection.at(i)->ID];
      for (int j = i; j < edge->_intersection.size(); j++) {
        int n1_j = node_mapsto_index.at(n1)[edge->_intersection.at(j)->ID];
        int n2_j = node_mapsto_index.at(n2)[edge->_intersection.at(j)->ID];
        Constraint Linking_constraint("Linking_bags_" + to_string(n1) + "," +
                                      to_string(n2) + ",(" + to_string(i) +
                                      "," + to_string(j) + ")");
        Linking_constraint = X_subs[n1](n1_i, n1_j) - X_subs[n2](n2_i, n2_j);
        _model.add_constraint(Linking_constraint = 0);
      }
    }
  }
  int i = 0, j = 0;
  func_ obj;
  for (auto a : _graph->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;

    auto key = pair<int, int>(i, j);
    int index_bag = index_maps2_bag[key].front();
    int ni = node_mapsto_index.at(index_bag)[a->src->ID];
    int nj = node_mapsto_index.at(index_bag)[a->dest->ID];
    if (ni < nj)
      obj += a->weight * ((_K - 1) * X_subs[index_bag](ni, nj) + 1) / _K;
    else
      obj += a->weight * ((_K - 1) * X_subs[index_bag](nj, ni) + 1) / _K;
  }
  _model.set_objective(min(obj));

  for (i = 0; i < _graph->nodes.size(); i++) {
    Constraint diag("(" + to_string(i) + ")");
    auto key = pair<int, int>(i, i);
    int index_bag = index_maps2_bag[key].front();
    int ni = node_mapsto_index.at(index_bag)[_graph->nodes.at(i)->ID];
    diag = X_subs[index_bag](ni, ni) - 1;
    _model.add_constraint(diag = 0);
  }

  for (auto a : _graph->_chordalextension->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    auto key = pair<int, int>(i, j);
    int index_bag = index_maps2_bag[key].front();
    int ni = node_mapsto_index.at(index_bag)[a->src->ID];
    int nj = node_mapsto_index.at(index_bag)[a->dest->ID];

    Constraint bound("(" + to_string(i) + "," + to_string(j) + ")");
    if (ni < nj) {
      bound = X_subs[index_bag](ni, nj) + 1 / (_K - 1);
    } else {
      bound = X_subs[index_bag](nj, ni) + 1 / (_K - 1);
    }
    _model.add_constraint(bound >= 0);
  }
}

void Minkmodel::add_triangle() {
  auto n = _graph->nodes.size();
  auto zij = (*(var<bool>*)(_model.get_var("zij")));

  for (auto i = 0; i < n; i++)
    for (auto h = i + 1; h < n; h++)
      for (auto j = h + 1; j < n; j++) {
        Constraint Triangle1("Triangle1(" + to_string(i) + "," + to_string(h) +
                             "," + to_string(j) + ")");
        Triangle1 = zij(i, h) + zij(h, j) - zij(i, j);
        Constraint Triangle2("Triangle2(" + to_string(i) + "," + to_string(h) +
                             "," + to_string(j) + ")");
        Triangle2 = zij(i, h) + zij(i, j) - zij(h, j);
        Constraint Triangle3("Triangle3(" + to_string(i) + "," + to_string(h) +
                             "," + to_string(j) + ")");
        Triangle3 = zij(i, j) + zij(h, j) - zij(i, h);
        _model.add_constraint(Triangle1 <= 1);
        _model.add_constraint(Triangle2 <= 1);
        _model.add_constraint(Triangle3 <= 1);
      }
}

void Minkmodel::add_triangle_lifted() {
  auto n = _graph->nodes.size();
  auto X = (*(var<bool>*)(_model.get_var("X")));

  for (auto i = 0; i < n; i++)
    for (auto h = i + 1; h < n; h++)
      for (auto j = h + 1; j < n; j++) {
        Constraint Triangle1("Triangle1(" + to_string(i) + "," + to_string(h) +
                             "," + to_string(j) + ")");
        Triangle1 = X(i, h) + X(h, j) - X(i, j);
        Constraint Triangle2("Triangle2(" + to_string(i) + "," + to_string(h) +
                             "," + to_string(j) + ")");
        Triangle2 = X(i, h) + X(i, j) - X(h, j);
        Constraint Triangle3("Triangle3(" + to_string(i) + "," + to_string(h) +
                             "," + to_string(j) + ")");
        Triangle3 = X(i, j) + X(h, j) - X(i, h);
        _model.add_constraint(Triangle1 <= 1);
        _model.add_constraint(Triangle2 <= 1);
        _model.add_constraint(Triangle3 <= 1);
      }
}

void Minkmodel::add_clique() {
  auto n = _graph->nodes.size();
  auto zij = (*(var<bool>*)(_model.get_var("zij")));
  if (_K > 2) {
    for (auto i = 0; i < n; i++)
      for (auto h = i + 1; h < n; h++)
        for (auto j = h + 1; j < n; j++)
          for (auto l = j + 1; l < n; l++) {
            Constraint Clique("ZClique(" + to_string(i) + "," + to_string(h) +
                              "," + to_string(j) + ", " + to_string(l) + ")");
            Clique = zij(i, h) + zij(i, j) + zij(i, l) + zij(h, j) + zij(h, l) +
                     zij(j, l);
            _model.add_constraint(Clique >= 1);
          }
  } else {
    for (auto i = 0; i < n - 1; i++)
      for (auto h = i + 1; h < n; h++)
        for (auto j = h + 1; j < n; j++) {
          Constraint Clique("Clique(" + to_string(i) + "," + to_string(h) +
                            "," + to_string(j) + ")");
          Clique = zij(i, h) + zij(i, j) + zij(h, j);
          _model.add_constraint(Clique >= 1);
        }
  }
}

void Minkmodel::add_clique_lifted() {
  auto n = _graph->nodes.size();
  auto X = (*(var<bool>*)(_model.get_var("X")));
  if (_K > 2) {
    for (auto i = 0; i < n; i++)
      for (auto h = i + 1; h < n; h++)
        for (auto j = h + 1; j < n; j++)
          for (auto l = j + 1; l < n; l++) {
            Constraint Clique("Clique(" + to_string(i) + "," + to_string(h) +
                              "," + to_string(j) + ", " + to_string(l) + ")");
            Clique = X(i, h) + X(i, j) + X(i, l) + X(h, j) + X(h, l) + X(j, l);
            _model.add_constraint(Clique >= -0.5 * _K);
          }
  } else {
    for (auto i = 0; i < n - 1; i++)
      for (auto h = i + 1; h < n; h++)
        for (auto j = h + 1; j < n; j++) {
          Constraint Clique("Clique(" + to_string(i) + "," + to_string(h) +
                            "," + to_string(j) + ")");
          Clique = X(i, h) + X(h, j) + X(i, j);
          _model.add_constraint(Clique >= -0.5 * _K);
        }
  }
}

// construct a temporary container;
vector<unsigned> temp;

// we generate indices of clique constraints recursively.
// this is based on the idea that we firstly choose one index and then choose
// k-1 indices from the remaining indices.
void Minkmodel::nchoosek(int bag_id, int offset, int K) {
  // if K is 0, we have a K tuple index set at hand.
  if (K == 0) {
    std::string key;
    for (unsigned i = 0; i < temp.size(); ++i) {
      key += to_string(temp[i]);
      if (i < temp.size() - 1) {
        key += ",";
      }
    }

    auto iter = _cliqueid->find(key);
    if (iter == _cliqueid->end()) {
      _cliqueid->insert(make_pair<>(key, temp));
    }
    return;
  }

  // Otherwise we continue to generate the indices.
  auto bag = _graph->_bags.at(bag_id);
  for (int i = offset; i <= bag.size() - K; i++) {
    temp.push_back(bag[i]->ID);
    nchoosek(bag_id, i + 1, K - 1);
    temp.pop_back();
  }
}

void Minkmodel::cliquetree_decompose() {
  int i1, i2, i3;
  for (int i = 0; i < _graph->_bags.size(); i++) {
    auto bag = _graph->_bags.at(i);
    if (bag.size() < 3) {
      continue;
    }
    for (int j = 0; j < bag.size() - 2; j++)
      for (int h = j + 1; h < bag.size() - 1; h++)
        for (int l = h + 1; l < bag.size(); l++) {
          i1 = bag[j]->ID;  // zero index.
          i2 = bag[h]->ID;
          i3 = bag[l]->ID;
          // cout << "(i1, i2, i3) " << i1 << " " << i2 << " " << i3
          // << endl;
          if (_ids.count(make_tuple(i1, i2, i3)) == 0) {
            _ids.insert(make_tuple(i1, i2, i3));
          } else {
            continue;
          }
        }
    if (bag.size() > _K && _K > 1) {
      nchoosek(i, 0, _K + 1);
    }
  }
  cout << "size of triangle inequalties: " << _ids.size() << endl;
  cout << "size of clique inequalties: " << _cliqueid->size() << endl;
}

void Minkmodel::add_3Dcuts() {
  auto X = (*(var<bool>*)(_model.get_var("X")));
  int i1, i2, i3;
  for (auto it : _ids) {
    i1 = get<0>(it);
    i2 = get<1>(it);
    i3 = get<2>(it);
    Constraint SDP3("SDP3(" + to_string(i1) + "," + to_string(i2) + "," +
                    to_string(i3) + ")");
    SDP3 = -2 * X(i1, i2) * X(i2, i3) * X(i1, i3);
    SDP3 -= 1;
    SDP3 += power(X(i1, i2), 2);
    SDP3 += power(X(i1, i3), 2);
    SDP3 += power(X(i2, i3), 2);
    _model.add_constraint(SDP3);
  }
};

void Minkmodel::add_triangle_tree() {
  int i1, i2, i3;
  auto zij = (*(var<bool>*)(_model.get_var("zij")));
  for (auto it : _ids) {
    i1 = get<0>(it);
    i2 = get<1>(it);
    i3 = get<2>(it);
    // cout << "(i1, i2, i3): " << i1 << ", " << i2 << ", " << i3 << endl;
    Constraint Triangle1("ZTriangle1(" + to_string(i1) + "," + to_string(i2) +
                         "," + to_string(i3) + ")");
    Triangle1 = zij(i1, i2) + zij(i1, i3) - zij(i2, i3);
    Constraint Triangle2("ZTriangle2(" + to_string(i1) + "," + to_string(i2) +
                         "," + to_string(i3) + ")");
    Triangle2 = zij(i1, i3) + zij(i2, i3) - zij(i1, i2);
    Constraint Triangle3("ZTriangle3(" + to_string(i1) + "," + to_string(i2) +
                         "," + to_string(i3) + ")");
    Triangle3 = zij(i2, i3) + zij(i1, i2) - zij(i1, i3);
    _model.add_constraint(Triangle1 <= 1);
    _model.add_constraint(Triangle2 <= 1);
    _model.add_constraint(Triangle3 <= 1);
  }
}

void Minkmodel::add_triangle_lifted_tree() {
  int i1, i2, i3;
  auto X = (*(var<bool>*)(_model.get_var("X")));
  for (auto it : _ids) {
    i1 = get<0>(it);
    i2 = get<1>(it);
    i3 = get<2>(it);
    // cout << "(i1, i2, i3): " << i1 << ", " << i2 <<", " << i3 <<endl;
    Constraint Triangle1("Triangle1(" + to_string(i1) + "," + to_string(i2) +
                         "," + to_string(i3) + ")");
    Triangle1 = X(i1, i2) + X(i1, i3) - X(i2, i3);
    Constraint Triangle2("Triangle2(" + to_string(i1) + "," + to_string(i2) +
                         "," + to_string(i3) + ")");
    Triangle2 = X(i1, i3) + X(i2, i3) - X(i1, i2);
    Constraint Triangle3("Triangle3(" + to_string(i1) + "," + to_string(i2) +
                         "," + to_string(i3) + ")");
    Triangle3 = X(i2, i3) + X(i1, i2) - X(i1, i3);
    _model.add_constraint(Triangle1 <= 1);
    _model.add_constraint(Triangle2 <= 1);
    _model.add_constraint(Triangle3 <= 1);
  }
}

void Minkmodel::add_clique_tree() {
  //    int i1,i2,i3,i4;
  auto zij = (*(var<bool>*)(_model.get_var("zij")));
  if (_K > 2) {
    for (auto it : (*_cliqueid)) {
      auto key = it.first;
      auto value = it.second;
      Constraint Clique("ZClique[" + key + "]");
      for (int i = 0; i < value.size() - 1; i++) {
        auto id1 = value[i];
        for (int j = i + 1; j < value.size(); j++) {
          auto id2 = value[j];
          if (id1 <= id2)
            Clique += zij(id1, id2);
          else
            Clique += zij(id2, id1);
        }
      }
      _model.add_constraint(Clique >= 1);
      // Clique.print();
    }
  } else {
    for (auto it : _ids) {
      auto i1 = get<0>(it);
      auto i2 = get<1>(it);
      auto i3 = get<2>(it);
      // cout << "(i1, i2, i3): " << i1 << ", " << i2 << ", " << i3<<
      // endl;
      Constraint Clique("ZClique(" + to_string(i1) + "," + to_string(i2) + "," +
                        to_string(i3) + ")");
      Clique = zij(i1, i2) + zij(i1, i3) + zij(i2, i3);
      _model.add_constraint(Clique >= 1);
    }
  }
}

void Minkmodel::add_general_clique() {
  //  t=1, q=2,...,k-1
  auto zij = (*(var<bool>*)(_model.get_var("zij")));
  for (int i = 0; i < _graph->_bags.size(); i++) {
    auto bag = _graph->_bags.at(i);
    if (bag.size() > _K + 1) {
      Constraint GClique("General_Clique" + to_string(i));
      int t = floor(bag.size() / _K);
      int q = (bag.size() - t * _K);
      for (int j = 0; j < bag.size() - 1; j++) {
        int idj = bag[j]->ID;
        for (int k = j + 1; k < bag.size(); k++) {
          int idk = bag[k]->ID;
          if (idj < idk)
            GClique += zij(idj, idk);
          else
            GClique += zij(idk, idj);
        }
      }
      _model.add_constraint(GClique >= 0.5 * t * (t - 1) * (t - q) +
                                           0.5 * t * (t + 1) * q);
    }
  }
}

void Minkmodel::node_edge_formulation() {
  var<bool> x("x");
  _model.add_var(x ^ ((int)(_K * _graph->nodes.size())));
  var<bool> y("y");
  // the number of arcs in the chordal extension
  _model.add_var(y ^ (_graph->arcs.size()));

  func_ obj_node_edge;
  int i = 0, j = 0;
  for (auto a : _graph->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    if (i <= j)
      obj_node_edge += (a->weight) * y(i, j);
    else
      obj_node_edge += (a->weight) * y(j, i);
  }
  _model.set_objective(min(obj_node_edge));

  // add assignment constraints

  for (i = 0; i < _graph->nodes.size(); i++) {
    Constraint Assign("Assignment" + to_string(i));
    for (int c = 0; c < _K; c++) Assign += x(i, c);
    _model.add_constraint(Assign = 1);
  }

  // add consistency constraints
  for (int c = 0; c < _K; c++) {
    for (auto a : _graph->arcs) {
      i = (a->src)->ID;
      j = (a->dest)->ID;
      Constraint Consistency1("Consistency1[" + to_string(i) + "," +
                              to_string(j) + ", " + to_string(c) + "]");
      Constraint Consistency2("Consistency2[" + to_string(i) + "," +
                              to_string(j) + ", " + to_string(c) + "]");
      Constraint Consistency3("Consistency3[" + to_string(i) + "," +
                              to_string(j) + ", " + to_string(c) + "]");
      if (i <= j) {
        Consistency1 = x(i, c) + x(j, c) - y(i, j);
        Consistency2 = x(i, c) - x(j, c) + y(i, j);
        Consistency3 = -1 * x(i, c) + x(j, c) + y(i, j);
      } else {
        Consistency1 = x(i, c) + x(j, c) - y(j, i);
        Consistency2 = x(i, c) - x(j, c) + y(j, i);
        Consistency3 = x(j, c) - x(i, c) + y(j, i);
      }
      _model.add_constraint(Consistency1 <= 1);
      _model.add_constraint(Consistency2 <= 1);
      _model.add_constraint(Consistency3 <= 1);
    }
  }
}

int Minkmodel::solve(int output, bool relax) {
  solver s(_model, _solver);
  cout << "Running the relaxation model \n";
  s.run(output, relax);
  return 1;
}

void Minkmodel::construct_fsol() {
  // construct a feasible solution from z
  auto zij = (*(var<bool>*)(_model.get_var("zij")));
  param<int> sol("sol");
  int i = 0, j = 0;
  if (_type == MIP) {
    cout << "The MIP solution is: " << endl;
    for (i = 0; i < _graph->nodes.size() - 1; i++) {
      for (j = i + 1; j < _graph->nodes.size(); j++) {
        sol(i, j) = zij(i, j).getvalue();
        cout << sol(i, j).to_str() << endl;
      }
    }
  } else {
    cout << "The constructed solution is: " << endl;
    for (auto a : _graph->_chordalextension->arcs) {
      i = (a->src)->ID;
      j = (a->dest)->ID;
      if (i <= j) {
        sol(i, j) = zij(i, j).getvalue();
      } else {
        // generally, will never enter
        cerr << "something wrong with lables of chordal extension graph";
        exit(1);
        sol(j, i) = zij(j, i).getvalue();
      }
    }

    Node* n = nullptr;
    Node* nn = nullptr;
    Arc* arc_chordal = nullptr;
    bool allzeros = true;
    double temp = 0;

    for (i = 0; i < _graph->_chordalextension->nodes.size() - 1; i++) {
      n = (_graph->_chordalextension->get_node(to_string(i + 1)));
      for (j = i + 1; j < _graph->_chordalextension->nodes.size(); j++) {
        nn = _graph->_chordalextension->get_node(to_string(j + 1));
        // cout<< "(n, nn): " << n->ID << ", " << nn->ID << endl;
        if (n->is_connected(nn)) {
          cout << sol(i, j).to_str() << endl;
        } else {
          string name = to_string(_graph->_chordalextension->arcs.size() + 1);
          arc_chordal = new Arc(name);
          arc_chordal->id = _graph->_chordalextension->arcs.size();
          arc_chordal->src = n;
          arc_chordal->dest = nn;
          arc_chordal->weight = 0;
          arc_chordal->connect();
          // now find the maximal clique containing edge (i,j)

          // for neigbour i intersect neighbour j
          set<int> idi;
          set<int> idj;
          Debug("neighbours of " << i << ": ");

          for (auto a : n->branches) {
            idi.insert(a->neighbour(n)->ID);
            Debug(a->neighbour(n)->ID << ", ");
          }
          Debug(endl << "neighbours of " << j << ": ");

          for (auto a : nn->branches) {
            idj.insert(a->neighbour(nn)->ID);
            Debug(a->neighbour(nn)->ID << ", ");
          }
          Debug(endl);

          // intersection of idi and idj.
          std::vector<int> inter;
          set_intersection(idi.begin(), idi.end(), idj.begin(), idj.end(),
                           std::back_inserter(inter));

          // how to get values of zij
          if (inter.size() == 0) {
            sol(i, j) = 1;
          } else {
            allzeros = true;
            for (auto h : inter) {
              temp = 0;
              Debug("(i, j, h): " << i << " " << j << " " << h << endl);
              if (h <= i) {
                temp += sol(h, i).getvalue();
                // cout << sol(h,i).getvalue() << endl;
              } else {
                temp += sol(i, h).getvalue();
                // cout << sol(i, h).getvalue() << endl;
              }
              if (h <= j) {
                temp += sol(h, j).getvalue();
                // cout << sol(h,j).getvalue() << endl;
              } else {
                temp += sol(j, h).getvalue();
                // cout << sol(j,h).getvalue() << endl;
              }

              // cout << "temp: " << temp << endl;
              // DebugOn("xhi + xhj = " << temp);

              if (temp == 2) {
                sol(i, j) = 1;
                allzeros = false;
                break;
              } else if (temp == 1) {
                sol(i, j) = 0;
                allzeros = false;
                break;
              } else
                continue;
            }
            if (allzeros) {
              sol(i, j) = 1;
            }
          }
          cout << sol(i, j).to_str() << endl;
          (_graph->_chordalextension)->add_arc(arc_chordal);
        }
      }
    }
  }
}
