//
//  Node.cpp
//  Cycle_Basis_PF
//  adapt from power tools
//  Created by Sumiran on 17/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#include <gravity/Arc.h>
#include <gravity/Node.h>
#include <algorithm>
#include <limits.h>
#include <iostream>

using namespace std;

Node::Node(){};

Node::Node(int id) : _name(""), ID(id), fill_in(0){};
Node::Node(string name, int id) : _name(name), ID(id), fill_in(0){};
Node::~Node(){};

Node* Node::clone() {
  Node* copy = new Node();
  copy->ID = ID;
  copy->_name = _name;
  copy->fill_in = 0;
  return copy;
};

/*
 @brief Adds a to the list of incident arcs
 */
void Node::addArc(Arc* a) { branches.push_back(a); }

/*
 @brief Find and remove incident arc from list of branches
 @return 0 if a was found and removed, -1 oterwise
 */
int Node::removeArc(Arc* a) {
  vector<Arc*>::iterator it = branches.begin();
  while (it != branches.end()) {
    if ((*it) == a) {
      it = branches.erase(it);
      return 0;
    }
    it++;
  }
  return -1;
}

bool Node::is_connected(Node* n) {
  for (auto a : branches) {
    if (n->ID == a->neighbour(this)->ID) {
      return true;
    }
  }
  for (auto a : n->branches) {
    if (ID == a->neighbour(n)->ID) {
      return true;
    }
  }
  return false;
}

void Node::update_fill_in(Node* n) {
  Node* nn = nullptr;
  for (auto a : branches) {
    nn = a->neighbour(this);  // this node.
    // if nn is null
    if (nn->ID == n->ID) {
      continue;  // self connect
    }
    if (!n->is_connected(nn)) {
      fill_in++;  // if this node is connected to node.
    }
  }
}

std::vector<Node*> Node::get_neighbours() {
  vector<Node*> res;
  set<Node*> temp;
  for (auto a : branches) {
    // if(a->_dest->_id=_id && std::find(res.begin(),res.end(), a->_src)==
    // res.end()){
    if (a->dest->ID == ID) {
      temp.insert(a->src);
    }

    // if(a->_src->_id==_id && std::find(res.begin(),res.end(), a->_dest)==
    // res.end() ){ if(a->_src->_id==_id && std::find(res.begin(),res.end(),
    // a->_dest)== res.end() ){
    if (a->src->ID == ID) {
      temp.insert(a->dest);
    }
  }
  res.assign(temp.begin(), temp.end());
  std::sort(res.begin(), res.end(),
            [](const Node* a, const Node* b) -> bool { return a->ID < b->ID; });
  return res;
}
