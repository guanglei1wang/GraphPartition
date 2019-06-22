//
//  Node.h
//

#ifndef Node_h
#define Node_h
#include <vector>
#include <string>
#include <set>

class Arc;

class Node{
    
public:
    
    std::string _name;
    int ID;
    std::vector<Arc*> branches; // do we really need this? as a node's branch is dependent on graph.
    Node* predecessor; // for shortest path algorithm
    /* the number of edges needed to make the subgraph formed by adjacent nodes a clique */
  int fill_in = 0;
  //int fill_in;

    
    // constructions
    Node();
    Node(int id);
    Node(std::string name, int id);
    ~Node();
    Node* clone();
    
    void addArc(Arc* a);
    int removeArc(Arc* a);
    void update_fill_in(Node* n);
    
    /*
     @brief Returns true if n is an adjacent node.
     */
    bool is_connected(Node* n);
    
    /* return its neighbours */
    std::vector<Node*> get_neighbours();
    
};

#endif
