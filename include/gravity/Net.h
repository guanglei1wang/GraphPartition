//  Net.h
//  Network
//  Copyright (c) All rights reserved
//

#ifndef Cycle_Basis_PF_Net_h
#define Cycle_Basis_PF_Net_h
#include <map>
#include <set>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <gravity/Node.h>
#include <gravity/Arc.h>
#include <gravity/Path.h>

/** remove methods are based on making object nullptr **/
class Net{
public:
    std::string _name;
    /** Set of nodes */
    std::vector<Node*> nodes;
    
    /** Set of arcs */
    std::vector<Arc*> arcs;
    
    /** Mapping the arcs to their source-destination */
    std::map<std::string, std::set<Arc*>*> lineID;
    
    /** Mapping the node name to its position in the vector, key = node id */
    std::map<std::string, Node*> nodeID;
    
    /** Vector of cycles forming a cycle basis */
    std::vector<Path*> cycle_basis;
    
    
    bool duplicate(int n1, int n2, int id1);
    
    /** Clone network */
    Net* _clone;
    Net* _chordalextension;
    
    
    /** Cloning */
    Net* clone();
    
    
    /* Constructors */
    Net();
    
    /* Destructors */
    ~Net();
   
    /** Tree decomposition bags */
    std::vector<std::vector<Node*>> _bags;
    
    /** Modifiers */

    void add_node(Node* n);
    bool add_arc(Arc* a);
    
    /* Accessors */
    
    Node* get_node(std::string id);
    /** returns the arc formed by node ids n1 and n2 */
    Arc* get_arc(Node* n1, Node* n2);
    
    bool has_arc(unsigned n1, unsigned n2) {
        return get_arc(n1,n2)!=nullptr;
    }
    /** returns the arc formed by node ids n1 and n2 */
    Arc* get_arc(int n1, int n2);
    
    void print_nodes();
    char* readline(FILE *input);
    void exit_input_error(int line_num);
    void readFile(std::string fname);
    void topology(std::string fname,bool complement);
    void readrudy(std::string fname);
    Net get_complement(std::string fname);//  (i, j): i<j \E
    
    /**  @brief Remove node and all incident arcs from the network
     @note Does not remove the incident arcs from the list of arcs in the network!
     @return the id of the node removed
     */
    std::string remove_end_node();
    
    /**  @brief Remove node and all incident arcs from the network
     @note Does not remove the incident arcs from the list of arcs in the network!
     @return the id of the node removed
     */
    void remove_arc(Arc* a);
    
    void remove_node(Node* n); 
    
    /** Compute the tree decomposition bags **/
    void get_tree_decomp_bags(bool print_bags = false);

    /** Compute the tree decomposition bags **/
    void get_cliquebags (bool print_clique = false);
  
    /** compute the minimum spanning tree **/ 
    std::vector<Arc*> minimal_spanning_tree_kruskal();
    Net* get_clique_tree_kruskal();
    
    /* a depth-first visit (from root) of a tree to lable predecessors*/
    void depth_first_vist(Node* root);
   // int test(std::string fname);
    int test();

    
};
#endif
