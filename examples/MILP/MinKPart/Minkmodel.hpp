//
//  model.hpp
//  Gravity
//
//  Created by Guanglei Wang on 19/6/17.
//
//
#ifndef Minkmodel_hpp
#define Minkmodel_hpp

#include <stdio.h>
#include <stdlib.h>

#include <gravity/Net.h>
#include <gravity/model.h>
#include <gravity/solver.h>

typedef enum { MIP, MIP_tree, SDP, SDP_tree,SDP_tree_global, Node_edge } ModelType;

class Minkmodel {
   public:
    ModelType _type;
    Model _model;
    SolverType _solver;
    double _K;
    Net* _graph;

    set<tuple<int, int, int>>
        _ids;  // indices for generating triangle inequalities.
    shared_ptr<map<std::string, vector<unsigned>>>
        _cliqueid;  // indices for generating clique constraints where vector
                    // stores the list of indices and string stores the key.
    set<tuple<int, int, int, int>> _ids4;  //

    Minkmodel();
    Minkmodel(ModelType type, Net* graph, double K);
    Minkmodel(ModelType type, Net* graph, double K, SolverType solver);

    ~Minkmodel();
    void reset();
    void build();

    /** Variables */
    void add_vars_origin();
    void add_vars_origin_tree();
    void add_vars_lifted();
    void add_vars_lifted_tree();
    void add_vars_lifted_tree_compact();

    /** Constraints */
    // different formulations
    // void add_obj();  included in add_vars_*
    // void add_obj_lifted();
    void add_triangle();
    void add_clique();
    void add_clique_tree();
    void add_triangle_lifted();
    void add_triangle_tree();
    void add_triangle_lifted_tree();
    void add_clique_lifted();
    void add_clique_lifted_tree();
    void add_general_clique();
    void add_wheel();
    void add_bicycle();
    void add_3Dcuts();
    void cliquetree_decompose();  // generate ids of
    void nchoosek(int, int, int);  // this function generates indices for clique constraints.
    void node_edge_formulation();
    //  post root node relaxation
    bool check_eigenvalues();
    void add_eigcut();
    void construct_fsol();

    /** Presolve */
    // void propagate_bounds();
    /** Solve */
    int solve(int output, bool relax);
    void print();
};
#endif /* model_hpp */
