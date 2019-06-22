//
//  model.hpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//

#ifndef model_hpp
#define model_hpp

#include <stdio.h>
#include <gravity/constraint.h>
#include <map>
#include <unordered_set>
#include <math.h>
#include <vector>
#include <thread>
#ifdef USE_IPOPT
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>
#endif
#ifdef USE_GUROBI
#include <gurobi_c++.h>
#endif
#ifdef USE_BONMIN
#include <coin/BonTMINLP.hpp>
#endif

using namespace std;

class Model {
    
protected:
    string                          _name;
    set<pair<size_t,size_t>>        _hess; /* A set representing pairs of variables linked in the hessian */
    vector<shared_ptr<func_>>       _functions;
    void add_var(param_* v);        //Add variables without reallocating memory
    void add_param(param_* v);      //Add variables without reallocating memory
    
public:
    size_t                          _nb_vars = 0;
    size_t                          _nb_params = 0;
    size_t                          _nb_cons = 0;
    size_t                          _nnz_g = 0; /* Number of non zeros in the Jacobian */
    size_t                          _nnz_h = 0; /* Number of non zeros in the Hessian */
    map<unsigned,param_*>           _params; /**< Sorted map pointing to all parameters contained in this model */
    map<unsigned,param_*>           _vars; /**< Sorted map pointing to all variables contained in this model. Note that a variable is a parameter with a bounds attribute. */
    map<string,param_*>             _params_name; /**< Sorted map pointing to all parameters contained in this model */
    map<string,param_*>             _vars_name; /**< Sorted map pointing to all variables contained in this model. Note that a variable is a parameter with a bounds attribute. */
    map<unsigned,Constraint*>        _cons; /**< Sorted map (increasing index) pointing to all constraints contained in this model */
    map<string,Constraint*>         _cons_name; /**< Sorted map (increasing index) pointing to all constraints contained in this model */
    map<unsigned, set<Constraint*>>        _v_in_cons; /**< Set of constraints where each variable appears, indexed by variable ids */
    map<pair<unsigned,unsigned>, map<pair<unsigned,unsigned>,set<pair<int,int>>>>            _hess_link; /**< Set of variables linked to one another in the hessian, indexed by pairs of variable ids, a pair contains the metavar id and the instanciated var id. The last set contains pairs of indices of constraints where both variables appear, the pair stores the metaconstraint id and the instanciated constraint id */
    func_                           _obj; /** Objective function */
    ObjectiveType                   _objt; /** Minimize or maximize */
    double                          _obj_val = 0;/** Objective function value */
    /** Constructor */
    //@{
    Model();
    Model(const string& name){
        _name = name;
    };
    //@}
    
    /* Destructor */
    ~Model();
    
    /* Accessors */
    
    size_t get_nb_vars() const;
    
    size_t get_nb_cons() const;
    
    size_t get_nb_nnz_g() const;
    
    size_t get_nb_nnz_h() const;
    
    param_* get_var(const string& vname) const;
    
    Constraint* get_constraint(const string& name) const;
    
    bool has_var(const param_& v) const{
        return (_vars.count(v.get_vec_id())!=0);
    };

    
    
    /* Modifiers */
    
    void set_x(const double* x); // Assign values to all variables based on array x.
    
    
    void add_var(param_& v); //Add variables by copying variable
    void del_var(const param_& v);
    
    
    void add_param(param_& v); //Add variables by copying variable
    void del_param(const param_& v);
    
    
    
    void add_constraint(const Constraint& c);
    
    void embed(func_& f);/**<  Transfer all variables and parameters to the model, useful for a centralized memory management. */
    
    void embed(expr& e);/**<  Transfer all variables and parameters to the model, useful for a centralized memory management. */
    
    void del_constraint(const Constraint& c);
    void set_objective(const func_& f, ObjectiveType t = minimize);
    void set_objective(pair<func_*, ObjectiveType> p);
    void set_objective_type(ObjectiveType);
    void init_indices();// Initialize the indices of all variables involved in the model
    void check_feasible(const double* x);
    void fill_in_maps();/*< Fill the _hess and _v_in_ maps to link variables with their constraints and compute the Jacobian & Hessian matrices */
    void fill_in_var_bounds(double* x_l ,double* x_u);
    void fill_in_var_init(double* x);
    void fill_in_cstr_bounds(double* g_l ,double* g_u);
    void fill_in_obj(const double* x , double& res, bool new_x);
    void fill_in_grad_obj(const double* x , double* res, bool new_x);
    void fill_in_cstr(const double* x , double* res, bool new_x);
    void fill_in_jac(const double* x , double* res, bool new_x);
    void fill_in_jac_nnz(int* iRow , int* jCol);
    void fill_in_hess(const double* x , double obj_factor, const double* lambda, double* res, bool new_x);
    void eval_funcs_parallel(const double* x , int start, int end);
    void fill_in_hess_nnz(int* iRow , int* jCol);
#ifdef USE_IPOPT
    void fill_in_var_linearity(Ipopt::TNLP::LinearityType* param_types);
    void fill_in_cstr_linearity(Ipopt::TNLP::LinearityType* const_types);
#endif
    
#ifdef USE_BONMIN
    void fill_in_param_types(Bonmin::TMINLP::VariableType* param_types);
#endif    
    
    void add_on_off(const Constraint& c, var<bool>& on);
    void add_on_off(var<>& v, var<bool>& on);
    
    void add_on_off_McCormick(string name, var<>& v, var<>& v1, var<>& v2, var<bool>& on);
    void add_McCormick(string name, var<>& v, var<>& v1, var<>& v2);

    
    friend std::vector<int> bounds(int parts, int mem);
    
    /* Operators */
    
    
    
    /* Output */
    void print_functions() const;
    void print() const;
    void print_solution() const;
    void print_constraints() const;
    
};


pair<func_*, ObjectiveType> max(const func_& f);
pair<func_*, ObjectiveType> min(const func_& f);


#endif /* model_hpp */
