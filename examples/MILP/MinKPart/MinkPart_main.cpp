//
//  MinKpartition.cpp
//  Gravity
//
//  Created by Guanglei Wang on 13/6/17.
//
//
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include "Minkmodel.hpp"

using namespace std;
#define EPS 0.00001
#define DebugOn(x) cout << x
#define DebugOff(x)
//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time() {
  LARGE_INTEGER time, freq;
  if (!QueryPerformanceFrequency(&freq)) {
    //  Handle error
    return 0;
  }
  if (!QueryPerformanceCounter(&time)) {
    //  Handle error
    return 0;
  }
  return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time() {
  FILETIME a, b, c, d;
  if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0) {
    //  Returns total user time.
    //  Can be tweaked to include kernel times as well.
    return (double)(d.dwLowDateTime |
                    ((unsigned long long)d.dwHighDateTime << 32)) *
           0.0000001;
  } else {
    //  Handle error
    return 0;
  }
}

//  Posix/Linux
#else
#include <sys/time.h>
#include <time.h>
double get_wall_time() {
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time() { return (double)clock() / CLOCKS_PER_SEC; }
#endif

#ifdef USE_MOSEK
double mosek_reduce(const Net* _graph, const double _K) {
  mosek::fusion::Model::t M = new mosek::fusion::Model("mink_reduce");
  auto _M = monty::finally([&]() { M->dispose(); });
  M->setLogHandler(
      [=](const std::string& msg) { std::cout << msg << std::flush; });

  // delare variables Y and restrict it in PSD cone of dimension graph->nodes.
  // int n = _graph->nodes.size();
  mosek::fusion::Variable::t Y = M->variable(
      "Y", _graph->nodes.size(), mosek::fusion::Domain::inPSDCone());

  int i = 0, j = 0;
  M->constraint(Y->diag(), mosek::fusion::Domain::equalsTo(1.0));
  for (auto a : _graph->_chordalextension->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    M->constraint(Y->index(i, j),
                  mosek::fusion::Domain::greaterThan(-1 / (_K - 1)));
  }

  monty::rc_ptr<::mosek::fusion::Expression> expr =
      mosek::fusion::Expr::constTerm(0.0);
  for (auto a : _graph->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    // note that i never equals to j.
    if (i < j) {
      expr = mosek::fusion::Expr::add(
          expr,
          mosek::fusion::Expr::mul(a->weight * (_K - 1) / _K, Y->index(i, j)));
    } else {
      expr = mosek::fusion::Expr::add(
          expr,
          mosek::fusion::Expr::mul(a->weight * (_K - 1) / _K, Y->index(j, i)));
    }
    // constant part wij/k
    expr = mosek::fusion::Expr::add(expr, a->weight / _K);
  }

  M->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
  M->writeTask("bigcone.task.gz");
  M->solve();
  std::cout << "Cost = " << M->primalObjValue() << std::endl;
  return M->primalObjValue();
}

// convert indices of the upper triangular part of a matrix to indices of a
// vector.
int vector_index(int row, int col, size_t n) {
  int start = (2 * n - row + 1) * row / 2;  // sum of the first (row-1) rows
  assert(col >= row);
  return start + col - row;
}

/** chordal decomposition with additional variables */
double mosek_reduce2(Net* _graph, double _K) {
  mosek::fusion::Model::t M = new mosek::fusion::Model("mink_reduce2");
  auto _M = monty::finally([&]() { M->dispose(); });
  M->setLogHandler(
      [=](const std::string& msg) { std::cout << msg << std::flush; });
  mosek::fusion::Variable::t Y = M->variable(
      "", _graph->nodes.size() + _graph->_chordalextension->arcs.size(),
      mosek::fusion::Domain::unbounded());

  int i = 0, j = 0;
  for (int i = 0; i < _graph->nodes.size(); i++) {
    M->constraint(Y->index(i), mosek::fusion::Domain::equalsTo(1.0));
  }

  for (auto a : _graph->_chordalextension->arcs) {
    int index = a->id + _graph->nodes.size();
    M->constraint("", Y->index(index),
                  mosek::fusion::Domain::greaterThan(-1 / (_K - 1)));
  }

  for (auto clique : _graph->_bags) {
    // for each clique, select the trucated submatrix associated with the
    // indices in the clique...
    monty::rc_ptr<::mosek::fusion::Expression> expr =
        mosek::fusion::Expr::constTerm(0.0);
    mosek::fusion::Variable::t sub =
        M->variable("sub", clique.size(), mosek::fusion::Domain::inPSDCone());
    // linking sub and Y
    for (int i = 0; i < clique.size(); i++)
      for (int j = i; j < clique.size(); j++) {
        int index1 = clique[i]->ID + 1;
        int index2 = clique[j]->ID + 1;
        if (index1 != index2) {
          auto arc = _graph->_chordalextension->get_arc(index1, index2);
          expr = mosek::fusion::Expr::add(
              expr, Y->index(_graph->nodes.size() + arc->id));
        } else {
          expr = mosek::fusion::Expr::add(expr, Y->index(clique[i]->ID));
        }
        expr = mosek::fusion::Expr::add(
            expr, mosek::fusion::Expr::mul(-1, sub->index(i, j)));
        M->constraint("", expr, mosek::fusion::Domain::equalsTo(0.0));
      }
  }

  monty::rc_ptr<::mosek::fusion::Expression> expr =
      mosek::fusion::Expr::constTerm(0.0);
  for (auto a : _graph->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    int index = _graph->nodes.size() + a->id;
    expr = mosek::fusion::Expr::add(
        expr,
        mosek::fusion::Expr::mul(a->weight * (_K - 1) / _K, Y->index(index)));
    expr = mosek::fusion::Expr::add(expr, a->weight / _K);
  }

  M->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
  M->solve();
  std::cout << "Cost = " << M->primalObjValue() << std::endl;
  return M->primalObjValue();
}

shared_ptr<monty::ndarray<int, 2>> bag2index(const vector<Node*> clique) {
  size_t N = clique.size() * clique.size();
  auto index = make_shared<monty::ndarray<int, 2>>(monty::shape(N, 2));

  for (int i = 0; i < N; i++) {
    int row = i / clique.size();
    int col = i % clique.size();
    (*index)(i, 0) = clique[row]->ID;  // row index
    (*index)(i, 1) = clique[col]->ID;  // column index
    // cout << "(rowID, colID): " << clique[row]->ID << "," <<
    // clique[col]->ID << endl;
  }
  return index;
}

/** chordal decompositon without introducing new variables */
double mosek_reduce3(Net* _graph, double _K) {
  mosek::fusion::Model::t M = new mosek::fusion::Model("mink_reduce3");
  auto _M = monty::finally([&]() { M->dispose(); });
  M->setLogHandler(
      [=](const std::string& msg) { std::cout << msg << std::flush; });
  int n = _graph->nodes.size();
  // define an nxn dimentional symmetric variable;
  auto X = M->variable(
      n, mosek::fusion::Domain::symmetric(mosek::fusion::Domain::unbounded()));
  // cout << X->toString() << endl;

  M->constraint(X->diag(), mosek::fusion::Domain::equalsTo(1.0));
  //
  for (auto clique : _graph->_bags) {
    shared_ptr<monty::ndarray<int, 2>> index = bag2index(clique);
    auto c = mosek::fusion::Var::reshape(X->pick(index), clique.size(),
                                         clique.size());
    auto con = M->constraint(c, mosek::fusion::Domain::inPSDCone());
  }

  int i = 0, j = 0;
  for (auto a : _graph->_chordalextension->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    if (i < j)
      M->constraint("", X->index(i, j),
                    mosek::fusion::Domain::greaterThan(-1 / (_K - 1)));
    else
      M->constraint("", X->index(j, i),
                    mosek::fusion::Domain::greaterThan(-1 / (_K - 1)));
  }
  monty::rc_ptr<::mosek::fusion::Expression> expr =
      mosek::fusion::Expr::constTerm(0.0);
  for (int l = 0; l < _graph->arcs.size(); l++) {
    auto a = _graph->arcs[l];
    i = (a->src)->ID;
    j = (a->dest)->ID;
    if (i < j) {
      expr = mosek::fusion::Expr::add(
          expr,
          mosek::fusion::Expr::mul(a->weight * (_K - 1) / _K, X->index(i, j)));
    } else {
      expr = mosek::fusion::Expr::add(
          expr,
          mosek::fusion::Expr::mul(a->weight * (_K - 1) / _K, X->index(j, i)));
    }
    // constant part wij/k
    expr = mosek::fusion::Expr::add(expr, a->weight / _K);
  }
  M->objective("objective", mosek::fusion::ObjectiveSense::Minimize, expr);
  M->writeTask("smallcone.task.gz");
  M->solve();
  std::cout << "Cost = " << M->primalObjValue() << std::endl;
  return M->primalObjValue();
}

double mosekcode(Net* _graph, double _K) {
  mosek::fusion::Model::t M = new mosek::fusion::Model("mink");
  auto _M = monty::finally([&]() { M->dispose(); });
  M->setLogHandler(
      [=](const std::string& msg) { std::cout << msg << std::flush; });

  mosek::fusion::Variable::t Y =
      M->variable("Y", mosek::fusion::Domain::inPSDCone(_graph->nodes.size()));

  int i = 0, j = 0;
  M->constraint(Y->diag(), mosek::fusion::Domain::equalsTo(1.0));

  for (i = 0; i < _graph->nodes.size(); i++)
    for (j = i + 1; j < _graph->nodes.size(); j++) {
      M->constraint("", Y->index(i, j),
                    mosek::fusion::Domain::greaterThan(-1 / (_K - 1)));
    }

  monty::rc_ptr<::mosek::fusion::Expression> expr =
      mosek::fusion::Expr::constTerm(0.0);
  // expr is a pointer to the Expression.
  for (auto a : _graph->arcs) {
    i = (a->src)->ID;
    j = (a->dest)->ID;
    if (i <= j) {
      expr = mosek::fusion::Expr::add(
          expr,
          mosek::fusion::Expr::mul(a->weight * (_K - 1) / _K, Y->index(i, j)));
    } else {
      expr = mosek::fusion::Expr::add(
          expr,
          mosek::fusion::Expr::mul(a->weight * (_K - 1) / _K, Y->index(j, i)));
    }
    expr = mosek::fusion::Expr::add(expr, a->weight / _K);
  }

  M->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
  M->solve();

  std::cout << "Cost = " << M->primalObjValue() << std::endl;
  return M->primalObjValue();
}
#endif

int main(int argc, const char* argv[]) {
  double k = 3;
  ModelType mt = MIP;
  bool relax = false;
  int output = 0;
  SolverType solver = cplex;
  const char* fname;
  const char* type;
  const char* relaxation;

  if (argc > 3) {
    fname = argv[1];       // problem instance data
    k = atoi(argv[2]);     // The min_k-partition parameter k.
    type = argv[3];        // Model type parameter
    relaxation = argv[4];  // true: LP relaxation, false: IP

    cout << "Model type: " << type << endl;

    if (strcmp(type, "MIP") == 0) {
      mt = MIP;
      solver = cplex;
      if (k > 3) {
        cout << "WARNING: we suggest k be an integer <= 3 for the MIP "
                "model."
             << endl;
        exit(1);
      }
    } else if (strcmp(type, "SDP") == 0) {
      mt = SDP;
      solver = mosek_;
    } else if (strcmp(type, "MIP_tree") == 0) {
      mt = MIP_tree;
      solver = cplex;
    } else if (strcmp(type, "SDP_tree") == 0) {
      mt = SDP_tree;
      solver = mosek_;
    } else if (strcmp(type, "SDP_tree_global") == 0) {
      mt = SDP_tree_global;
      solver = mosek_;
    } else if (strcmp(type, "Node_edge") == 0) {
      mt = Node_edge;
      solver = cplex;
    } else {
      cerr << "invalid input on type selection, please enter: 'MIP', "
              "'SDP', 'MIP_tree', 'SDP_tree' or 'Node_edge' "
           << endl;
      exit(1);
    }

    if (strcmp(relaxation, "true") == 0)
      relax = true;
    else {
      relax = false;
    }
  } else {
    fname = "../../data_sets/Minkcut/random10_100.txt";
    // fname = "../../data_sets/Minkcut/band100_3.txt";
    // fname = "../../data_sets/Minkcut/spinglass2g_1111.txt";
    // fname = "../../data_sets/Minkcut/toy_chopra.txt";
    // fname = "../../data_sets/Minkcut/toy.txt";
    k = 3;
    relax = true;
    mt = SDP_tree;
    mt = SDP_tree_global;
    solver= cplex;
    solver = mosek_;
  }
  // Graph reading and clique-tree generation.
  Net* graph = new Net();
  graph->readrudy(fname);
  // get cliques if using new models.
  if (mt == SDP_tree || mt == MIP_tree || mt== SDP_tree_global) {
    graph->get_cliquebags();
  }
  // Buid model
  //mosek_reduce(graph,k); // using mosek native language.
  Minkmodel mymodel(mt, graph, k, solver);
  mymodel.build();

  double wall0 = get_wall_time();
  double cpu0 = get_cpu_time();
  mymodel.solve(output, relax);

  double wall1 = get_wall_time();
  double cpu1 = get_cpu_time();

  cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
  cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";
  // get the optimal solution as follows.
  ofstream outfile("MkP_result.txt", ios_base::app);
  if (!outfile) {
    cout << "Oops! Uable to save session data! \n";
    exit(1);
  } else {
    outfile << k << "," << graph->nodes.size() << "," << (cpu1 - cpu0) << ", "
            << mymodel._model._obj_val << endl;
  }  
  return    0;
}
