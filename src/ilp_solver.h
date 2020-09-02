#ifndef ILP_SOLVER_H
#define ILP_SOLVER_H
#include "darknet.h"

#ifndef GUROBI_ILP
#define GUROBI_ILP
#endif

#ifdef GUROBI_ILP
#include "gurobi_c.h"
int solve_checkpoint_ilp(ilp_solver *ilp, float const budget);
#endif
void init_ilp_solver(ilp_solver *ilp, int n);
void print_ilp_matrices(ilp_solver *ilp);
#endif

struct edge {
    int u;
    int v;
};
typedef struct edge edge;