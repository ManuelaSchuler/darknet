#ifndef ILP_SOLVER_H
#define ILP_SOLVER_H
#include "darknet.h"

void init_ilp_solver(ilp_solver *ilp, int n);
#ifdef GLPK_ILP
int solve_example_ilp();
int solve_checkpoint_ilp();
#endif
void print_ilp_matrix(matrix m);
void print_ilp_matrices(ilp_solver *ilp);
#endif

