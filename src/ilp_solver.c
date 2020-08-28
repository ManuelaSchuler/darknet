#include "ilp_solver.h"

void init_ilp_solver(ilp_solver* ilp, int n){
    ilp->r = make_matrix(n,n);
    ilp->s = make_matrix(n,n);
    fill_matrix_eye(ilp->r);
    fill_matrix(ilp->s, 0.0);
}

void print_int_matrix(matrix m)
{
    int i, j;
    printf("%d X %d Matrix:\n",m.rows, m.cols);
    for(i = 0; i < m.rows; ++i){
        printf("|  ");
        for(j = 0; j < m.cols; ++j){
            printf("%i ", (int)(m.vals[i][j]+0.5f));
        }
        printf(" |\n");
    }
}

void print_ilp_matrices(ilp_solver *ilp){
    print_int_matrix(ilp->r);
    print_int_matrix(ilp->s);
}
