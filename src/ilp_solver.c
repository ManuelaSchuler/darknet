#include "ilp_solver.h"

#ifdef GLPK_ILP
int solve_checkpoint_ilp(){
    glp_prob *mip = glp_create_prob();
    glp_set_prob_name(mip, "checkpoint_mip");
    glp_set_obj_dir(mip, GLP_MIN);
    // TODO: Problem statement
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    int err = glp_intopt(mip, &parm);
    printf("err: %i", err);

    glp_delete_prob(mip);
    return 0;
}

int solve_example_ilp(){
/*
    https://gist.github.com/msakai/2450935
    Maximize
        obj: x1 + 2 x2 + 3 x3 + x4
    Subject To
        c1: - x1 + x2 + x3 + 10 x4 <= 20
        c2: x1 - 3 x2 + x3 <= 30
        c3: x2 - 3.5 x4 = 0
    Bounds
        0 <= x1 <= 40
        2 <= x4 <= 3
    General
        x4
    End
*/
    glp_prob *mip = glp_create_prob();
    glp_set_prob_name(mip, "example_ilp");
    glp_set_obj_dir(mip, GLP_MAX);

    glp_add_rows(mip, 3);
    glp_set_row_name(mip, 1, "c1");
    glp_set_row_bnds(mip, 1, GLP_UP, 0.0, 20.0);
    glp_set_row_name(mip, 2, "c2");
    glp_set_row_bnds(mip, 2, GLP_UP, 0.0, 30.0);
    glp_set_row_name(mip, 3, "c3");
    glp_set_row_bnds(mip, 3, GLP_FX, 0.0, 0);

    glp_add_cols(mip, 4);
    glp_set_col_name(mip, 1, "x1");
    glp_set_col_bnds(mip, 1, GLP_DB, 0.0, 40.0);
    glp_set_obj_coef(mip, 1, 1.0);
    glp_set_col_name(mip, 2, "x2");
    glp_set_col_bnds(mip, 2, GLP_LO, 0.0, 0.0);
    glp_set_obj_coef(mip, 2, 2.0);
    glp_set_col_name(mip, 3, "x3");
    glp_set_col_bnds(mip, 3, GLP_LO, 0.0, 0.0);
    glp_set_obj_coef(mip, 3, 3.0);
    glp_set_col_name(mip, 4, "x4");
    glp_set_col_bnds(mip, 4, GLP_DB, 2.0, 3.0);
    glp_set_obj_coef(mip, 4, 1.0);
    glp_set_col_kind(mip, 4, GLP_IV);

    int ia[1+9], ja[1+9];
    double ar[1+9];
    ia[1]=1,ja[1]=1,ar[1]=-1;   // a[1,1] = -1
    ia[2]=1,ja[2]=2,ar[2]=1;    // a[1,2] = 1
    ia[3]=1,ja[3]=3,ar[3]=1;    // a[1,3] = 1
    ia[4]=1,ja[4]=4,ar[4]=10;   // a[1,4] = 10
    ia[5]=2,ja[5]=1,ar[5]=1;    // a[2,1] = 1
    ia[6]=2,ja[6]=2,ar[6]=-3;   // a[2,2] = -3
    ia[7]=2,ja[7]=3,ar[7]=1;    // a[2,3] = 1
    ia[8]=3,ja[8]=2,ar[8]=1;    // a[3,2] = 1
    ia[9]=3,ja[9]=4,ar[9]=-3.5; // a[3,4] = -3.5
    glp_load_matrix(mip, 9, ia, ja, ar);

    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    int err = glp_intopt(mip, &parm);
    printf("err: %i", err);

    double z = glp_mip_obj_val(mip);
    double x1 = glp_mip_col_val(mip, 1);
    double x2 = glp_mip_col_val(mip, 2);
    double x3 = glp_mip_col_val(mip, 3);
    double x4 = glp_mip_col_val(mip, 4);
    printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g, x4 = %g\n", z, x1, x2, x3, x4);
    // z = 122.5; x1 = 40; x2 = 10.5; x3 = 19.5, x4 = 3

    glp_delete_prob(mip);
    return 0;
}
#endif

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
