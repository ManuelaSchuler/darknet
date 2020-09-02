#include "ilp_solver.h"

int getLengthOfInt(int i){
    return snprintf( NULL, 0, "%d", i );
}

int getLengthForVarName(int t, int i){
    // format: %c [ %d , %d ]
    // extra chars for:
    //       type [ , ]
    int extra_chars = 4;
    return getLengthOfInt(t) + getLengthOfInt(i) + extra_chars + 1;
}

float getCosts(int i){
    // TODO: return costs of node i
    return 1.0;
}
#ifdef GUROBI_ILP
int solve_checkpoint_ilp(ilp_solver *ilp, float const budget){
    GRBenv   *env   = NULL;
    GRBmodel *model = NULL;
    int       error = 0;
    int       optimstatus;
    double    objval;

    int t,i,idx;
    int T = ilp->r.cols;
    float ram = 4e6;
    double gcd = budget/ram;
    int E = T;

    double    obj[T*T];
    char      vtype[T*T];
    char*     names[T*T];
    double    sol[T*T];
    double    lb[T*T];
    double    ub[T*T];
    int       beg[T*T];
    int       ind[T*T];
    double    val[T*T];

    error = GRBemptyenv(&env);
    if (error) goto QUIT;
    error = GRBsetstrparam(env, "LogFile", "checkpoint_ilp.log");
    if (error) goto QUIT;
    error = GRBstartenv(env);
    if (error) goto QUIT;
    error = GRBnewmodel(env, &model, "checkpoint_ilp", 0, NULL, NULL, NULL, NULL, NULL);
    if (error) goto QUIT;
    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
    if (error) goto QUIT;

    /* R */
    for (t = 0 ; t < T; ++t){
        for (i = 0 ; i < T; ++i){
            idx = (t*T) +i;
            vtype[idx] = GRB_BINARY;
            int l = getLengthForVarName(t, i);
            names[idx] = malloc(l);
            snprintf( names[idx], l, "%c[%d,%d]", 'R', t, i );
            obj[idx] = getCosts(i);
        }
    }
    error = GRBaddvars(model, T*T, 0, NULL, NULL, NULL, obj, NULL, NULL, vtype, names);
    if (error) goto QUIT;
    /* S */
    for (i = 0 ; i < T*T; ++i){
        names[i][0] = 'S';
    }
    error = GRBaddvars(model, T*T, 0, NULL, NULL, NULL, NULL, NULL, NULL, vtype, names);
    if (error) goto QUIT;
    /* F */
    for (i = 0 ; i < T*E; ++i){
        names[i][0] = 'F';
    }
    error = GRBaddvars(model, T*E, 0, NULL, NULL, NULL, NULL, NULL, NULL, vtype, names);
    if (error) goto QUIT;
    /* U */
    for (i = 0 ; i < T*T; ++i){
        vtype[i] = GRB_CONTINUOUS;
        lb[i] = 0.0;
        ub[i] = gcd;
        names[i][0] = 'U';
    }
    error = GRBaddvars(model, T*T, 0, NULL, NULL, NULL, NULL, lb, ub, vtype, names);
    if (error) goto QUIT;



    error = GRBoptimize(model);
    if (error) goto QUIT;
    error = GRBwrite(model, "checkpoint_ilp.lp");
    if (error) goto QUIT;
    // Capture solution information
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
    if (error) goto QUIT;
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if (error) goto QUIT;
    printf("\nOptimization complete\n");
    if (optimstatus == GRB_OPTIMAL) {
        printf("Optimal objective: %.4e\n", objval);
        printf("\n");
    } else if (optimstatus == GRB_INF_OR_UNBD) {
        printf("Model is infeasible or unbounded\n");
    } else {
        printf("Optimization was stopped early\n");
    }
QUIT:
    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }
    for (i = 0 ; i < T*T; ++i){
        free(names[i]);
    }
    GRBfreemodel(model);
    GRBfreeenv(env);
    return 0;
}
#endif

void init_ilp_solver(ilp_solver* ilp, int n){
    ilp->r = make_matrix(n,n);
    ilp->s = make_matrix(n,n);
    fill_matrix_eye(ilp->r);
    fill_matrix(ilp->s, 0.0);
}

void print_int_matrix(matrix const m)
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
