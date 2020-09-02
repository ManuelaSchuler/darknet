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

edge** get_edge_list(int T, int E){
    // TODO: This assumes a linear graph
    // + additional connection from first node to last node:
    // load_data:labels -> accuracy:labels
    int t,i;
    edge** edges = calloc(E, sizeof(struct edge));
    for (t = 0 ; t < T-1; ++t){
        edge* e = calloc(1, sizeof(struct edge));
        e->u = t;
        e->v = t+1;
        edges[t] = e;
    }
    edge* e = calloc(1, sizeof(struct edge));
    e->u = 0;
    e->v = T-1;
    edges[E-1] = e;
    for (i = 0 ; i < E; i++){
        edge *e = edges[i];
        printf("edge %i -> %i\n", e->u, e->v);
    }
    return edges;
}

double get_compute_costs(int i){
    double costs[] = {0.00014611584477390574, 0.0014940974900433343, 0.0001300606158239152, 2.7387530341769212e-05};
    return costs[i];
}

float get_memory_costs(int i){
    double costs[] = {0.01, 0.31, 0.01, 0.00};
    return costs[i];
}

#ifdef GUROBI_ILP
int solve_checkpoint_ilp(ilp_solver *ilp, float const budget){
    GRBenv   *env   = NULL;
    GRBmodel *model = NULL;
    int       error = 0;
    int       optimstatus;
    int       sol_count;
    double    objval;

    int t,i,idx,offset,counter;
    int T = ilp->r.cols;
    float ram = 1;
    double gcd = budget/ram;
    int E = T;
    struct edge **edges = get_edge_list(T, E);

    double    obj[T*T];
    char      vtype[T*T];
    char*     names[T*T];
    double    r_sol[T*T];
    double    s_sol[T*T];
    double    lb[T*T];
    double    ub[T*T];
    int       beg[T*T];
    int       ind[T*T];
    int       r_ind[T*T];
    int       s_ind[T*T];
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
            obj[idx] = get_compute_costs(i); // only dependent on node index!
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
        ind[i] = i;
        val[i] = 1.0;
        names[i][0] = 'U';
    }
    error = GRBaddvars(model, T*T, 0, NULL, NULL, NULL, NULL, lb, ub, vtype, names);
    if (error) goto QUIT;

    // add Constraints
    // R sum constraint == T (exactly one new evaluation per timestep)
    for (t = 0 ; t < T; ++t){ ind[t] = (t*T)+t; }
    GRBaddconstr(model, T, ind, val, GRB_EQUAL, T, NULL);
    // R area above diagonal == 0
    counter = 0;
    for (t = 0 ; t < T; ++t)
        for (i = t+1 ; i < T; ++i)
            r_ind[counter++] = (t*T) +i;
    GRBaddconstr(model, counter, r_ind, val, GRB_EQUAL, 0, NULL);
    // U range constraints (0, ..., gcd)
    offset = 2*(T*T)+(T*E);
    for (t = 0 ; t < T; ++t){
        for (i = 0 ; i < T; ++i){
            idx = (t*T) +i;
            ind[idx] = offset+idx;
            GRBaddconstr(model, 1, ind+idx, val+idx, GRB_GREATER_EQUAL, 0.0, NULL);
            GRBaddconstr(model, 1, ind+idx, val+idx, GRB_LESS_EQUAL, gcd, NULL);
        }
    }
    // S area above diagonal == 0
    counter = 0;
    offset = 1*(T*T);
    for (t = 0 ; t < T; ++t){
        for (i = t ; i < T; ++i){
            s_ind[counter++] = offset +(t*T)+i;
        }
    }
    GRBaddconstr(model, counter, s_ind, val, GRB_EQUAL, 0, NULL);

    // ensure all checkpoints are in memory
    offset = 1*(T*T);   // offset of S
    counter = 0;
    for (t = 0 ; t < T-1; ++t){
        for (i = 0 ; i < T; ++i){
            s_ind[0] = (t*T)+i;               // R[t,i]
            s_ind[1] = offset + (t*T)+i;      // S[t,i]
            s_ind[2] = offset + ((t+1)*T)+i;  // S[t+1,i]
            val[0] = -1.0;
            val[1] = -1.0;
            val[2] = +1.0;
            GRBaddconstr(model, 3, s_ind, val, GRB_LESS_EQUAL, 0, NULL);
        }
    }

    int u,v;
    for (i = 0 ; i < E; ++i){
        edge* e = edges[i];
        for (t = 0 ; t < T; ++t){
            s_ind[0] = (t*T)+e->v;           // R[t,v]
            s_ind[1] = (t*T)+e->u;           // R[t,u]
            s_ind[2] = offset + (t*T)+e->u;  // S[t,u]
            val[0] = +1.0;
            val[1] = -1.0;
            val[2] = -1.0;
            GRBaddconstr(model, 3, s_ind, val, GRB_LESS_EQUAL, 0, NULL);
        }
    }

    error = GRBoptimize(model);
    if (error) goto QUIT;
    error = GRBwrite(model, "checkpoint_ilp.lp");
    if (error) goto QUIT;
    // Capture solution information
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
    if (error) goto QUIT;
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if (error) goto QUIT;
    error = GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &sol_count);
    if (error) goto QUIT;
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, T*T, r_sol);
    if (error) goto QUIT;
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, T*T, T*T, s_sol);
    if (error) goto QUIT;
    printf("\nOptimization complete\n");
    if (optimstatus == GRB_OPTIMAL) {
        printf("Optimal objective: %.4e\n", objval);
        for (t = 0 ; t < T; ++t){
            for (i = 0 ; i < T; ++i){
                ilp->r.vals[t][i] = (int)r_sol[(t*T)+i];
                ilp->s.vals[t][i] = (int)s_sol[(t*T)+i];
            }
        }
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
    for (i = 0 ; i < T*T; ++i)
        free(names[i]);
    for (i = 0 ; i < E; ++i)
        free(edges[i]);
    free(edges);
    GRBfreemodel(model);
    GRBfreeenv(env);
    return 0;
}
#endif

void init_ilp_solver(ilp_solver* ilp, int n){
    ilp->r = make_matrix(n,n);
    ilp->s = make_matrix(n,n);
}

void print_ilp_matrix(matrix const m)
{
    int i, j, value;
    printf("%d X %d Matrix:\n",m.rows, m.cols);
    for(i = 0; i < m.rows; ++i){
        printf("|  ");
        for(j = 0; j < m.cols; ++j){
            value = (int)(m.vals[i][j]+0.5f);
            if (value == 1)
                printf("\u2588");
            else
                printf(" ");
        }
        printf(" |\n");
    }
}

void print_ilp_matrices(ilp_solver *ilp){
    print_ilp_matrix(ilp->r);
    print_ilp_matrix(ilp->s);
}
