//#include "mex.h"
#include <bits/stdc++.h>
#include <math.h>
#include "Eigen/Core"
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "OptimLib/optim.hpp"
#include "constraints.h"
#include "testdata.h"

static void printvec(const Eigen::VectorXd& v) {
    for (int i = 0; i < v.size(); i++) {
        printf("%.2lf ", v[i]);
    }
    printf("\n");
}
/*
    TODO: is it a good idea to take X as a flat array?
          certainly makes translation easier
*/

// calculates value of objective function at X
static double objfun(const Eigen::VectorXd& X, DistConstraints& E, DistConstraints& L,
                     DistConstraints& U, double *w, double *f) {
    
    double ret = 0.0, temp;

    // dimensions and size of X
    int d = 3, n = X.size()/d;

    // number of bounds etc.
    int ne = E.size(), nl = L.size(), nu = U.size();

    // weights
    double we = w[0], wl = w[1], wu = w[2], wr = w[3];

    // idk what these are
    double f_hb  = f[0], f_tau = f[1], f_tal = f[2], f_vdw = f[3];

    /* equality constraints */
    for (int i = 0 ; i < ne ; i++) {
        int ti = (E[i].i - 1) * d, tj = (E[i].j - 1) * d;
        double target = E[i].val, dist = 0;

        for (int j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        dist = sqrt(dist);
        temp = target - dist;
        ret += we * temp * temp;
    }
    
    
    /* lower bounds */
    for (int i = 0 ; i < nl ; i++) {
        int ti = (L[i].i - 1) * d, tj = (L[i].j - 1) * d, type = L[i].type;
        double target = L[i].val, dist = 0;

        for (int j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        dist = sqrt(dist);
        temp = target - dist;
        if (temp > 0){
            if (type == 1){
                ret += f_vdw * wl * temp * temp;
            }else if (type == -2){
                ret += f_tal*wl * temp * temp;
            }
        }
    }
    
    
    /* upper bounds */
    for (int i = 0 ; i < nu ; i++) {
        int ti = (U[i].i - 1) * d, tj = (U[i].j - 1) * d, type = U[i].type;
        double target = U[i].val, dist = 0;

        for (int j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        dist = sqrt(dist);
        temp = target - dist;
        if (temp < 0){
            if (type >= 0){
                ret += wu * temp * temp;
            }else if (type == -1){
                ret += f_hb*wu * temp * temp;
            }else if (type == -2){
                ret += f_tau*wu * temp * temp;
            }
        }
    }
    
    
    /* regularization */
    for (int i = 0 ; i < n ; i++) {
        int ti = i*d;
        double dist = 0;
        for (int j = 0; j < d ; j++) {
            temp = X[ti+j];
            dist += temp * temp;
        }
        ret += wr * dist;
    }

    return ret;
}

// calculates value of gradient of objective function at X
static void gradfun(const Eigen::VectorXd& X, DistConstraints& E, DistConstraints& L,
                    DistConstraints& U, double* w, double* f,
                    Eigen::VectorXd* G) {
    
    double temp;

    // dimensions and size of X
    int d = 3, n = X.size()/d;

    // number of bounds etc.
    int ne = E.size(), nl = L.size(), nu = U.size();

    // weights
    double we = 2*w[0], wl = 2*w[1], wu = 2*w[2], wr = 2*w[3];

    // idk what these are
    double f_hb  = f[0], f_tau = f[1], f_tal = f[2], f_vdw = f[3];

    /*double *X, *E, *L, *U, *w, *G, *f;
    int i, j, ti, tj, type;
    double we, wl, wu, wr, dist, sdist, temp, target, tempG;
    double f_hb, f_tau, f_tal, f_vdw;

    int ne, nl, nu, n , d;
    
    X  = mxGetPr(prhs[0]);
    d  = mxGetM(prhs[0]);
    n  = mxGetN(prhs[0]);
    
    
    E = mxGetPr(prhs[1]);
    ne = mxGetM(prhs[1]);
    
    
    L = mxGetPr(prhs[2]);
    nl = mxGetM(prhs[2]);
    
    U = mxGetPr(prhs[3]);
    nu = mxGetM(prhs[3]);
    
    w  = mxGetPr(prhs[4]);
    we = 2*w[0];
    wl = 2*w[1];
    wu = 2*w[2];
    wr = 2*w[3];
    
    f = mxGetPr(prhs[5]);
    
    f_hb  = f[0];
    f_tau = f[1];
    f_tal = f[2];
    f_vdw = f[3];*/
    
    /* mxArray for output */
    //plhs[0] = mxCreateDoubleMatrix(d, n, mxREAL);
    /* pointer to output data */
    //G = mxGetPr(plhs[0]);
    
    
    /* equality constraints */
    for (int i = 0 ; i < ne ; i++) {
        int ti = (E[i].i - 1) * d, tj = (E[i].j - 1) * d;
        double target = E[i].val, dist = 0;

        for (int j = 0; j < d; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        double sdist = sqrt(dist);
        double tempG = (1 - target/sdist);
        for (int j = 0; j < d ; j++) {
            (*G)[ti+j] += we*(X[ti+j] - X[tj+j])*tempG;
            (*G)[tj+j] -= we*(X[ti+j] - X[tj+j])*tempG;
        }
    }
    
    /* lower bounds */
    for (int i = 0 ; i < nl ; i++) {
        int ti = (L[i].i - 1) * d, tj = (L[i].j - 1) * d, type = L[i].type;
        double target = L[i].val, dist   = 0;

        for (int j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        double sdist = sqrt(dist);
        double tempG = (1 - target/sdist);
        if (sdist < target) {
            for (int j = 0; j < d ; j++) {
                if (type == 1){
                    (*G)[ti+j] += f_vdw*wl*(X[ti+j] - X[tj+j])*tempG;
                    (*G)[tj+j] -= f_vdw*wl*(X[ti+j] - X[tj+j])*tempG;
                } else if (type == -2){
                    (*G)[ti+j] += f_tal*wl*(X[ti+j] - X[tj+j])*tempG;
                    (*G)[tj+j] -= f_tal*wl*(X[ti+j] - X[tj+j])*tempG;
                }else{
                    printf("error: unknown lower bound!\n");
                    exit(-1);
                }
            }
        }
    }
    
    
    /* upper bounds */
    for (int i = 0 ; i < nu ; i++) {
        int ti = (U[i].i - 1) * d, tj = (U[i].j - 1) * d, type = U[i].type;
        double target = U[i].val, dist = 0;

        for (int j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        double sdist = sqrt(dist);
        double tempG = (1 - target/sdist);
        if (sdist > target) {
            for (int j = 0; j < d ; j++) {
                if (type >= 0){
                    (*G)[ti+j] += wu*(X[ti+j] - X[tj+j])*tempG;
                    (*G)[tj+j] -= wu*(X[ti+j] - X[tj+j])*tempG;
                } else if (type == -1){
                    (*G)[ti+j] += f_hb*wu*(X[ti+j] - X[tj+j])*tempG;
                    (*G)[tj+j] -= f_hb*wu*(X[ti+j] - X[tj+j])*tempG;
                } else if (type == -2){
                    (*G)[ti+j] += f_tau*wu*(X[ti+j] - X[tj+j])*tempG;
                    (*G)[tj+j] -= f_tau*wu*(X[ti+j] - X[tj+j])*tempG;
                }else{
                    printf("error: unknown upper bound!\n");
                    exit(-1);
                }
            }
        }
    }
    
    
    /* regularization */
    for (int i = 0 ; i < n ; i++) {
        int ti = i*d;
        for (int j = 0; j < d ; j++) {
            (*G)[ti+j] += wr*X[ti+j];
        }
    }
}

class FGData {
public:
    DistConstraints *E, *L, *U;
    double *w, *f;
};

/*class FGFunction {
public:
    DistConstraints E, L, U;
    double W[4], F[4];
    FGFunction(DistConstraints& e, DistConstraints& l, DistConstraints& u,
               double w[4], double f[4]) {

        E = e; L = l; U = u;
        for (int i = 0; i < 4; i++) {
            W[i] = w[i];
            F[i] = f[i];
        }
    }
    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
        gradfun(x, E, L, U, W, F, grad);
        double val = objfun(x, E, L, U, W, F);
        //printvec(x);
        //printf("(%lf)\n\n", val);
        return val;
    }
};*/

double fgcalc(const Eigen::VectorXd& x, Eigen::VectorXd* grad, void* dat) {
    FGData* fd = (FGData*)dat;
    gradfun(x, *(fd->E), *(fd->L), *(fd->U), fd->w, fd->f, grad);
    double val = objfun(x, *(fd->E), *(fd->L), *(fd->U), fd->w, fd->f);
    //printvec(x);
    //printf("(%lf)\n\n", val);
    return val;
}

int main(void) {
    test_w[3] = -0.0030;
    printf("running...\n");
    printf("objfun = %lf\n", objfun(test_rawX, test_equality_cons, test_lo_bounds,
                           test_up_bounds, test_w, test_f));
    printf("gradfun = ");
    Eigen::VectorXd g(test_rawX.size());
    g.setZero();
    gradfun(test_rawX, test_equality_cons, test_lo_bounds, test_up_bounds,
            test_w, test_f, &g);
    for (int i = 0; i < g.size(); i++) {
        printf("%lf ", g[i]);
    }
    g.setZero();
    printf("\n\n");

    // try doing the bfgs thingy
    //FGFunction f(test_equality_cons, test_lo_bounds, test_up_bounds, test_w, test_f);
    FGData data = {
        .E = &test_equality_cons,
        .L = &test_lo_bounds,
        .U = &test_up_bounds,
        .w = test_w,
        .f = test_f
    };
    optim::algo_settings_t settings = {
        .print_level = 1,
    };
    bool suc = optim::bfgs(test_rawX, fgcalc, &data, settings);
    if (suc)
        printf("success\n");
    else
        printf("failed\n");
    printvec(test_rawX);

    /*double out;
    LBFGSpp::LBFGSParam<double> p;
    p.max_iterations = 250;
    p.ftol = 1.0e-6;
    p.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_WOLFE;

    LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchBacktracking> solver(p);
    solver.minimize(f, test_rawX, out);
    printf("minimisation = %lf\n", out);
    for (int i = 0; i < test_rawX.size(); i++) {
        printf("%lf ", test_rawX[i]);
    }*/

    printf("\ndone\n");
    return 0;
}
