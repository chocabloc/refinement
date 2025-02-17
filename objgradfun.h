//#include "mex.h"
#include <bits/stdc++.h>
#include <cmath>
#include "fgdata.h"
#include "Eigen/Core"

static void printvec(const Eigen::VectorXd& v) {
    for (double i : v)
        printf("%.4lf ", i);
    printf("\n");
}

// calculates value of objective function at X
static double objfun(const Eigen::VectorXd& X, FGData* dat) {
    
    double ret = 0.0, temp;

    auto& E = *(dat->E);
    auto& L= *(dat->L);
    auto& U = *(dat->U);
    auto w = dat->w, f = dat->f;

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
static void gradfun(const Eigen::VectorXd& X, FGData* dat, Eigen::VectorXd* G) {
    
    double temp;

    G->setZero();
    auto& E = *(dat->E);
    auto& L= *(dat->L);
    auto& U = *(dat->U);
    auto w = dat->w, f = dat->f;

    // dimensions and size of X
    int d = 3, n = X.size()/d;

    // number of bounds etc.
    int ne = E.size(), nl = L.size(), nu = U.size();

    // weights
    double we = 2*w[0], wl = 2*w[1], wu = 2*w[2], wr = 2*w[3];

    // idk what these are
    double f_hb  = f[0], f_tau = f[1], f_tal = f[2], f_vdw = f[3];
    
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

static double fgcalc(const Eigen::VectorXd& x, Eigen::VectorXd* grad, void* dat) {
    auto* fd = (FGData*)dat;
    if (grad)
        gradfun(x, fd, grad);
    double val = objfun(x, fd);

    /*
    printf("x = ");
    printvec(x);

    if (grad) {
        printf("g = ");
        printvec(*grad);
    }
    printf("(%lf)\n\n", val);
    */

    return val;
}
