#pragma once
#include "fgdata.h"
#include "objgradfun.h"
#include "Eigen/Core"
#include "LBFGSpp/LBFGS.h"

class FGFunction {
public:
    FGData* dat;
    explicit FGFunction(Eigen::VectorXd& init_x, FGData* d) {
        dat = d;

        /*
            re-calculate w[3] how it is done in the matlab code:
            w(4) = ow(4)*obj/(25*trace(X0'*X0));
        */
        double s = 0.0, ow3 = dat->w[3];
        for (auto a: init_x)
            s += a*a;
        dat->w[3] = 0.0;
        double obval = objfun(init_x, dat),
               k = (ow3*obval/25.0)*(1.0/s);
        dat->w[3] = k;
    }

    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) const
    {
        gradfun(x, dat, &grad);
        double val = objfun(x, dat);
        //printvec(x);
        //printf("(%lf)\n\n", val);
        return val;
    }

    double runBFGS(Eigen::VectorXd& init_x) {
        double out;

        // parameters for LBFGS
        // TODO: maybe we need to tweak these?
        LBFGSpp::LBFGSParam p;
        p.max_iterations = 500;//1000;
        //p.ftol = 1.0e-14;
        p.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
        //p.epsilon = 1.0e-14;

        // run the solver
        LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchBacktracking> solver(p);
        solver.minimize(*this, init_x, out);

        return out;
    }
};

