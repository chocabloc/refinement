#include <bits/stdc++.h>
#include "LBFGSpp/LBFGS.h"

/*namespace BFGS {
    using namespace Eigen;
    using Pars = LBFGSpp::LBFGSParam<double>;
    using FGFunc = std::function<double(std::vector<double>& x,
                                        std::vector<double>& g)>;
    using VecX = std::vector<double>;

    class Solver {
    private:
        size_t dim;
        FGFunc fg;
        VectorXd X;
        Pars pars;

        double fgwrapper(VectorXd& x, VectorXd& g) {
            return fg(x, )
        }
    public:
        Solver(int nvar, FGFunc func, VecX& x, Pars pars) {
            dim = nvar; fg = func; this->pars = pars;
            X = Eigen::Map<VectorXd>(x.data(), nvar);
        }

        Solver(int nvar, FGFunc func, VecX& x): Solver(nvar, func, x, Pars()) {} 

        int solve() {
            LBFGSpp::LBFGSSolver<double> solver(pars);

        }
    };
}*/

double f(Eigen::VectorXd& x, Eigen::VectorXd& g) {
    return 0.0;
}

int main(void) {
        using namespace LBFGSpp;
        Eigen::VectorXd x {{0, 0, 0}};
        double out;
        LBFGSParam<double> p;
        LBFGSSolver<double> solver(p);
        solver.minimize(f, x, out);
}