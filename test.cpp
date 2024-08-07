#include "fgfunction.h"
#include "testdata.h"
#include <chrono>
using namespace std;

int main() {
    //printf("running LBFGS on test data...\n");
    FGData data = {
        .E = &test_equality_cons,
        .L = &test_lo_bounds,
        .U = &test_up_bounds,
        .w = test_w,
        .f = test_f
    };

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    FGFunction f(test_rawX, &data);
    double out = f.runBFGS(test_rawX);
    end = chrono::system_clock::now();
    chrono::duration<double> time = end - start;
    printf("%lf", time.count());

    /*
    printf("minimisation = %lf\n", out);
    for (double i : test_rawX)
        printf("%lf ", i);
    printf("\ndone\n");
*/

    // OptimLib stuff, couldn't get it to work well
    /*
    printf("objfun = %lf\n", objfun(test_rawX, &data));
    printf("gradfun = ");
    Eigen::VectorXd g(test_rawX.size());
    g.setZero();
    gradfun(test_rawX, &data, &g);
    for (int i = 0; i < g.size(); i++) {
        printf("%lf ", g[i]);
    }
    g.setZero();
    printf("\n\n");
*/

    // try doing the bfgs thingy
    /*optim::algo_settings_t settings = {
        //.print_level = 1,
        .iter_max = 250,
        .bfgs_settings = {
            .wolfe_cons_1 = 1e-4,
            .wolfe_cons_2 = 0.5
        }
    };
    bool suc = optim::bfgs(test_rawX, fgcalc, &data, settings);
    if (suc)
        printf("success\n");
    else
        printf("failed\n");
    printvec(test_rawX);*/


    return 0;
}