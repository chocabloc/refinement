## Refinement Stuff

To run on test data (which is in `testdata.h`), do `g++ -I. -O3 test.cpp -o test`, and then `./test`.

To run on your own data, include `fgfunction.h`, put your data in a `FGData` structure,
create a `FGFunction` object, and call `runBFGS()` on that. The dependencies `LBFGSpp`, `fgdata.h`, and `objgradfun.h` 
should be present in the same directory, and `Eigen` should be in your include path.

Example usage:
```c++
#include "fgfunction.cpp"

int main() {
    // this is your data
    FGData data = {
        .E = &equality_constraints,
        .L = &lower_bounds,
        .U = &upper_bounds,
        .w = w,
        .f = f
    };
    
    // init_x is your starting point for the minimsation
    // (which is X0 in postProcessingMe.m)
    FGFunction f(init_x, &data);
    
    // run the BFGS
    // note that passing init_x will overwrite it with the result,
    // so pass a copy if you want to preserve it
    double out = f.runBFGS(init_x);

    // print out the minimal objective function value and the result
    printf("minimisation = %lf\n", out);
    for (double i : test_rawX)
        printf("%lf ", i);
    printf("\ndone\n");
}
```