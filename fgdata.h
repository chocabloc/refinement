#pragma once
#include <bits/stdc++.h>

// class for Distance Constraints
class DistConstraint {
public:
    int i, j, type;
    double val;
    DistConstraint(int I, int J, int T, double V): i(I), j(J), type(T), val(V) {} 
};
using DistConstraints = std::vector<DistConstraint>;

// class to store bounds and other optimisation data
class FGData {
public:
    DistConstraints *E, *L, *U;
    double *w, *f;
};
