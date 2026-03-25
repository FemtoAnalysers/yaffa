/* ROOT wrappers for functions to be used in TF1s */

#ifndef ROOTFUNCTIONS_HXX
#define ROOTFUNCTIONS_HXX

#include "Functions.hxx"

// Source function for 3 identical particles
double SourceAAA(double* x, double* p) {
    // Variables
    double hyperRadius = x[0];

    // Parameters
    double rho0 = p[0];

    return _SourceAAA(hyperRadius, rho0);
}

// Source function for 2 particles
double SourceGauss(double* x, double* p) {
    // Variables
    double rStar = x[0];
    
    // Parameters
    double r0 = p[0];

    return _SourceGauss(rStar, r0);
}

#endif
