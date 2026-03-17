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

#endif
