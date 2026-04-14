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

// Source function for 3 identical particles
double SourceCountsAAA(double* x, double* p) {
    // Variables
    double hyperRadius = x[0];

    // Parameters
    double norm = p[0];
    double rho0 = p[1];

    return norm * _SourceAAA(hyperRadius, rho0);
}


// Source function for 2 particles
double SourceGauss(double* x, double* p) {
    // Variables
    double rStar = x[0];
    
    // Parameters
    double r0 = p[0];

    return _SourceGauss(rStar, r0);
}

// Source function for 2 particles
double SourceCountsGauss(double* x, double* p) {
    // Variables
    double rStar = x[0];
    
    // Parameters
    double norm = p[0];
    double r0 = p[1];

    return norm * _SourceGauss(rStar, r0);
}

// Source function for 2 particles
double SourceAAAJC(double* x, double* p) {
    // Variables
    double r12 = x[0];
    double r312 = x[1];
    
    // Parameters
    double r0 = p[0];

    return _SourceAAAJC(r12, r312, r0);
}

// Source function for 2 particles
double SourceCountsAAAJC(double* x, double* p) {
    // Variables
    double r12 = x[0];
    double r312 = x[1];
    
    // Parameters
    double norm = p[0];
    double r0 = p[1];

    return norm * _SourceAAAJC(r12, r312, r0);
}
#endif
