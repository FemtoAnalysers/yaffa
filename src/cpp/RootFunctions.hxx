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

// Source function for 2 identical particles including the effect of resonances
double SourceCountsGaussResonances(double* x, double* p) {
    // Variables
    double rStar = x[0];

    // Parameters
    double norm = p[0];
    double f = p[1];
    double rp = p[2];
    double rs = p[3];

    double source_pp = f * f * _SourceGauss(rStar, rp);
    double source_ps = 2 * f * (1 - f) * _SourceGauss(rStar, std::sqrt((rp * rp + rs * rs) / 2));
    double source_ss = (1 - f) * (1 - f) * _SourceGauss(rStar, rs);

    return norm * (source_pp + source_ps + source_ss);
}

// Source function for 3 identical particles including the effect of resonances
double SourceCountsAAAGaussResonances(double* x, double* p) {
    // Variables
    double hyperRadius = x[0];

    // Parameters
    double norm = p[0];
    double f = p[1]; // Fraction of primordinal particles
    double rp = p[2]; // Single-particle radius of primordial particles
    double rs = p[3]; // Single-particle radius of secondary particles

    double source_ppp = pow(f, 3) * _SourceAAA(hyperRadius, rp);
    double source_pps = 3 * f * f * (1 - f) * _SourceAAApprAvg(hyperRadius, rp, rs);
    double source_pss = 3 * f * pow(1 - f, 2) * _SourceAAApprAvg(hyperRadius, rs, rp); // Same as ppr with rp <--> rs
    double source_sss = pow(1 - f, 3) * _SourceGauss(hyperRadius, rs);

    return norm * (source_ppp + source_pps + source_pss + source_sss);
}

double SourceCountsAAApprAvg(double *x, double *p) {
    double hyperRadius = x[0];

    double norm = p[0];
    double rp = p[1];
    double rs = p[2];

    return norm * _SourceAAApprAvg(hyperRadius, rp, rs);
}

double SourceCountsAAAprrAvg(double *x, double *p) {
    double hyperRadius = x[0];

    double norm = p[0];
    double rp = p[1];
    double rs = p[2];

    // The same as primary-primary-resonances but with r_prim and r_reso switched
    return norm * _SourceAAApprAvg(hyperRadius, rs, rp);
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
