/* ROOT wrappers for functions to be used in TF1s */

#ifndef ROOTFUNCTIONS_HXX
#define ROOTFUNCTIONS_HXX

#include "Functions.hxx"

/*
Source functions are of type:
    SourcePdf<X>     -> normalized to one by construction
    SourceCounts<X>  -> p[0] is arbitrary normalization
*/

// 2-body source functions ---------------------------------------------------------------------------------------------

// Gaussian source for 2 particles
double SourcePdfGauss(double* x, double* p) {
    // Variables
    double rStar = x[0];

    // Parameters
    double r0 = p[0];

    return _SourcePdfGauss(rStar, r0);
}

// Gaussian source for 2 identical particles including the effect of resonances
double SourcePdfGaussResonances(double* x, double* p) {
    // Variables
    double rStar = x[0];

    // Parameters
    double f = p[0];
    double rp = p[1];
    double delta = p[2];  // Delta radius: rs = rp + delta. p[2] must be limited > 0

    return _SourcePdfGaussResonances(rStar, f, rp, delta);
}

// 3-body source functions ---------------------------------------------------------------------------------------------

// Gaussian source for 3 identical particles, in the hyper-radius
double SourcePdfAAAHypRad(double* x, double* p) {
    // Variables
    double hypRad = x[0];

    // Parameters
    double rho0 = p[0];

    return _SourcePdfAAAHypRad(hypRad, rho0);
}

// Source for 3 identical particles where 2 are primary and the 3rd one originates from a resonance, in the hyper-radius
double SourcePdfAAApprHypRad(double* x, double* p) {
    // Variables
    double hypRad = x[0];

    // Parameters
    double rp = p[0];
    double rs = p[1];

    return _SourcePdfAAApprHypRad(hypRad, rp, rs);
}

// Source for 3 identical particles where 1 is primary and 2 originate from a resonance, in the hyper-radius
double SourcePdfAAAprrHypRad(double* x, double* p) {
    // Variables
    double hypRad = x[0];

    // Parameters
    double rp = p[0];
    double rs = p[1];

    // The same as primary-primary-resonances but with r_prim and r_reso switched
    return _SourcePdfAAApprHypRad(hypRad, rs, rp);
}

// Source for 3 identical particles including the effect of resonances, in the hyper-radius
double SourcePdfAAAGaussResonancesHypRad(double* x, double* p) {
    // Variables
    double hypRad = x[0];

    // Parameters
    double f = p[0];   // Fraction of primordial particles
    double rp = p[1];  // Single-particle radius of primordial particles
    double rs = p[2];  // Single-particle radius of secondary particles

    return _SourcePdfAAAGaussResonancesHypRad(hypRad, f, rp, rs);
}

// Hyper-angle distribution for 3 identical particles of the same kind. Since the source in this case is hypercentral,
// only the Jacobian survives. No parameters
double SourcePdfAAAHypAngle(double* x, double* p) {
    // Variables
    double hypAngle = x[0];

    return _SourcePdfAAAHypAngle(hypAngle);
}

// Hyper-angle distribution for 3 identical particles where 2 are primary and the 3rd one originates from a resonance
double SourcePdfAAApprHypAngle(double* x, double* p) {
    // Variables
    double hypAngle = x[0];

    // Parameters
    double rp = p[0];
    double rs = p[1];

    return _SourcePdfAAApprHypAngle(hypAngle, rp, rs);
}

// Hyper-angle distribution for 3 identical particles where 1 is primary and 2 originate from a resonance
double SourcePdfAAAprrHypAngle(double* x, double* p) {
    // Variables
    double hypAngle = x[0];

    // Parameters
    double rp = p[0];
    double rs = p[1];

    // The same as primary-primary-resonances but with r_prim and r_reso switched
    return _SourcePdfAAApprHypAngle(hypAngle, rs, rp);
}

// Hyper-angle distribution for 3 identical particles including the effect of resonances
double SourcePdfAAAGaussResonancesHypAngle(double* x, double* p) {
    // Variables
    double hypAngle = x[0];

    // Parameters
    double f = p[0];   // Fraction of primordial particles
    double rp = p[1];  // Single-particle radius of primordial particles
    double rs = p[2];  // Single-particle radius of secondary particles

    return _SourcePdfAAAGaussResonancesHypAngle(hypAngle, f, rp, rs);
}

// Gaussian source for 3 identical particles in Jacobi coordinates (r12, r3,12)
double SourcePdfAAAJC(double* x, double* p) {
    // Variables
    double r12 = x[0];
    double r312 = x[1];

    // Parameters
    double r0 = p[0];

    return _SourcePdfAAAJC(r12, r312, r0);
}

// Source for 3 identical particles where 2 are primary and the 3rd one originates from a resonance, in
// (hyper-radius, hyper-angle)
double SourcePdfAAAppr(double* x, double* p) {
    // Variables
    double hypRad = x[0];
    double hypAngle = x[1];

    // Parameters
    double rp = p[0];
    double rs = p[1];

    return _SourcePdfAAAppr(hypRad, hypAngle, rp, rs);
}

// Counts source functions ---------------------------------------------------------------------------------------------

double SourceCountsGauss(double* x, double* p) { return p[0] * SourcePdfGauss(x, p + 1); }
double SourceCountsGaussResonances(double* x, double* p) { return p[0] * SourcePdfGaussResonances(x, p + 1); }
double SourceCountsAAAHypRad(double* x, double* p) { return p[0] * SourcePdfAAAHypRad(x, p + 1); }
double SourceCountsAAApprHypRad(double* x, double* p) { return p[0] * SourcePdfAAApprHypRad(x, p + 1); }
double SourceCountsAAAprrHypRad(double* x, double* p) { return p[0] * SourcePdfAAAprrHypRad(x, p + 1); }
double SourceCountsAAAGaussResonancesHypRad(double* x, double* p) {
    return p[0] * SourcePdfAAAGaussResonancesHypRad(x, p + 1);
}
double SourceCountsAAAHypAngle(double* x, double* p) { return p[0] * SourcePdfAAAHypAngle(x, p + 1); }
double SourceCountsAAApprHypAngle(double* x, double* p) { return p[0] * SourcePdfAAApprHypAngle(x, p + 1); }
double SourceCountsAAAprrHypAngle(double* x, double* p) { return p[0] * SourcePdfAAAprrHypAngle(x, p + 1); }
double SourceCountsAAAGaussResonancesHypAngle(double* x, double* p) {
    return p[0] * SourcePdfAAAGaussResonancesHypAngle(x, p + 1);
}
double SourceCountsAAAJC(double* x, double* p) { return p[0] * SourcePdfAAAJC(x, p + 1); }
double SourceCountsAAAppr(double* x, double* p) { return p[0] * SourcePdfAAAppr(x, p + 1); }
#endif
