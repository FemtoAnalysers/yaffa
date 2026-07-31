/* Various mathematical functions to be used in fits etc. */
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_hyperg.h"

#ifndef FUNCTIONS_HXX
#define FUNCTIONS_HXX

#include <cmath>

// Source function for 3 identical particles
// Ref.: PRC 109, 034006 (2024) (Eq. 40, 41)
// DOI: https://doi.org/10.1103/PhysRevC.109.034006
double _SourceAAA(double hyperRadius,  // Hyper-radius defined as in 3B NOTES
                  double rho0          // source size
) {
    return exp(-hyperRadius * hyperRadius / rho0 / rho0) / pow(rho0, 6) * pow(hyperRadius, 5);
}

// Gaussian source for 2B
double _SourceGauss(double rStar, double r0) {
    return 4 * M_PI * rStar * rStar * exp(-rStar * rStar / 4 / r0 / r0) / pow(4 * M_PI * r0 * r0, 1.5);
}

// Gaussian source for 3 identical particles expressed in Jacobi coordinates. Based on Mathematica calculation
double _SourceAAAJC(double r12, double r312, double r0) {
    double arg = - (3 * r12 * r12 + 4 * r312 * r312) / 12 / r0 / r0;
    double norm = pow(2 * sqrt(3) * M_PI * r0 * r0, -3);
    double jacobianr12 = 4 * M_PI * r12 * r12;
    double jacobianr312 = 4 * M_PI * r312 * r312;

    return jacobianr12 * jacobianr312 * norm * exp(arg);
}

/*
Regularized confluent hypergeometric function = 0F1(a, z) / Gamma(a).
See https://reference.wolfram.com/language/ref/Hypergeometric0F1Regularized.html
*/
double Hypergeometric0F1Regularized(double a, double z) {
    return gsl_sf_hyperg_0F1(a, z) / gsl_sf_gamma(a);
}

// Source function for 3 identical particles where 2 are primary and the 3rd one originates from a resonance
double _SourceAAApprAvg(double hypRad, double rp, double rs) {
    double rp2 = rp * rp;
    double rs2 = rs * rs;

    double z = (std::pow(rp2 - rs2, 2) * std::pow(hypRad, 4)) / (64 * std::pow(rp, 4) * std::pow(rp2 + 2 * rs2, 2));
    double chgr = Hypergeometric0F1Regularized(2, z);

    double arg = -(((2 * rp2 + rs2) * std::pow(hypRad, 2)) / (4 * rp2 * (rp2 + 2 * rs2)));
    double norm = 3 * std::sqrt(3) / (64 * pow(rp, 3) * pow(rp2 + 2 * rs2, 3. / 2));

    return norm * std::exp(arg) * std::pow(hypRad, 5) * chgr;
}

#endif
