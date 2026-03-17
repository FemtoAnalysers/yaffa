/* Various mathematical functions to be used in fits etc. */

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

#endif
