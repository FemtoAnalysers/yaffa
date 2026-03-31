#ifndef CONSTANTS_H
#define CONSTANTS_H

// Usefull constants
const double PI_value = 3.141592653589793;

// Physics constants
namespace phys {
    const double c_value = 299792458; // m/s
    const double Alpha = 1.0 / 137.035999;
    const double hbar_c = 197.3269631;
}

// Particle Mass in MeV
namespace ParticleMass { 
    const double MassProton = 938.2720813; 
    const double MassLambda = 1115.683;
    const double MassPionPlusMinus = 139.57061;
    const double MassPionZero = 134.9770;
}

// LambdaParameters
namespace LambdaParameters_pp {
    const double LambdaPar_pp_Gen = 0.67; // fm
    const double LambdaPar_pp_lam = 0.203;
    const double LambdaPar_pp_Flat = 0.116;
    const double LambdaPar_pp_Fake = 0.011;
}

#endif  