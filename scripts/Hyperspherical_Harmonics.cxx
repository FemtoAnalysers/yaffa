#include <TApplication.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TTree.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
// #include
// </Users/raffaele/alice/sw/osx_arm64/boost/v1.83.0-alice1-local1/include/boost/math/special_functions/spherical_harmonic.hpp>
// #include </opt/homebrew/lib/python3.9/site-packages/pythran/boost/math/special_functions/bessel.hpp>
// #include <gsl/gsl_sf.h>
#include <TGenPhaseSpace.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <stdio.h>
#include <tgmath.h>

#include <cstdlib>
#include <ctime>
#include <random>

#include "../src/cpp/Functions.hxx"
#include "../src/cpp/RootFunctions.hxx"

using namespace std;

auto RandomGen = new TRandom();
const double Pi = TMath::Pi();

const double fermi = 1. / 197.3269631;
const double m1 = 938.2720813;
const double m2 = 938.2720813;
const double m3 = 1115.683;
const double mu12 = m1 * m2 / (m2 + m1);
const double mu312 = m3 * (m2 + m1) / (m1 + m2 + m3);

const double Mass = m1 / 2.;

// const double alpha = 4.*(1.+pow(m3,2)/pow(m2+m3,2)+pow(m3,2)/pow(m1+m3,2));
// const double gamma = 4.*pow(m1+m2+m3,2)/pow(m1+m2,2)*(pow(m2,2)/pow(m2+m3,2)+pow(m1,2)/pow(m1+m3,2));
// const double Mass = mu12*alpha;

const double rhoMin = 900;
const double rhoMax = 1000;

const int Nsteps = 25000;
const double h = 0.04;

double calculate_coefficient(const int j_a, const int j_b, const int j_c, const int m_a, const int m_b, const int m_c) {
    const double wigner_3j = gsl_sf_coupling_3j(j_a, j_b, j_c, m_a, m_b, -m_c);
    assert((j_a - j_b + m_c) % 2 == 0);
    const int j = (j_a - j_b + m_c) / 2;
    double result = std::sqrt(j_c + 1) * wigner_3j;
    result *= pow(-1, j);

    // logg[LResonances].debug("CG: ", result, " I1: ", j_a, " I2: ", j_b,
    //                         " IR: ", j_c, " iz1: ", m_a, " iz2: ", m_b,
    //                         " izR: ", m_c);
    //
    return result;
}

// Function to calculate Jacobi polynomials P_n^{(\alpha, \beta)}(x)
double jacobiPolynomial(int n, double alpha, double beta, double x) {
    if (n == 0) {
        return 1.0;
    }
    if (n == 1) {
        return 0.5 * ((2 * (alpha + 1)) + (alpha + beta + 2) * (x - 1));
    }

    double P0 = 1.0;
    double P1 = 0.5 * ((2 * (alpha + 1)) + (alpha + beta + 2) * (x - 1));
    double Pn;

    for (int k = 1; k < n; ++k) {
        double a1 = 2 * (k + alpha) * (k + beta) * (2 * k + alpha + beta);
        double a2 = (2 * k + alpha + beta - 1) * (alpha * alpha - beta * beta);
        double a3 = (2 * k + alpha + beta) * (2 * k + alpha + beta - 1) * (2 * k + alpha + beta - 2);
        double a4 = 2 * (k + 1) * (k + alpha + beta + 1) * (2 * k + alpha + beta);

        Pn = ((a2 + a3 * x) * P1 - a1 * P0) / a4;
        P0 = P1;
        P1 = Pn;
    }

    return Pn;
}

double NormHH(int n, double l, double nu) {
    return sqrt((2. * nu * tgamma(nu - n) * tgamma(n + 1)) / (tgamma(nu - n - l - 1. / 2.) * tgamma(n + l + 3. / 2.)));
}

// Function to compute the spherical harmonics Y_lm
std::complex<double> spherical_harmonic(int l, int m, double theta, double phi) {
    // Compute the associated Legendre polynomial
    double plm = gsl_sf_legendre_sphPlm(l, std::abs(m), std::cos(theta));

    // Calculate the normalization factor
    double norm_factor = 1.;  // sqrt((2.0 * l + 1.0));
                              // std::sqrt((2.0 * l + 1.0)*
                              //           tgamma(1+l - abs(m)) / tgamma(1+l + abs(m)));

    // Compute the real and imaginary parts
    double real_part = norm_factor * plm * std::cos(m * phi);
    double imag_part = norm_factor * plm * std::sin(m * phi);

    return std::complex<double>(real_part, imag_part);
}

double IntegrateSH(int l1, int m1, int l2, int m2) {
    int MCsteps = 10000000;
    std::complex<double> HCP(0, 0);

    for (int i = 0; i < MCsteps; i++) {
        double theta12 = RandomGen->Uniform(0, Pi);
        double phi12 = RandomGen->Uniform(0, 2. * Pi);

        HCP = HCP + conj(spherical_harmonic(l1, m1, theta12, phi12)) * spherical_harmonic(l2, m2, theta12, phi12) *
                        sin(theta12);
    }

    HCP = HCP / double(MCsteps) * pow(Pi, 2) * 2.;
    return abs(HCP);
}

std::complex<double> HypersphericalHarmonics(int l12, int l312, int nu, int L, int M, double phi12, double theta12,
                                             double phi312, double theta312, double phi) {
    std::complex<double> SH(0, 0);

    for (int M12 = -l12; M12 <= l12; M12++) {
        for (int M312 = -l312; M312 <= l312; M312++) {
            double cg_coefficient = calculate_coefficient(2 * l12, 2 * l312, 2 * L, 2 * M12, 2 * M312, 2 * M);
            SH = SH + cg_coefficient * spherical_harmonic(l12, M12, theta12, phi12) *
                          spherical_harmonic(l312, M312, theta312, phi312);
        }
    }

    return SH * NormHH(nu, l312, l12 + l312 + 2. * nu + 2) * pow(cos(phi), l312) * pow(sin(phi), l12) *
           jacobiPolynomial(nu, l12 + 1. / 2., l312 + 1. / 2., cos(2. * phi));
}

double IntegrateHH(int l12_1, int l312_1, int nu_1, int L_1, int M_1, int l12_2, int l312_2, int nu_2, int L_2,
                   int M_2) {
    int MCsteps = 1000000;
    std::complex<double> HCP(0, 0);

    for (int i = 0; i < MCsteps; i++) {
        double theta12 = RandomGen->Uniform(0, Pi);
        double theta312 = RandomGen->Uniform(0, Pi);
        double phi12 = RandomGen->Uniform(0, 2. * Pi);
        double phi312 = RandomGen->Uniform(0, 2. * Pi);
        double phi = RandomGen->Uniform(0, Pi / 2.);

        HCP =
            HCP + conj(HypersphericalHarmonics(l12_1, l312_1, nu_1, L_1, M_1, phi12, theta12, phi312, theta312, phi)) *
                      HypersphericalHarmonics(l12_2, l312_2, nu_2, L_2, M_2, phi12, theta12, phi312, theta312, phi) *
                      pow(sin(phi) * cos(phi), 2) * sin(theta12) * sin(theta312);
    }

    HCP = HCP / double(MCsteps) * pow(Pi, 5) * 2.;
    return abs(HCP);
}

std::complex<double> HypersphericalHarmonicsJJz(int K, int l12, int l312, double S, int L, double J, double Jz,
                                                double phi12, double theta12, double phi312, double theta312,
                                                double phi) {
    std::complex<double> SHJ(0, 0);

    double nu = (K - l12 - l312) / 2.;

    for (int M = -L; M <= L; M++) {
        for (int i = 0; i < int(2. * S + 1.); i++) {
            double MS = -S + i;
            double cg_coefficient = calculate_coefficient(2 * S, 2 * L, 2 * J, 2 * MS, 2 * M, 2 * Jz);
            SHJ = SHJ +
                  cg_coefficient * HypersphericalHarmonics(l12, l312, nu, L, M, phi12, theta12, phi312, theta312, phi);
        }
    }

    return SHJ;
}

double IntegrateHHJz(int K_1, int l12_1, int l312_1, double S_1, int L_1, double J_1, double Jz_1, int K_2, int l12_2,
                     int l312_2, double S_2, int L_2, double J_2, double Jz_2) {
    int MCsteps = 10000000;
    std::complex<double> HCP(0, 0);

    if (abs(J_1 - J_2) > 0.5 || abs(Jz_1 - Jz_2) > 0.5)
        return abs(HCP);

    for (int i = 0; i < MCsteps; i++) {
        double theta12 = RandomGen->Uniform(0, Pi);
        double theta312 = RandomGen->Uniform(0, Pi);
        double phi12 = RandomGen->Uniform(0, 2. * Pi);
        double phi312 = RandomGen->Uniform(0, 2. * Pi);
        double phi = RandomGen->Uniform(0, Pi / 2.);

        HCP = HCP + conj(HypersphericalHarmonicsJJz(K_1, l12_1, l312_1, S_1, L_1, J_1, Jz_1, phi12, theta12, phi312,
                                                    theta312, phi)) *
                        HypersphericalHarmonicsJJz(K_2, l12_2, l312_2, S_2, L_2, J_2, Jz_2, phi12, theta12, phi312,
                                                   theta312, phi) *
                        pow(sin(phi) * cos(phi), 2) * sin(theta12) * sin(theta312);
    }

    HCP = HCP / double(MCsteps) * pow(Pi, 5) * 2.;
    return abs(HCP);
}

void CalculateNumberOfStates(int K) {
    double NoS = pow(K + 2, 2) * (K + 1) * (K + 3) / 12.;
    cout << "The total number of states is " << NoS << endl;
    int Nnu = int(K / 2.) + 1;
    int count = 0;

    for (int nu = 0; nu < Nnu; nu++) {
        int L12 = K - 2 * nu;
        for (int j = 0; j < L12 + 1; j++) {
            int L1 = j;
            int L2 = L12 - j;
            for (int L = abs(L1 - L2); L <= (L1 + L2); L++) {
                for (int Lz = -L; Lz <= L; Lz++) {
                    cout << K << "  " << L1 << "    " << L2 << "   " << nu << "  " << L << " " << Lz << endl;
                    count++;
                }
            }
        }
    }
    cout << count << endl;
}

double DotSpin(double S1, int S23, double MS1, double S2, int S12, double MS2) {
    double result = 0;
    if (abs(S1 - S2) < 0.5 && abs(MS1 - MS2) < 0.5) {
        if (S12 == 0 && S23 == 0)
            result = sqrt(1. / 4.);
        if (S12 == 0 && S23 == 1)
            result = sqrt(3. / 4.);
        if (S12 == 1 && S23 == 0 && S1 == 1. / 2.)
            result = sqrt(3. / 4.);
        if (S12 == 1 && S23 == 1 && S1 == 1. / 2.)
            result = sqrt(1. / 4.);
        if (S12 == 1 && S23 == 0 && S1 == 3. / 2.)
            result = 1.;
        if (S12 == 1 && S23 == 1 && S1 == 3. / 2.)
            result = 1.;
    }

    return result;
}

double IntegrateHHCC(int K, int L, double ML, int l12_1, int l312_1, int l12_2, int l312_2) {
    int MCsteps = 1000000;
    std::complex<double> HCP(0, 0);

    for (int i = 0; i < MCsteps; i++) {
        double theta12 = RandomGen->Uniform(0, Pi);
        double theta312 = RandomGen->Uniform(0, Pi);
        double phi12 = RandomGen->Uniform(0, 2. * Pi);
        double phi312 = RandomGen->Uniform(0, 2. * Pi);
        double phi = RandomGen->Uniform(0, Pi / 2.);

        double r12 = cos(phi);
        double r12x = cos(phi) * sin(theta12) * cos(phi12);
        double r12y = cos(phi) * sin(theta12) * sin(phi12);
        double r12z = cos(phi) * cos(theta12);
        double r312 = sin(phi);
        double r312x = sin(phi) * sin(theta312) * cos(phi312);
        double r312y = sin(phi) * sin(theta312) * sin(phi312);
        double r312z = sin(phi) * cos(theta312);

        double a11 = -sqrt((m1 * m3) / ((m1 + m2) * (m2 + m3)));
        double a12 = -sqrt((m2 * (m1 + m2 + m3)) / ((m2 + m3) * (m1 + m2)));
        double r23x = a11 * r12x + a12 * r312x;
        double r23y = a11 * r12y + a12 * r312y;
        double r23z = a11 * r12z + a12 * r312z;
        double r23 = sqrt(pow(r23x, 2) + pow(r23y, 2) + pow(r23z, 2));
        double a21 = -a12;
        double a22 = a11;
        double r123x = a21 * r12x + a22 * r312x;
        double r123y = a21 * r12y + a22 * r312y;
        double r123z = a21 * r12z + a22 * r312z;
        double r123 = sqrt(pow(r123x, 2) + pow(r123y, 2) + pow(r123z, 2));
        double phi_2 = atan(r123 / r23);
        double theta123 = acos(r123z / r123);
        double phi123 = atan(r123y / r123x);
        double theta23 = acos(r23z / r23);
        double phi23 = atan(r23y / r23x);

        double b11 = -sqrt((m2 * m3) / ((m1 + m2) * (m1 + m3)));
        double b12 = sqrt((m1 * (m1 + m2 + m3)) / ((m1 + m3) * (m1 + m2)));
        double r31x = b11 * r12x + b12 * r312x;
        double r31y = b11 * r12y + b12 * r312y;
        double r31z = b11 * r12z + b12 * r312z;
        double r31 = sqrt(pow(r31x, 2) + pow(r31y, 2) + pow(r31z, 2));
        double b21 = -b12;
        double b22 = b11;
        double r231x = b21 * r12x + b22 * r312x;
        double r231y = b21 * r12y + b22 * r312y;
        double r231z = b21 * r12z + b22 * r312z;
        double r231 = sqrt(pow(r231x, 2) + pow(r231y, 2) + pow(r231z, 2));
        double phi_3 = atan(r231 / r31);
        double theta231 = acos(r231z / r231);
        double phi231 = atan(r231y / r231x);
        double theta31 = acos(r31z / r31);
        double phi31 = atan(r31y / r31x);

        int nu_1 = (K - l12_1 - l312_1) / 2.;
        int nu_2 = (K - l12_2 - l312_2) / 2.;
        HCP = HCP + conj(HypersphericalHarmonics(l12_1, l312_1, nu_1, L, ML, phi23, theta23, phi123, theta123, phi_2)) *
                        HypersphericalHarmonics(l12_2, l312_2, nu_2, L, ML, phi23, theta23, phi123, theta123, phi_2) *
                        pow(sin(phi) * cos(phi), 2) * sin(theta12) * sin(theta312);
    }

    HCP = HCP / double(MCsteps) * pow(Pi, 5) * 2.;
    return abs(HCP);
}

double aCoeff23(int K, int L, double S, double J, double Jz, int l12, int l312, int S12, int l23, int l123, int S23) {
    std::complex<double> SHJ(0, 0);

    for (int M = -L; M <= L; M++) {
        for (int i = 0; i < int(2. * S + 1.); i++) {
            double MS = -S + i;
            double cg_coefficient = calculate_coefficient(2 * S, 2 * L, 2 * J, 2 * MS, 2 * M, 2 * Jz);
            SHJ = SHJ + pow(cg_coefficient, 2) * IntegrateHHCC(K, L, M, l23, l123, l12, l312) *
                            DotSpin(S, S23, MS, S, S12, MS);
        }
    }

    return abs(SHJ);
}

double IntegrateHHCC31(int K, int L, double ML, int l12_1, int l312_1, int l12_2, int l312_2) {
    int MCsteps = 1000000;
    std::complex<double> HCP(0, 0);

    for (int i = 0; i < MCsteps; i++) {
        double theta12 = RandomGen->Uniform(0, Pi);
        double theta312 = RandomGen->Uniform(0, Pi);
        double phi12 = RandomGen->Uniform(0, 2. * Pi);
        double phi312 = RandomGen->Uniform(0, 2. * Pi);
        double phi = RandomGen->Uniform(0, Pi / 2.);

        double r12 = cos(phi);
        double r12x = cos(phi) * sin(theta12) * cos(phi12);
        double r12y = cos(phi) * sin(theta12) * sin(phi12);
        double r12z = cos(phi) * cos(theta12);
        double r312 = sin(phi);
        double r312x = sin(phi) * sin(theta312) * cos(phi312);
        double r312y = sin(phi) * sin(theta312) * sin(phi312);
        double r312z = sin(phi) * cos(theta312);

        double a11 = -sqrt((m1 * m3) / ((m1 + m2) * (m2 + m3)));
        double a12 = -sqrt((m2 * (m1 + m2 + m3)) / ((m2 + m3) * (m1 + m2)));
        double r23x = a11 * r12x + a12 * r312x;
        double r23y = a11 * r12y + a12 * r312y;
        double r23z = a11 * r12z + a12 * r312z;
        double r23 = sqrt(pow(r23x, 2) + pow(r23y, 2) + pow(r23z, 2));
        double a21 = -a12;
        double a22 = a11;
        double r123x = a21 * r12x + a22 * r312x;
        double r123y = a21 * r12y + a22 * r312y;
        double r123z = a21 * r12z + a22 * r312z;
        double r123 = sqrt(pow(r123x, 2) + pow(r123y, 2) + pow(r123z, 2));
        double phi_2 = atan(r123 / r23);
        double theta123 = acos(r123z / r123);
        double phi123 = atan(r123y / r123x);
        double theta23 = acos(r23z / r23);
        double phi23 = atan(r23y / r23x);

        double b11 = -sqrt((m2 * m3) / ((m1 + m2) * (m1 + m3)));
        double b12 = sqrt((m1 * (m1 + m2 + m3)) / ((m1 + m3) * (m1 + m2)));
        double r31x = b11 * r12x + b12 * r312x;
        double r31y = b11 * r12y + b12 * r312y;
        double r31z = b11 * r12z + b12 * r312z;
        double r31 = sqrt(pow(r31x, 2) + pow(r31y, 2) + pow(r31z, 2));
        double b21 = -b12;
        double b22 = b11;
        double r231x = b21 * r12x + b22 * r312x;
        double r231y = b21 * r12y + b22 * r312y;
        double r231z = b21 * r12z + b22 * r312z;
        double r231 = sqrt(pow(r231x, 2) + pow(r231y, 2) + pow(r231z, 2));
        double phi_3 = atan(r231 / r31);
        double theta231 = acos(r231z / r231);
        double phi231 = atan(r231y / r231x);
        double theta31 = acos(r31z / r31);
        double phi31 = atan(r31y / r31x);

        int nu_1 = (K - l12_1 - l312_1) / 2.;
        int nu_2 = (K - l12_2 - l312_2) / 2.;
        HCP = HCP + conj(HypersphericalHarmonics(l12_1, l312_1, nu_1, L, ML, phi31, theta31, phi231, theta231, phi_3)) *
                        HypersphericalHarmonics(l12_2, l312_2, nu_2, L, ML, phi31, theta31, phi231, theta231, phi_3) *
                        pow(sin(phi) * cos(phi), 2) * sin(theta12) * sin(theta312);
    }

    HCP = HCP / double(MCsteps) * pow(Pi, 5) * 2.;
    return abs(HCP);
}

double aCoeff31(int K, int L, double S, double J, double Jz, int l12, int l312, int S12, int l23, int l123, int S23) {
    std::complex<double> SHJ(0, 0);

    for (int M = -L; M <= L; M++) {
        for (int i = 0; i < int(2. * S + 1.); i++) {
            double MS = -S + i;
            double cg_coefficient = calculate_coefficient(2 * S, 2 * L, 2 * J, 2 * MS, 2 * M, 2 * Jz);
            SHJ = SHJ + pow(cg_coefficient, 2) * IntegrateHHCC31(K, L, M, l23, l123, l12, l312) *
                            DotSpin(S, S23, MS, S, S12, MS);
        }
    }

    return abs(SHJ);
}

double IntegrateHHJzCC(int K_1, int l12_1, int l312_1, double S_1, int L_1, double J_1, double Jz_1, int K_2, int l12_2,
                       int l312_2, double S_2, int L_2, double J_2, double Jz_2) {
    int MCsteps = 10000000;
    std::complex<double> HCP(0, 0);

    if (abs(J_1 - J_2) > 0.5 || abs(Jz_1 - Jz_2) > 0.5)
        return abs(HCP);

    for (int i = 0; i < MCsteps; i++) {
        double theta12 = RandomGen->Uniform(0, Pi);
        double theta312 = RandomGen->Uniform(0, Pi);
        double phi12 = RandomGen->Uniform(0, 2. * Pi);
        double phi312 = RandomGen->Uniform(0, 2. * Pi);
        double phi = RandomGen->Uniform(0, Pi / 2.);

        double CosTh = cos(theta12) * cos(theta312) + sin(theta12) * sin(theta312) * cos(phi12 - phi312);

        double r12 = cos(phi);
        double r12x = cos(phi) * sin(theta12) * cos(phi12);
        double r12y = cos(phi) * sin(theta12) * sin(phi12);
        double r12z = cos(phi) * cos(theta12);
        double r312 = sin(phi);
        double r312x = sin(phi) * sin(theta312) * cos(phi312);
        double r312y = sin(phi) * sin(theta312) * sin(phi312);
        double r312z = sin(phi) * cos(theta312);

        double a11 = -sqrt((m1 * m3) / ((m1 + m2) * (m2 + m3)));
        double a12 = -sqrt((m2 * (m1 + m2 + m3)) / ((m2 + m3) * (m1 + m2)));
        double r23x = a11 * r12x + a12 * r312x;
        double r23y = a11 * r12y + a12 * r312y;
        double r23z = a11 * r12z + a12 * r312z;
        double r23 = sqrt(pow(r23x, 2) + pow(r23y, 2) + pow(r23z, 2));
        double a21 = -a12;
        double a22 = a11;
        double r123x = a21 * r12x + a22 * r312x;
        double r123y = a21 * r12y + a22 * r312y;
        double r123z = a21 * r12z + a22 * r312z;
        double r123 = sqrt(pow(r123x, 2) + pow(r123y, 2) + pow(r123z, 2));
        double phi_2 = atan(r123 / r23);
        double theta123 = acos(r123z / r123);
        double phi123 = atan(r123y / r123x);
        double theta23 = acos(r23z / r23);
        double phi23 = atan(r23y / r23x);

        double b11 = -sqrt((m2 * m3) / ((m1 + m2) * (m1 + m3)));
        double b12 = sqrt((m1 * (m1 + m2 + m3)) / ((m1 + m3) * (m1 + m2)));
        double r31x = b11 * r12x + b12 * r312x;
        double r31y = b11 * r12y + b12 * r312y;
        double r31z = b11 * r12z + b12 * r312z;
        double r31 = sqrt(pow(r31x, 2) + pow(r31y, 2) + pow(r31z, 2));
        double b21 = -b12;
        double b22 = b11;
        double r231x = b21 * r12x + b22 * r312x;
        double r231y = b21 * r12y + b22 * r312y;
        double r231z = b21 * r12z + b22 * r312z;
        double r231 = sqrt(pow(r231x, 2) + pow(r231y, 2) + pow(r231z, 2));
        double phi_3 = atan(r231 / r31);
        double theta231 = acos(r231z / r231);
        double phi231 = atan(r231y / r231x);
        double theta31 = acos(r31z / r31);
        double phi31 = atan(r31y / r31x);

        HCP = HCP + conj(HypersphericalHarmonicsJJz(K_1, l12_1, l312_1, S_1, L_1, J_1, Jz_1, phi23, theta23, phi123,
                                                    theta123, phi_2)) *
                        HypersphericalHarmonicsJJz(K_2, l12_2, l312_2, S_2, L_2, J_2, Jz_2, phi12, theta12, phi312,
                                                   theta312, phi) *
                        pow(sin(phi) * cos(phi), 2) * sin(theta12) * sin(theta312);
    }

    HCP = HCP / double(MCsteps) * pow(Pi, 5) * 2.;
    return abs(HCP);
}

double Vpp(double r, int l12, double S12) {
    double V0 = -30.45;  // in MeV
    double r0 = 1.815;   // in fm
    double V1 = 0;       // in MeV
    double r1 = 1.;      // in fm

    double Vpotpp = V0 * exp(-pow(r / r0, 2)) * (1 - 1. / 2. * S12 * (S12 + 1)) +
                    V1 * exp(-pow(r / r1, 2)) * 1. / 2. * S12 * (S12 + 1);
    if (l12 != 0)
        Vpotpp = 0;

    return Vpotpp;
}

double VpL(double r, int S12, double Stot) {
    double V0 = -29.9;  // in MeV
    double r0 = 1.47;   // in fm
    double V1 = -20.6;  // in MeV
    double r1 = 1.547;  // in fm
    double S0factor = 0.;
    double S1factor = 0.;

    if (S12 == 0) {
        S0factor = 1. / 4.;
        S1factor = 3. / 4.;
    }
    if (S12 == 1) {
        if (Stot == 1. / 2.) {
            S0factor = 3. / 4.;
            S1factor = 1. / 4.;
        }
        if (Stot == 3. / 2.) {
            S0factor = 1.;
            S1factor = 1.;
        }
    }

    return V0 * exp(-pow(r / r0, 2)) * S0factor + V1 * exp(-pow(r / r1, 2)) * S1factor;
}

double CoulombScreenedPot(double x, double qq) {
    double alphaQED = 1. / 137.;
    double xScreen = 100.;
    return qq * alphaQED / (x * fermi) * expl(-powl(x / xScreen, 4));
}

double VppL(double x, double y) {
    double W3 = 11.795;  // in MeV
    double rho3 = 2.0;   // in fm

    return W3 * exp(-(pow(x, 2) + pow(y, 2)) / pow(rho3, 2));
}

double Vpot(double r12, double r23, double r31, int l12, int S12, double Stot) {
    return Vpp(r12, l12, S12) + CoulombScreenedPot(r12, 1.) + VpL(r23, S12, Stot) +
           VpL(r31, S12, Stot);  // + VppL(r23,r31);
    // return Vpp(r12,S12) + VpL(r23,S12,Stot) + VpL(r31,S12,Stot);
}

const int Nstates = 76;

void HypercentralPotentials(double Kin, double JtotIN) {
    std::ifstream inputFile;
    inputFile.open(Form("States12With_K%i_J%1.1f.txt", int(Kin), JtotIN));

    double K_1[Nstates], S_1[Nstates], J_1[Nstates], Jz_1[Nstates];
    int l12_1[Nstates], l312_1[Nstates], L_1[Nstates], S12_1[Nstates];

    double K_2[Nstates], S_2[Nstates], J_2[Nstates], Jz_2[Nstates];
    int l12_2[Nstates], l312_2[Nstates], L_2[Nstates], S12_2[Nstates];

    for (int i = 0; i < Nstates; i++) {
        inputFile >> K_1[i] >> l12_1[i] >> l312_1[i] >> S_1[i] >> L_1[i] >> S12_1[i] >> J_1[i] >> Jz_1[i];
        K_2[i] = K_1[i];
        l12_2[i] = l12_1[i];
        l312_2[i] = l312_1[i];
        S_2[i] = S_1[i];
        L_2[i] = L_1[i];
        S12_2[i] = S12_1[i];
        J_2[i] = J_1[i];
        Jz_2[i] = Jz_1[i];
    }

    TGraph* hPot[Nstates][Nstates];
    TString Name[Nstates][Nstates];

    for (int count = 0; count < Nstates; count++) {
        cout << count << "  " << K_1[count] << "  " << l12_1[count] << "  " << l312_1[count] << "  " << S_1[count]
             << "  " << L_1[count] << "  " << S12_1[count] << "  " << J_1[count] << "  " << Jz_1[count] << endl;
        for (int countprime = count; countprime < Nstates; countprime++) {
            hPot[count][countprime] = new TGraph();
            cout << count << "-" << countprime << endl;
            double rho = 0.;
            for (int u = 0; u < Nsteps; u++) {
                double uoffset = 0.001;
                double stepsize = 0.01;
                int MCsteps = 1000;  // 10000;
                if (rho <= 1.)
                    rho = u * stepsize + uoffset;
                if (rho > 1. && rho <= 4.) {
                    stepsize = 0.1;
                    rho = rho + stepsize;
                    MCsteps = 1000;  // 10000;
                }
                if (rho > 4. && rho <= 50.) {
                    stepsize = 0.2;
                    rho = rho + stepsize;
                    MCsteps = 100;  // 1000;
                }
                if (rho > 50.) {
                    stepsize = 1.;
                    rho = rho + stepsize;
                    MCsteps = 100;  // 1000;
                }
                if (rho > rhoMax)
                    continue;
                std::complex<double> HCP = 0;
                for (int i = 0; i < MCsteps; i++) {
                    if (abs(S_1[count] - S_2[countprime]) > 0.5)
                        continue;
                    if (abs(S12_1[count] - S12_2[countprime]) > 0.5)
                        continue;

                    double theta12 = RandomGen->Uniform(0, Pi);
                    double theta312 = RandomGen->Uniform(0, Pi);
                    double phi12 = RandomGen->Uniform(0, 2. * Pi);
                    double phi312 = RandomGen->Uniform(0, 2. * Pi);
                    double phi = RandomGen->Uniform(0, Pi / 2.);

                    double CosTh = cos(theta12) * cos(theta312) + sin(theta12) * sin(theta312) * cos(phi12 - phi312);

                    double r12 = sqrt(Mass / mu12) * rho * cos(phi);
                    double r312 = sqrt(Mass / mu312) * rho * sin(phi);
                    double r31 =
                        sqrt(pow(m2 / (m1 + m2) * r12, 2) + pow(r312, 2) - 2. * r312 * m2 / (m1 + m2) * r12 * CosTh);
                    double r23 =
                        sqrt(pow(m1 / (m1 + m2) * r12, 2) + pow(r312, 2) + 2. * r312 * m1 / (m1 + m2) * r12 * CosTh);

                    HCP =
                        HCP + conj(HypersphericalHarmonicsJJz(K_1[count], l12_1[count], l312_1[count], S_1[count],
                                                              L_1[count], J_1[count], Jz_1[count], phi12, theta12,
                                                              phi312, theta312, phi)) *
                                  Vpot(r12, r23, r31, l12_1[count], S12_1[count], S_1[count]) *
                                  HypersphericalHarmonicsJJz(K_2[countprime], l12_2[countprime], l312_2[countprime],
                                                             S_2[countprime], L_2[countprime], J_2[countprime],
                                                             Jz_2[countprime], phi12, theta12, phi312, theta312, phi) *
                                  pow(sin(phi) * cos(phi), 2) * sin(theta12) * sin(theta312);
                }
                HCP = HCP / double(MCsteps) * pow(Pi, 5) * 2.;
                hPot[count][countprime]->SetPoint(u, rho, real(HCP));
                Name[count][countprime] = Form("Potential_Matrix_%i_%i", count, countprime);
            }
        }
    }

    TFile* OutPut = new TFile(Form("OutputPot_K%i_J%1.1f.root", int(Kin), JtotIN), "RECREATE");
    OutPut->cd();
    for (int i = 0; i < Nstates; i++) {
        for (int j = i; j < Nstates; j++) {
            hPot[i][j]->Write(Name[i][j].Data());
        }
    }

    OutPut->Close();
}

void HypercentralPotentialsNEW(double Kin, double JtotIN) {
    std::ifstream inputFile;
    inputFile.open(Form("StatesWith_K%i_J%1.1f.txt", int(Kin), JtotIN));

    double K_1[Nstates], S_1[Nstates], J_1[Nstates], Jz_1[Nstates];
    int l12_1[Nstates], l312_1[Nstates], L_1[Nstates], S12_1[Nstates];

    double K_2[Nstates], S_2[Nstates], J_2[Nstates], Jz_2[Nstates];
    int l12_2[Nstates], l312_2[Nstates], L_2[Nstates], S12_2[Nstates];

    for (int i = 0; i < Nstates; i++) {
        inputFile >> K_1[i] >> l12_1[i] >> l312_1[i] >> S_1[i] >> L_1[i] >> S12_1[i] >> J_1[i] >> Jz_1[i];
        K_2[i] = K_1[i];
        l12_2[i] = l12_1[i];
        l312_2[i] = l312_1[i];
        S_2[i] = S_1[i];
        L_2[i] = L_1[i];
        S12_2[i] = S12_1[i];
        J_2[i] = J_1[i];
        Jz_2[i] = Jz_1[i];
    }

    TGraph* hPot[Nstates][Nstates];
    TString Name[Nstates][Nstates];
    std::complex<double> HCP[Nstates][Nstates];

    double rho = 0.;
    for (int u = 0; u < Nsteps; u++) {
        double uoffset = 0.001;
        double stepsize = 0.01;
        int MCsteps = 500;
        if (rho <= 1.)
            rho = u * stepsize + uoffset;
        if (rho > 1. && rho <= 4.) {
            stepsize = 0.1;
            rho = rho + stepsize;
            // MCsteps = 500;
        }
        if (rho > 4. && rho <= 50.) {
            stepsize = 0.2;
            rho = rho + stepsize;
            // MCsteps = 100;
        }
        if (rho > 50.) {
            stepsize = 1.;
            rho = rho + stepsize;
            // MCsteps = 100;
        }
        if (rho > rhoMax)
            continue;
        cout << rho << " fm" << endl;
        for (int i = 0; i < MCsteps; i++) {
            double theta12 = RandomGen->Uniform(0, Pi);
            double theta312 = RandomGen->Uniform(0, Pi);
            double phi12 = RandomGen->Uniform(0, 2. * Pi);
            double phi312 = RandomGen->Uniform(0, 2. * Pi);
            double phi = RandomGen->Uniform(0, Pi / 2.);

            double CosTh = cos(theta12) * cos(theta312) + sin(theta12) * sin(theta312) * cos(phi12 - phi312);

            double r12 = sqrt(Mass / mu12) * rho * cos(phi);
            double r312 = sqrt(Mass / mu312) * rho * sin(phi);
            double r31 = sqrt(pow(m2 / (m1 + m2) * r12, 2) + pow(r312, 2) - 2. * r312 * m2 / (m1 + m2) * r12 * CosTh);
            double r23 = sqrt(pow(m1 / (m1 + m2) * r12, 2) + pow(r312, 2) + 2. * r312 * m1 / (m1 + m2) * r12 * CosTh);
            for (int count = 0; count < Nstates; count++) {
                // cout<<count<<"  "<<K_1[count]<<"  "<<l12_1[count]<<"  "<<l312_1[count]<<"  "<<S_1[count]<<"
                // "<<L_1[count]<<"  "<<S12_1[count]<<"  "<<J_1[count]<<"  "<<Jz_1[count]<<endl;
                for (int countprime = count; countprime < Nstates; countprime++) {
                    if (abs(S_1[count] - S_2[countprime]) > 0.5)
                        continue;
                    if (abs(S12_1[count] - S12_2[countprime]) > 0.5)
                        continue;
                    if (u == 0 && i == 0) {
                        hPot[count][countprime] = new TGraph();
                    }
                    if (i == 0)
                        HCP[count][countprime] = 0.;
                    HCP[count][countprime] =
                        HCP[count][countprime] +
                        conj(HypersphericalHarmonicsJJz(K_1[count], l12_1[count], l312_1[count], S_1[count], L_1[count],
                                                        J_1[count], Jz_1[count], phi12, theta12, phi312, theta312,
                                                        phi)) *
                            Vpot(r12, r23, r31, l12_1[count], S12_1[count], S_1[count]) *
                            HypersphericalHarmonicsJJz(K_2[countprime], l12_2[countprime], l312_2[countprime],
                                                       S_2[countprime], L_2[countprime], J_2[countprime],
                                                       Jz_2[countprime], phi12, theta12, phi312, theta312, phi) *
                            pow(sin(phi) * cos(phi), 2) * sin(theta12) * sin(theta312);
                    if (i == MCsteps - 1) {
                        HCP[count][countprime] = HCP[count][countprime] / double(MCsteps) * pow(Pi, 5) * 2.;
                        hPot[count][countprime]->SetPoint(u, rho, real(HCP[count][countprime]));
                        Name[count][countprime] = Form("Potential_Matrix_%i_%i", count, countprime);
                    }
                }
            }
        }
    }

    TFile* OutPut = new TFile(Form("OutputPot_K%i_J%1.1f_test.root", int(Kin), JtotIN), "RECREATE");
    OutPut->cd();
    for (int i = 0; i < Nstates; i++) {
        for (int j = i; j < Nstates; j++) {
            hPot[i][j]->Write(Name[i][j].Data());
        }
    }

    OutPut->Close();
}

void ExtractMatrixTEST(double Kin, double JtotIN, double rho) {
    std::ifstream inputFile;
    inputFile.open(Form("States12With_K%i_J%1.1f.txt", int(Kin), JtotIN));

    double K[Nstates], S[Nstates], J[Nstates], Jz[Nstates];
    int l12[Nstates], l312[Nstates], L[Nstates], S12[Nstates];
    TGraph* hPot[Nstates][Nstates];
    TString Name[Nstates][Nstates];

    TFile* InPut = new TFile(Form("OutputPot_K%i_J%1.1f.root", int(Kin), JtotIN), "READ");
    InPut->cd();

    for (int i = 0; i < Nstates; i++) {
        inputFile >> K[i] >> l12[i] >> l312[i] >> S[i] >> L[i] >> S12[i] >> J[i] >> Jz[i];
        for (int j = i; j < Nstates; j++) {
            Name[i][j] = Form("Potential_Matrix_%i_%i", i, j);
            hPot[i][j] = (TGraph*)InPut->Get(Name[i][j].Data());
        }
    }
    inputFile.close();
    InPut->Close();

    // cout<<"CIAO"<<endl;

    double data[Nstates * Nstates];
    double aMatrix[Nstates][Nstates];

    // double rho = 0.;
    int count = 0;
    for (int i = 0; i < Nstates; i++) {
        for (int j = i; j < Nstates; j++) {
            if (i == j) {
                aMatrix[i][j] = K[i] * (K[i] + 4) + (2. * Mass) * pow(rho * fermi, 2) * hPot[i][j]->Eval(rho);
                count++;
                data[i * Nstates + i] = aMatrix[i][j];
            } else {
                aMatrix[i][j] = (2. * Mass) * pow(rho * fermi, 2) * hPot[i][j]->Eval(rho);
                aMatrix[j][i] = aMatrix[i][j];
                count++;
                data[i * Nstates + j] = aMatrix[i][j];
                count++;
                data[j * Nstates + i] = aMatrix[i][j];
            }
            // data[count] = aMatrix[i][j];
            // count++;
        }
    }
    for (int i = 0; i < count; i++) {
        // if(data[i]<1.e-3) data[i] = 0.;
        cout << data[i] << "    ";
        if ((i + 1) % Nstates == 0)
            cout << endl << endl;
    }

    gsl_matrix_view m = gsl_matrix_view_array(data, Nstates, Nstates);
    gsl_vector* eval = gsl_vector_alloc(Nstates);
    gsl_matrix* evec = gsl_matrix_alloc(Nstates, Nstates);
    gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(Nstates);

    gsl_eigen_symmv(&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

    for (int i = 0; i < Nstates; i++) {
        double eval_i = gsl_vector_get(eval, i);
        gsl_vector_view evec_i = gsl_matrix_column(evec, i);
        cout << rho << "    " << i << "    " << eval_i << endl;
    }
}

void ExtractMatrix(double Kin, double JtotIN) {
    std::ifstream inputFile;
    inputFile.open(Form("States12With_K%i_J%1.1f.txt", int(Kin), JtotIN));

    double K[Nstates], S[Nstates], J[Nstates], Jz[Nstates];
    int l12[Nstates], l312[Nstates], L[Nstates], S12[Nstates];
    TGraph* hPot[Nstates][Nstates];
    TString Name[Nstates][Nstates];

    TFile* InPut = new TFile(Form("OutputPot_K%i_J%1.1f.root", int(Kin), JtotIN), "READ");
    InPut->cd();

    for (int i = 0; i < Nstates; i++) {
        inputFile >> K[i] >> l12[i] >> l312[i] >> S[i] >> L[i] >> S12[i] >> J[i] >> Jz[i];
        for (int j = i; j < Nstates; j++) {
            Name[i][j] = Form("Potential_Matrix_%i_%i", i, j);
            hPot[i][j] = (TGraph*)InPut->Get(Name[i][j].Data());
        }
    }
    inputFile.close();
    InPut->Close();

    // cout<<"CIAO"<<endl;

    double data[Nstates * Nstates];
    double aMatrix[Nstates][Nstates];

    TGraph* hU[Nstates];
    TGraph* hEV[Nstates][Nstates];

    double rho = 0.;
    for (int u = 0; u < Nsteps; u++) {
        double uoffset = 0.001;
        double stepsize = 0.01;
        if (rho <= 1.)
            rho = u * stepsize + uoffset;
        if (rho > 1. && rho < 4.) {
            stepsize = 0.1;
            rho = rho + stepsize;
            // MCsteps = 500;
        }
        if (rho >= 4.) {
            stepsize = 0.2;
            rho = rho + stepsize;
            // MCsteps = 100;
        }
        if (rho >= rhoMax)
            continue;
        // double rho = r*0.1;
        int count = 0;
        for (int i = 0; i < Nstates; i++) {
            for (int j = i; j < Nstates; j++) {
                if (i == j) {
                    aMatrix[i][j] = K[i] * (K[i] + 4) + (2. * Mass) * pow(rho * fermi, 2) * hPot[i][j]->Eval(rho);
                    count++;
                    data[i * Nstates + i] = aMatrix[i][j];
                } else {
                    aMatrix[i][j] = (2. * Mass) * pow(rho * fermi, 2) * hPot[i][j]->Eval(rho);
                    aMatrix[j][i] = aMatrix[i][j];
                    count++;
                    data[i * Nstates + j] = aMatrix[i][j];
                    count++;
                    data[j * Nstates + i] = aMatrix[i][j];
                }
                // data[count] = aMatrix[i][j];
                // count++;
            }
        }

        gsl_matrix_view m = gsl_matrix_view_array(data, Nstates, Nstates);
        gsl_vector* eval = gsl_vector_alloc(Nstates);
        gsl_matrix* evec = gsl_matrix_alloc(Nstates, Nstates);
        gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(Nstates);

        gsl_eigen_symmv(&m.matrix, eval, evec, w);
        gsl_eigen_symmv_free(w);
        gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

        for (int i = 0; i < Nstates; i++) {
            if (u == 0)
                hU[i] = new TGraph();
            double eval_i = gsl_vector_get(eval, i);
            gsl_vector_view evec_i = gsl_matrix_column(evec, i);
            hU[i]->SetPoint(u, rho, eval_i / pow(rho * fermi, 2) / (2. * Mass));
            // if(i==0){
            for (int j = 0; j < Nstates; j++) {
                if (u == 0)
                    hEV[i][j] = new TGraph();
                hEV[i][j]->SetPoint(u, rho, gsl_matrix_get(evec, j, i));
            }
            //}

            // if(i==0) cout<<rho<<"   "<<eval_i<<"    "<<gsl_matrix_get(evec,0,i)<<"    "<<gsl_matrix_get(evec,1,i)<<"
            // "<<gsl_matrix_get(evec,2,i)<<"    "<<gsl_matrix_get(evec,3,i)<<endl;
        }
        // Free allocated memory
        // gsl_matrix_free(m);
        // gsl_vector_free(eval);
        // gsl_matrix_free(evec);
        // gsl_eigen_symmv_free(w);
    }

    TFile* OutputPot = new TFile("UnPot.root", "RECREATE");
    OutputPot->cd();
    for (int i = 0; i < Nstates; i++) {
        hU[i]->RemovePoint(0);
        hU[i]->Write(Form("AdiabatiC_ch%i", i));
    }
    for (int i = 0; i < Nstates; i++) {
        for (int j = 0; j < Nstates; j++) {
            hEV[i][j]->Write(Form("Ch%i_%i", i, j));
            if (j >= i)
                hPot[i][j]->Write(Form("Potential_Matrix_%i_%i", i, j));
        }
    }
    OutputPot->Close();
}

std::complex<double> AdiabaticBasis(int n0, double Kin, double JtotIN, double rho, double phi12, double theta12,
                                    double phi312, double theta312, double phi) {
    std::ifstream inputFile;
    inputFile.open(Form("StatesWith_K%i_J%1.1f.txt", int(Kin), JtotIN));

    double K, S, J, Jz;
    int l12, l312, L, S12;

    std::complex<double> SHJ(0, 0);

    TGraph* hAd[Nstates];
    TFile* InputPot = new TFile("UnPot.root", "READ");
    InputPot->cd();
    for (int j = 0; j < Nstates; j++) {
        inputFile >> K >> l12 >> l312 >> S >> L >> S12 >> J >> Jz;
        hAd[j] = (TGraph*)InputPot->Get(Form("Ch%i_%i", n0, j));
        SHJ = SHJ + hAd[j]->Eval(rho) *
                        HypersphericalHarmonicsJJz(K, l12, l312, S, L, J, Jz, phi12, theta12, phi312, theta312, phi);
    }
    InputPot->Close();
    inputFile.close();

    return SHJ;
}
TGraph* IntegrateAB(int n0, int n1, double Kin, double JtotIN) {
    std::ifstream inputFile;
    inputFile.open(Form("StatesWith_K%i_J%1.1f.txt", int(Kin), JtotIN));

    double K[Nstates], S[Nstates], J[Nstates], Jz[Nstates];
    int l12[Nstates], l312[Nstates], L[Nstates], S12[Nstates];

    TGraph *hAd0[Nstates], *hAd1[Nstates];
    TFile* InputPot = new TFile("UnPot.root", "READ");
    InputPot->cd();
    for (int i = 0; i < Nstates; i++) {
        hAd0[i] = (TGraph*)InputPot->Get(Form("Ch%i_%i", n0, i));
        hAd1[i] = (TGraph*)InputPot->Get(Form("Ch%i_%i", n1, i));
        inputFile >> K[i] >> l12[i] >> l312[i] >> S[i] >> L[i] >> S12[i] >> J[i] >> Jz[i];
    }
    inputFile.close();
    InputPot->Close();

    int MCsteps = 500;
    TGraph* hHCP = new TGraph();

    double rho = 0.;
    for (int u = 0; u < Nsteps; u++) {
        double uoffset = 0.001;
        double stepsize = 0.1;
        if (rho <= 1.)
            rho = u * stepsize + uoffset;
        if (rho > 1. && rho < 4.) {
            stepsize = 0.1;
            rho = rho + stepsize;
            // MCsteps = 500;
        }
        if (rho >= 4.) {
            stepsize = 0.2;
            rho = rho + stepsize;
            // MCsteps = 100;
        }
        if (rho >= rhoMax)
            continue;

        std::complex<double> HCP(0, 0);

        for (int i = 0; i < MCsteps; i++) {
            double theta12 = RandomGen->Uniform(0, Pi);
            double theta312 = RandomGen->Uniform(0, Pi);
            double phi12 = RandomGen->Uniform(0, 2. * Pi);
            double phi312 = RandomGen->Uniform(0, 2. * Pi);
            double phi = RandomGen->Uniform(0, Pi / 2.);
            for (int j = 0; j < Nstates; j++) {
                for (int u = 0; u < Nstates; u++) {
                    if (abs(J[u] - J[j]) > 0.5 || abs(Jz[u] - Jz[j]) > 0.5)
                        continue;
                    if (abs(S[u] - S[j]) > 0.5 || abs(S12[u] - S12[j]) > 0.5)
                        continue;
                    HCP = HCP + conj(hAd1[u]->Eval(rho) * HypersphericalHarmonicsJJz(K[u], l12[u], l312[u], S[u], L[u],
                                                                                     J[u], Jz[u], phi12, theta12,
                                                                                     phi312, theta312, phi)) *
                                    hAd0[j]->Eval(rho) *
                                    HypersphericalHarmonicsJJz(K[j], l12[j], l312[j], S[j], L[j], J[j], Jz[j], phi12,
                                                               theta12, phi312, theta312, phi) *
                                    pow(sin(phi) * cos(phi), 2) * sin(theta12) * sin(theta312);
                }
            }
            ///    HCP = HCP + conj(AdiabaticBasis(n0, Kin, JtotIN, rho, phi12, theta12, phi312, theta312,
            ///    phi))*AdiabaticBasis(n1, Kin, JtotIN, rho, phi12, theta12, phi312, theta312,
            ///    phi)*pow(sin(phi)*cos(phi),2)*sin(theta12)*sin(theta312);
        }

        HCP = HCP / double(MCsteps) * pow(Pi, 5) * 2.;
        cout << rho << "    " << abs(HCP) << endl;
        hHCP->SetPoint(u, rho, abs(HCP));
    }
    return hHCP;
}

void SolveSEQ(int n0, double Q, TGraph*& hu) {
    // READ THE POT MATRIX
    TGraph *hAdS[Nstates], *hAdP;

    TFile* InputPot = new TFile("UnPot.root", "READ");
    InputPot->cd();
    hAdP = (TGraph*)InputPot->Get(Form("AdiabatiC_ch%i", n0));
    // hAdP=(TGraph*)InputPot->Get(Form("Potential_Matrix_%i_%i",n0,n0));
    for (int i = 0; i < Nstates; i++) {
        hAdS[i] = (TGraph*)InputPot->Get(Form("Ch%i_%i", n0, i));
    }
    InputPot->Close();

    cout << "Potential limit at rho = 40 fm: " << hAdP->Eval(40.) * (2. * Mass) * pow(40. * fermi, 2)
         << endl;  //+CoulombScreenedPot(50., 1.)<<endl;

    double u[Nsteps], f[Nsteps], rho;

    hu = new TGraph();
    double rho0 = 0.;
    u[0] = 0.;
    f[0] = 0.;
    hu->SetPoint(0, rho0, u[0]);

    double rho1 = 0.0001;
    u[1] = 1.e-5;
    f[1] = 15. / 4. / pow(rho1 * fermi, 2) + 2. * Mass * hAdP->Eval(rho1) - pow(Q, 2);
    hu->SetPoint(1, rho1, u[1]);
    // cout<<rho0<<"    "<<u[0]<<endl;
    // cout<<rho1<<"    "<<u[1]<<endl;

    for (int i = 2; i < Nsteps; i++) {
        rho = i * h;
        // if(rho>350) continue;
        // cout<<rho<<endl;
        f[i] = 15. / 4. / pow(rho * fermi, 2) + 2. * Mass * hAdP->Eval(rho) - pow(Q, 2);
        u[i] = 1. / (1 - pow(h * fermi, 2) / 12. * f[i]) *
               (2. * u[i - 1] * (1 + 5. * pow(h * fermi, 2) / 12. * f[i - 1]) -
                u[i - 2] * (1 - pow(h * fermi, 2) / 12. * f[i - 2]));
        hu->SetPoint(i, rho, u[i]);

        // cout<<rho<<"    "<<u[i]<<endl;
    }
}

/*
void diagonalise ()
{

double data[Nstates*Nstates];
int count = 0;
for (int i = 0; i < Nstates; i++)
{
    for (int j = i; j < Nstates; j++)
    {
        data[count] = aMatrix[i][j];
        count++;
    }

}

gsl_matrix_view m = gsl_matrix_view_array (data, Nstates, Nstates);
gsl_vector *eval = gsl_vector_alloc (Nstates);
gsl_matrix *evec = gsl_matrix_alloc (Nstates, Nstates);
gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (Nstates);

gsl_eigen_symmv (&m.matrix, eval, evec, w);
gsl_eigen_symmv_free (w);
//gsl_eigen_jacobi(&m.matrix, eval, evec);
//gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);




for (int i = 0; i < Nstates; i++)
{
    double eval_i = gsl_vector_get (eval, i);
    gsl_vector_view evec_i = gsl_matrix_column (evec, i);

    cout<<"eigenvalue = "<<eval_i<<endl;
    cout<<"eigenvector = "<<gsl_matrix_get(evec,0,i)<<" "<<gsl_matrix_get(evec,1,i)<<"
"<<gsl_matrix_get(evec,2,i)<<endl;

    printf ("eigenvalue = %g\n", eval_i);
    printf ("eigenvector = \n");
    gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
}


//return 0;
}
*/

void SaveStates12(double Kin, double JtotIN) {
    int count = 0;
    double Stot = 0.;

    std::ofstream outputFile(Form("States12With_K%i_J%1.1f.txt", int(Kin), JtotIN), std::ofstream::out);

    for (int K = Kin; K <= 150; K += 2) {
        int NoS = pow(K + 2, 2) * (K + 1) * (K + 3) / 12.;
        cout << "K = " << K << endl;
        int Nnu = int(K / 2.) + 1;

        cout << "Nst.  l12   l312  Stot   L   S12  |  Jtot    MJ |  PL    PS  | Keep? " << endl;

        for (int nu = 0; nu < Nnu; nu++) {
            int Ltot12 = K - 2 * nu;
            for (int j = 0; j < Ltot12 + 1; j++) {
                int l12 = j;
                if (l12 > 0)
                    continue;
                int l312 = Ltot12 - j;
                for (int L = abs(l12 - l312); L <= (l12 + l312); L++) {
                    for (int S12 = 0; S12 <= 1; S12++) {
                        // if(S12!=0) continue;
                        for (int cStot = 0; cStot <= int((S12 + 1. / 2.) - abs(S12 - 1. / 2.)); cStot++) {
                            Stot = abs(S12 - 1. / 2.) + cStot;
                            // if(Stot>0.5) continue;
                            for (int cJtot = 0; cJtot <= int(Stot + L - abs(Stot - L)); cJtot++) {
                                double Jtot = abs(Stot - L) + cJtot;
                                if (Jtot != JtotIN)
                                    continue;
                                // for (int cMJ = 0; cMJ < int(2.*Jtot+1); cMJ++)
                                //{
                                double MJ = -Jtot;  // + cMJ;
                                double ParityL = pow(-1, l12);
                                double ParityS = pow(-1, S12 + 1);
                                TString Result = "YES";
                                if (ParityL * ParityS > 0)
                                    Result = "NO";
                                // if(S12==1) Result = "NO";
                                if (Result == "NO")
                                    continue;
                                // if(Parity>0&&S12==1) Result = "NO";
                                outputFile << K << " " << l12 << " " << l312 << " " << Stot << " " << L << " " << S12
                                           << " " << Jtot << " " << MJ << endl;
                                cout << count << "       " << l12 << "   " << l312 << "    " << Stot << "    " << L
                                     << "    " << S12 << "  |  " << Jtot << "    " << MJ << " | " << ParityL << "  "
                                     << ParityS << " | " << Result.Data() << endl;

                                count++;

                                //}
                            }
                        }
                    }
                }
            }
        }
    }

    outputFile.close();
}

void SaveStates23(double Kin, double JtotIN) {
    int count = 0;
    double Stot = 0.;

    std::ofstream outputFile(Form("States23With_K%i_J%1.1f.txt", int(Kin), JtotIN), std::ofstream::out);

    for (int K = Kin; K <= 3; K += 2) {
        int NoS = pow(K + 2, 2) * (K + 1) * (K + 3) / 12.;
        cout << "K = " << K << endl;
        int Nnu = int(K / 2.) + 1;

        cout << "Nst.  l12   l312  Stot   L   S12  |  Jtot    MJ |  PL    PS  | Keep? " << endl;

        for (int nu = 0; nu < Nnu; nu++) {
            int Ltot12 = K - 2 * nu;
            for (int j = 0; j < Ltot12 + 1; j++) {
                int l12 = j;
                int l312 = Ltot12 - j;
                for (int L = abs(l12 - l312); L <= (l12 + l312); L++) {
                    // if(L!=0) continue;
                    for (int S12 = 0; S12 <= 1; S12++) {
                        // if(S12!=0) continue;
                        for (int cStot = 0; cStot <= int((S12 + 1. / 2.) - abs(S12 - 1. / 2.)); cStot++) {
                            Stot = abs(S12 - 1. / 2.) + cStot;
                            // if(Stot>0.5) continue;
                            for (int cJtot = 0; cJtot <= int(Stot + L - abs(Stot - L)); cJtot++) {
                                double Jtot = abs(Stot - L) + cJtot;
                                if (Jtot != JtotIN)
                                    continue;
                                // for (int cMJ = 0; cMJ < int(2.*Jtot+1); cMJ++)
                                //{
                                double MJ = -Jtot;  // + cMJ;
                                double ParityL = pow(-1, l12);
                                double ParityS = pow(-1, S12 + 1);
                                TString Result = "YES";
                                if (ParityL * ParityS > 0)
                                    Result = "NO";
                                // if(S12==1) Result = "NO";
                                // if(Result=="NO") continue;
                                // if(Parity>0&&S12==1) Result = "NO";
                                outputFile << K << " " << l12 << " " << l312 << " " << Stot << " " << L << " " << S12
                                           << " " << Jtot << " " << MJ << endl;
                                cout << count << "       " << l12 << "   " << l312 << "    " << Stot << "    " << L
                                     << "    " << S12 << "  |  " << Jtot << "    " << MJ << " | " << ParityL << "  "
                                     << ParityS << " | " << Result.Data() << endl;

                                count++;

                                //}
                            }
                        }
                    }
                }
            }
        }
    }

    outputFile.close();
}

int ComputeStates(double K) {
    int count = 0;
    int countTOT = 0;
    double Stot = 0.;

    int NoS = pow(K + 2, 2) * (K + 1) * (K + 3) / 12.;
    // cout<<"K = "<<K<<endl;
    int Nnu = int(K / 2.) + 1;

    // cout<<"Nst.  l12   l312  Stot   L   S12  |  Jtot    MJ |  PL    PS  | Keep? "<<endl;

    for (int nu = 0; nu < Nnu; nu++) {
        int Ltot12 = K - 2 * nu;
        for (int j = 0; j < Ltot12 + 1; j++) {
            int l12 = j;
            int l312 = Ltot12 - j;
            for (int L = abs(l12 - l312); L <= (l12 + l312); L++) {
                for (int S12 = 0; S12 <= 1; S12++) {
                    for (int cStot = 0; cStot <= int((S12 + 1. / 2.) - abs(S12 - 1. / 2.)); cStot++) {
                        Stot = abs(S12 - 1. / 2.) + cStot;
                        for (int cJtot = 0; cJtot <= int(Stot + L - abs(Stot - L)); cJtot++) {
                            double Jtot = abs(Stot - L) + cJtot;
                            for (int cMJ = 0; cMJ < int(2. * Jtot + 1); cMJ++) {
                                double MJ = -Jtot + cMJ;
                                double ParityL = pow(-1, l12);
                                double ParityS = pow(-1, S12 + 1);
                                TString Result = "YES";
                                if (ParityL * ParityS > 0)
                                    Result = "NO";
                                if (Result == "YES")
                                    count++;
                                //                                    cout<<countTOT<<"       "<<l12<<"   "<<l312<<"
                                //                                    "<<Stot<<"    "<<L<<"    "<<S12<<"  |  "<<Jtot<<"
                                //                                    "<<MJ<<" | "<<ParityL<<"  "<<ParityS<<" |
                                //                                    "<<Result.Data()<<endl;

                                countTOT++;
                            }
                        }
                    }
                }
            }
        }
    }

    //    cout<<"Total expected states: "<<NoS*8.<<endl;
    //    cout<<"Correct symmetry states: "<<count<<endl;

    return count;
}
/*
int PrintNStates(double K, double JtotIN){

    int count = 0;
    int countJ = 0;
    int countTOT = 0;
    double Stot = 0.;



        int NoS = pow(K+2,2)*(K+1)*(K+3)/12.;
        //cout<<"K = "<<K<<endl;
        int Nnu = int(K/2.) + 1;

        //cout<<"Nst.  l12   l312  Stot   L   S12  |  Jtot    MJ |  PL    PS  | Keep? "<<endl;

        for (int nu = 0; nu < Nnu; nu++)
        {
            int Ltot12 = K - 2*nu;
            for (int j = 0; j < Ltot12+1; j++)
            {
                int l12 = j;
                int l312 = Ltot12 - j;
                for (int L = abs(l12-l312); L <= (l12+l312); L++)
                {
                    for (int S12 = 0; S12 <= 1; S12++)
                    {
                        for (int cStot = 0; cStot <= int((S12+1./2.)-abs(S12-1./2.)); cStot++)
                        {
                            Stot = abs(S12-1./2.) + cStot;
                            for (int cJtot = 0; cJtot <= int(Stot+L -abs(Stot-L)); cJtot++)
                            {
                                double Jtot = abs(Stot-L) + cJtot;
                                for (int cMJ = 0; cMJ < int(2.*Jtot+1); cMJ++)
                                {
                                    double MJ = -Jtot + cMJ;
                                    double ParityL = pow(-1,l12);
                                    double ParityS = pow(-1,S12+1);
                                    TString Result = "YES";
                                    if(ParityL*ParityS>0) Result = "NO";
                                    if(Result=="YES") count++;
                                    if(Jtot==JtotIN&&MJ==-JtotIN&&Result=="YES") countJ++;
                                    //cout<<countTOT<<"       "<<l12<<"   "<<l312<<"    "<<Stot<<"    "<<L<<" "<<S12<<"
|  "<<Jtot<<"    "<<MJ<<" | "<<ParityL<<"  "<<ParityS<<" | "<<Result.Data()<<endl;

                                    countTOT++;
                                }

                            }

                        }
                    }

                }

            }

        }


    //cout<<"Total expected states: "<<NoS*8.<<endl;
    //cout<<"Correct symmetry states: "<<count<<endl;
    //cout<<"Total number of states with J="<<JtotIN<<": "<<countJ<<endl;

    return countJ;

}


TGraph* UfreeExact(int K, double Q){

    double u[Nsteps],rho;

    TGraph *hu = new TGraph();
//pow(i,L) --> L=0 is 1, L=1 is i, L=2 is -1, L=3 is -i, L=4 is 1, L=5 is i
// if(L is pari) sign=Cos(L*Pi/2)
// if(L is dispari) sign = Sin(L*Pi/2)
    int sign = int(cos(K*TMath::Pi()/2.));
    if(sign==0) sign = int(sin(K*TMath::Pi()/2.));
    for (int i = 0; i < Nsteps; i++)
    {
        rho = i*h;
        u[i] = sign*sqrt(Q*rho*fermi)*gsl_sf_bessel_Jn(K+2, Q*rho*fermi);
        hu -> SetPoint(i, rho, u[i]);
    }

    return hu;

}

TGraph *hNULL = new TGraph();

void NullHisto(){

    hNULL -> SetPoint(0,0,0);
    hNULL -> SetPoint(1,100,0);

}
*/
void RescaleGraph(TGraph*& graph, double xScale, double yScale) {
    int nPoints = graph->GetN();

    for (int i = 0; i < nPoints; ++i) {
        double x, y;
        graph->GetPoint(i, x, y);

        // Rescale x and y coordinates
        x *= xScale;
        y *= yScale;

        // Set the modified point back to the graph
        graph->SetPoint(i, x, y);
    }
}

/*
const int NcoupledEq = 1;
TGraph *hPot[NcoupledEq][NcoupledEq];

void OpenPotMatrix(int Kin, double Jin){

    NullHisto();

    TFile *InputPot = new TFile(Form("OutputPot_K%i_J%1.1f.root",Kin,Jin),"READ");
    for (int i = 0; i < NcoupledEq; i++)
    {
        for (int j = i; j < NcoupledEq; j++)
        {
//            hPot[i][j] = new TGraph();
//            hPot[j][i] = new TGraph();
//            hPot[i][j] -> SetPoint(0,0,0);
//            hPot[i][j] -> SetPoint(1,100,0);
//            hPot[j][i] -> SetPoint(0,0,0);
//            hPot[j][i] -> SetPoint(1,100,0);

            hPot[i][j] = (TGraph*)InputPot->Get(Form("Potential_Matrix_%i_%i",i,j));
            //hPot[i][j] ->SetPoint(hPot[i][j] -> GetN(),20,0.);
            //hPot[i][j] ->SetPoint(hPot[i][j] -> GetN()+1,500,0.);
            hPot[j][i] = (TGraph*)InputPot->Get(Form("Potential_Matrix_%i_%i",i,j));
            //hPot[j][i] ->SetPoint(hPot[j][i] -> GetN(),20,0.);
            //hPot[j][i] ->SetPoint(hPot[j][i] -> GetN()+1,500,0.);
        //    if(i!=j) RescaleGraph(hPot[i][j],1.,-1.);
        //    if(i!=j) RescaleGraph(hPot[j][i],1.,-1.);
        }
    }
    InputPot->Close();
}
*/
// Function to perform numerical integration using the trapezoidal rule
long double IntegrateGraph(TGraph* graph, double x_min, double x_max) {
    int nPoints = graph->GetN();
    double* x = graph->GetX();
    double* y = graph->GetY();

    // Find the range within the interval [x_min, x_max]
    std::vector<double> x_in_range;
    std::vector<double> y_in_range;

    for (int i = 0; i < nPoints; ++i) {
        if (x[i] >= x_min && x[i] <= x_max) {
            x_in_range.push_back(x[i]);
            y_in_range.push_back(y[i]);
        }
    }

    // Perform trapezoidal integration
    long double integral = 0.0;
    for (size_t i = 0; i < x_in_range.size() - 1; ++i) {
        double dx = x_in_range[i + 1] - x_in_range[i];
        double avg_y = 0.5 * (y_in_range[i] + y_in_range[i + 1]);
        integral += dx * avg_y;
    }

    return integral;
}

/*
void SolveSEQ(int Kin, double JtotIN, double Q, TGraph *hu[]){

    // READ THE POT MATRIX
    OpenPotMatrix(Kin, JtotIN);


    double K[NcoupledEq],Stot[NcoupledEq],Jtot[NcoupledEq],MJ[NcoupledEq];
    int l12[NcoupledEq],l312[NcoupledEq],L[NcoupledEq],S12[NcoupledEq];


    double u[Nsteps][NcoupledEq],f[Nsteps][NcoupledEq],sterm[Nsteps][NcoupledEq],rho, rho0, rho1;
    double u0[Nsteps][NcoupledEq];

    TGraph *hu0[NcoupledEq];
    TGraph *hUFreeExact[NcoupledEq];
    double MaxFree[NcoupledEq];

    std::ifstream inputFile;
    inputFile.open(Form("StatesWith_K%i_J%1.1f.txt",int(Kin),JtotIN));
    double NoS[NcoupledEq];

    for (int i = 0; i < NcoupledEq; i++)
    {
        hu[i] = new TGraph();
        inputFile>>K[i]>>l12[i]>>l312[i]>>Stot[i]>>L[i]>>S12[i]>>Jtot[i]>>MJ[i];

        hu0[i] = UfreeExact(K[i], Q);

        rho0 = 0.;
        u[0][i] = 0.;
        f[0][i] = 0.;
        sterm[0][i] = 0.;
        hu[i] -> SetPoint(0, rho0, u[0][i]);

        //NoS[i] = pow(K[i]+2,2)*(K[i]+1)*(K[i]+3)/12.;
        NoS[i] = PrintNStates(K[i],1./2.);//pow(K[i]+2,2)*(K[i]+1)*(K[i]+3)/12.;

        rho1 = 0.0001;
        u[1][i] = hu0[i] -> Eval(rho1);
        f[1][i] = (K[i]+3./2.)*(K[i]+5./2.)/pow(rho1*fermi,2) + 2.*Mass*hPot[i][i]->Eval(rho1) - pow(Q,2);

    }


    int sign0[NcoupledEq],sign1[NcoupledEq];
    double Parity = 0;
    for (int i = 0; i < NcoupledEq; i++)
    {
        sterm[1][i] = 0.;
        Parity = pow(-1,K[i]);
        if(Parity>0){
            sign0[i] = pow(-1,K[i]/2.);
            for (int j = 0; j < NcoupledEq; j ++)
            {
                sign1[j] = pow(-1,K[j]/2.);
                if(j==i) continue;
                sterm[1][i] = sterm[1][i] + 2.*NoS[j]/NoS[i]*Mass*hPot[i][j]->Eval(rho1)*u[1][j];
            }
        }
        if(Parity<0){
            sign0[i] = pow(-1,(K[i]-1.)/2.);
            for (int j = 0; j < NcoupledEq; j ++)
            {
                sign1[j] = pow(-1,(K[j]-1.)/2.);
                if(j==i) continue;
                sterm[1][i] = sterm[1][i] + 2.*NoS[j]/NoS[i]*Mass*hPot[i][j]->Eval(rho1)*u[1][j];
            }
        }
        hu[i] -> SetPoint(1, rho1, u[1][i]);
    }

    bool condition[NcoupledEq];
    for (int i = 2; i < Nsteps; i++)
    {
        rho = i*h;
        //if(rho>10) continue;
        //cout<<rho<<endl;
        for (int l = 0; l < NcoupledEq; l++)
        {
            f[i][l] = (K[l]+3./2.)*(K[l]+5./2.)/pow(rho*fermi,2) + 2.*Mass*hPot[l][l]->Eval(rho) - pow(Q,2);
            u0[i][l] = hu0[l] ->
Eval(rho);// 1./(1-pow(h*fermi,2)/12.*f[i][l])*(2.*u[i-1][l]*(1+5.*pow(h*fermi,2)/12.*f[i-1][l]) -
u[i-2][l]*(1-pow(h*fermi,2)/12.*f[i-2][l])); condition[l] = true;
        }
        bool TotalCondition = true;
        int lN =0;
        do{
            for (int l = 0; l < NcoupledEq; l++)
            {
                sterm[i][l] = 0.;
                for (int j = 0; j < NcoupledEq; j++)
                {
                    if(j>0) u[i][j] = hu0[j]->Eval(rho);
                    if(j==l) continue;
                    if(j<l){
                        sterm[i][l] = sterm[i][l] + 2.*NoS[j]/NoS[l]*Mass*hPot[l][j]->Eval(rho)*u[i][j]; //(NoS[j]*8.);
                    }else if(j>l){
                        sterm[i][l] = sterm[i][l] + 2.*NoS[j]/NoS[l]*Mass*hPot[l][j]->Eval(rho)*u[i][j]; //(NoS[j]*8.);
                    }
                }

                if(l==0) u[i][l] = 1./(1-pow(h*fermi,2)/12.*f[i][l])*(2.*u[i-1][l]*(1+5.*pow(h*fermi,2)/12.*f[i-1][l]) -
u[i-2][l]*(1-pow(h*fermi,2)/12.*f[i-2][l])+ pow(h*fermi,2)/12.*(sterm[i][l]+10.*sterm[i-1][l]+sterm[i-2][l]));

                if(abs((u[i][l] - u0[i][l])/u[i][l])<1.e-10&&condition[l]){
                    condition[l] = false;
                    TotalCondition = false;
                }
                u0[i][l] = u[i][l];
            }
            for (int l = 0; l < NcoupledEq; l++)
            {
                if(condition[l]) TotalCondition = true;
            }
            for (int l = 0; l < NcoupledEq; l++)
            {
                if(!TotalCondition) hu[l] -> SetPoint(i, rho, u[i][l]);
            }
        }while(TotalCondition);

    }

}


void SolveCoupledSEQ(int Kin, double JtotIN, double Q, TGraph *hu[]){

    // READ THE POT MATRIX
    OpenPotMatrix(Kin, JtotIN);


    double K[NcoupledEq],Stot[NcoupledEq],Jtot[NcoupledEq],MJ[NcoupledEq];
    int l12[NcoupledEq],l312[NcoupledEq],L[NcoupledEq],S12[NcoupledEq];


    double u[Nsteps][NcoupledEq],f[Nsteps][NcoupledEq],sterm[Nsteps][NcoupledEq],rho, rho0, rho1;
    double u0[Nsteps][NcoupledEq];

    TGraph *hu0[NcoupledEq];
    TGraph *hUFreeExact[NcoupledEq];
    double MaxFree[NcoupledEq];

    std::ifstream inputFile;
    inputFile.open(Form("StatesWith_K%i_J%1.1f.txt",int(Kin),JtotIN));
    double NoS[NcoupledEq];

    SolveSEQ(0, 1./2., Q, hu0);


    for (int i = 0; i < NcoupledEq; i++)
    {
        hu[i] = new TGraph();
        inputFile>>K[i]>>l12[i]>>l312[i]>>Stot[i]>>L[i]>>S12[i]>>Jtot[i]>>MJ[i];

        //hu0[i] = UfreeExact(K[i], Q);

        rho0 = 0.;
        u[0][i] = 0.;
        f[0][i] = 0.;
        sterm[0][i] = 0.;
        hu[i] -> SetPoint(0, rho0, u[0][i]);

        NoS[i] = PrintNStates(K[i],1./2.);//pow(K[i]+2,2)*(K[i]+1)*(K[i]+3)/12.;

        rho1 = 0.0001;
        u[1][i] = hu0[i] -> Eval(rho1);
        f[1][i] = (K[i]+3./2.)*(K[i]+5./2.)/pow(rho1*fermi,2) + 2.*Mass*hPot[i][i]->Eval(rho1) - pow(Q,2);

    }


    int sign0[NcoupledEq],sign1[NcoupledEq];
    double Parity = 0;
    for (int i = 0; i < NcoupledEq; i++)
    {
        sterm[1][i] = 0.;
        Parity = pow(-1,K[i]);
        if(Parity>0){
            sign0[i] = pow(-1,K[i]/2.);
            for (int j = 0; j < NcoupledEq; j ++)
            {
                sign1[j] = pow(-1,K[j]/2.);
                if(j==i) continue;
                sterm[1][i] = sterm[1][i] + 2.*NoS[j]/NoS[i]*Mass*hPot[i][j]->Eval(rho1)*u[1][j];
            }
        }
        if(Parity<0){
            sign0[i] = pow(-1,(K[i]-1.)/2.);
            for (int j = 0; j < NcoupledEq; j ++)
            {
                sign1[j] = pow(-1,(K[j]-1.)/2.);
                if(j==i) continue;
                sterm[1][i] = sterm[1][i] + 2.*NoS[j]/NoS[i]*Mass*hPot[i][j]->Eval(rho1)*u[1][j];
            }
        }
        hu[i] -> SetPoint(1, rho1, u[1][i]);
    }

    bool condition[NcoupledEq];
    for (int i = 2; i < Nsteps; i++)
    {
        rho = i*h;
        //if(rho>10) continue;
        //cout<<rho<<endl;
        for (int l = 0; l < NcoupledEq; l++)
        {
            f[i][l] = (K[l]+3./2.)*(K[l]+5./2.)/pow(rho*fermi,2) + 2.*Mass*hPot[l][l]->Eval(rho) - pow(Q,2);
            u0[i][l] = hu0[l] ->
Eval(rho);// 1./(1-pow(h*fermi,2)/12.*f[i][l])*(2.*u[i-1][l]*(1+5.*pow(h*fermi,2)/12.*f[i-1][l]) -
u[i-2][l]*(1-pow(h*fermi,2)/12.*f[i-2][l])); condition[l] = true;
        }
        bool TotalCondition = true;
        int lN =0;
        do{
            for (int l = 0; l < NcoupledEq; l++)
            {
                sterm[i][l] = 0.;
                for (int j = 0; j < NcoupledEq; j++)
                {
                    if(j<l){
                        sterm[i][l] = sterm[i][l] + 2.*NoS[j]/NoS[l]*Mass*hPot[l][j]->Eval(rho)*u[i][j]; //(NoS[j]*8.);
                    }else if(j>l){
                        sterm[i][l] = sterm[i][l] + 2.*NoS[j]/NoS[l]*Mass*hPot[l][j]->Eval(rho)*u[i-1][j];
//(NoS[j]*8.);
                    }
                }

                u[i][l] = 1./(1-pow(h*fermi,2)/12.*f[i][l])*(2.*u[i-1][l]*(1+5.*pow(h*fermi,2)/12.*f[i-1][l]) -
u[i-2][l]*(1-pow(h*fermi,2)/12.*f[i-2][l])+ pow(h*fermi,2)/12.*(sterm[i][l]+10.*sterm[i-1][l]+sterm[i-2][l]));

                if(abs((u[i][l] - u0[i][l])/u[i][l])<1.e-10&&condition[l]){
                    condition[l] = false;
                    TotalCondition = false;
                }
                u0[i][l] = u[i][l];
            }
            for (int l = 0; l < NcoupledEq; l++)
            {
                if(condition[l]) TotalCondition = true;
            }
            for (int l = 0; l < NcoupledEq; l++)
            {
                if(!TotalCondition)	hu[l] -> SetPoint(i, rho, u[i][l]);
            }
        }while(TotalCondition);

    }
}
*/
TGraph* AntisymmetricFreeWF() {
    TGraph* hPsi2FreeAS;

    hPsi2FreeAS = new TGraph();

    RandomGen->SetSeed(4);

    int Nsteps2 = 20000;

    double thetar1, thetar2;
    double phir1, phir2;
    double phir;
    double thetaq1, thetaq2;
    double phiq1, phiq2;
    double phiq;
    double CosTh, CosAlpha;
    double Test1;
    double Norm;
    double differential, differential1, differential2;

    int Nmax = 5000;
    double xOffset = 0.001;

    for (int j = 0; j < Nmax; j++) {
        Test1 = 0.;
        Norm = 0.;
        double x = 0.1 * j + xOffset;
        for (int i = 0; i < Nsteps2; i++) {
            thetar1 = RandomGen->Uniform(0, Pi);
            thetar2 = RandomGen->Uniform(0, Pi);
            phir = RandomGen->Uniform(0, Pi / 2.);
            thetaq1 = RandomGen->Uniform(0, Pi);
            thetaq2 = RandomGen->Uniform(0, Pi);
            phiq = RandomGen->Uniform(0, Pi / 2.);

            differential = pow(sin(phir) * cos(phir), 2) * sin(thetar1) * sin(thetar2) * pow(sin(phiq) * cos(phiq), 2) *
                           sin(thetaq1) * sin(thetaq2);
            Test1 = Test1 - 1. / 2. * cos(2. * x * fermi * cos(phir) * cos(phiq) * (cos(thetar1) * cos(thetaq1))) *
                                differential;

            Norm = Norm + differential;
        }
        Test1 = Test1 / Norm;

        hPsi2FreeAS->SetPoint(j, x, 1 + Test1);
    }

    return hPsi2FreeAS;
}

TGraph* HigherPW(int K0, int Kmax) {
    TGraph* HPW = new TGraph();
    int Nmax = 10000;
    double xOffset = 0.0001;
    double* NK = new double[Kmax];

    for (int i = 0; i < Nmax; i++) {
        double x = 0.1 * i + xOffset;
        double Sum = 0.;
        for (int K = K0; K < Kmax; K++) {
            if (i == 0)
                NK[K] = ComputeStates(K);
            Sum = Sum + pow(gsl_sf_bessel_Jn(K + 2, x), 2) * NK[K];
        }
        Sum = Sum * 1. / 4. * pow(2., 6) / pow(x, 4);
        HPW->SetPoint(i, x, Sum);
    }

    delete[] NK;
    return HPW;
}

/*
TGraph* WaveFunction(double Q3){

    TGraph *hWF = new TGraph();
    TGraph *hWF00 = new TGraph();
    TGraph *HPW = HigherPW(1, 30);
    TGraph *hWF0 = HigherPW(0, 1);
    int Nmax = 100;
    double rhoOffset = 0.001;

    TGraph *u0[NcoupledEq];
    SolveCoupledSEQ(0,1./2.,Q3,u0);


    double WF00 = 0.;
    for (int j = 0; j < Nsteps; j++)
    {
        double rho = j*h + rhoOffset;
        WF00 = 1./4.*pow(2.,6)/pow(Q3*rho*fermi,4)*2.*pow(u0[0]->Eval(rho),2)/(Q3*rho*fermi);

        hWF00 -> SetPoint(j,rho,WF00);
    }

    double MaxFree = IntegrateGraph(hWF0,30,50);
    double MaxNotFree = IntegrateGraph(hWF00,30,50);
    RescaleGraph(hWF00, 1., MaxFree/MaxNotFree);


    double WF = 0.;
    for (int j = 0; j < Nsteps; j++)
    {
        double rho = j*h + rhoOffset;
//        WF = WF + 1./4.*pow(2.,6)/(pow(Q3,4)*pow(rho0*fermi,6))*rho*h*fermi*fermi*exp(-pow(rho/rho0,2))*
        WF = hWF00->Eval(rho) + HPW->Eval(Q3*rho*fermi);
        //WF0 = 1./4.*pow(2.,6)/pow(Q3*rho*fermi,4)*2.*pow(u0[0]->Eval(rho),2)/(Q3*rho*fermi) + HPW->Eval(Q3*rho*fermi);

        hWF -> SetPoint(j,rho,WF);
    }

    return hWF;

}
*/
TGraph* WaveFunction2(double Q3, TGraph* HPW, TGraph* hWF0) {
    TGraph* hWF = new TGraph();
    TGraph* hWF00 = new TGraph();
    TGraph* hWF0test = new TGraph();
    int Nmax = 100;
    double rhoOffset = 0.0001;

    TGraph *u0, *u[Nstates];
    cout << "Start solving SEQ" << endl;
    //    SolveCoupledSEQ(0,1./2.,Q3,u0);
    //    for (int i = 0; i < Nstates; i++)
    for (int i = 0; i < 1; i++) {
        SolveSEQ(i, Q3, u[i]);
    }

    TGraph* hAdS[Nstates];

    TFile* InputPot = new TFile("UnPot.root", "READ");
    InputPot->cd();
    for (int i = 0; i < Nstates; i++) {
        hAdS[i] = (TGraph*)InputPot->Get(Form("Ch%i_%i", i, i));
    }
    InputPot->Close();

    // SolveSEQ(1,Q3,u0);
    // RescaleGraph(u0, 1., 1.e-282);

    double WF00 = 0.;
    double WF0 = 0.;
    for (int j = 0; j < Nsteps; j++) {
        double rho = j * h + rhoOffset;
        if (rho > rhoMax)
            continue;
        double Utot = 0.;
        for (int i = 0; i < 1; i++) {
            Utot = Utot + pow(u[i]->Eval(rho), 2);
        }
        // WF00 = 1./4.*pow(2.,6)/pow(Q3*rho*fermi,4)*2.*pow(u0->Eval(rho),2)/(Q3*rho*fermi);
        WF00 = 1. / 4. * pow(2., 6) / pow(Q3 * rho * fermi, 4) * 2. * Utot / (Q3 * rho * fermi);
        hWF00->SetPoint(j, rho, WF00);
        hWF0test->SetPoint(j, rho, hWF0->Eval(rho * Q3 * fermi));
    }

    double MaxFree = IntegrateGraph(hWF0test, rhoMin, rhoMax);
    double MaxNotFree = IntegrateGraph(hWF00, rhoMin, rhoMax);
    RescaleGraph(hWF00, 1., MaxFree / MaxNotFree);
    cout << "Look here: " << MaxFree << " has to be equal to " << MaxNotFree << endl;

    double WF = 0.;
    for (int j = 0; j < Nsteps; j++) {
        double rho = j * h + rhoOffset;
        if (rho > rhoMax)
            continue;
        //        WF = WF + 1./4.*pow(2.,6)/(pow(Q3,4)*pow(rho0*fermi,6))*rho*h*fermi*fermi*exp(-pow(rho/rho0,2))*
        WF = hWF00->Eval(rho);  //+ HPW->Eval(Q3*rho*fermi);
        // WF0 = 1./4.*pow(2.,6)/pow(Q3*rho*fermi,4)*2.*pow(u0[0]->Eval(rho),2)/(Q3*rho*fermi) +
        // HPW->Eval(Q3*rho*fermi);

        hWF->SetPoint(j, rho, WF);
    }

    return hWF;
}

TGraph* CorrelationFunction(double rho0) {
    double Nmax = 50;
    double rhoOffset = 0.001;
    double Qoffset = 0.;

    TGraph* HPW = HigherPW(1, 40);
    TGraph* hWF0 = HigherPW(0, 1);
    //    TGraph *hAB = IntegrateAB(0, 0, 0, 1./2.);

    TGraph* hCF = new TGraph();
    TGraph* hCF0 = new TGraph();

    TFile* OutWF = new TFile("WF.root", "RECREATE");
    for (int i = 1; i < Nmax; i++) {
        double Q3 = i * 4 + Qoffset;
        double CF = 0.;
        double CF0 = 0.;
        cout << "Q = " << Q3 << " MeV/c" << endl;
        TGraph* WF = WaveFunction2(Q3, HPW, hWF0);
        TGraph* hWF0test = new TGraph();
        OutWF->cd();
        WF->Write(Form("Q3%i", int(Q3)));
        for (int j = 0; j < Nsteps; j++) {
            double rho = j * h + rhoOffset;
            if (rho > rhoMax)
                continue;
            CF = CF + WF->Eval(rho) * pow(rho, 5) * h / pow(rho0, 6) * exp(-pow(rho / rho0, 2));
            hWF0test->SetPoint(j, rho, hWF0->Eval(rho * Q3 * fermi));
            CF0 = CF0 + hWF0->Eval(rho * Q3 * fermi) * pow(rho, 5) * h / pow(rho0, 6) * exp(-pow(rho / rho0, 2));
        }
        hWF0test->Write(Form("Q300%i", int(Q3)));
        hCF->SetPoint(i, Q3 * sqrt(6.35) / 1000., CF);
        hCF0->SetPoint(i, Q3 * sqrt(6.35) / 1000., CF0);
    }
    hCF->Write("CF00");
    hCF0->Write("CF0");
    OutWF->Close();
    return hCF;
}


double rp = 1;
double rs = 2;
int nRho = 8;
double maxRho = 8;
int nIter = 15;

double SandwitchSource(double hypRad, int l12, int l312, int nu, int L, int M, int nIter) {

    // int l12, int l312, int nu, int L, int M, double phi12, double theta12,
    //                                          double phi312, double theta312, double phi

    complex<double> sum = 0;

    double dOmega = (std::numbers::pi / nIter)         // dtheta12
                    * (std::numbers::pi / nIter)       // dtheta312
                    * (2 * std::numbers::pi / nIter)   // dphi12
                    * (2 * std::numbers::pi / nIter)   // dphi312
                    * (std::numbers::pi / 2 / nIter);  // dphi

    // dOmega = (std::numbers::pi / 2 / nIter);  // dphi
                    
    for (int iTheta312 = 0; iTheta312 < nIter; iTheta312++) {
        double theta312 = std::numbers::pi * double(iTheta312) / nIter;

        for (int iTheta12 = 0; iTheta12 < nIter; iTheta12++) {
            double theta12 = std::numbers::pi * double(iTheta12) / nIter;

            for (int iPhi12 = 0; iPhi12 < nIter; iPhi12++) {
                double phi12 = 2 * std::numbers::pi * double(iPhi12) / nIter;

                for (int iPhi312 = 0; iPhi312 < nIter; iPhi312++) {
                    double phi312 = 2 * std::numbers::pi * double(iPhi312) / nIter;

                    for (int iPhi = 0; iPhi < nIter; iPhi++) {
                        double phi = std::numbers::pi / 2 * double(iPhi) / nIter;

                        complex<double> hh1 =
                            HypersphericalHarmonics(l12, l312, nu, L, M, phi12, theta12, phi312, theta312, phi);
                        complex<double> hh2 =
                            HypersphericalHarmonics(l12, l312, nu, L, M, phi12, theta12, phi312, theta312, phi);

                        double measure = sin(theta12) * sin(theta312) * pow(cos(phi), 2) * pow(sin(phi), 2) * std::pow(hypRad, 5);

                        sum += conj(hh1) * hh2 * measure * dOmega * _SourceAAAppr(hypRad, phi, rp, rs) * pow(numbers::pi, 3);
                    }
                }
            }
        }
    }


    std::cout << "sum" << sum << std::endl;
    return sum.real();
}

void Sandwitch(int l12, int l312, int nu, int L, int M, int nIter) {
    // int l12, int l312, int nu, int L, int M, double phi12, double theta12,
    //                                          double phi312, double theta312, double phi

    complex<double> sum = 0;

    double dOmega = (std::numbers::pi / nIter)         // dtheta12
                    * (std::numbers::pi / nIter)       // dtheta312
                    * (2 * std::numbers::pi / nIter)   // dphi12
                    * (2 * std::numbers::pi / nIter)   // dphi312
                    * (std::numbers::pi / 2 / nIter);  // dphi

    for (int iTheta312 = 0; iTheta312 < nIter; iTheta312++) {
        double theta312 = std::numbers::pi * double(iTheta312) / nIter;

        for (int iTheta12 = 0; iTheta12 < nIter; iTheta12++) {
            double theta12 = std::numbers::pi * double(iTheta12) / nIter;

            for (int iPhi12 = 0; iPhi12 < nIter; iPhi12++) {
                double phi12 = 2 * std::numbers::pi * double(iPhi12) / nIter;

                for (int iPhi312 = 0; iPhi312 < nIter; iPhi312++) {
                    double phi312 = 2 * std::numbers::pi * double(iPhi312) / nIter;

                    for (int iPhiIter = 0; iPhiIter < nIter; iPhiIter++) {
                        double phi = std::numbers::pi / 2 * double(iPhiIter) / nIter;

                        complex<double> hh1 =
                            HypersphericalHarmonics(l12, l312, nu, L, M, phi12, theta12, phi312, theta312, phi);
                        complex<double> hh2 =
                            HypersphericalHarmonics(l12, l312, nu, L, M, phi12, theta12, phi312, theta312, phi);

                        double measure = sin(theta12) * sin(theta312) * pow(cos(phi), 2) * pow(sin(phi), 2);

                        sum += conj(hh1) * hh2 * measure * dOmega;
                    }
                }
            }
        }
    }
    std::cout << "sum" << sum << std::endl;
}


void Hyperspherical_Harmonics() {
    // Sandwitch(0, 0, 0, 0, 0, 30);

    double *rhos = new double[nRho];
    double *cf = new double[nRho];

    for (int iRho = 0; iRho < nRho; iRho++) {   
        double rho = maxRho * double(iRho) / nRho;
        
        rhos[iRho] = rho;
        cf[iRho] = SandwitchSource(rho, 0, 0, 0, 0, 0, nIter);
    }
    
    TGraph *gCF = new TGraph(nRho, rhos, cf);
    
    TF1 *fSourceAAAppp = new TF1("fSourceAAAppp", SourceCountsAAA, 0, maxRho, 2);
    fSourceAAAppp->SetParameter(0, 1);
    fSourceAAAppp->SetParameter(1, 2 * rp);

    TF1 *fSourceAAAppr = new TF1("fSourceAAAppr", SourceCountsAAApprAvg, 0, maxRho, 3);
    fSourceAAAppr->SetParameter(0, 1);
    fSourceAAAppr->SetParameter(1, rp);
    fSourceAAAppr->SetParameter(2, rs);

    fSourceAAAppp->Draw("");
    fSourceAAAppr->Draw("same");
    gCF->Draw("lp");

}
