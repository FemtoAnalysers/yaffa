#ifndef YAFFA_WAVEFUNCTION_H
#define YAFFA_WAVEFUNCTION_H

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

// Tabulated squared modulus of a femtoscopic wave function:
//   |psi(k*, r*)|^2     for a two-body system, or
//   |Psi(Q3, rho)|^2    for a three-body system,
// plus the Koonin-Pratt integral that turns it into a correlation function.
//
// The table is rectangular: one value for every (momentum, radius) pair.
// Both axes are stored explicitly and may be non-uniformly spaced.
// Indexing convention, in memory and on disk: row = momentum, column = radius.
class WaveFunction {
   public:
    // momentum : k* (nBody == 2) or Q3 (nBody == 3), one entry per row
    // radius   : r* (nBody == 2) or hyper-radius rho (nBody == 3), one per column
    // values   : row-major, values[iMom * radius.size() + iRad]
    // nBody    : number of particles (2 or 3); sets the Koonin-Pratt Jacobian
    // system      : free-form label ("pp", "ppp", ...)
    // description : free-form notes (reference, assumptions, ...); may be
    //               multi-line. Written into the file header as comment lines for
    //               the human reader; not parsed back when a file is loaded.
    WaveFunction(std::vector<double> momentum, std::vector<double> radius,
                 std::vector<double> values, int nBody, std::string system = "",
                 std::string description = "");

    // Read the standardized text format written by Save().
    explicit WaveFunction(const std::string& filename);

    void Save(const std::string& filename) const;
    void Print() const;

    const std::vector<double>& Momentum() const { return fMomentum; }
    const std::vector<double>& Radius() const { return fRadius; }
    const std::vector<double>& Values() const { return fValues; }  // row-major, see At()
    int GetNBody() const { return fNBody; }
    const std::string& System() const { return fSystem; }
    const std::string& Description() const { return fDescription; }

    // |psi|^2 at grid node (iMom, iRad).
    double At(std::size_t iMom, std::size_t iRad) const {
        return fValues[iMom * fRadius.size() + iRad];
    }

    // Koonin-Pratt correlation function, evaluated at every stored momentum:
    //
    //   C(mom) = Int J(r) S(r) |psi(mom, r)|^2 dr  /  Int J(r) S(r) dr
    //
    // with J(r) = r^(3*nBody - 4), i.e. r^2 for two bodies and r^5 for three
    // (the hyper-radius volume element). S is the emission source; it need not
    // be normalized. The integral is a trapezoid sum over the stored radius
    // axis, so non-uniform spacing is handled correctly.
    //
    // Overload 1: S given as a callable S(r).
    // Overload 2: S given as its values already sampled on Radius().
    std::vector<double> CorrelationFunction(const std::function<double(double)>& source) const;
    std::vector<double> CorrelationFunction(const std::vector<double>& sourceOnRadiusAxis) const;

   private:
    void CheckShape() const;
    double Jacobian(double r) const { return std::pow(r, 3 * fNBody - 4); }

    std::vector<double> fMomentum;  // rows
    std::vector<double> fRadius;    // columns
    std::vector<double> fValues;    // row-major, size fMomentum.size() * fRadius.size()
    int fNBody;
    std::string fSystem;  // free-form label: "pp", "ppp", ...
    std::string fDescription;  // free-form notes; may be multi-line
};

#endif  // YAFFA_WAVEFUNCTION_H
