#include "WaveFunction.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Text format written by Save() / read by the filename constructor
// ---------------------------------------------------------------------------
//
//   # yaffa.WaveFunction 1           <- format magic + version (must be line 1)
//   #
//   # ... free-text comment lines: the quantity, the two axes (name, unit,
//   #     point count, range, spacing), the table order, and any notes
//   #     passed as `description`, so the file is self-explanatory. Ignored on read.
//   #
//   # system = pp                    <- parsed: free-form system label
//   # potential = AV18               <- parsed: interaction potential
//   # waves = spd                    <- parsed: partial waves included
//   # nbody = 2                      <- parsed: 2 or 3 (sets axis meaning + KP Jacobian)
//   #
//   mom\radius  <rad_0> <rad_1> ... <rad_{N-1}>   <- radius axis; the corner label
//                                                    names the two axes, ignored on read
//   <mom_0>  v00 v01 ... v0{N-1}                   <- one row per momentum value
//   <mom_1>  v10 v11 ... v1{N-1}
//   ...
//
// Every '#' line is skipped on read except the keys 'system', 'potential',
// 'waves' and 'nbody', so those lines must carry the bare value only (no trailing text); all other
// comments are for the human reader. Blank lines are skipped. Everything is
// whitespace-separated; axes may be non-uniform. Values are row-major:
// values[iMom * nRadius + iRad].

namespace {

std::string Trim(const std::string& s) {
    const std::size_t a = s.find_first_not_of(" \t");
    if (a == std::string::npos) return "";
    const std::size_t b = s.find_last_not_of(" \t");
    return s.substr(a, b - a + 1);
}

}  // namespace

WaveFunction::WaveFunction(std::vector<double> momentum, std::vector<double> radius,
                           std::vector<double> values, int nBody, std::string system,
                           std::string potential, std::string waves,
                           std::string description)
    : fMomentum(std::move(momentum)),
      fRadius(std::move(radius)),
      fValues(std::move(values)),
      fNBody(nBody),
      fSystem(std::move(system)),
      fPotential(std::move(potential)),
      fWaves(std::move(waves)),
      fDescription(std::move(description)) {
    if (fNBody != 2 && fNBody != 3)
        throw std::runtime_error("WaveFunction: nBody must be 2 or 3");
    CheckShape();
}

WaveFunction::WaveFunction(const std::string& filename) : fNBody(0) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("WaveFunction: cannot open " + filename);

    bool haveNBody = false;
    bool haveHeader = false;
    std::string line;
    while (std::getline(in, line)) {
        while (!line.empty() &&
               (line.back() == ' ' || line.back() == '\t' || line.back() == '\r'))
            line.pop_back();
        if (line.empty()) continue;

        if (line[0] == '#') {  // metadata "key = value" or free comment
            const std::string body = line.substr(1);
            const std::size_t eq = body.find('=');
            if (eq == std::string::npos) continue;
            const std::string key = Trim(body.substr(0, eq));
            const std::string val = Trim(body.substr(eq + 1));
            if (key == "system") {
                fSystem = val;
            } else if (key == "potential") {
                fPotential = val;
            } else if (key == "waves") {
                fWaves = val;
            } else if (key == "nbody") {
                fNBody = std::atoi(val.c_str());
                if (fNBody != 2 && fNBody != 3)
                    throw std::runtime_error("WaveFunction: nbody must be 2 or 3");
                haveNBody = true;
            }
            continue;
        }

        std::istringstream iss(line);
        if (!haveHeader) {  // first data line: "mom\radius" corner label + radius axis
            std::string corner;
            iss >> corner;
            for (double r; iss >> r;) fRadius.push_back(r);
            haveHeader = true;
        } else {  // momentum value + one table row
            double mom;
            iss >> mom;
            fMomentum.push_back(mom);
            std::size_t count = 0;
            for (double v; iss >> v; count++) fValues.push_back(v);
            if (count != fRadius.size())
                throw std::runtime_error(
                    "WaveFunction: a data row has the wrong number of columns");
        }
    }

    if (!haveNBody) throw std::runtime_error("WaveFunction: file is missing '# nbody = ...'");
    if (!haveHeader) throw std::runtime_error("WaveFunction: file has no data");
    CheckShape();
}

void WaveFunction::CheckShape() const {
    if (fValues.size() != fMomentum.size() * fRadius.size())
        throw std::runtime_error(
            "WaveFunction: values.size() != momentum.size() * radius.size()");
}

// Human-readable one-line summary of an axis: count, endpoints and spacing.
static std::string AxisSummary(const std::vector<double>& a) {
    std::ostringstream s;
    s << std::setprecision(12);
    s << a.size() << " points";
    if (a.empty()) return s.str();
    s << ", range [" << a.front() << ", " << a.back() << "]";
    if (a.size() >= 2) {
        const double step = (a.back() - a.front()) / static_cast<double>(a.size() - 1);
        bool uniform = true;
        for (std::size_t i = 1; i < a.size() && uniform; i++)
            if (std::fabs((a[i] - a[i - 1]) - step) > 1e-6 * std::fabs(step)) uniform = false;
        if (uniform)
            s << ", uniform step " << step;
        else
            s << ", non-uniform spacing";
    }
    return s.str();
}

void WaveFunction::Save(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("WaveFunction: cannot open " + filename);
    out << std::setprecision(17);

    const bool three = (fNBody == 3);
    const char* psiName = three ? "|Psi(Q3, rho)|^2" : "|psi(k*, r*)|^2";
    const char* momName = three ? "Q3 (three-body momentum transfer)" : "k* (pair relative momentum)";
    const char* radName = three ? "rho (three-body hyper-radius)" : "r* (pair relative separation)";

    out << "# yaffa.WaveFunction 1\n";
    out << "#\n";
    out << "# Tabulated squared modulus of a femtoscopic wave function, " << psiName << ",\n";
    out << "# sampled on a rectangular (momentum, radius) grid. It is stored one row\n";
    out << "# per momentum value; integrating it against an emission source over the\n";
    out << "# radius axis (the Koonin-Pratt formula) gives the correlation function.\n";
    out << "#\n";
    out << "# quantity      : " << psiName << "\n";
    out << "# momentum axis : " << momName << ", MeV/c -- one value per row\n";
    out << "#                 " << AxisSummary(fMomentum) << "\n";
    out << "# radius axis   : " << radName << ", fm -- one value per column\n";
    out << "#                 " << AxisSummary(fRadius) << "\n";
    out << "# table order   : row-major, values[iMom * nRadius + iRad]\n";
    out << "#\n";
    out << "# data block below:\n";
    out << "#   row 1   : mom\\radius  r_0 r_1 ... r_{nRadius-1}   (the radius axis)\n";
    out << "#   row i+1 : m_i  v(i,0) v(i,1) ... v(i,nRadius-1)    where v(i,j) holds "
        << psiName << " at (m_i, r_j)\n";
    out << "#\n";
    if (!fDescription.empty()) {  // free-form notes, for the human reader only
        std::istringstream cs(fDescription);
        std::string cl;
        while (std::getline(cs, cl)) out << "# " << cl << "\n";
        out << "#\n";
    }
    out << "# system = " << fSystem << "\n";
    out << "# potential = " << fPotential << "\n";
    out << "# waves = " << fWaves << "\n";
    out << "# nbody = " << static_cast<int>(fNBody) << "\n";
    out << "#\n";

    out << "mom\\radius";  // corner label: column 1 is the momentum axis, the rest
                           // of this row is the radius axis. Ignored on read.
    for (std::size_t j = 0; j < fRadius.size(); j++) out << ' ' << fRadius[j];
    out << '\n';

    for (std::size_t i = 0; i < fMomentum.size(); i++) {
        out << fMomentum[i];
        for (std::size_t j = 0; j < fRadius.size(); j++) out << ' ' << At(i, j);
        out << '\n';
    }
}

void WaveFunction::Print() const {
    std::cout << "WaveFunction[" << fSystem << "] " << fNBody << "-body\n";
    std::cout << "  momentum: " << fMomentum.size() << " points";
    if (!fMomentum.empty())
        std::cout << " in [" << fMomentum.front() << ", " << fMomentum.back() << "]";
    std::cout << "\n  radius:   " << fRadius.size() << " points";
    if (!fRadius.empty())
        std::cout << " in [" << fRadius.front() << ", " << fRadius.back() << "]";
    std::cout << std::endl;
}

std::vector<double> WaveFunction::CorrelationFunction(
    const std::vector<double>& sourceOnRadiusAxis) const {
    const std::size_t nMom = fMomentum.size();
    const std::size_t nRad = fRadius.size();
    if (sourceOnRadiusAxis.size() != nRad)
        throw std::runtime_error("WaveFunction::CorrelationFunction: source has wrong size");

    // Radius-dependent weight J(r) * S(r); identical for every momentum.
    std::vector<double> weight(nRad);
    for (std::size_t j = 0; j < nRad; j++)
        weight[j] = Jacobian(fRadius[j]) * sourceOnRadiusAxis[j];

    std::vector<double> cf(nMom, 0.);
    for (std::size_t i = 0; i < nMom; i++) {
        double num = 0., den = 0.;
        for (std::size_t j = 0; j + 1 < nRad; j++) {  // trapezoid over the radius axis
            const double dr = fRadius[j + 1] - fRadius[j];
            num += 0.5 * dr * (weight[j] * At(i, j) + weight[j + 1] * At(i, j + 1));
            den += 0.5 * dr * (weight[j] + weight[j + 1]);
        }
        cf[i] = (den != 0.) ? num / den : 0.;
    }
    return cf;
}

std::vector<double> WaveFunction::CorrelationFunction(
    const std::function<double(double)>& source) const {
    std::vector<double> sampled(fRadius.size());
    for (std::size_t j = 0; j < fRadius.size(); j++) sampled[j] = source(fRadius[j]);
    return CorrelationFunction(sampled);
}
