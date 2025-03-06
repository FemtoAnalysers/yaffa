/*
Script to compute the same- and mixed event distributions from pythia events.
IMPORTANT: this script must be run with Pthia 8.310 as other versions were found to fail.
    * 8.304 sometimes fails to force the decay channels of some hadrons like Delta- and Sigma(1385)-
    * 8.312 was observed to crash during color reconnection mode runs with exit code 11 (segmentation violation)
*/

// C++ libraries
#include <array>
#include <deque>
#include <map>
#include <string>
#include <vector>
#include <cxxabi.h>  // For demangling

// 3rd party libraries
#include "yaml-cpp/yaml.h"

// ROOT libraries
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include "Math/Vector4D.h"
#include "Math/Boost.h"
#include "TRandom3.h"

// ALICE libraries
#define PYTHIA_V 8312

#if PYTHIA_V == 8312
#include "Pythia8/Pythia.h"
#include "Pythia8/ParticleData.h"
#else
#include "TPythia.h"
#include "TParticleData.h"
#endif

#if false
#define DEBUG(msg, ...) do { printf(msg, ##__VA_ARGS__); } while(0)
#define DEBUG_VAR(var) do { std::cerr << #var << ": " << var << std::endl; } while (0)
#else
#define DEBUG(msg, ...)
#define DEBUG_VAR(var)
#endif

const double HBARC =  197.3269804; // um=MeV*fm

static YAML::Node cfgMother;
static YAML::Node cfgPart0;
static YAML::Node cfgPart1;

// Find Bin
template<typename T>
int FindBin(T value, std::vector<T> list) {
    for (size_t i = 0; i < list.size() - 1; i++) {
        if (list[i] < value && value <= list[i + 1]) {
            return i;
        }
    }

    return -1;
}


void SetProcess(Pythia8::Pythia &pythia, std::string process) {
    if (process == "SoftQCD") {
        pythia.readString("SoftQCD:all = on");
        return;
    }
    
    if (process == "HardQCD") {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
        return;
    }
    
    if (process == "NonDiffractive") {
        pythia.readString("SoftQCD:nonDiffractive = on");
        return;
    }
    
    throw std::invalid_argument("Process not implemented. Exit!");
}

bool haveCommonAncestor(Pythia8::Particle p1, Pythia8::Particle p2) {
    int status1 = std::abs(p1.status());
    int status2 = std::abs(p2.status());

    if (!(81 <= status1 && status1 <= 89)) {
        return false;
        // std::cerr << "Error: the 'haveCommonAncestor' method received a non-primary hadron. Status: " << status1 << std::endl;
        // throw std::invalid_argument("abs(status) of p1 must be in [81, 89]");
    }
    if (!(81 <= status2 && status2 <= 89)) {
        return false;
        // std::cerr << "Error: the 'haveCommonAncestor' method received a non-primary hadron" << std::endl;
        // throw std::invalid_argument("abs(status) of p2 must be in [81, 89]");
    }

    return p1.mother1() == p2.mother1() && p1.mother2() == p2.mother2() && p1.mother1() < p1.mother2() && p1.mother1() > 0 && p2.mother1() > 0;
}

void SetTune(Pythia8::Pythia &pythia, std::string tune) {
    if (tune == "Monash") {
        pythia.readString(Form("Tune:pp = 14"));
        return;
    }
    
    if (tune == "CRMode0") {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 2.9");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.43");
        pythia.readString("ColourReconnection:timeDilationMode = 0");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
        return;
    }
    
    if (tune == "CRMode2") {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.20");
        pythia.readString("ColourReconnection:timeDilationMode = 2");
        pythia.readString("ColourReconnection:timeDilationPar = 0.18");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
        return;
    }
    
    if (tune == "CRMode3") {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.15");
        pythia.readString("ColourReconnection:timeDilationMode = 3");
        pythia.readString("ColourReconnection:timeDilationPar = 0.073");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
        return;
    }
    
    throw std::invalid_argument("Tune not implemented. Exit!");
}

/*
Integral for Blast Wave. Taken from AliPWGFunc.
    r = x[0]
    mass = par[0]
    beta = par[1]
    temp = par[2]
    n = par[3]
    norm = par[4]
*/
double BlasWaveIntegrand(const double *x, const double *par) {
    double r = x[0];

    double mass = par[0];
    double pT = par[1];
    double beta_s = par[2];
    double temp = par[3];
    double n = par[4];

    // Keep beta within reasonable limits
    double beta = beta_s * TMath::Power(r, n);
    if (beta > 0.9999999999999999) beta = 0.9999999999999999;

    double mT = TMath::Sqrt(mass * mass + pT * pT);

    double rho0 = TMath::ATanH(beta);
    double arg00 = pT * TMath::SinH(rho0) / temp;
    if (arg00 > 700.) arg00 = 700.;  // Avoid Floating Point Exception
    double arg01 = mT * TMath::CosH(rho0) / temp;
    double f0 = r * mT * TMath::BesselI0(arg00) * TMath::BesselK1(arg01);

    return f0;
}

/*
Blast Wave. Taken from AliPWGFunc.
    pT = x[0]

    mass = par[0]
    beta = par[1]
    temp = par[2]
    n = par[3]
    norm = par[4]
*/
double BlastWave(const double *x, const double *par) {
    double pT = x[0];

    double mass = par[0];
    double beta = par[1];
    double temp = par[2];
    double n = par[3];
    double norm = par[4];

    static TF1 *fIntBG = new TF1("fIntBG", BlasWaveIntegrand, 0, 1, 5);
    fIntBG->SetNpx(1000);
    fIntBG->SetParameters(mass, pT, beta, temp, n);
    double result = fIntBG->Integral(0, 1);

    return x[0] * result * norm;
}

/*
Non-relativistic Breit-Wigner distribution with threhsold. Not normalized to 1.
---

parameters:
---
    0) mass of the resonance
    1) width of the resonance
    2) threshold
*/
double BreitWigner(double *x, double *par) {
    double energy = x[0];

    double mass = par[0];
    double width = par[1];
    double threshold = par[2];

    if (mass < threshold) {
        throw std::runtime_error("BreitWigner: invalid parameters. Mass must be larger than the threshold.");
    }   

    if (energy < threshold) {
        return 0;
    }

    return width / 2 / std::numbers::pi / (pow(energy - mass, 2) + 0.25 * width * width);
}

/*
Non-relativistic Breit-Wigner distribution with threhsold normalized to 1.
---
Reference:
    F. Giacosa and V. Shastry
    Sill distribution: Genesis and salient features
    https://doi.org/10.1393/ncc/i2024-24186-8

parameters:
---
    0) mass of the resonance
    1) width of the resonance
    2) threshold
*/
// double BreitWigner(double *x, double *par) {
//     double energy = x[0];

//     double mass = par[0];
//     double width = par[1];
//     double threshold = par[2];

//     if (mass < threshold) {
//         throw std::runtime_error("BreitWigner: invalid parameters. Mass must be larger than the threshold.");
//     }   

//     if (energy < threshold) {
//         return 0;
//     }

//     return width / 2 / std::numbers::pi / (pow(energy - mass, 2) + 0.25 * gamma * gamma);
// }


float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2) {
    ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
    ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
    ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
    ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);

    ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;
    float kStar = 0.5 * trackRelK.P();
    return kStar;
}

/*
Non-relativistic Sill distribution normalized to one.
---
Reference:
    F. Giacosa and V. Shastry
    Sill distribution: Genesis and salient features
    https://doi.org/10.1393/ncc/i2024-24186-8

Functional form (LaTeX):
d^{\text{nrSill}}(E) = \frac{\gamma \sqrt{E - E_{\text{th}}}}{2\pi} \left[ (E - M)^2 + \frac{1}{4} \left( \gamma \sqrt{E - E_{\text{th}}} \right)^2 \right]^{-1} \theta(E - E_{\text{th}})

with:
\gamma = \Gamma / \sqrt{M - E_{\text{th}}

parameters:
---
    0) mass of the resonance
    1) width of the resonance
    2) threshold
*/
double Sill(double *x, double *par) {
    double energy = x[0];

    double mass = par[0];
    double width = par[1];
    double threshold = par[2];

    if (mass < threshold) {
        throw std::runtime_error("Sill: invalid parameters. Mass must be larger than the threshold.");
    }   

    if (energy < threshold) {
        return 0;
    }

    double gamma = width / std::sqrt(mass - threshold);
    double norm = gamma * std::sqrt(energy - threshold) / 2 / std::numbers::pi;

    return norm / (pow(energy - mass, 2) + 0.25 * gamma * gamma * (energy - threshold));
}

float ComputeKstar(Pythia8::Particle part1, Pythia8::Particle part2) {
    ROOT::Math::PxPyPzMVector p1(part1.px(), part1.py(), part1.pz(), part1.m());
    ROOT::Math::PxPyPzMVector p2(part2.px(), part2.py(), part2.pz(), part2.m());

    return ComputeKstar(p1, p2);
}

// Compute transverse mass (m_T) of the pair
double ComputeMt(Pythia8::Particle part1, Pythia8::Particle part2) {
    ROOT::Math::PxPyPzMVector p1(part1.px(), part1.py(), part1.pz(), part1.m());
    ROOT::Math::PxPyPzMVector p2(part2.px(), part2.py(), part2.pz(), part2.m());

    double kT = 0.5 * (p1 + p2).Pt();
    double m = 0.5 * (p1.M() + p2.M());
    return pow(kT * kT + m * m, 0.5);
}

// Return true if the PDG code corresponds to a charged and long-lived particle
bool IsDetectable(const int &pdg) {
    int absPdg = std::abs(pdg);

    return absPdg == 11 ||   // electrons
           absPdg == 13 ||   // muons
           absPdg == 211 ||  // pions
           absPdg == 321 ||  // kaons
           absPdg == 2212;   // protons
}

bool TriggerHM(const std::vector<Pythia8::Particle> *particles) {
    // evaluate multiplicity at forward rapidity
    int nChForward = 0;
    for (size_t iPart = 3; iPart < particles->size(); iPart++) {
        auto &part = (*particles)[iPart];

        int pdg = std::abs(part.id());
        float eta = part.eta();

        if (IsDetectable(pdg) && ((-3.7 < eta && eta < -1.7) || (2.8 < eta && eta < 5.1))) {  // V0A and V0C acceptance
            nChForward++;
        }
    }
    return nChForward > 130;
}

bool Trigger(const std::vector<Pythia8::Particle> *particles, std::string trigger) {
    if (trigger == "MB") {
        return true;
    }

    if (trigger == "HM") {
        return TriggerHM(particles);
    }

    throw std::invalid_argument("Invalid trigger");
}

// Compute the multiplicity of charged particles in the TPC acceptance
int ComputeMultTPC(const std::vector<Pythia8::Particle> *particles) {
    int mult = 0;
    for (size_t iPart = 3; iPart < particles->size(); iPart++) {
        Pythia8::Particle part = (*particles)[iPart];

        if (part.isFinal() && std::abs(part.eta()) < 0.8 && IsDetectable(part.id())) mult++;
    }
    return mult;
}

// Set the default selections
template<typename T1, typename T2>
void SetDefault(YAML::Node &&node, T1 min, T2 max) {
    if (!node.IsDefined() || node.IsNull()) {
        node.push_back(min);
        node.push_back(max);
    }
}

void SetDefaults(YAML::Node &cfg) {
    SetDefault(cfg["status"], -300, 300);
    SetDefault(cfg["pt"], 0, 100);
    SetDefault(cfg["eta"], -100, 100);
    SetDefault(cfg["y"], -10, 10);
    SetDefault(cfg["prodvtx"], -1., 1000);

    // Apply defaults to the daughters
    if (cfg["daus"].IsDefined() && cfg["daus"].IsSequence() && cfg["daus"].size() > 0) {
        for (auto &&node : cfg["daus"]) {
            SetDefaults(node);
        }
    }
}

template<typename T>
bool IsSelected(T value, const YAML::Node node, bool includeExtremes = false) {
    std::array<T, 2> range = node.as<std::array<T, 2>>();
    if (includeExtremes) return range[0] <= value && value <= range[1];
    return range[0] < value && value < range[1];
}

bool IsSelected(const std::vector<Pythia8::Particle> *particles, int iPart, const YAML::Node &cfgSelections) {
    auto& part = (*particles)[iPart];

    if (std::abs(part.id()) != cfgSelections["pdg"].as<int>()) return false;
    if (!IsSelected(part.status(), cfgSelections["status"], true)) return false;
    if (!IsSelected(part.pT(), cfgSelections["pt"])) return false;
    if (!IsSelected(part.eta(), cfgSelections["eta"])) return false;
    if (!IsSelected(part.y(), cfgSelections["y"])) return false;
    double prodvtx = pow(part.xProd() * part.xProd() + part.yProd() * part.yProd() + part.zProd() * part.zProd(), 0.5);
    if (!IsSelected(prodvtx, cfgSelections["prodvtx"])) return false;

    if (cfgSelections["daus"].IsDefined() && cfgSelections["daus"].IsSequence() && cfgSelections["daus"].size() > 0) {
        auto dauIdx = part.daughterList();
        for (long unsigned int iDau = 0; iDau < dauIdx.size(); iDau++) {
            if (!IsSelected(particles, dauIdx[iDau], cfgSelections["daus"][iDau])) return false;
        }
    }

    return true;
}

std::string GetDaughters(YAML::Node cfg) {
    // Check if the "daus" key exists and is a sequence
    if (!cfg["daus"].IsDefined() || !cfg["daus"].IsSequence() || cfg["daus"].size() == 0) return "";

    // Check if "pdg" is valid
    if (!cfg["pdg"].IsScalar()) return "";

    // Gather the daughter particles into a string
    std::string daus;
    for (const auto &dau : cfg["daus"]) {
        daus += " " + dau["pdg"].as<std::string>();
    }

    return daus;
}

int GetParticlesInDecayChain(const std::vector<Pythia8::Particle> *particles, int iPart, YAML::Node cfgMom, std::vector<int> &part0, std::vector<int> &part1) {
    const auto& mom = (*particles)[iPart];

    std::multiset<int> daughters;
    for (int iDau = mom.daughter1(); iDau <= mom.daughter2(); iDau++) {
        auto dau = (*particles)[iDau];
        daughters.insert(std::abs(dau.id()));
    }
    std::multiset<int> daughtersTarget;
    for (size_t iCfg = 0; iCfg < cfgMom.size(); iCfg++) {
        daughtersTarget.insert(cfgMom[iCfg]["pdg"].as<int>());
    }

    if (daughters != daughtersTarget) return 0;

    DEBUG("Start analyzing the decay tree of pdg=%d idx=%d\n", mom.id(), iPart);
    DEBUG("Start analyzing daus:\n");
    for (int iDau = mom.daughter1(); iDau <= mom.daughter2(); iDau++) {
        auto dau = (*particles)[iDau];

        DEBUG("    Checking now daughter with pdg=%d, idx=%d\n", dau.id(), iDau);
        // find the corresponding particle in the cfg file
        int dauCfgIdx = -1;
        for (size_t iCfg = 0; iCfg < cfgMom.size(); iCfg++) {
            DEBUG("        Check if compatible with config n. %lu/%lu\n", iCfg+1, cfgMom.size());
            if (std::abs(dau.id()) == cfgMom[iCfg]["pdg"].as<int>()) {
                DEBUG("        --> Match in configuration!\n");
                dauCfgIdx = iCfg;
                break;
            }
        }

        if (dauCfgIdx == -1) {
            DEBUG("    Daughter not found. Stop here.\n");
            return 0;
        }

        const YAML::Node& dauCfg = cfgMom[dauCfgIdx];

        if (dauCfg["daus"].IsNull()) {
            DEBUG("Daus are null. daupdg = %d cfg0=%d cfg1=%d\n", dau.id(), cfgPart0["pdg"].as<int>(), cfgPart1["pdg"].as<int>());
            if (!dauCfg["use"].as<bool>()) continue;

            if (std::abs(dau.id()) == cfgPart0["pdg"].as<int>()) {
                DEBUG("push-back in part0\n");
                part0.push_back(iDau);
            } else if  (std::abs(dau.id()) == cfgPart1["pdg"].as<int>()) {
                DEBUG("push-back in part1\n");
                part1.push_back(iDau);
            }
        } else if (dauCfg["daus"].IsSequence()) {
            DEBUG("Recursively looking into the daughters\n");
            if (!GetParticlesInDecayChain(particles, iDau, dauCfg["daus"], part0, part1)) return 0;
        }
    }

    return 1;
}

std::tuple<int, double, double> ComputeBinning(double xMin, double xMax, int precision=3) {
    double mul = std::pow(10, precision);
    int xMinMul = (int) (xMin * mul);
    int xMaxMul = (int) (xMax * mul);

    int nBins = xMaxMul - xMinMul;

    return {nBins, xMinMul / mul, xMaxMul / mul};
}

/*
Returns the name of the type as a string.
*/
template <typename T>
std::string getType() {
    int status;
    std::string name = abi::__cxa_demangle(typeid(T).name(), 0, 0, &status);

    if (status == 0) {
        return name;
    } else {
        return typeid(T).name();
    }
}

/*
Load a variable from a YAML node and check its validity.
*/
template<typename T>
T load(YAML::Node node, std::string name) {
    if (!node) {
        throw YAML::InvalidNode("Invalid node!");
    }

    if (!node[name].IsDefined()) {
        throw std::runtime_error("Key '" + name + "' is not defined!");
    }

    try {
        return node[name].as<T>();
    } catch (const YAML::BadConversion&) {
        throw std::runtime_error("Key '" + name + "' must be of type " + getType<T>() + "!");
    }
}

void MakeDistr(
    std::string inFileName = "",
    std::string oFileName = "Distr.root",
    std::string cfgFile = "cfg_makedistr_example.yml",
    int seed = 31
    ) {

    bool doGenerate = inFileName == "";

    // Load PDG
    TDatabasePDG *PDG = TDatabasePDG::Instance();

    // Load simulation settings
    YAML::Node cfg = YAML::LoadFile(cfgFile.data());
    size_t nEvents = cfg["nevts"].as<unsigned int>();
    unsigned int md = cfg["mixdepth"].as<int>();
    bool rejevtwopairs = cfg["rejevtwopairs"].as<bool>();
    std::vector<double> mTBins = load<std::vector<double>>(cfg, "mTBins");
    std::vector<double> mTMins(mTBins);
    mTMins.pop_back();
    std::vector<double> mTMaxs(mTBins);
    mTMaxs.erase(mTMaxs.begin());

    // Load selections for mother particle
    cfgMother = YAML::Clone(cfg["decaychain"]);
    cfgMother["daus"] = YAML::Null; // Don't check the daughters. // todo: understand why if this line is removed all hists are empty

    // Load selections for part 0
    cfgPart0 = cfg["part0"];

    // Load selections for part 1. If they don't exist, use the same as part 0 (same-particle pairs e.g. pp)
    cfgPart1 = cfg["part1"].IsNull() ? cfgPart0 : cfg["part1"];

    std::cout << "\033[34mParticle selections before defaults (Mother)\033[0m" << std::endl;
    std::cout << cfgMother<< std::endl;

    std::cout << "\033[34mParticle selections before defaults (part0)\033[0m" << std::endl;
    std::cout << cfgPart0<< std::endl;

    std::cout << "\033[34mParticle selections before defaults (part1)\033[0m" << std::endl;
    std::cout << cfgPart1<< std::endl;

    SetDefaults(cfgMother);
    SetDefaults(cfgPart0);
    SetDefaults(cfgPart1);

    std::cout << "\033[34mParticle selections after defaults (Mother)\033[0m" << std::endl;
    std::cout << cfgMother<< std::endl;

    std::cout << "\033[34mParticle selections after defaults (part0)\033[0m" << std::endl;
    std::cout << cfgPart0<< std::endl;

    std::cout << "\033[34mParticle selections after defaults (part1)\033[0m" << std::endl;
    std::cout << cfgPart1<< std::endl;

    auto pdg0 = cfgPart0["pdg"].as<int>();
    auto pdg1 = cfgPart1["pdg"].as<int>();

    Pythia8::Pythia pythia;
    pythia.readString("Next:numberShowEvent = 0");
    SetProcess(pythia, cfg["process"].as<std::string>());
    SetTune(pythia, cfg["tune"].as<std::string>());

    std::string trigger = cfg["trigger"].as<std::string>();

    // Set decay channel for part0
    if (std::string daus = GetDaughters(cfgPart0); daus != ""){
        pythia.readString(std::to_string(pdg0) + ":onMode = off");
        pythia.readString(std::to_string(pdg0) + ":onIfMatch =" + daus);
    }

    // Set decay channel for part1
    if (std::string daus = GetDaughters(cfgPart1); daus != ""){
        pythia.readString(std::to_string(pdg1) + ":onMode = off");
        pythia.readString(std::to_string(pdg1) + ":onIfMatch =" + daus);
    }

    if (cfg["injection"].size() > 1) {
        printf("Error. The BW mass limit is not implemented for more than 1 injection. Exit!");
        exit(1);
    }

    for (const auto &part : cfg["injection"]) {
        int myPdg = part["pdg"].as<int>();
        std::string name = part["name"].as<std::string>();
        std::string antiname = part["antiname"].as<std::string>();
        int spin = part["spin"].as<int>();
        int charge = part["charge"].as<int>();
        int color = 0;
        int mass = part["mass"].as<double>();
        double width = part["width"].as<double>();
        double tau0 = HBARC / width * 1.e-12; // Conversion fm -> mm
        double mMin = mass * 0.5;
        double mMax = 0; // If mMax < mMin then no upper limit is imposed

        pythia.particleData.addParticle(myPdg, name, antiname, spin, charge, color, mass, width, mMin, mMax, tau0);
        int meMode = part["meMode"].as<int>();
        pythia.particleData.readString(Form("%d:addChannel = 1 1 %d %s", myPdg, meMode, part["daus"].as<std::string>().data()));
    }

    std::cout << "Applying the following customization to pythia:" << std::endl;
    for (const auto &line : cfg["customization"]) {
        std::string lineStr = line.as<std::string>();
        std::cout << "   * " << lineStr << std::endl;
        pythia.readString(lineStr.data());
    }
    std::cout << "End of customization." << std::endl;

    // Compute the minimum mass that a resonance modelled with a Breit-Wigner can have
    double threshold=1.e12;
    if (cfg["injection"].size() > 0) {
        auto particleEntry = pythia.particleData.particleDataEntryPtr(cfg["injection"][0]["pdg"].as<int>());
        for (int iDecChn = 0; iDecChn < particleEntry->sizeChannels(); iDecChn++) {
            auto channel = particleEntry->channel(iDecChn);

            // Loop over the daughters in the decay channel
            double sum = 0;
            for (int iDau = 0; iDau < channel.multiplicity(); iDau++) {
                int dauPdg = channel.product(iDau);

                sum += pythia.particleData.particleDataEntryPtr(dauPdg).get()->m0();
            }
            if (sum < threshold) threshold = sum;
        }
        printf("Mass threshold: %.6f GeV\n", threshold);
    }

    // Load Pt and y distributions
    int pdgMother = cfg["decaychain"]["pdg"].as<int>();

    std::string kinemFile = load<std::string>(cfg["decaychain"], "kinemfile");
    std::string efficiency0 = load<std::string>(cfg["part0"], "efficiency");
    std::string efficiency1 = load<std::string>(cfg["part1"], "efficiency");
    std::string ptshape = load<std::string>(cfg["decaychain"], "ptshape");
    std::string yshape = load<std::string>(cfg["decaychain"], "yshape");
    std::string etashape = load<std::string>(cfg["decaychain"], "etashape");
    bool corr = load<bool>(cfg["decaychain"], "corr");

    if (kinemFile != "" && ptshape != "") {
        printf("\033[31mError: you can't set the kinematics both from 'ptshape' and 'kinemfile'. Kinematics from file will be used.\033[0m\n");
    }

    if (kinemFile != "" && yshape != "") {
        printf("\033[31mError: you can't set the kinematics both from 'yshape' and 'kinemfile'. Kinematics from file will be used.\033[0m\n");
    }

    if (cfg["injection"].size() == 0 && ptshape != "") {
        printf("\033[33mWarning: you are specifying a pt-shape but the mother particle is not injected. Set ptshape: '' to silence this warning.\033[0m\n");
    }

    if (cfg["injection"].size() == 0 && yshape != "") {
        printf("\033[33mWarning: you are specifying a y-shape but the mother particle is not injected. Set yshape: '' to silence this warning.\033[0m\n");
    }

    if (etashape != "") {
        printf("\033[33mWarning: you should never use a 1D eta-shape because it's correlated with pT. Be careful!\033[0m\n");
    }

    if (etashape != "" && yshape != "") {
        printf("\031[31mError: decide if you want a y-shape or eta-shape. Exit!\033[0m\n");
        exit(1);
    }

    // Load 2D histogram for kinematics
    TH2D *hYvsPt = nullptr;
    TH1D *hEff0 = nullptr;
    TH1D *hEff1 = nullptr;
    TH1D *hPt = nullptr;
    TH1D *hY = nullptr;
    TH1D *hEta = nullptr;
    TF1 *fEff0 = nullptr;
    TF1 *fEff1 = nullptr;
    TF1 *fPt = nullptr;
    TF1 *fY = nullptr;
    TF1 *fEta = nullptr;

    if (kinemFile != "") {
        TFile *fKinem = TFile::Open(kinemFile.data());
        if (!fKinem) exit(1); // TFile::Open already prints an error message

        hYvsPt = (TH2D *) fKinem->Get("hYvsPt");
        if (!hYvsPt) {
            printf("Error: cannot get object 'hYvsPt' from file '%s'. Exit!\n", kinemFile.data());
            printf("Available objects: \n");
            fKinem->ls();
            exit(1);
        }
        hYvsPt->SetDirectory(0);
        fKinem->Close();

        if (!corr) {
            hPt = (TH1D *) hYvsPt->ProjectionX("hPt", 1, -1);
            hY = (TH1D *) hYvsPt->ProjectionY("hY", 1, -1);

            hPt->SetDirectory(0);
            hY->SetDirectory(0);

            delete hYvsPt;
            hYvsPt = nullptr;
        }
    }

    // Load Efficiency
    if (size_t pos = efficiency0.find(':'); pos != std::string::npos) {
        TFile *effFile = TFile::Open(efficiency0.substr(0, pos).data());
        if (!effFile) exit(1);
        hEff0 = (TH1D *) effFile->Get(efficiency0.substr(pos + 1).data());
        hEff0->SetDirectory(0);
        effFile->Close();
    } else {
        fEff0 = new TF1("fEff0", efficiency0.data(), 0, 10);
    }

    // Load Efficiency
    if (size_t pos = efficiency1.find(':'); pos != std::string::npos) {
        TFile *effFile = TFile::Open(efficiency1.substr(0, pos).data());
        if (!effFile) exit(1);
        hEff1 = (TH1D *) effFile->Get(efficiency1.substr(pos + 1).data());
        hEff1->SetDirectory(0);
        effFile->Close();
    } else {
        fEff1 = new TF1("fEff1", efficiency1.data(), 0, 10);
    }

    // Set pT shape
    size_t pos = ptshape.find(':');
    if (pos != std::string::npos) {
        std::string first = ptshape.substr(0, pos);
        std::string second = ptshape.substr(pos + 1);

        if (first == "blastwave") {
            // first = "blastwave"
            // second = version of blastwave parameters

            // format: Journal, volume number, year, page range
            std::map<std::string, std::array<double, 3>> blastwavePars = {
                {"PLB_728_2014_2538", {0.26, 0.166, 3.9}}, // peripheral pPb collisions. See Tab 5. //! This is in <beta_T> not beta_s!!!
                {"EPJC_80_2020_693_extrap", {0.529863, 0.1420758, 1.126957}}, // From fit and extrapolation of the hep data //! This is in <beta_T> not beta_s!!!
                {"EPJC_80_2020_693_extrap_betas", {0.828429408, 0.1420758, 1.126957}}, // From fit and extrapolation of the hep data. beta_s computed as <beta_T> * (n+2)/2
                {"EPJC_80_2020_693_extrap_betas_betaTUpper", {0.8521, 0.1420758, 1.126957}}, // From fit and extrapolation of the hep data. beta_s computed as <beta_T> * (n+2)/2; syst variation
                {"EPJC_80_2020_693_extrap_betas_betaTLower", {0.8052, 0.1420758, 1.126957}}, // From fit and extrapolation of the hep data. beta_s computed as <beta_T> * (n+2)/2; syst variation
                {"EPJC_80_2020_693_extrap_betas_TKinUpper", {0.828429408, 0.156, 1.126957}}, // From fit and extrapolation of the hep data. beta_s computed as <beta_T> * (n+2)/2; syst variation
                {"EPJC_80_2020_693_extrap_betas_TKinLower", {0.828429408, 0.128, 1.126957}}, // From fit and extrapolation of the hep data. beta_s computed as <beta_T> * (n+2)/2; syst variation
                {"EPJC_80_2020_693_extrap_betas_NUpper", {0.8440, 0.1420758, 1.185}}, // From fit and extrapolation of the hep data. beta_s computed as <beta_T> * (n+2)/2; syst variation
                {"EPJC_80_2020_693_extrap_betas_NLower", {0.8133, 0.1420758, 1.069}}, // From fit and extrapolation of the hep data. beta_s computed as <beta_T> * (n+2)/2; syst variation
            };
            
            printf("Using Blast-Wave function for pt distribution\n");
            double mass;
            if (cfg["injection"].size() > 0) {
                mass = cfg["injection"][0]["mass"].as<double>();
            } else {
                mass = pythia.particleData.particleDataEntryPtr(pdgMother).get()->m0();
            }

            auto [beta, Tkin, n] = blastwavePars[second];
            fPt = new TF1("fPt", BlastWave, 0, 10, 5);
            fPt->SetParameter(0, mass);
            fPt->SetParameter(1, beta);
            fPt->SetParameter(2, Tkin);
            fPt->SetParameter(3, n);
            fPt->SetParameter(4, 1);
        } else {
            // first = file name
            // second = path of histogram inside root file
            TFile *ptFile = TFile::Open(first.data());
            if (!ptFile) exit(1);
            hPt = (TH1D *) ptFile->Get(second.data());
            hPt->SetDirectory(0);
            ptFile->Close();
        }
    } else {
        fPt = new TF1("fPt", ptshape.data(), 0, 10);
    }

    // Set rapidity shape
    if (size_t pos = yshape.find(':'); pos != std::string::npos) {
        TFile *yFile = TFile::Open(yshape.substr(0, pos).data());
        if (!yFile) exit(1);
        hY = (TH1D *) yFile->Get(yshape.substr(pos + 1).data());
        hY->SetDirectory(0);
        yFile->Close();
    } else if (yshape != "") {
        fY = new TF1("fY", cfg["decaychain"]["yshape"].as<std::string>().data(), -10, 10);
    }

    // Set pseudo-rapidity shape
    if (size_t pos = etashape.find(':'); pos != std::string::npos) {
        std::string fileName = etashape.substr(0, pos);
        std::string objName = etashape.substr(pos + 1);

        TFile *etaFile = TFile::Open(fileName.data());
        if (!etaFile) exit(1);
        if (!etaFile) {
            printf("\033[31mError: file '%s' could not be loaded. Exit!\033[0m\n", fileName.data());
            exit(1);
        }

        hEta = (TH1D *) etaFile->Get(objName.data());
        if (!hEta) {
            printf("\033[31mError: etashape '%s' could not be loaded. Exit!\033[0m\n", objName.data());
            exit(1);
        }

        hEta->SetDirectory(0);
        etaFile->Close();
    } else if (etashape != ""){
        printf("Using formula");
        fEta = new TF1("fEta", cfg["decaychain"]["etashape"].as<std::string>().data(), -10, 10);
    }

    std::string lineShape = load<std::string>(cfg["injection"][0], "lineshape");
    TF1 *fLineShape;
    if (lineShape == "breitwigner") {
        fLineShape = new TF1("fLineShape", BreitWigner, 0, 20, 3);
        fLineShape->SetParameter(0, cfg["injection"][0]["mass"].as<double>()); // Value in GeV
        fLineShape->SetParameter(1, cfg["injection"][0]["width"].as<double>() / 1000); // Value in MeV
        fLineShape->SetParameter(2, threshold);
    } else if (lineShape == "sill") {
        fLineShape = new TF1("fLineShape", Sill, 0, 20, 3);
        fLineShape->SetParameter(0, cfg["injection"][0]["mass"].as<double>()); // Value in GeV
        fLineShape->SetParameter(1, cfg["injection"][0]["width"].as<double>() / 1000); // Value in MeV
        fLineShape->SetParameter(2, threshold);
    } else {
        throw std::invalid_argument("Lineshape not implemented!");
    }
    fLineShape->SetTitle(";#it{M} (GeV/#it{c}^2);Probability");
    fLineShape->SetNpx(100000);

    // Setting the seed here is not sufficient to ensure reproducibility, setting the seed of gRandom is necessary
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", cfg["sqrts"].as<double>() * 1000); // from TeV to GeV
    pythia.init();

    gRandom->SetSeed(seed); // Set the seed to ensure reproducibility of the evevents generated by pythia

    // QA histograms
    TH1D * hEvt = new TH1D("hEvt", ";;Counts", 1, 0.5, 1.5);
    hEvt->GetXaxis()->SetBinLabel(1, "Events");

    TH1D * hEvtMult = new TH1D("hEvtMult", ";#it{N}_{ch}|_{|#eta|<0.8};Counts", 100, 0., 100);

    // Pairs
    std::map<std::tuple<int, int, int>, std::map<std::string, TH1D *>> hSE;
    std::map<std::tuple<int, int, int>, TH1D *> hME;
    
    // Pair QA
    std::map<std::tuple<int, int, int>, TH2D *> hPairMultSE;
    std::map<std::tuple<int, int, int>, TH2D *> hPtMotherVsKstar;

    int nPart0 = PDG->GetParticle(pdg0)->AntiParticle() && pdg0 != pdg1 ? 2 : 1;
    int nPart1 = PDG->GetParticle(pdg1)->AntiParticle() ? 2 : 1;
    for (int iPart0 = 0; iPart0 < nPart0; iPart0++) {
        for (int iPart1 = 0; iPart1 < nPart1; iPart1++) {
            // mT = -1 corresponds to mT integrated
            std::tuple<int, int, int> pair = {iPart0, iPart1, -1};
            DEBUG("Inserting histograms for pair (%d, %d)\n", iPart0, iPart1);
            hSE.insert({pair, {}});
            hSE[pair].insert(std::pair<std::string, TH1D *>({"Common",  new TH1D(Form("hSE%d%dCommon", iPart0, iPart1), ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)}));
            hSE[pair].insert(std::pair<std::string, TH1D *>({"NonCommon", new TH1D(Form("hSE%d%dNonCommon", iPart0, iPart1), ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)}));
            hME.insert({pair, new TH1D(Form("hME%d%d", iPart0, iPart1), ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)});
            hPairMultSE.insert({pair, new TH2D(Form("hPairMultSE%d%d", iPart0, iPart1), ";#it{N}_{0};#it{N}_{1};Counts", 51, -0.5, 50.5, 31, -0.5, 50.5)});
            hPtMotherVsKstar.insert({pair, new TH2D(Form("hPtMotherVsKstar%d%d", iPart0, iPart1), ";#it{k}* (GeV/#it{c});#it{p}_{T}^{Mother} (GeV/#it{c});Counts", 2000, 0, 2, 1000, 0, 10)});
            
            for (size_t iMt = 0; iMt < mTMins.size(); iMt++) {
                std::tuple<int, int, int> pair = {iPart0, iPart1, iMt};
                DEBUG("Inserting histograms for pair (%d, %d) iMt=%d\n", iPart0, iPart1, iMt);
                hSE.insert({pair, {}});
                hSE[pair].insert(std::pair<std::string, TH1D *>({"Common",  new TH1D(Form("hSE%d%d_mT%zu_Common", iPart0, iPart1, iMt), ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)}));
                hSE[pair].insert(std::pair<std::string, TH1D *>({"NonCommon", new TH1D(Form("hSE%d%d_mT%zu_NonCommon", iPart0, iPart1, iMt), ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)}));
                hME.insert({pair, new TH1D(Form("hME%d%d_mT%zu", iPart0, iPart1, iMt), ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)});
                hPairMultSE.insert({pair, new TH2D(Form("hPairMultSE%d%d_mT%zu", iPart0, iPart1, iMt), ";#it{N}_{0};#it{N}_{1};Counts", 51, -0.5, 50.5, 31, -0.5, 50.5)});
                hPtMotherVsKstar.insert({pair, new TH2D(Form("hPtMotherVsKstar%d%d_mT%zu", iPart0, iPart1, iMt), ";#it{k}* (GeV/#it{c});#it{p}_{T}^{Mother} (GeV/#it{c});Counts", 2000, 0, 2, 1000, 0, 10)});
            }
        }
    }

    // Single-particle QA for part 0
    double mass = pythia.particleData.particleDataEntryPtr(pdg0).get()->m0();
    auto [nBins, massMin, massMax] = ComputeBinning(mass*0.7, mass*1.3);
    std::map<std::string, TH1*> hQA0 = {
        {"mass", new TH1D("hMass0", ";#it{M} (GeV/#it{c}^{2});Counts", nBins, massMin, massMax)},
        {"pt", new TH1D("hPt0", ";#it{p}_{T} (GeV/#it{c});Counts", 1000, 0, 10)},
        {"y", new TH1D("hY0", ";#it{y};Counts", 200, -10, 10)},
        {"eta", new TH1D("hEta0", ";#eta;Counts", 200, -10, 10)},
    };

    // Single-particle QA for part 1
    mass = pythia.particleData.particleDataEntryPtr(pdg1).get()->m0();
    std::tie(nBins, massMin, massMax) = ComputeBinning(mass*0.7, mass*1.3);
    std::map<std::string, TH1*> hQA1 = {
        {"mass", new TH1D("hMass1", ";#it{M} (GeV/#it{c}^{2});Counts", nBins, massMin, massMax)},
        {"pt", new TH1D("hPt1", ";#it{p}_{T} (GeV/#it{c});Counts", 1000, 0, 10)},
        {"y", new TH1D("hY1", ";#it{y};Counts", 200, -10, 10)},
        {"eta", new TH1D("hEta1", ";#eta;Counts", 200, -10, 10)},
    };

    // Single-particle QA for Mother particle
    if (cfg["injection"].size() > 0) {
        mass = cfg["injection"][0]["mass"].as<double>();
    } else {
        mass = pythia.particleData.particleDataEntryPtr(pdgMother).get()->m0();
    }
    std::tie(nBins, massMin, massMax) = ComputeBinning(mass*0.7, mass*1.3);
    std::map<std::string, TH1*> hQAMother = {
        {"mass", new TH1D("hMassMother", ";#it{M} (GeV/#it{c}^{2});Counts", nBins, massMin, massMax)},
        {"pt", new TH1D("hPtMother", ";#it{p}_{T} (GeV/#it{c});Counts", 1000, 0, 10)},
        {"y", new TH1D("hYMother", ";#it{y};Counts", 200, -10, 10)},
        {"eta", new TH1D("hEtaMother", ";#eta;Counts", 200, -10, 10)},
    };

    std::vector<int> part0{};
    std::vector<int> part1{};
    std::deque<std::vector<Pythia8::Particle>> partBuffer{};

    for (size_t iEvent = 0; iEvent < nEvents; iEvent++) {
        part0.clear();
        part1.clear();

        DEBUG("\n\nGenerating a new event\n");
        if (cfg["injection"].IsDefined() && cfg["injection"].IsSequence() && cfg["injection"].size() == 0) {
            pythia.next();
        } else if (cfg["injection"].IsDefined() && cfg["injection"].IsSequence() && cfg["injection"].size() > 0) {
            pythia.event.reset();

            for (const auto& inj : cfg["injection"]) {
                DEBUG("\n\nInjecting a new particle\n");

                int myPdg = inj["pdg"].as<int>();
                double mass = fLineShape->GetRandom();;
                double pt;
                double y = std::nan("");
                double eta = std::nan("");

                if (hYvsPt) {
                    hYvsPt->GetRandom2(pt, y);
                } else if ((hPt || fPt) && (hY || fY || hEta || fEta)) {
                    if (hPt) {
                        pt = hPt->GetRandom();
                    } else if (fPt) {
                        pt = fPt->GetRandom();
                    } else {
                        printf("Error: undefined pt distribution for injected particle. Exit!\n");
                        exit(1);
                    }

                    if (hY) {
                        y = hY->GetRandom();
                    } else if (fY) {
                        y = fY->GetRandom();
                    } else if (hEta) {
                        eta = hEta->GetRandom();
                    } else if (fEta) {
                        eta = fEta->GetRandom();
                    } else {
                        printf("Error: both rapidity and pseudo-rapidity distribution for injected particle are undefined. Exit!\n");
                        exit(1);
                    }
                } else {
                    printf("Error: pt and/or y distributions are not properly defined. Exit!\n");
                    exit(1);
                }
                double phi = gRandom->Uniform(2 * TMath::Pi());
                double tau = gRandom->Exp(1);
                double mt = TMath::Sqrt(mass * mass + pt * pt);
                double pz;
                if (y == y) {
                    pz = mt * TMath::SinH(y);
                } else if (eta == eta) {
                    pz = pt * TMath::SinH(eta);
                } else {
                    printf("Error: pt and/or y distributions are not properly defined. Exit!\n");
                    exit(1);
                }

                Pythia8::Particle myPart;
                myPart.id(myPdg);
                myPart.status(81);
                myPart.m(mass);
                myPart.xProd(0.);
                myPart.yProd(0.);
                myPart.zProd(0.);
                myPart.tProd(0.);
                myPart.e(TMath::Sqrt(mt * mt + pz * pz));
                myPart.px(pt * TMath::Cos(phi));
                myPart.py(pt * TMath::Sin(phi));
                myPart.pz(pz);
                myPart.tau(tau);

                // remove all particles generated in the event and append the Lambda(1520)
                pythia.event.append(myPart);
                pythia.particleData.mayDecay(myPdg, true);
            }
            
            // force the decay of the Lambda
            pythia.moreDecays();
        } else {
            std::cerr << "Error in injection configuration. Exit!" << std::endl;
            exit(1);
        }

        const std::vector<Pythia8::Particle> *particles;
        if (pythia.event.size() > 0) {
            particles = pythia.event.particles();
        }

        if (!Trigger(particles, trigger)) {
            continue;
        }

        hEvt->Fill(1);

        // Part 0 is the event, 1 and 2 the beams. In case the hadrons are injected there are no beam particles
        for (size_t iPart = 1; iPart < particles->size(); iPart++) {
            auto &part = (*particles)[iPart];

            if (cfg["decaychain"]["enable"].as<bool>()) {
                if (!IsSelected(particles, iPart, cfgMother)) continue;

                DEBUG("\n\n==========================================================================================================\n");
                if (GetParticlesInDecayChain(particles, iPart, cfg["decaychain"]["daus"], part0, part1)) {
                part0.erase(std::remove_if(part0.begin(), part0.end(), [&particles](int iPart){return !IsSelected(particles, iPart, cfgPart0);}), part0.end());
                part1.erase(std::remove_if(part1.begin(), part1.end(), [&particles](int iPart){return !IsSelected(particles, iPart, cfgPart1);}), part1.end());

                // Fill QA
                hQAMother["mass"]->Fill(part.m());
                hQAMother["pt"]->Fill(part.pT());
                hQAMother["y"]->Fill(part.y());
                hQAMother["eta"]->Fill(part.eta());

                DEBUG("size after loading particles: %zu %zu\n", part0.size(), part1.size());
                break;
                } else {
                    part0.clear();
                    part1.clear();
                }
            } else {
                if (IsSelected(particles, iPart, cfgPart0)) {
                    part0.push_back(iPart);
                } else if (IsSelected(particles, iPart, cfgPart1)) {
                    part1.push_back(iPart);
                }
            }
        }

        // Fill QA for Part0
        for (const int& i0 : part0) {
            hQA0["mass"]->Fill((*particles)[i0].m());
            hQA0["pt"]->Fill((*particles)[i0].pT());
            hQA0["y"]->Fill((*particles)[i0].y());
            hQA0["eta"]->Fill((*particles)[i0].eta());
        }

        // Fill QA for Part1
        for (const int& i1 : part1) {
            hQA1["mass"]->Fill((*particles)[i1].m());
            hQA1["pt"]->Fill((*particles)[i1].pT());
            hQA1["y"]->Fill((*particles)[i1].y());
            hQA1["eta"]->Fill((*particles)[i1].eta());
        }

        // Skip events without pairs
        if (rejevtwopairs && (part0.size() == 0 || part1.size() == 0)) continue;

        int mult = ComputeMultTPC(particles);
        hEvtMult->Fill(mult);

        int mult0plus = std::count_if(part0.begin(), part0.end(), [&particles](int iPart) { return (*particles)[iPart].id() > 0; });
        int mult1plus = std::count_if(part1.begin(), part1.end(), [&particles](int iPart) { return (*particles)[iPart].id() > 0; });
        hPairMultSE[std::tuple<int, int, int>({0, 0, -1})]->Fill(mult0plus, mult1plus);
        if (nPart0>1) hPairMultSE[std::tuple<int, int, int>({1, 0, -1})]->Fill(part0.size() - mult0plus, mult1plus);
        if (nPart1>1) hPairMultSE[std::tuple<int, int, int>({0, 1, -1})]->Fill(mult0plus, part1.size() - mult1plus);
        if (nPart0 > 1 && nPart1 > 1) hPairMultSE[std::tuple<int, int, int>({1, 1, -1})]->Fill(part0.size() - mult0plus, part1.size() - mult1plus);

        DEBUG("Particle multiplicities in this event: n(%d)=%zu, n(%d)=%zu\n", pdg0, part0.size(), pdg1, part1.size());

        //! WARNING: this value makes sense ONLY for events with a SINGLE injected resonance! Don't use for realistic events
        double ptMother = (*particles)[1].pT(); // The injected particle is always at index=1.

        // Same event
        DEBUG("Start same-event pairing\n");
        for (size_t i0 = 0; i0 < part0.size(); i0++) {
            const auto p0 = (*particles)[part0[i0]];

            double eff0 = 1;
            if (fEff0) {
                eff0 = fEff0->Eval(p0.pT());
            } else if (hEff0) {
                eff0 = hEff0->GetBinContent(hEff0->FindBin(p0.pT()));
            }

            // don't pair twice in case of same-part femto
            int start = pdg0 == pdg1 ? i0 + 1 : 0;
            auto buffer = pdg0 == pdg1 ? part0 : part1;
            for (size_t i1 = start; i1 < buffer.size(); i1++) {
                const auto p1 = (*particles)[buffer[i1]];
                double kStar = ComputeKstar(p0, p1);
                double mT = ComputeMt(p0, p1);
                int iMt = FindBin(mT, mTBins);
                
                std::tuple<int, int, int> pair = {pdg0 == pdg1 ? 0 : p0.id() < 0, pdg0 == pdg1 ? p0.id() * p1.id() < 0 : p1.id() < 0, iMt};
                DEBUG("    SE(idx=%zu, idx=%zu): pdg0=%d pdg1=%d  --->  (%d, %d)\n", i0, i1, p0.id(), p1.id(), pair.first, pair.second);

                double eff1 = 1;
                if (fEff1) {
                    eff1 = fEff1->Eval(p1.pT());
                } else if (hEff1) {
                    eff1 = hEff1->GetBinContent(hEff1->FindBin(p1.pT()));
                }

                std::string ancestor = haveCommonAncestor(p0, p1) ? "Common" : "NonCommon";
                // std::cout << std::get<0>(pair) << "  " << std::get<1>(pair) << "  "<< std::get<2>(pair) << "  " << ancestor <<std::endl;
                hSE[pair][ancestor]->Fill(kStar, eff0 * eff1);
                hPtMotherVsKstar[pair]->Fill(kStar, ptMother, eff0 * eff1);
            }
        }

        // Mixed event
        DEBUG("Start mixed-event pairing\n");
        for (size_t i0 = 0; i0 < part0.size(); i0++) {
            const auto p0 = pythia.event[part0[i0]];

            double eff0 = 1;
            if (fEff0) {
                eff0 = fEff0->Eval(p0.pT());
            } else if (hEff0) {
                eff0 = hEff0->GetBinContent(hEff0->FindBin(p0.pT()));
            }

            for (size_t iME = 0; iME < partBuffer.size(); iME++) {
                for (size_t i1 = 0; i1 < partBuffer[iME].size(); i1++) {
                    const auto p1 = partBuffer[iME][i1];
                    double kStar = ComputeKstar(p0, p1);
                    double mT = ComputeMt(p0, p1);
                    int iMt = FindBin(mT, mTBins);

                    std::tuple<int, int, int> pair = {pdg0 == pdg1 ? 0 : p0.id() < 0, pdg0 == pdg1 ? p0.id() * p1.id() < 0 : p1.id() < 0, iMt};
                    DEBUG("    ME(idx0=%zu, iMix=%zu, idx1=%zu): pdg0=%d pdg1=%d   (%d, %d)\n", i0, iME, i1, p0.id(), p1.id(), pair.first, pair.second);

                    double eff1 = 1;
                    if (fEff1) {
                        eff1 = fEff1->Eval(p1.pT());
                    } else if (hEff1) {
                        eff1 = hEff1->GetBinContent(hEff1->FindBin(p1.pT()));
                    }

                    hME[pair]->Fill(kStar, eff0 * eff1);
                }
            }
        }

        // todo: check if the buffer type(0 or 1) changes anything
        auto partX = pdg0 == pdg1 ? part0 : part1;
        partBuffer.push_back({});
        for (const auto &iPart : partX) partBuffer[partBuffer.size() - 1].push_back((*particles)[iPart]);
        if (partBuffer.size() > md) partBuffer.pop_front();
    }

    TFile *oFile = TFile::Open(oFileName.data(), "recreate");
    if (!oFile) exit(1);

    hEvtMult->Write();
    hEvt->Write();

    oFile->mkdir("qa0");
    oFile->cd("qa0");
    for (const auto &[key, hist] : hQA0) {
        hist->Write();
    }
    oFile->mkdir("qa1");
    oFile->cd("qa1");
    for (const auto &[key, hist] : hQA1) {
        hist->Write();
    }
    oFile->mkdir("qaMother");
    oFile->cd("qaMother");
    fLineShape->Write();
    for (const auto &[key, hist] : hQAMother) {
        hist->Write();
    }

    for (int iPart0 = 0; iPart0 < nPart0; iPart0++) {
        for (int iPart1 = 0; iPart1 < nPart1; iPart1++) {
            std::tuple<int, int, int> pair = {iPart0, iPart1, -1};
            std::string pairName = Form("p%d%d", iPart0, iPart1);
            oFile->mkdir(pairName.data());
            oFile->cd(pairName.data());
            TH2D* hSETot = (TH2D *) hSE[pair]["Common"]->Clone("hSE");
            hSETot->Add(hSE[pair]["NonCommon"]);
            hSETot->Write("hSE");
            hME[pair]->Write("hME");

            hSE[pair]["Common"]->Write("hSECommon");
            hSE[pair]["NonCommon"]->Write("hSENonCommon");
            hPairMultSE[pair]->Write("hPairMult");
            hPtMotherVsKstar[pair]->Write("hPtMotherVsKstar");

            for (size_t iMt = 0; iMt < mTMins.size(); iMt++) {
                std::tuple<int, int, int> pair = {iPart0, iPart1, iMt};
                std::string pairName = Form("p%d%d/mT%zu", iPart0, iPart1, iMt);
                oFile->mkdir(pairName.data());
                oFile->cd(pairName.data());

                TH2D* hSETot = (TH2D *) hSE[pair]["Common"]->Clone("hSE");
                hSETot->Add(hSE[pair]["NonCommon"]);
                hSETot->Write("hSE");
                hME[pair]->Write("hME");

                hSE[pair]["Common"]->Write("hSECommon");
                hSE[pair]["NonCommon"]->Write("hSENonCommon");
                hPairMultSE[pair]->Write("hPairMult");
                hPtMotherVsKstar[pair]->Write("hPtMotherVsKstar");
            }
        }
    }

    oFile->Close();
    std::cout << "Output saved in " << oFileName << std::endl;
}


int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Expected 3 parameters" << std::endl;
        return 1;
    }
    std::string inFileName(argv[1]);
    std::string oFileName(argv[2]);
    std::string cfg(argv[3]);

    MakeDistr(inFileName, oFileName, cfg, -1);
    return 0;
}
