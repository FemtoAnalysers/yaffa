/*
Script to compute the same- and mixed event distributions from pythia events.
*/

// C++ libraries
#include <array>
#include <deque>
#include <map>
#include <string>
#include <vector>

// 3rd party libraries
#include "yaml-cpp/yaml.h"

// ROOT libraries
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TMath.h"
#include "Math/Vector4D.h"
#include "Math/Boost.h"
#include "TRandom3.h"

// ALICE libraries
#include "Pythia8/Pythia.h"
#include "Pythia8/ParticleData.h"

#if false
#define DEBUG(msg, ...) do { printf(msg, ##__VA_ARGS__); } while(0)
#define DEBUG_VAR(var) do { std::cerr << #var << ": " << var << std::endl; } while (0)
#else
#define DEBUG(msg, ...)
#define DEBUG_VAR(var)
#endif

const double HBARC =  197.3269804; // um=MeV*fm

static YAML::Node cfgPart0;
static YAML::Node cfgPart1;

enum processes { kSoftQCD = 0, kNonDiffractive, kHardQCD };
std::map<processes, const char*> processNames = {
    {kSoftQCD, "SoftQCD"},
    {kNonDiffractive, "NonDiffractive"},
    {kHardQCD, "HardQCD"},
};

// todo: not implemented
enum triggers { kMB = 0, kHM };
std::map<triggers, const char*> triggerNames = {
    {kMB, "MB"},
    {kHM, "HM"},
};

float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2) {
    ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
    ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
    ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
    ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);

    ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;
    float kStar = 0.5 * trackRelK.P();
    return kStar;
}

float ComputeKstar(Pythia8::Particle part1, Pythia8::Particle part2) {
    ROOT::Math::PxPyPzMVector p1(part1.px(), part1.py(), part1.pz(), part1.m());
    ROOT::Math::PxPyPzMVector p2(part2.px(), part2.py(), part2.pz(), part2.m());

    return ComputeKstar(p1, p2);
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

// Compute the multiplicity of charged particles in the TPC acceptance
int ComputeMultTPC(const Pythia8::Pythia &pythia) {
    int mult = 0;
    for (int iPart = 3; iPart < pythia.event.size(); iPart++) {
        Pythia8::Particle part = pythia.event[iPart];

        if (part.isFinal() && std::abs(part.eta()) < 0.8 && IsDetectable(part.id())) mult++;
    }
    return mult;
}

// Set the default selections
template<typename T1, typename T2>
void SetDefault(YAML::Node &&node, T1 min, T2 max) {
    if (node.IsNull()) {
        node.push_back(min);
        node.push_back(max);
    }
}

void SetDefaults(YAML::Node &cfg) {
    SetDefault(cfg["status"], -300, 300);
    SetDefault(cfg["pt"], 0, 100);
    SetDefault(cfg["eta"], -10, 10);
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

bool IsSelected(const Pythia8::Pythia &pythia, int iPart, const YAML::Node &cfgSelections) {
    auto& part = pythia.event[iPart];

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
            if (!IsSelected(pythia, dauIdx[iDau], cfgSelections["daus"][iDau])) return false;
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

void GetParticlesInDecayChain(const Pythia8::Pythia &pythia, int iPart, YAML::Node cfgMom, std::vector<int> &part0, std::vector<int> &part1) {
    const auto& mom = pythia.event[iPart];

    DEBUG("Start analyzing the decay tree of pdg=%d idx=%d\n", mom.id(), iPart);
    DEBUG("Start analyzing daus:\n");
    for (int iDau = mom.daughter1(); iDau <= mom.daughter2(); iDau++) {
        auto dau = pythia.event[iDau];

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
            return;
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
            GetParticlesInDecayChain(pythia, iDau, dauCfg["daus"], part0, part1);
        }
    }
}

void MakeDistr(
    std::string oFileName = "Distr.root",
    std::string cfgFile = "cfg_makedistr_example.yml",
    int seed = 31
    ) {
    // Load PDG
    auto PDG = TDatabasePDG::Instance();

    // Load simulation settings
    YAML::Node cfg = YAML::LoadFile(cfgFile.data());
    size_t nEvents = cfg["nevts"].as<unsigned int>();
    unsigned int md = cfg["mixdepth"].as<int>();
    bool rejevtwopairs = cfg["rejevtwopairs"].as<bool>();

    // Load selections for part 0
    cfgPart0 = cfg["part0"];
    // Load selections for part 1. If they don't exist, use the same as part 0 (same-particle pairs e.g. pp)
    cfgPart1 = cfg["part1"].IsNull() ? cfgPart0 : cfg["part1"];

    std::cout << "\033[34mParticle selections before defaults (part0)\033[0m" << std::endl;
    std::cout << cfgPart0<< std::endl;

    std::cout << "\033[34mParticle selections before defaults (part1)\033[0m" << std::endl;
    std::cout << cfgPart1<< std::endl;

    SetDefaults(cfgPart0);
    SetDefaults(cfgPart1);

    std::cout << "\033[34mParticle selections after defaults (part0)\033[0m" << std::endl;
    std::cout << cfgPart0<< std::endl;

    std::cout << "\033[34mParticle selections after defaults (part1)\033[0m" << std::endl;
    std::cout << cfgPart1<< std::endl;

    auto pdg0 = cfgPart0["pdg"].as<int>();
    auto pdg1 = cfgPart1["pdg"].as<int>();

    Pythia8::Pythia pythia;
    pythia.readString("Next:numberShowEvent = 0");

    // Set process
    if (cfg["process"].as<std::string>() == kSoftQCD) {
        pythia.readString("SoftQCD:all = on");
    } else if (cfg["process"].as<std::string>() == "HardQCD") {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    } else if (cfg["process"].as<std::string>() == "NonDiffractive") {
        pythia.readString("SoftQCD:nonDiffractive = on");
    } else {
        std::cout << "\033[31mError: Process not implemented. Exit!\033[0m" << std::endl;
        exit(1);
    }

    // Set tune
    if (cfg["tune"].as<std::string>() == "Monash") {
        pythia.readString(Form("Tune:pp = 14"));
    } else if (cfg["tune"].as<std::string>() == "CRMode0") {
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
    } else if (cfg["tune"].as<std::string>() == "CRMode2") {
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
    } else if (cfg["tune"].as<std::string>() == "CRMode3") {
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
    } else {
        std::cout << "\033[31mError: Tune not implemented. Exit!\033[0m" << std::endl;
        exit(1);
    }

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
        pythia.particleData.readString(Form("%d:addChannel = 1 1 0 %s", myPdg, part["daus"].as<std::string>().data()));
    }

    std::cout << "Applying the following customization to pythia:" << std::endl;
    for (const auto &line : cfg["customization"]) {
        std::string lineStr = line.as<std::string>();
        std::cout << "   * " << lineStr << std::endl;
        pythia.readString(lineStr.data());
    }
    std::cout << "End of customization." << std::endl;

    // Compute the minimum mass that a resonance modelled with a Breit-Wigner can have
    double minBWMass=1.e12;
    auto particleEntry = pythia.particleData.particleDataEntryPtr(cfg["injection"][0]["pdg"].as<int>());
    for (int iDecChn = 0; iDecChn < particleEntry->sizeChannels(); iDecChn++) {
        auto channel = particleEntry->channel(iDecChn);

        // Loop over the daughters in the decay channel
        double sum = 0;
        for (int iDau = 0; iDau < channel.multiplicity(); iDau++) {
            int dauPdg = channel.product(iDau);
            sum += PDG->GetParticle(dauPdg)->Mass();
        }
        if (sum < minBWMass) minBWMass = sum;   
    }
    printf("Mass limit:. %.3f", minBWMass);

    // Setting the seed here is not sufficient to ensure reproducibility, setting the seed of gRandom is necessary
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", cfg["sqrts"].as<double>() * 1000); // from TeV to GeV
    pythia.init();

    gRandom->SetSeed(seed); // Set the seed to ensure reproducibility of the evevents generated by pythia

    // QA histograms
    TH1D * hEvtMult = new TH1D("hEvtMult", ";#it{N}_{ch}|_{|#eta|<0.8};Counts", 100, 0., 100);

    // Pairs
    std::map<std::pair<int, int>, TH1D *> hSE, hME;
    
    // Pair QA
    std::map<std::pair<int, int>, TH2D *> hPairMultSE;

    int nPart0 = PDG->GetParticle(pdg0)->AntiParticle() && pdg0 != pdg1 ? 2 : 1;
    int nPart1 = PDG->GetParticle(pdg1)->AntiParticle() ? 2 : 1;
    for (int iPart0 = 0; iPart0 < nPart0; iPart0++) {
        for (int iPart1 = 0; iPart1 < nPart1; iPart1++) {
            std::pair<int, int> pair = {iPart0, iPart1};
            DEBUG("Inserting histograms for pair (%d, %d)\n", iPart0, iPart1);
            hSE.insert({pair, new TH1D(Form("hSE%d%d", iPart0, iPart1), ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)});
            hME.insert({pair, new TH1D(Form("hME%d%d", iPart0, iPart1), ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)});
            hPairMultSE.insert({pair, new TH2D(Form("hPairMultSE%d%d", iPart0, iPart1), ";#it{N}_{0};#it{N}_{1};Counts", 51, -0.5, 50.5, 31, -0.5, 50.5)});
        }
    }

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
                double mass;
                do {
                    mass = gRandom->BreitWigner(inj["mass"].as<double>(), inj["width"].as<double>()/1000);
                } while (mass < minBWMass*1.001); // Correction factor needed for precision issues

                double pt = gRandom->Exp(1);
                double y = gRandom->Gaus();
                double phi = gRandom->Uniform(2 * TMath::Pi());
                double tau = gRandom->Exp(1);
                double mt = TMath::Sqrt(mass * mass + pt * pt);
                double pz = TMath::SinH(y) * mt;

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
            cerr << "Error in injection configuration. Exit!" << std::endl;
            exit(1);
        }

        // Part 0 is the event, 1 and 2 the beams. In case the hadrons are injected there are no beam particles
        for (int iPart = 1; iPart < pythia.event.size(); iPart++) {
            auto &part = pythia.event[iPart];

            int pdg = part.id();
            int absPdg = std::abs(pdg);

            if (cfg["decaychain"]["enable"].as<bool>()) {
                int pdgMother = cfg["decaychain"]["pdg"].as<int>();
                if (absPdg != std::abs(pdgMother)) continue;

                DEBUG("\n\n==========================================================================================================\n");
                GetParticlesInDecayChain(pythia, iPart, cfg["decaychain"]["daus"], part0, part1);
                part0.erase(std::remove_if(part0.begin(), part0.end(), [&pythia](int iPart){return !IsSelected(pythia, iPart, cfgPart0);}), part0.end());
                part1.erase(std::remove_if(part1.begin(), part1.end(), [&pythia](int iPart){return !IsSelected(pythia, iPart, cfgPart1);}), part1.end());

                DEBUG("size after loading particles: %zu %zu\n", part0.size(), part1.size());

                break;
            } else {
                if (IsSelected(pythia, iPart, cfgPart0)) {
                    part0.push_back(iPart);
                } else if (IsSelected(pythia, iPart, cfgPart1)) {
                    part1.push_back(iPart);
                }
            }
        }

        // Skip events without pairs
        if (cfg["rejevtwopairs"] && (part0.size() == 0 || part1.size() == 0)) continue;

        int mult = ComputeMultTPC(pythia);
        hEvtMult->Fill(mult);

        int mult0plus = std::count_if(part0.begin(), part0.end(), [&pythia](int iPart) { return pythia.event[iPart].id() > 0; });
        int mult1plus = std::count_if(part1.begin(), part1.end(), [&pythia](int iPart) { return pythia.event[iPart].id() > 0; });
        hPairMultSE[std::pair<int, int>({0, 0})]->Fill(mult0plus, mult1plus);
        if (nPart0>1) hPairMultSE[std::pair<int, int>({1, 0})]->Fill(part0.size() - mult0plus, mult1plus);
        if (nPart1>1) hPairMultSE[std::pair<int, int>({0, 1})]->Fill(mult0plus, part1.size() - mult1plus);
        if (nPart0 > 1 && nPart1 > 1) hPairMultSE[std::pair<int, int>({1, 1})]->Fill(part0.size() - mult0plus, part1.size() - mult1plus);

        DEBUG("Particle multiplicities in this event: n(%d)=%zu, n(%d)=%zu\n", pdg0, part0.size(), pdg1, part1.size());

        // Same event
        DEBUG("Start same-event pairing\n");
        for (size_t i0 = 0; i0 < part0.size(); i0++) {
            const auto p0 = pythia.event[part0[i0]];

            // don't pair twice in case of same-part femto
            int start = pdg0 == pdg1 ? i0 + 1 : 0;
            auto buffer = pdg0 == pdg1 ? part0 : part1;
            for (size_t i1 = start; i1 < buffer.size(); i1++) {
                const auto p1 = pythia.event[buffer[i1]];
                double kStar = ComputeKstar(p0, p1);
                std::pair<int, int> pair = {pdg0 == pdg1 ? 0 : p0.id() < 0, pdg0 == pdg1 ? p0.id() * p1.id() < 0 : p1.id() < 0};
                DEBUG("    SE(idx=%zu, idx=%zu): pdg0=%d pdg1=%d  --->  (%d, %d)\n", i0, i1, p0.id(), p1.id(), pair.first, pair.second);

                hSE[pair]->Fill(kStar);
            }
        }

        // Mixed event
        DEBUG("Start mixed-event pairing\n");
        for (size_t i0 = 0; i0 < part0.size(); i0++) {
            const auto p0 = pythia.event[part0[i0]];

            for (size_t iME = 0; iME < partBuffer.size(); iME++) {
                for (size_t i1 = 0; i1 < partBuffer[iME].size(); i1++) {
                    const auto p1 = partBuffer[iME][i1];
                    double kStar = ComputeKstar(p0, p1);
                    std::pair<int, int> pair = {pdg0 == pdg1 ? 0 : p0.id() < 0, pdg0 == pdg1 ? p0.id() * p1.id() < 0 : p1.id() < 0};
                    DEBUG("    ME(idx0=%zu, iMix=%zu, idx1=%zu): pdg0=%d pdg1=%d   (%d, %d)\n", i0, iME, i1, p0.id(), p1.id(), pair.first, pair.second);

                    hME[pair]->Fill(kStar);
                }
            }
        }

        // todo: check if the buffer type(0 or 1) changes anything
        auto partX = pdg0 == pdg1 ? part0 : part1;
        partBuffer.push_back({});
        for (const auto &iPart : partX) partBuffer[partBuffer.size() - 1].push_back(pythia.event[iPart]);
        if (partBuffer.size() > md) partBuffer.pop_front();
    }

    TFile *oFile = new TFile(oFileName.data(), "recreate");
    hEvtMult->Write();

    for (int iPart0 = 0; iPart0 < nPart0; iPart0++) {
        for (int iPart1 = 0; iPart1 < nPart1; iPart1++) {
            std::pair<int, int> pair = {iPart0, iPart1};
            std::string pairName = Form("p%d%d", iPart0, iPart1);
            oFile->mkdir(pairName.data());
            oFile->cd(pairName.data());
            hSE[pair]->Write("hSE");
            hME[pair]->Write("hME");
            hPairMultSE[pair]->Write("hPairMult");
        }
    }

    oFile->Close();
    std::cout << "Output saved in " << oFileName << std::endl;
}
