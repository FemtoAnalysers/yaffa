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
#include "Math/Vector4D.h"
#include "Math/Boost.h"
#include "TRandom3.h"

// ALICE libraries
#include "Pythia8/Pythia.h"

#if false
#define DEBUG(msg, ...) do { printf(msg, ##__VA_ARGS__); } while(0)
#define DEBUG_VAR(var) do { std::cerr << #var << ": " << var << std::endl; } while (0)
#else
#define DEBUG(msg, ...)
#define DEBUG_VAR(var)
#endif

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
void SetDefault(YAML::Node node, T1 min, T2 max) {
    if (node.IsNull()) {
        node.push_back(min);
        node.push_back(max);
    }
}


YAML::Node ApplyDefaults(const YAML::Node &cfg) {
    YAML::Node cfgNew = cfg;

    SetDefault(cfgNew["status"], -300, 300);
    SetDefault(cfgNew["pt"], 0, 100);
    SetDefault(cfgNew["eta"], -10, 10);
    SetDefault(cfgNew["y"], -10, 10);
    SetDefault(cfgNew["prodvtx"], -1., 1000);

    return cfgNew;
}


template<typename T>
bool IsSelected(T value, const YAML::Node node, bool includeExtremes = false) {
    std::array<T, 2> range = node.as<std::array<T, 2>>();
    if (includeExtremes) return range[0] <= value && value <= range[1];
    return range[0] < value && value < range[1];
}

bool IsSelected(const Pythia8::Pythia &pythia, int iPart, const YAML::Node &cfgSelections) {
    auto& part = pythia.event[iPart];

    if (!IsSelected(part.status(), cfgSelections["status"], true)) return false;
    if (!IsSelected(part.pT(), cfgSelections["pt"])) return false;
    if (!IsSelected(part.eta(), cfgSelections["eta"])) return false;
    if (!IsSelected(part.y(), cfgSelections["y"])) return false;
    double prodvtx = pow(part.xProd() * part.xProd() + part.yProd() * part.yProd() + part.zProd() * part.zProd(), 0.5);
    if (!IsSelected(prodvtx, cfgSelections["prodvtx"])) return false;

    return true;
}


void MakeDistr(
    std::string oFileName = "Distr.root",
    std::string cfgFile = "cfg_makedistr_example.yml",
    int seed = 31
    ) {
    // Load simulation settings

    YAML::Node cfg = YAML::LoadFile(cfgFile.data());
    size_t nEvents = cfg["nevts"].as<unsigned int>();
    unsigned int md = cfg["mixdepth"].as<int>();
    bool rejevtwopairs = cfg["rejevtwopairs"].as<bool>();

    // Load selections for part 0
    cfgPart0 = ApplyDefaults(cfg["part0"]);
    auto pdg0 = cfgPart0["pdg"].as<int>();

    // Load selections for part 1. If they don't exist, use the same as part 0 (same-particle pairs e.g. pp)
    cfgPart1 = cfg["part1"].IsNull() ? cfgPart0 : ApplyDefaults(cfg["part1"]);
    auto pdg1 = cfgPart1["pdg"].as<int>();

    Pythia8::Pythia pythia;

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

    auto PDG = TDatabasePDG::Instance();
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

    std::vector<Pythia8::Particle> part0{};
    std::vector<Pythia8::Particle> part1{};
    std::deque<std::vector<Pythia8::Particle>> partBuffer{};

    for (size_t iEvent = 0; iEvent < nEvents; iEvent++) {
        part0.clear();
        part1.clear();

        DEBUG("\n\nGenerating a new event\n");
        pythia.next();

        // Part 0 is the event, 1 and 2 the beams ==> start from 3
        for (int iPart = 3; iPart < pythia.event.size(); iPart++) {
            auto &part = pythia.event[iPart];

            int pdg = part.id();
            int absPdg = std::abs(pdg);

            if (absPdg == pdg0 && IsSelected(pythia, iPart, cfgPart0)) {
                part0.push_back(pythia.event[iPart]);
            } else if (absPdg == pdg1 && IsSelected(pythia, iPart, cfgPart1)) {
                part1.push_back(pythia.event[iPart]);
            }
        }

        DEBUG("Particles found for pdg0=%d: %zu\n", pdg0, part0.size());
        for (size_t iPart = 0; iPart < part0.size(); iPart++) {
            auto part = part0[iPart];
            DEBUG("    %zu: pdg=%d  px=%.4f\n", iPart, part.id(), part.px());
        }
        DEBUG("Particles found for pdg1=%d: %zu\n", pdg1, part1.size());
        for (size_t iPart = 0; iPart < part1.size(); iPart++) {
            auto part = part1[iPart];
            DEBUG("    %zu: pdg=%d  px=%.4f\n", iPart, part.id(), part.px());
        }

        // Skip events without pairs // todo: reimplement
        // if (cfg["rejevtwopairs"] && (part0.size() == 0 || part1.size() == 0)) continue;

        int mult = ComputeMultTPC(pythia);
        hEvtMult->Fill(mult);

        int mult0plus = std::count_if(part0.begin(), part0.end(), [](auto p) { return p.id() > 0; });
        int mult1plus = std::count_if(part1.begin(), part1.end(), [](auto p) { return p.id() > 0; });
        hPairMultSE[std::pair<int, int>({0, 0})]->Fill(mult0plus, mult1plus);
        if (nPart0>1) hPairMultSE[std::pair<int, int>({1, 0})]->Fill(part0.size() - mult0plus, mult1plus);
        if (nPart1>1) hPairMultSE[std::pair<int, int>({0, 1})]->Fill(mult0plus, part1.size() - mult1plus);
        if (nPart0 > 1 && nPart1 > 1) hPairMultSE[std::pair<int, int>({1, 1})]->Fill(part0.size() - mult0plus, part1.size() - mult1plus);

        DEBUG("Particle multiplicities in this event: n(%d)=%zu, n(%d)=%zu\n", pdg0, part0.size(), pdg1, part1.size());

        // Same event
        DEBUG("Start same-event pairing\n");
        for (size_t i0 = 0; i0 < part0.size(); i0++) {
            const auto p0 = part0[i0];

            // don't pair twice in case of same-part femto
            int start = pdg0 == pdg1 ? i0 + 1 : 0;
            auto buffer = pdg0 == pdg1 ? part0 : part1;
            for (size_t i1 = start; i1 < buffer.size(); i1++) {
                const auto p1 = buffer[i1];
                double kStar = ComputeKstar(p0, p1);
                std::pair<int, int> pair = {pdg0 == pdg1 ? 0 : p0.id() < 0, pdg0 == pdg1 ? p0.id() * p1.id() < 0 : p1.id() < 0};
                DEBUG("    SE(idx=%zu, idx=%zu): pdg0=%d pdg1=%d  --->  (%d, %d)\n", i0, i1, p0.id(), p1.id(), pair.first, pair.second);

                hSE[pair]->Fill(kStar);
            }
        }

        // Mixed event
        DEBUG("Start mixed-event pairing\n");
        for (size_t i0 = 0; i0 < part0.size(); i0++) {
            const auto p0 = part0[i0];

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
        partBuffer.push_back(pdg0 == pdg1 ? part0 : part1);
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
