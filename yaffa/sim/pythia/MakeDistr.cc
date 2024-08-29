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

    // Load selections for part 0
    YAML::Node cfgPart0 = cfg["part0"];
    auto pdg0 = cfgPart0["pdg"].as<int>();
    auto statusrange0 = cfgPart0["status"].IsNull() ? std::array<int, 2>({-300, 300}) : cfgPart0["status"].as<std::array<int, 2>>();
    auto ptrange0 =     cfgPart0["pt"].IsNull() ? std::array<double, 2>({0, 9.e9}) : cfgPart0["pt"].as<std::array<double, 2>>();
    auto etarange0 =    cfgPart0["eta"].IsNull() ? std::array<double, 2>({-9.e9, 9.e9}) : cfgPart0["eta"].as<std::array<double, 2>>();
    auto yrange0 =      cfgPart0["y"].IsNull() ? std::array<double, 2>({-9.e9, 9.e9}) : cfgPart0["y"].as<std::array<double, 2>>();
    auto prodvtxrange0 = cfgPart0["prodvtx"].IsNull() ? std::array<double, 2>({-1., 9.e9}) : cfgPart0["prodvtx"].as<std::array<double, 2>>();

    // Load selections for part 1. If they don't exist, use the same as part 0 (same-particle pairs e.g. pp)
    YAML::Node cfgPart1 = cfg["part1"].IsNull() ? cfg["part0"] : cfg["part1"];
    auto pdg1 = cfgPart1["pdg"].as<int>();
    auto statusrange1 = cfgPart1["status"].IsNull() ? std::array<int, 2>({-300, 300}) : cfgPart1["status"].as<std::array<int, 2>>();
    auto ptrange1 = cfgPart1["pt"].IsNull() ? std::array<double, 2>({0, 9.e9}) : cfgPart1["pt"].as<std::array<double, 2>>();
    auto etarange1 = cfgPart1["eta"].IsNull() ? std::array<double, 2>({-9.e9, 9.e9}) : cfgPart1["eta"].as<std::array<double, 2>>();
    auto yrange1 = cfgPart1["y"].IsNull() ? std::array<double, 2>({-9.e9, 9.e9}) : cfgPart1["y"].as<std::array<double, 2>>();
    auto prodvtxrange1 = cfgPart1["prodvtx"].IsNull() ? std::array<double, 2>({-1., 9.e9}) : cfgPart1["prodvtx"].as<std::array<double, 2>>();

    // QA histograms
    TH1D * hEvtMult = new TH1D("hEvtMult", ";#it{N}_{ch}|_{|#eta|<0.8};Counts", 100, 0., 100);

    // The number of necessary histogram depends on the nature of the pairs. Possible scenarios:
    //  1) One of the 2 particles is the charge-conjugate of itself (e.g. p-Phi, JPsi-Phi) ==> Only 1 histogram is needed
    //  2) Both particles are different from their charge-conjugate (e.g. p-Lambda, p-p) ==> 2 histograms are needed
    std::map<std::string, TH1D *> hSE = {{"p00", new TH1D("hSE00", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)}};
    std::map<std::string, TH1D *> hME = {{"p00", new TH1D("hME00", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)}};

    std::map<std::string, TH2D *> hPairMultSE = {{"p00", new TH2D("hPairMultSE00", ";#it{N}_{0};#it{N}_{1};Counts", 31, -0.5, 30.5, 31, -0.5, 30.5)}};

    // Check if there is at least one particle that is the charge-conjugate of itself
    auto PDG = TDatabasePDG::Instance();
    bool isSelfCC = !(PDG->GetParticle(pdg0)->AntiParticle() && PDG->GetParticle(pdg1)->AntiParticle());

    if (!isSelfCC) {
        hSE.insert({"p01", new TH1D("hSE01", ";#it{k}* (GeV/#it{c});Counts", 2000, 0., 2.)});
        hME.insert({"p01", new TH1D("hME01", ";#it{k}* (GeV/#it{c});Counts", 2000, 0., 2.)});

        hPairMultSE.insert({"p01", new TH2D("hPairMultSE01", ";#it{N}_{0};#it{N}_{1};Counts", 31, -0.5, 30.5, 31, -0.5, 30.5)});
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

            if (absPdg == pdg0) {
                // Selections for part 0
                if (!(statusrange0[0] <= part.status() && part.status() <= statusrange0[1])) continue;
                if (!(ptrange0[0] < part.pT() && part.pT() < ptrange0[1])) continue;
                if (!(etarange0[0] < part.eta() && part.eta() < etarange0[1])) continue;
                if (!(yrange0[0] < part.y() && part.y() < yrange0[1])) continue;
                double prodvtx = pow(part.xProd() * part.xProd() + 
                                     part.yProd() * part.yProd() + 
                                     part.zProd() * part.zProd(), 0.5);
                if (!(prodvtxrange0[0] < prodvtx && prodvtx < prodvtxrange0[1])) continue;

                part0.push_back(pythia.event[iPart]);
            } else if (absPdg == pdg1) {
                // Selections for Part 1
                if (!(statusrange1[0] <= part.status() && part.status() <= statusrange1[1])) continue;
                if (!(ptrange1[0] < part.pT() && part.pT() < ptrange1[1])) continue;
                if (!(etarange1[0] < part.eta() && part.eta() < etarange1[1])) continue;
                if (!(yrange1[0] < part.y() && part.y() < yrange1[1])) continue;
                double prodvtx = pow(part.xProd() * part.xProd() + 
                                     part.yProd() * part.yProd() + 
                                     part.zProd() * part.zProd(), 0.5);
                if (!(prodvtxrange1[0] < prodvtx && prodvtx < prodvtxrange1[1])) continue;

                part1.push_back(pythia.event[iPart]);
            }
        }

        DEBUG("n(%d)=%zu, n(%d)=%zu\n", pdg0, part0.size(), pdg1, part1.size());

        // Skip events without pairs
        if (cfg["rejevtwopairs"] && (part0.size() == 0 || part1.size() == 0)) continue;

        int mult = ComputeMultTPC(pythia);
        hEvtMult->Fill(mult);

        // The single-particle multiplicities are computed differently depending on the type of pair
        if (isSelfCC) { // At least one of the particles is its charge-conjugate e.g. phi-p, phi-Pi0
            hPairMultSE["p00"]->Fill(part0.size(), part1.size());
        } else {
            int mult0plus = std::count_if(part0.begin(), part0.end(), [](auto p) { return p.id() > 0; });
            int mult0minus = part0.size() - mult0plus;

            if (pdg0 == pdg1) { // Same particle femto e.g. p-p or Lambda-Lambda
                hPairMultSE["p00"]->Fill(mult0plus, mult0plus);
                hPairMultSE["p00"]->Fill(mult0minus, mult0minus);

                hPairMultSE["p01"]->Fill(mult0plus, mult0minus);
            } else { // Different particle femto e.g. p-Pi, Lambda-K
                int mult1plus = std::count_if(part1.begin(), part1.end(), [](auto p) { return p.id() > 0; });
                int mult1minus = part1.size() - mult1plus;

                hPairMultSE["p00"]->Fill(mult0plus, mult1plus);
                hPairMultSE["p00"]->Fill(mult0minus, mult1minus);
                hPairMultSE["p01"]->Fill(mult0plus, mult1minus);
                hPairMultSE["p01"]->Fill(mult0minus, mult1plus);
            }
        }

        // Same event
        for (size_t i0 = 0; i0 < part0.size(); i0++) {
            const auto p0 = part0[i0];

            // don't pair twice in case of same-part femto
            int start = pdg0 == pdg1 ? i0 + 1 : 0;
            auto buffer = pdg0 == pdg1 ? part0 : part1;
            for (size_t i1 = start; i1 < buffer.size(); i1++) {
                const auto p1 = buffer[i1];

                double kStar = ComputeKstar(p0, p1);
                std::string pair = p0.id() * p1.id() > 0 || isSelfCC ? "p00" : "p01";
                DEBUG("SE(%zu, %zu): pdg0: %d pdg1: %d   %s\n", i0, i1, p0.id(), p1.id(), pair.data());

                hSE[pair]->Fill(kStar);
            }
        }

        // Mixed event
        for (size_t i0 = 0; i0 < part0.size(); i0++) {
            const auto p0 = part0[i0];

            for (size_t iME = 0; iME < partBuffer.size(); iME++) {
                for (size_t i1 = 0; i1 < partBuffer[iME].size(); i1++) {
                    const auto p1 = partBuffer[iME][i1];

                    double kStar = ComputeKstar(p0, p1);
                    std::string pair = p0.id() * p1.id() > 0 || isSelfCC ? "p00" : "p01";

                    DEBUG("ME(%zu, %zu, %zu): pdg0: %d pdg1: %d   %s\n", i0, iME, i1, p0.id(), p1.id(), pair.data());
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
    oFile->mkdir("p00");
    oFile->cd("p00");
    hSE["p00"]->Write("hSE");
    hME["p00"]->Write("hME");
    hPairMultSE["p00"]->Write("hPairMult");

    if (!isSelfCC) {
        oFile->mkdir("p01");
        oFile->cd("p01");
        hSE["p01"]->Write("hSE");
        hME["p01"]->Write("hME");
        hPairMultSE["p01"]->Write("hPairMult");
    }

    oFile->Close();
}
