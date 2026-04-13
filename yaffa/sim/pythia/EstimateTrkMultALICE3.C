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


// todo: not implemented
enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };
std::map<tunes, const char*> tuneNames = {
    {kMonash, "Monash"},
    {kCRMode0, "CRMode0"},
    {kCRMode2, "CRMode2"},
    {kCRMode3, "CRMode3"},
};

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
int ComputeMult(const Pythia8::Pythia &pythia, double maxEtaTrk) {
    int mult = 0;
    for (int iPart = 3; iPart < pythia.event.size(); iPart++) {
        Pythia8::Particle part = pythia.event[iPart];

        if (part.isFinal() && std::abs(part.eta()) < maxEtaTrk && IsDetectable(part.id())) mult++;
    }
    return mult;
}

bool IsInAcc(Pythia8::Particle part, const Pythia8::Pythia &pythia, double maxEtaTrk) {
    if (std::abs(part.id()) == 421 || std::abs(part.id()) == 411) { // D0 and D+
        for (const int &dau : part.daughterList()) if (std::abs(pythia.event[dau].eta()) > maxEtaTrk) return false;
        return true;
    } else if (std::abs(part.id()) == 413) { // D*+(2010)
        Pythia8::Particle d1, d2, pion, Dzero;
        d1 = pythia.event[part.daughter1()];
        d2 = pythia.event[part.daughter2()];
        if (std::abs(d1.id()) == 421) {
            Dzero=d1;
            pion=d2;
        } else {
            Dzero=d2;
            pion=d1;
        }

        return std::abs(pion.eta()) < maxEtaTrk && IsInAcc(Dzero, pythia, maxEtaTrk);
    } else {
        printf("Decay of particle %d is not implemented. Exit!\n", part.id());
        exit(1);
    }
    return true;
}



void EstimateTrkMultALICE3(
    long unsigned int nEvents = 10000,
    std::string oFileName = "../../secrets/alice3/TrkMultALICE3.root",
    int seed = 31
    ) {
    Pythia8::Pythia pythia;

    // Set process
    pythia.readString("SoftQCD:nonDiffractive = on");

    // Set tune (CR2)
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

    // Setting the seed here is not sufficient to ensure reproducibility, setting the seed of gRandom is necessary
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", 14000); // from TeV to GeV
    pythia.init();

    gRandom->SetSeed(seed); // Set the seed to ensure reproducibility of the evevents generated by pythia

    std::vector<double> etas = {4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.9};
    std::map<int, TH1D*> hTrkMult = {};
    for (const auto &eta : etas) {
        hTrkMult.insert({std::round(eta*10), new TH1D(Form("hTrkMult_EtaLT%d", (int)std::round(eta*10)), Form(";#it{N}_{ch}_{|#eta_{trk}|<%.1f};Counts", eta), 501, -0.5, 500.5)});
    }

    for (size_t iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia.next();

        for (const auto &eta :etas) {
            hTrkMult[(int)std::round(eta*10)]->Fill(ComputeMult(pythia, eta));
        }
    }

    TFile *oFile = new TFile(oFileName.data(), "recreate");
    for (const auto &eta : etas) {
        hTrkMult[(int)std::round(eta*10)]->Write();
    }
    oFile->Close();
}
