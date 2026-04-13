#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TSystem.h>

#include <array>
#include <deque>
#include <map>
#include <string>
#include <vector>

// ALICE libraries
#include "Pythia8/Pythia.h"

enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };

enum processes { kSoftQCD = 0, kNonDiffractive, kHardQCD };
std::map<processes, const char*> processNames = {
    {kSoftQCD, "SoftQCD"},
    {kNonDiffractive, "NonDiffractive"},
    {kHardQCD, "HardQCD"},
};
void TuneHMTrigger(int nEvents, tunes tune, processes process, std::string oFileName, int seed);
bool IsDetectable(int absPdg);

bool IsDetectable(int absPdg) {
    if (absPdg == 11 ||   // electrons
        absPdg == 13 ||   // muons
        absPdg == 211 ||  // pions
        absPdg == 321 ||  // kaons
        absPdg == 2212    // protons
    )
        return true;

    return false;
}

void TuneHMTrigger(int nEvents, tunes tune, processes process, std::string oFileName, int seed) {
    int energy = 13000;

    Pythia8::Pythia pythia;

    // Set process
    if (process == kSoftQCD) {
        pythia.readString("SoftQCD:all = on");
    } else if (process == kHardQCD) {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    } else if (process == kNonDiffractive) {
        pythia.readString("SoftQCD:nonDiffractive = on");
    } else {
        std::cout << "\033[31mError: Process not implemented. Exit!\033[0m" << std::endl;
        exit(1);
    }

    // set tune
    if (tune == kMonash) {
        pythia.readString(Form("Tune:pp = 14"));
    } else if (tune == kCRMode0) {
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
    } else if (tune == kCRMode2) {
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
    } else if (tune == kCRMode3) {
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
    }

    // Setting the seed here is not sufficient to ensure reproducibility, setting the seed of gRandom is necessary
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", 13000); // um = GeV
    pythia.init();

    gRandom->SetSeed(seed); // Set the seed to ensure reproducibility of the evevents generated by pythia

    // Define output
    TH2I *hMult = new TH2I("hMult", "hMult;particles in V0A or V0C;Partices in |#eta|<0.5;Entries", 201, -0.5, 200.5,
                           201, -0.5, 200.5);

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia.next();


        // evaluate multiplicity at forward rapidity
        int nChForward = 0;
        int nChMid = 0;
        for (auto iPart = 3; iPart < pythia.event.size(); ++iPart) {
            auto &part = pythia.event[iPart];
    
            if (!part.isFinal() || !IsDetectable(std::abs(part.id()))) continue;
            float eta = part.eta();

            if (std::abs(eta) < 0.5) {
                nChMid++;
            } else if ((-3.7 < eta && eta < -1.7) || (2.8 < eta && eta < 5.1)) {  // V0A and V0C acceptance
                nChForward++;
            }
        }
        hMult->Fill(nChForward, nChMid);
    }

    TFile outFile(oFileName.data(), "recreate");
    hMult->Write();
    outFile.Close();
}
