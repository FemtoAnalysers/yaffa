#include <TFile.h>
#include <TH1.h>

#include "Pythia8/Pythia.h"
#include "TStopwatch.h"

#endif

using namespace Pythia8;

namespace {
enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };
enum processes { kSoftQCD = 0, kNonDiffractive, kHardQCD };
enum triggers { kMB = 0, kHM };
}  // namespace

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

void CreateDataset(int nEvents = 100000, double maxRunTime = 0.5,
                   int energy = 13.6,  // TeV
                   int tune = kCRMode2, int process = kSoftQCD, int trigger = kHM, int seed = 42,
                   std::string oFileName = "auto") {
    float energy = 13000;
    TStopwatch timer;
    timer.Start();

    Pythia pythia;

    if (process == kSoftQCD) {
        pythia.readString("SoftQCD:all = on");
    } else if (process == kHardQCD) {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    }

    // Set tune
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

    // Init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", 0));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", energy);
    pythia.init();

    if (oFileName == "auto") {
        oFileName = Form("PythiaProd_%spp%dTeV_%s_%s.root", trigger, energy, tune, process);
    }
    auto oFile = new TFile(oFileName.data(), "recreate");
    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    TTree* tEvents = new TTree("tEvents", "tEvents");
    tEvents->Branch("particles", "TClonesArray", &particles);

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        if (timer.RealTime() > maxRunTime * 3600) {
            int time = int(timer.RealTime());  // in seconds
            int nHours = int(time / 3600);
            int nMins = int((time % 3600) / 60);

            printf("Reached max run time: %d:%d\n", nHours, nMins);
            break;
        } else {
            timer.Continue();
        }

        pythia.GenerateEvent();
        pythia.ImportParticles(particles, "All");

        if (trigger == kHM) {
            // evaluate multiplicity at forward rapidity
            int nChForward = 0;
            for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
                TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
                int pdg = std::abs(particle->GetPdgCode());
                int status = std::abs(particle->GetStatusCode());
                float eta = particle->Eta();

                if (IsDetectable(pdg) && ((-3.7 < eta && eta < -1.7) || (2.8 < eta && eta < 5.1)) &&
                    status == 1) {  // V0A and V0C acceptance
                    nChForward++;
                }
            }
            if (nChForward > 130) {
                tEvents->Fill();
            }
        } else {
            std::cout << "Trigger not implemented. Exit!\n" << std::endl;
        }
    }
    oFile->Close();
}
