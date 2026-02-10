#include "TFile.h"
#include "TH1.h"
#include "TPythia8.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TParticle.h"

using namespace Pythia8;

enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };
enum processes { kSoftQCD = 0, kNonDiffractive, kHardQCD };
enum triggers { kMB = 0, kHM };

std::map<tunes, const char*> tuneNames = {
    {kMonash, "Monash"},
    {kCRMode0, "CRMode0"},
    {kCRMode2, "CRMode2"},
    {kCRMode3, "CRMode3"},
};

std::map<processes, const char*> processNames = {
    {kSoftQCD, "SoftQCD"},
    {kNonDiffractive, "NonDiffractive"},
    {kHardQCD, "HardQCD"},
};

std::map<triggers, const char*> triggerNames = {
    {kMB, "MB"},
    {kHM, "HM"},
};
    
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
                   double energy = 13.6,  // TeV
                   tunes tune = kCRMode2, processes process = kSoftQCD, triggers trigger = kHM, int seed = 42,
                   std::string oFileName = "auto") {
    TStopwatch timer;
    timer.Start();

    TPythia8 *pythia = TPythia8::Instance();
    
    // Set process
    if (process == kSoftQCD) {
        pythia->ReadString("SoftQCD:all = on");
    } else if (process == kHardQCD) {
        pythia->ReadString("HardQCD:hardccbar = on");
        pythia->ReadString("HardQCD:hardbbbar = on");
    } else if (process == kNonDiffractive) {
        pythia->ReadString("SoftQCD:nonDiffractive = on");
    }

    // Set tune
    if (tune == kMonash) {
        pythia->ReadString(Form("Tune:pp = 14"));
    } else if (tune == kCRMode0) {
        pythia->ReadString(Form("Tune:pp = 14"));
        pythia->ReadString("ColourReconnection:mode = 1");
        pythia->ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia->ReadString("ColourReconnection:m0 = 2.9");
        pythia->ReadString("ColourReconnection:allowJunctions = on");
        pythia->ReadString("ColourReconnection:junctionCorrection = 1.43");
        pythia->ReadString("ColourReconnection:timeDilationMode = 0");
        pythia->ReadString("StringPT:sigma = 0.335");
        pythia->ReadString("StringZ:aLund = 0.36");
        pythia->ReadString("StringZ:bLund = 0.56");
        pythia->ReadString("StringFlav:probQQtoQ = 0.078");
        pythia->ReadString("StringFlav:ProbStoUD = 0.2");
        pythia->ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia->ReadString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia->ReadString("BeamRemnants:remnantMode = 1");
        pythia->ReadString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode2) {
        pythia->ReadString(Form("Tune:pp = 14"));
        pythia->ReadString("ColourReconnection:mode = 1");
        pythia->ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia->ReadString("ColourReconnection:m0 = 0.3");
        pythia->ReadString("ColourReconnection:allowJunctions = on");
        pythia->ReadString("ColourReconnection:junctionCorrection = 1.20");
        pythia->ReadString("ColourReconnection:timeDilationMode = 2");
        pythia->ReadString("ColourReconnection:timeDilationPar = 0.18");
        pythia->ReadString("StringPT:sigma = 0.335");
        pythia->ReadString("StringZ:aLund = 0.36");
        pythia->ReadString("StringZ:bLund = 0.56");
        pythia->ReadString("StringFlav:probQQtoQ = 0.078");
        pythia->ReadString("StringFlav:ProbStoUD = 0.2");
        pythia->ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia->ReadString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia->ReadString("BeamRemnants:remnantMode = 1");
        pythia->ReadString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode3) {
        pythia->ReadString(Form("Tune:pp = 14"));
        pythia->ReadString("ColourReconnection:mode = 1");
        pythia->ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia->ReadString("ColourReconnection:m0 = 0.3");
        pythia->ReadString("ColourReconnection:allowJunctions = on");
        pythia->ReadString("ColourReconnection:junctionCorrection = 1.15");
        pythia->ReadString("ColourReconnection:timeDilationMode = 3");
        pythia->ReadString("ColourReconnection:timeDilationPar = 0.073");
        pythia->ReadString("StringPT:sigma = 0.335");
        pythia->ReadString("StringZ:aLund = 0.36");
        pythia->ReadString("StringZ:bLund = 0.56");
        pythia->ReadString("StringFlav:probQQtoQ = 0.078");
        pythia->ReadString("StringFlav:ProbStoUD = 0.2");
        pythia->ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia->ReadString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia->ReadString("BeamRemnants:remnantMode = 1");
        pythia->ReadString("BeamRemnants:saturation = 5");
    }

    // Init
    pythia->ReadString(Form("Random:seed %d", seed));
    pythia->Initialize(2212 /* p */, 2212 /* p */,  energy * 1000. /* TeV */);

    TFile *oFile = new TFile(oFileName.data(), "recreate");
    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    TTree *tEvents = new TTree("tEvents","tEvents");
    tEvents->Branch("events", &particles);

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        if (timer.RealTime() > maxRunTime * 3600) {
            int time = int(timer.RealTime());  // Time in seconds
            int nHours = int(time / 3600);
            int nMins = int((time % 3600) / 60);

            printf("Reached max run time: %d:%d\n", nHours, nMins);
            break;
        } else {
            timer.Continue();
        }

        pythia->GenerateEvent();
        pythia->ImportParticles(particles,"All");

        if (trigger == kMB) {
            for (auto iPart = 3; iPart < particles->GetEntriesFast(); ++iPart) {
                auto part = (TParticle*) particles->At(iPart);

                int pdg = std::abs(part->GetPdgCode());
                float eta = part->Eta();

                if (IsDetectable(pdg) && ((-0.5 < eta && eta < 0.5))) {  // V0A and V0C acceptance
                    tEvents->Fill();
                    break;
                }
            }
        }
        else if (trigger == kHM) {
            // evaluate multiplicity at forward rapidity
            int nChForward = 0;
            for (auto iPart = 3; iPart < particles->GetEntriesFast(); ++iPart) {
                auto part = (TParticle*) particles->At(iPart);

                int pdg = std::abs(part->GetPdgCode());
                float eta = part->Eta();

                if (IsDetectable(pdg) && ((-3.7 < eta && eta < -1.7) || (2.8 < eta && eta < 5.1))) {  // V0A and V0C acceptance
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
    // oFile->cd();
    tEvents->Write();
    oFile->Close();
}
