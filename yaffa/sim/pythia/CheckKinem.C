#include <TFile.h>
#include <TH1.h>

#include "Pythia8/Pythia.h"
#include "TStopwatch.h"

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

bool Trigger(Pythia8::Event event, triggers trigger) {
    if (trigger == kMB) {
        for (auto iPart = 3; iPart < event.size(); ++iPart) {
            auto& part = event[iPart];

            int pdg = std::abs(part.id());
            float eta = part.eta();

            // Require at least 1 charged track in the TPC acceptance
            if (IsDetectable(pdg) && ((-0.5 < eta && eta < 0.5))) return true;
        }
    } else if (trigger == kHM) {
        // evaluate multiplicity at forward rapidity
        int nChForward = 0;
        for (auto iPart = 3; iPart < event.size(); ++iPart) {
            auto& part = event[iPart];

            int pdg = std::abs(part.id());
            float eta = part.eta();

            if (IsDetectable(pdg) &&
                ((-3.7 < eta && eta < -1.7) || (2.8 < eta && eta < 5.1))) {  // V0A and V0C acceptance
                nChForward++;
            }
        }
        if (nChForward > 130) return true;
    } else {
        std::cout << "Trigger not implemented. Exit!\n" << std::endl;
    }

    return false;
}

void CheckKinem(int nEvents = 1000,           // Number of events
                double maxRunTime = 0.5,       // Max run time (um = h)
                double energy = 13,            // um = TeV
                tunes tune = kCRMode2,         // Tune
                processes process = kSoftQCD,  // Process
                triggers trigger = kMB,        // Trigger
                int pdg = 211,                 // pdg code of the particle
                int pdgMother = 3224,          // pdg code of the mother
                int seed = 42,                 // seed of the random number generator
                std::string oFileName = "auto") {
    TStopwatch timer;
    timer.Start();

    Pythia8::Pythia pythia;

    if (process == kSoftQCD) {
        pythia.readString("SoftQCD:all = on");
    } else if (process == kHardQCD) {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    } else if (process == kNonDiffractive) {
        pythia.readString("SoftQCD:nonDiffractive = on");
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
    pythia.settings.parm("Beams:eCM", energy * 1000);  // From TeV to GeV
    pythia.init();

    if (oFileName == "auto") {
        oFileName = Form("KinemDau_%spp%.1fTeV_%s_%s.root", triggerNames[trigger], energy, tuneNames[tune],
                         processNames[process]);
    }
    auto oFile = new TFile(oFileName.data(), "recreate");
    TH1D* hPtAll = new TH1D("hPtAll", ";#it{p}_{T} (GeV/c);Counts", 200, 0, 6);
    TH1D* hPt = new TH1D("hPt", ";#it{p}_{T} (GeV/c);Counts", 200, 0, 6);

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

        pythia.next();

        if (!Trigger(pythia.event, trigger)) continue;

        for (auto iPart = 3; iPart < pythia.event.size(); ++iPart) {
            auto& part = pythia.event[iPart];

            float eta = part.eta();
            float pt = part.pT();


            if (std::abs(part.id()) == pdg) {
                hPtAll->Fill(pt);

                if (std::abs(pythia.event[part.mother1()].id()) == pdgMother) {
                    hPtAll->Fill(pt);
                }
            }
        }
    }

    hPtAll->Write();
    hPt->Write();
    oFile->Close();
}
