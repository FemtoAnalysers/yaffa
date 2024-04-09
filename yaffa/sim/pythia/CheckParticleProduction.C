#include <array>
#include <map>
#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "Pythia8/Pythia.h"

using namespace Pythia8;

enum tunes { kMonash, kCRMode0, kCRMode2, kCRMode3 };

void CheckParticleProduction(int nEvents = 20000,
                             std::vector<int> pdgs = {102212, 102214, 102216, 112214, 122212, 202212, 202216, 212212,
                                                      212214},
                             tunes tune = tunes::kCRMode2, int seed = 52) {
    Pythia pythia;

    // Set the processes
    pythia.readString("SoftQCD:nonDiffractive = on");

    // Set the tune
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

    pythia.readString(Form("Random:seed %d", seed));
    pythia.readString("Random:setSeed = on");
    pythia.init();

    // Output
    TH1F* hPt = new TH1F("hPt", ";#it{k*};Counts", 1500, 0, 6);

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;

        for (int iPart = 3; iPart < pythia.event.size(); iPart++) {
            Particle part = pythia.event.at(iPart);
            if (std::find(pdgs.begin(), pdgs.end(), part.id()) != pdgs.end()) {
                hPt->Fill(part.pT());
            }
        }
    }

    TFile oFile(Form("/scratch5/ge86rim/an/LPi/sim/pythia/Nstar/Analysis_Results_%d.root", seed), "recreate");
    hPt->Write();
    oFile.Close();
}
