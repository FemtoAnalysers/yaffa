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

void RapidityVsMass(std::string process = "NonDiffractive", std::string tune = "Monash", std::string oFileName = "Distr.root", int seed = 1) {
    Pythia8::Pythia pythia;
    pythia.readString("Next:numberShowEvent = 0");
    // Set process
    if (process == "SoftQCD") {
        pythia.readString("SoftQCD:all = on");
    } else if (process == "HardQCD") {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    } else if (process == "NonDiffractive") {
        pythia.readString("SoftQCD:nonDiffractive = on");
    } else {
        std::cout << "\033[31mError: Process not implemented. Exit!\033[0m" << std::endl;
        exit(1);
    }

    // Set tune
    if (tune == "Monash") {
        pythia.readString(Form("Tune:pp = 14"));
    } else if (tune == "CRMode0") {
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
    } else if (tune == "CRMode2") {
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
    } else if (tune == "CRMode3") {
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
    pythia.settings.parm("Beams:eCM", 13000); // from TeV to GeV
    pythia.init();

    gRandom->SetSeed(seed); // Set the seed to ensure reproducibility of the evevents generated by pythia

    // QA histograms
    TH2D * hY = new TH2D("hY", ";;#sigma_{y}", 6, 0.5, 6.5, 200, -20, 20);
    hY->GetXaxis()->SetBinLabel(1, "#pi");
    hY->GetXaxis()->SetBinLabel(2, "K");
    hY->GetXaxis()->SetBinLabel(3, "p");
    hY->GetXaxis()->SetBinLabel(4, "#Lambda");
    hY->GetXaxis()->SetBinLabel(5, "D^{+}");
    hY->GetXaxis()->SetBinLabel(6, "#Sigma(1385)^{+}");

    for (size_t iEvent = 0; iEvent < 50000; iEvent++) {
        pythia.next();

        // Part 0 is the event, 1 and 2 the beams. In case the hadrons are injected there are no beam particles
        for (int iPart = 1; iPart < pythia.event.size(); iPart++) {
            auto &part = pythia.event[iPart];

            int status = std::abs(part.status());
            if (81 <= status && status <= 89) {
                if (int pdg = std::abs(part.id()); pdg == 211) hY->Fill(1, part.y());
                else if (pdg == 321) hY->Fill(2, part.y());
                else if (pdg == 2212) hY->Fill(3, part.y());
                else if (pdg == 3122) hY->Fill(4, part.y());
                else if (pdg == 411) hY->Fill(5, part.y());
                else if (pdg == 3224) hY->Fill(6, part.y());
            }
        }
    }

    auto PDG = TDatabasePDG::Instance();

    TGraphErrors *gSigmaVsMass = new TGraphErrors(1);
    for (int iBin = 0; iBin < hY->GetNbinsX(); iBin++) {
        auto hProj = hY->ProjectionY("", iBin+1, iBin+1);

        auto c = new TCanvas();
        hProj->Fit("gaus", "M", "", -6.5, 6.5);
        hProj->Draw();
        c->SaveAs(Form("cY_%d.pdf", iBin+1));
        
        auto fGauss = hProj->GetFunction("gaus");
        
        double mass;
        if (iBin == 0) mass = PDG->GetParticle(211)->Mass();
        else if (iBin == 1) mass = PDG->GetParticle(321)->Mass();
        else if (iBin == 2) mass = PDG->GetParticle(2212)->Mass();
        else if (iBin == 3) mass = PDG->GetParticle(3122)->Mass();
        else if (iBin == 4) mass = PDG->GetParticle(411)->Mass();
        else if (iBin == 5) mass = PDG->GetParticle(3224)->Mass();

        gSigmaVsMass->SetPoint(iBin, mass, fGauss->GetParameter(2));
        gSigmaVsMass->SetPointError(iBin, 0, fGauss->GetParError(2));
    }
    TFile *oFile = new TFile(oFileName.data(), "recreate");
    hY->Write();
    gSigmaVsMass->Write();

    oFile->Close();
    std::cout << "Output saved in " << oFileName << std::endl;
}
