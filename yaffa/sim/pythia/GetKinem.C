/*
Script to obtain the (pT, y) distribution from pythia events at generator level.
*/

// C++ libraries
#include <map>
#include <string>
#include <vector>

// ROOT libraries
#include "Riostream.h"
#include "TFile.h"
#include "TH2D.h"
#include "TRandom3.h"

// ALICE libraries
#include "Pythia8/Pythia.h"

void GetKinem(
    long unsigned int nEvents = 1000,
    int pdg = 2212,
    std::string tune = "Monash",
    std::string process = "NonDiffractive",
    int seed = 1
    ) {

    // Make output file name
    std::string oFileName = "Kinem";
    if (pdg == 211) oFileName += "_pi";
    else if (pdg == 321) oFileName += "_K";
    else if (pdg == 2212) oFileName += "_p";
    else if (pdg == 3324) oFileName += "_Xi1530zero";
    else if (pdg == 3314) oFileName += "_Xi1530";
    else {
        printf("name for pdg %d not implemented\n", pdg);
        exit(1);
    }
    oFileName += "_" + tune;
    oFileName += "_" + process;
    oFileName += "_" + std::to_string(seed);
    oFileName += ".root";

    // Create the generator
    Pythia8::Pythia pythia;

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

    // Setting the seed here is not sufficient to ensure reproducibility, setting the seed of gRandom is necessary
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", 14000); // from TeV to GeV
    pythia.init();

    gRandom->SetSeed(seed); // Set the seed to ensure reproducibility of the evevents generated by pythia

    TH2D *hYvsPt = new TH2D("hYvsPt", ";#it{p}_{T} (GeV/#it{c});#it{y};Counts", 200, 0, 10, 200, -10, 10);
    TH2D *hEtavsPt = new TH2D("hEtavsPt", ";#it{p}_{T} (GeV/#it{c});#eta;Counts", 200, 0, 10, 200, -10, 10);
    TH1D *hY = new TH1D("hY", ";#it{y};Counts", 200, -10, 10);
    TH1D *hEta = new TH1D("hEta", ";#eta;Counts", 200, -10, 10);
    TH1D *hPt = new TH1D("hPt", ";#it{p}_{T} (GeV/#it{c});Counts", 200, 0, 10);
    
    for (size_t iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia.next();

        // Part 0 is the event, 1 and 2 the beams ==> start from 3
        for (int iPart = 3; iPart < pythia.event.size(); iPart++) {
            auto &part = pythia.event[iPart];

            if (std::abs(part.id()) == pdg) {
                hYvsPt->Fill(part.pT(), part.y());
                hEtavsPt->Fill(part.pT(), part.eta());
                hY->Fill(part.y());
                hEta->Fill(part.eta());
                hPt->Fill(part.pT());
            }
        }
    }
    
    TFile *oFile = new TFile(oFileName.data(), "recreate");
    hYvsPt->Write();
    hEtavsPt->Write();
    hY->Write();
    hEta->Write();
    hPt->Write();
    oFile->Close();

    printf("Output saved in %s\n", oFileName.data());
}
