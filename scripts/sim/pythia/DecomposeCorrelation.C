#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include "Pythia8/Pythia.h"

#include "utils.hxx"


void DecomposeCorrelation(int nEvents=10000, int pdg1=3122, int pdg2=211, int seed=1, std::string outFileName= "AnalysisResults.root") {
    Pythia8::Pythia pythia;
    
    // set seed for simulation
    pythia.readString(Form("Random:seed %d", seed));
    pythia.readString("Random:setSeed = on");
    
    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:nonDiffractive = on");

    // Force the decays
    pythia.readString("3122:onMode = off"); // Lambda
    pythia.readString("3122:onIfMatch = 2212 211");
    
    pythia.readString("3224:onMode = off"); // Sigma(1385)
    pythia.readString("3224:onIfMatch = 3122 211");

    // init
    pythia.init();
    printf("Ehi\n");

    TFile oFile(outFileName.data(), "recreate");
    std::map<std::string, TH1F *> hSE = {
        {"tot", new TH1F("hSETot", ";#it{k*} (MeV/#it{c});Counts", 1500, 0., 3.)},
        {"prim1", new TH1F("hSEPrim1", ";#it{k*} (MeV/#it{c});Counts", 1500, 0., 3.)},
        {"prim2", new TH1F("hSEPrim1", ";#it{k*} (MeV/#it{c});Counts", 1500, 0., 3.)},
        {"rest", new TH1F("hSERest", ";#it{k*} (MeV/#it{c});Counts", 1500, 0., 3.)},
    };

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia.next();

        std::vector<int> parts1 = {};
        std::vector<int> parts2 = {};

        // Start from iPart = 3 since 0 is the event and 1 and 2 are the colliding particles
        for(int iPart=3; iPart<pythia.event.size(); iPart++) {
            int abspdg = std::abs(pythia.event[iPart].id());
            
            if (abspdg == pdg1 || abspdg == pdg2) {
                if (!IsInAcc(pythia.event, iPart)) continue;
                
                auto part = pythia.event[iPart];
                if (part.xProd() * part.xProd() + part.yProd() * part.yProd() + part.zProd() * part.zProd() > 0.1) continue;

                if (abspdg == pdg1) {
                    parts1.push_back(iPart);
                } else {
                    parts2.push_back(iPart);
                }
            }
        }

        // Same event
        for (const auto &iP1 : parts1) {
            auto p1 = pythia.event[iP1];
            for (const auto &iP2 : parts2) {
                auto p2 = pythia.event[iP2];

                if (p1.daughter1() <= iP2 && iP2 <= p1.daughter2()) continue; // part2 is direct dau of part1
                if (p2.daughter1() <= iP1 && iP1 <= p2.daughter2()) continue; // part1 is direct dau of part2

                auto kstar = RelativePairMomentum(p1, p2);

                hSE["tot"]->Fill(kstar);

                if (IsPrimary(pythia.event, iP1)) hSE["prim1"]->Fill(kstar);
                else hSE["rest"]->Fill(kstar);
            }
        }
    }

    for (const auto & [_, h] : hSE) {
        h->Write();
    }

    oFile.Close();
}
