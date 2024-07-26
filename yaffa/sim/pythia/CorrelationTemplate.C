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


void CorrelationTemplate(int nEvents=10000, int pdgM=3224, int pdgd1=3122, int pdgd2=211, int seed=1, std::string outFileName= "AnalysisResults.root") {
    Pythia8::Pythia pythia;
    
    // set seed for simulation
    pythia.readString(Form("Random:seed %d", seed));
    pythia.readString("Random:setSeed = on");
    
    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:all = on");

    // decay mode of interest, from PYTHIA config file

    // <particle id="3124" name="Lambda(1520)0" antiName="Lambda(1520)bar0" spinType="4" chargeType="0" colType="0" 
    //           m0="1.51950" mWidth="0.01560" mMin="1.40000" mMax="1.65000">
    //  <channel onMode="1" bRatio="0.0667000" products="3122 211 -211"/>
    //  <channel onMode="1" bRatio="0.0020000" products="3212 211 -211"/>
    // </particle>

    // Force the decays
    pythia.readString("3122:onMode = off"); // Lambda
    pythia.readString("3122:onIfMatch = 2212 211");
    
    pythia.readString("3224:onMode = off"); // Sigma(1385)
    pythia.readString("3224:onIfMatch = 3122 211");
    // pythia.readString("211:onMode = off");
    // pythia.readString("3224:onIfMatch = 3122 211 -211");
    // pythia.readString("3124:onIfMatch = 3212 211 -211");

    // init
    pythia.init();

    TFile oFile(outFileName.data(), "recreate");
    TH1F * hSE = new TH1F("hSE", ";#it{k*} (MeV/#it{c});Counts", 1500, 0., 3.);

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia.next();
        for(int iPart=3; iPart<pythia.event.size(); iPart++) {
            if(std::abs(pythia.event[iPart].id()) != pdgM) continue;

            std::set<int> parts1 = SelectInDaughterTree(pythia.event, iPart, {pdgd1});
            std::set<int> parts2 = SelectInDaughterTree(pythia.event, iPart, {pdgd2});

            std::set<int>::iterator it = parts1.begin();

            for (const auto &iP1 : parts1) {
                if (!IsInAcc(pythia.event, iP1)) continue;

                for (const auto &iP2 : parts2) {
                    if (!IsInAcc(pythia.event, iP2)) continue;

                    auto p1 = pythia.event[iP1];
                    auto p2 = pythia.event[iP2];

                    if (p1.daughter1() <= iP2 && iP2 <= p1.daughter2()) continue; // part2 is dau of part1
                    if (p2.daughter1() <= iP1 && iP1 <= p2.daughter2()) continue; // part1 is dau of part2

                    // std::copy_if(oldSet.begin(), oldSet.end(), std::inserter(newSet, newSet.end()), [](const T & value){/*predicate here*/});

                    auto kstar = RelativePairMomentum(pythia.event[iP1], pythia.event[iP2]);

                    hSE->Fill(kstar);
                }
            }
        }
    }

    hSE->Write();
    oFile.Close();
}
