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
using namespace Pythia8;

float RelativePairMomentum(TLorentzVector &PartOne, TLorentzVector &PartTwo) {
  TLorentzVector trackSum = PartOne + PartTwo;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartOneCMS = PartOne;
  TLorentzVector PartTwoCMS = PartTwo;

  PartOneCMS.Boost(-betax, -betay, -betaz);
  PartTwoCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;

  return 0.5 * trackRelK.P();
}

double BreitWigner(double *x, double *par) {

    double yield = par[0];
    double mean = par[1];
    double gamma = par[2];

    return yield * gamma / TMath::Pi() / (gamma * gamma + (x[0] - mean) * (x[0] - mean));
}

void Lambda1520ThreeBody(int nEvents=10000, int seed=42, 
                std::string kinemFilePath="/home/ktas/go98pog/hello_anton/SimLambda1520Kinem_MERGED.root", // kinematic information
                std::string outFilePath= "./Lambda1520ThreeBodyPtScan") {
    
    Pythia8::Pythia pythia;
    
    // set seed for simulation
    pythia.readString(Form("Random:seed %d", seed));
    pythia.readString("Random:setSeed = on");
    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:all = on");

    // switch on Lambda(1520) in p pi- pi0
    pythia.readString("102134:onMode = off");
    pythia.readString("211:onMode = off");
    pythia.readString("111:onMode = off");
    pythia.readString("102134:oneChannel = 1 1 1 3122 211 -211");
    // pythia.readString("102134:oneChannel = 1 1 0 3122 211 -211");
    int pdgLambda1520 = 102134;

    // init
    pythia.init();
  
    // Read histograms from Xi(1530) to simulate acceptance and pT of the particle
    TFile *kinemXi1530 = TFile::Open(kinemFilePath.data(), "r");
    TH1F * hAccXi1530 = (TH1F*)kinemXi1530->Get("hAcc");
    TH1F * hPtXi1530 = (TH1F*)kinemXi1530->Get("hPt");
    TH1F * hSampledBreitWigner = new TH1F("hSampledBreitWigner", ";#it{k*};Counts", 12500, 0, 3.);
    TH1F * hSampledPt = new TH1F("hSampledPt", ";#it{p}_{T} (GeV/c);Counts", 1000, 0, 10);
    TH1F * hSampledY = new TH1F("hSampledY", ";y;Counts", 1000, -10, 10);
    TH1F * hSampledEta = new TH1F("hSampledEta", ";#eta;Counts", 1000, -10, 10);
    TH1F * hSampledMass = new TH1F("hSampledMass", ";M;Counts", 1000, 1.4, 1.7);
    TH1F * hSELambda1520Pt015 = new TH1F("hSELambda1520Pt015", ";#it{k*};Counts", 1500, 0., 6.);
    TH1F * hSELambda1520Pt017 = new TH1F("hSELambda1520Pt017", ";#it{k*};Counts", 1500, 0., 6.);
    TH1F * hSELambda1520Pt019 = new TH1F("hSELambda1520Pt019", ";#it{k*};Counts", 1500, 0., 6.);
    TH1F * hSELambda1520NoAccPt015 = new TH1F("hSELambda1520NoAccPt015", ";#it{k*};Counts", 1500, 0., 6.);
    TH1F * hSELambda1520NoAccPt017 = new TH1F("hSELambda1520NoAccPt017", ";#it{k*};Counts", 1500, 0., 6.);
    TH1F * hSELambda1520NoAccPt019 = new TH1F("hSELambda1520NoAccPt019", ";#it{k*};Counts", 1500, 0., 6.);

    double LambdaPiPiThreshold = 2*TDatabasePDG::Instance()->GetParticle(211)->Mass() + 
                                 TDatabasePDG::Instance()->GetParticle(3122)->Mass();

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

        if(iEvent%5000 == 0)
          std::cout << "Processing event " << iEvent << endl;

        double massLambda1520 = 0.;
        while(massLambda1520 < LambdaPiPiThreshold) {
            massLambda1520 = gRandom->BreitWigner(1.519, 0.016);
            massLambda1520 = gRandom->BreitWigner(1.519, 0);
        }
        hSampledBreitWigner->Fill(massLambda1520);

        // initialize particle properties of Lambda(1520)
        double phiLambda1520 = gRandom->Rndm() * 2 * TMath::Pi();
        double etaLambda1520 = hAccXi1530->GetRandom();
        double thetaLambda1520 = 2 * TMath::ATan( TMath::Exp(-etaLambda1520) );
        double ptLambda1520 = hPtXi1530->GetRandom();
        double pxLambda1520 = ptLambda1520 * TMath::Cos(phiLambda1520);
        double pyLambda1520 = ptLambda1520 * TMath::Sin(phiLambda1520);
        double pzLambda1520 = ptLambda1520 * TMath::SinH(etaLambda1520);
        double pLambda1520 = TMath::Sqrt(ptLambda1520 * ptLambda1520 + pzLambda1520 * pzLambda1520);
        double ELambda1520 = TMath::Sqrt(massLambda1520 * massLambda1520 + pLambda1520 * pLambda1520);
        hSampledPt->Fill(ptLambda1520);
        
        Particle lambda1520;
        lambda1520.id(pdgLambda1520);
        lambda1520.status(81);
        lambda1520.m(massLambda1520);
        lambda1520.xProd(0.);
        lambda1520.yProd(0.);
        lambda1520.zProd(0.);
        lambda1520.tProd(0.);
        lambda1520.e(ELambda1520);
        lambda1520.px(pxLambda1520);
        lambda1520.py(pyLambda1520);
        lambda1520.pz(pzLambda1520);

        hSampledY->Fill(lambda1520.y());
        hSampledEta->Fill(lambda1520.eta());
        hSampledMass->Fill(lambda1520.m());
        
        // remove all particles generated in the event and append the Lambda(1520)
        pythia.event.reset();
        pythia.event.append(lambda1520);
        int idPartLambda = pythia.event[1].id();
        pythia.particleData.mayDecay(idPartLambda, true);
        
        // force the decay of the Lambda
        pythia.moreDecays();

        Particle lambda;
        Particle piPlus;
        Particle piMinus;

        for(int iPart=1; iPart<pythia.event.size(); iPart++) {
            if(pythia.event[iPart].id()==3122) {
                lambda = pythia.event[iPart];
            }
            if(pythia.event[iPart].id()==211 && pythia.event.at(pythia.event[iPart].mother1()).id() != 3122) {
                piPlus = pythia.event[iPart];
            }
            if(pythia.event[iPart].id()==-211 && pythia.event.at(pythia.event[iPart].mother1()).id() != 3122) {
                piMinus = pythia.event[iPart];
            }
        }
        TLorentzVector momLambda(lambda.px(), lambda.py(), lambda.pz(), lambda.e());
        TLorentzVector momPiMinus(piMinus.px(), piMinus.py(), piMinus.pz(), piMinus.e());
        TLorentzVector momPiPlus(piPlus.px(), piPlus.py(), piPlus.pz(), piPlus.e());
 
        float kStarLPiPlus = RelativePairMomentum(momLambda, momPiPlus);

        if(piPlus.pT() > 0.15 && piPlus.pT() <= 2. && lambda.pT()>0.3) {
            hSELambda1520NoAccPt015->Fill(kStarLPiPlus);
            if(abs(piPlus.eta()) <= 0.8) {
                hSELambda1520Pt015->Fill(kStarLPiPlus);
            }
        }

        if(piPlus.pT() > 0.17 && piPlus.pT() <= 2. && lambda.pT()>0.3) {
            hSELambda1520NoAccPt017->Fill(kStarLPiPlus);
            if(abs(piPlus.eta()) <= 0.8) {
                hSELambda1520Pt017->Fill(kStarLPiPlus);
            }
        }

        if(piPlus.pT() > 0.19 && piPlus.pT() <= 2. && lambda.pT()>0.3) {
            hSELambda1520NoAccPt019->Fill(kStarLPiPlus);
            if(abs(piPlus.eta()) <= 0.8) {
                hSELambda1520Pt019->Fill(kStarLPiPlus);
            }
        }

    }

    // Write histograms to file
    std::string outFileName = outFilePath + "_" + std::to_string(seed) + ".root";  
    TFile oFile(outFileName.data(), "recreate");
    oFile.cd();
    hAccXi1530->Write();
    hPtXi1530->Write();
    hSELambda1520Pt015->Write("hSE_015");
    hSELambda1520Pt017->Write("hSE_017");
    hSELambda1520Pt019->Write("hSE_019");
    hSELambda1520NoAccPt015->Write("hSE_NoAcceptanceCuts_015");
    hSELambda1520NoAccPt017->Write("hSE_NoAcceptanceCuts_015");
    hSELambda1520NoAccPt019->Write("hSE_NoAcceptanceCuts_015");
    hSampledBreitWigner->Write(); 
    hSampledPt->Write(); 
    hSampledY->Write(); 
    hSampledEta->Write(); 
    hSampledMass->Write(); 
    oFile.Close();
 
}