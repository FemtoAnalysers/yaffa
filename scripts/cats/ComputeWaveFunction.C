// The environment variable 'CATS' must be defined and must point to the DLM
// Cats reposotory. It is recommended to define it in your .env file
R__LOAD_LIBRARY($CATS/install/lib/libCATS.so)

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TMinuit.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CATS.h"
#include "CATSconstants.h"
#include "CommonAnaFunctions.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomp.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_Fitters.h"
#include "DLM_Histo.h"
#include "DLM_Potentials.h"
#include "DLM_ResponseMatrix.h"
#include "DLM_RootFit.h"
#include "DLM_RootWrapper.h"
#include "DLM_Source.h"
#include "TRandom.h"
#include "TString.h"
#include "TSystem.h"

const double RADIUS_STEP = 0.01; // um = fm
const double RADIUS_MAX = 6.25;
const double KSTAR_STEP = 0.4; // um = MeV
const double KSTAR_MAX = 400; // um = MeV

std::vector<std::vector<double>> Sum(std::vector<std::vector<double>> matrix1, std::vector<std::vector<double>> matrix2) {
    std::vector<std::vector<double>> sum;
    for (int i = 0; i  < matrix1.size(); i++){
        sum.push_back({});

        for (int j = 0; j  < matrix1.size(); j++){
            sum[sum.size() - 1].push_back(matrix1[i][j] + matrix2[i][j]);
        }
    }
    
    return sum;
}

TH2D* TableToTH2D(std::vector<std::vector<double>> matrix) {
    TH2D *hWF = new TH2D("hWF", "", std::round(RADIUS_MAX / RADIUS_STEP), 0, RADIUS_MAX, std::round(KSTAR_MAX / KSTAR_STEP), 0, KSTAR_MAX);
    for (int i = 0; i  < matrix.size(); i++){
        for (int j = 0; j  < matrix[i].size(); j++){
            hWF->SetBinContent(j+1, i+1, matrix[i][j]);
        }
    }
    return hWF;
}
void print(std::vector<std::vector<double>> matrix){

    printf("-->> %d %d\n", matrix.size(), matrix[0].size());
    for (int i = 0; i  < matrix.size(); i++){
        for (int j = 0; j  < matrix[i].size(); j++){
            printf("%.3f   ", matrix[i][j]);
        }
        printf("\n");
    }
}

void ToFile(std::string filename, std::vector<std::vector<double>> matrix) {
    std::ofstream out(filename);
    
    out << std::fixed << std::setprecision(6);

    out << "radius";
    for (int iKStar = 0; iKStar < matrix.size(); iKStar++) {
        out << "\t" << KSTAR_STEP / 2 + KSTAR_STEP * iKStar;
    }
    out << "\n";

    for (int iRadius = 0; iRadius < matrix[0].size(); iRadius++) {
        out << RADIUS_STEP / 2 + RADIUS_STEP * iRadius;
        for (int iKStar = 0; iKStar < matrix.size(); iKStar++) {
            out << " " << matrix[iKStar][iRadius] << " ";
        }
        out << "\n";
    }
}

std::vector<std::vector<double>> GetWaveFunction(CATS *cats, double weight = 1, int channel = -1) {
    std::vector<std::vector<double>> wf;

    if (channel < 0) {
        int nChannels = cats->GetNumChannels();
        double weightList[4];

        if (nChannels == 2) {  // without p-wave
            weightList[0] = 0.25;
            weightList[1] = 0.75;
        } else if (nChannels == 4) {  // p-wave
            weightList[0] = 3. / 12.;
            weightList[1] = 1. / 12.;
            weightList[2] = 3. / 12.;
            weightList[3] = 5. / 12.;
        }

        wf = GetWaveFunction(cats, weightList[0], 0);

        for (int iChn = 1; iChn < cats->GetNumChannels(); iChn++) {
            wf = Sum(wf, GetWaveFunction(cats, weightList[iChn], iChn));
        }
        return wf;
    }

    for (unsigned short iKStar = 0; iKStar < cats->GetNumMomBins(); iKStar++) {
        double kstar = cats->GetMomentum(iKStar);
        wf.push_back({});
        for (unsigned short iRadius = 0; RADIUS_STEP / 2 + RADIUS_STEP * iRadius < RADIUS_MAX; iRadius++) {
            double radius = RADIUS_STEP / 2 + RADIUS_STEP * iRadius;
            double WFvalue = cats->EvalWaveFun2(iKStar, radius, channel);
            wf[wf.size() - 1].push_back(WFvalue * weight);
        }
    }
    
    return wf;
}

void ComputeWaveFunction(int pdg1 = 2212, int pdg2 = 2212, double r0 = 1.25) {
    auto pdg = TDatabasePDG::Instance();

    // const double m1 = pdg->GetParticle(pdg1)->Mass() * 1000; // conversion GeV -> MeV
    // const double m2 = pdg->GetParticle(pdg2)->Mass() * 1000; // conversion GeV -> MeV
    const int nKStarBins = std::round(KSTAR_MAX / KSTAR_STEP);
    // const double KSTAR_MAX = 500.;
    // const int nRadBins = 2000;
    // double rmax = 20;

    // Roberta
    const double m1 = 938.2720813;  // conversion GeV -> MeV
    const double m2 = 938.2720813;  // conversion GeV -> MeV
    const double SourceRadValue = 1.25;
    double rmax = SourceRadValue * 5.;

    const double mu = m1 * m2 / (m1 + m2);

    CATS *cats = new CATS();
    cats->SetMomBins(nKStarBins, 0, KSTAR_MAX);
    cats->SetQ1Q2(1);                  // 1 same charge; 0 neutral charge; -1
    cats->SetQuantumStatistics(true);  // set it to true if you want quantum statistics
    cats->SetRedMass(mu);
    
    // Set the potential
    DLM_CommonAnaFunctions can;
    if (pdg1 == 2212 && pdg2 == 2212) {
        can.SetUpCats_pp(*cats, "AV18", "Gauss", 0, 0);
    } else {
        std::cout << "System not implemented. Exit!" << std::endl;
        exit(1);
    }

    cats->SetAnaSource(0, r0);
    cats->KillTheCat();

    TGraph* gCF = new TGraph();
    for (int i = 0; i < nKStarBins; i++) {
        double CorrFuncValue = cats->GetCorrFun(i);
        double kstar = cats->GetMomentum(i);
        gCF->SetPoint(i, kstar, CorrFuncValue);
    }

    // obtain the total wave function
    const int nChannels = cats->GetNumChannels();
    double weighList[nChannels];
    if (nChannels == 2) {  // without p-wave
        weighList[0] = 0.25;
        weighList[1] = 0.75;
    } else if (nChannels == 4) {  // p-wave
        weighList[0] = 3. / 12.;
        weighList[1] = 1. / 12.;
        weighList[2] = 3. / 12.;
        weighList[3] = 5. / 12.;
    } else {
        std::cerr << "Error: unexpected number of channels: " << nChannels << std::endl;
        return;
    }

    auto wf = GetWaveFunction(cats);
    ToFile("wf.dat", wf);
    // for (unsigned short iChn = 0; iChn < nChannels; iChn++) {
    //     break;
    // }
    TH2D *hWF = TableToTH2D(wf);


    TFile* OutputFile = new TFile("ppCF_WF_Source.root", "recreate");
    OutputFile->cd();
    gCF->Write("gCF");
    hWF->SetTitle("pp, AV18, |#psi|^{2};r (fm);k* (MeV/c);|#psi|^{2}");
    hWF->Write("hWF");
    OutputFile->Close();
}
