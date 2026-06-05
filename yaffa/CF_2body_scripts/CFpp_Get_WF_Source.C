R__LOAD_LIBRARY($HOME/DLM/install/lib/libCATS.so) // good for sunrise setting. To change in case of different installation path


#include "/home/feriorob/DLM/install/include/CATS.h"
#include "/home/feriorob/DLM/install/include/CATSconstants.h"
#include "/home/feriorob/DLM/install/include/DLM_Histo.h"
#include "/home/feriorob/DLM/install/include/DLM_Ck.h"
#include "/home/feriorob/DLM/install/include/DLM_CkDecomposition.h"
#include "/home/feriorob/DLM/install/include/DLM_CkDecomp.h"
#include "/home/feriorob/DLM/install/include/DLM_CkModels.h"
#include "/home/feriorob/DLM/install/include/DLM_Fitters.h"
#include "/home/feriorob/DLM/install/include/DLM_Potentials.h"
#include "/home/feriorob/DLM/install/include/DLM_Source.h"
#include "/home/feriorob/DLM/install/include/DLM_RootFit.h"
#include "/home/feriorob/DLM/install/include/DLM_ResponseMatrix.h"
#include "/home/feriorob/DLM/install/include/DLM_RootWrapper.h"
#include "/home/feriorob/DLM/install/include/CommonAnaFunctions.h"

#include "/home/feriorob/PhD_analysis/analysis/constants.h" // contains the constants used in the analysis, such as the mass of the proton, lambdapar etc.
//#include "/home/feriorob/PhD_analysis/analysis/utils.h" // contains some utility functions, such as SetGraphStyle, GetColor etc.
#include "/home/feriorob/PhD_analysis/analysis/utils.cpp"
#include "/home/feriorob/PhD_analysis/analysis/utils_physics.cpp"

#include <iostream>
#include <complex>
#include "TString.h"
#include "TSystem.h"
#include <fstream>
#include <string>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TChain.h>
#include <TCutG.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <sstream>
#include <TMinuit.h>
#include <TCanvas.h>
#include "TRandom.h"
#include <TGraphErrors.h>
#include <vector>

void CFpp_Get_WF_Source()
{
    const double massProton = ParticleMass::MassProton;
    const double mu12=massProton/2.; // m1*m2/(m1+m2)
    const double SourceRadValue=1.249;
    const int nBins = 200;
    const double kmin = 0.;
    const double kmax = 400.;
    const int nRadBins = 400;
    double rmax = SourceRadValue*5.;

    CATS Kitty;
    Kitty.SetMomBins(nBins, kmin, kmax);
    CATSparameters* cParSource = new CATSparameters(CATSparameters::tSource, 1, true);
    cParSource->SetParameter(0, SourceRadValue); 
    Kitty.SetAnaSource(GaussSource, *cParSource); // Set a Gaussian source with parameters saved in cParSource
    Kitty.SetQ1Q2(1); // 1 same charge; 0 neutral charge; -1
    Kitty.SetQuantumStatistics(true); // set it to true if you want quantum statistics
    Kitty.SetRedMass(mu12);
    Use_P_D_waves(Kitty, true, true); // true true means p- and d-waves included
    Kitty.KillTheCat(); 

    TGraph *gCFpp = new TGraph(); //more simple possible ppCF
    for (int i=0; i<nBins; i++){
        double CorrFuncValue=Kitty.GetCorrFun(i);
        double kstar = Kitty.GetMomentum(i);
        gCFpp->SetPoint(i, kstar, CorrFuncValue);
    }

    // obtain the total wave function
    unsigned short NumChannels=Kitty.GetNumChannels();
    double weighList[NumChannels];
    if(NumChannels==2){
        weighList[0] = 0.25;
        weighList[1] = 0.75; // weights for the two channels, in case of no p-waves
    }
    if(NumChannels==4){
        weighList[0] = 3./12.;
        weighList[1] = 1./12.;
        weighList[2] = 3./12.;
        weighList[3] = 5./12.; // weights for the four channels, in case of p-waves
    }
    if(NumChannels!=2 && NumChannels!=4){
        std::cerr << "Error: unexpected number of channels: " << NumChannels << std::endl;
        return;
    }

    TH2D *hWFpp = new TH2D("hWFpp","hWFpp",nRadBins,0.,rmax,nBins,kmin,kmax);
    for (unsigned short iMom=0; iMom<nBins; iMom++){
        for (unsigned short iRad=0; iRad<nRadBins; iRad++){
            double radius=hWFpp->GetXaxis()->GetBinCenter(iRad+1);
            double WFvalue=0.;
            for (unsigned short usCh=0; usCh<NumChannels; usCh++){
                WFvalue+=Kitty.EvalWaveFun2(iMom,radius,usCh)*weighList[usCh];
                }
                hWFpp->SetBinContent(iRad+1,iMom+1,WFvalue);
            }
        }

    TFile *OutputFile = new TFile("/home/feriorob/PhD_analysis/analysis/output_files/CF_pp_outputs/CF_pp_mid_steps/ppCF_WF_Source.root","recreate");
    OutputFile->cd();
    gCFpp->Write("gCFpp");
    hWFpp->GetXaxis()->SetTitle("r (fm)");
    hWFpp->GetYaxis()->SetTitle("k* (MeV/c)");
    hWFpp->SetTitle("pp, AV18, |#psi|^{2} ");
    hWFpp->Write("hWFpp");
    OutputFile->Close();
}