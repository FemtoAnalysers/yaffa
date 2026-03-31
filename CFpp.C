R__LOAD_LIBRARY($HOME/DLM/install/lib/libCATS.so) // good for sunrise setting. To change in case of different installation path


#include "/home/feriorob/DLM/install/include/CATS.h"
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

using namespace std;

void CATS_pp_AV181_SwaveOnly(CATS & Kitty, const bool & pwaves, const bool & dwaves);

//---------------------------

void ppCF() {

    const double massProton = ParticleMass::MassProton;
    const double mu12=massProton/2.; // m1*m2/(m1+m2) //mass of the proton in MeV from PDG
    double hbarc = 197.3269631;

    const int nBins = 250;
    const double kmin = 0.;
    const double kmax = 1000.;
    const double SourceRadValue = 1.249; // Effective Source for pp in Run 2 data (pp collisions)

    cout<<"Solving for AV18"<<endl;

    CATS KittyAV18;
    KittyAV18.SetMomBins(nBins, kmin, kmax);
    CATSparameters* cParSource = new CATSparameters(CATSparameters::tSource, 1, true);
    cParSource->SetParameter(0, SourceRadValue); 
    KittyAV18.SetAnaSource(GaussSource, *cParSource); // Set a Gaussian source with parameters saved in cParSource
    KittyAV18.SetQ1Q2(1); // 1 same charge; 0 neutral charge; -1
    KittyAV18.SetQuantumStatistics(true); // set it to true if you want quantum statistics
    KittyAV18.SetRedMass(mu12);
    CATS_pp_AV181_SwaveOnly(KittyAV18, true, true); // true true means p- and d-wabes included
    KittyAV18.KillTheCat(); 

    TGraph *gCAV18 = new TGraph();

  	for (size_t i = 0; i < nBins; i++) {
        double CFav18 = KittyAV18.GetCorrFun(i);
        double kstarav18 = KittyAV18.GetMomentum(i);
		    gCAV18 -> SetPoint(i,kstarav18,CFav18);

	  }

    TFile *file = new TFile("ppCF_results.root","recreate");

    file->cd();
        gCAV18->Write("gCAV18");
    file->Close();

}

void CATS_pp_AV181_SwaveOnly(CATS & Kitty, const bool & pwaves, const bool & dwaves) {
  //the 4 channels for pp are:
  //s=0: 1S0 + 3D1
  //s=1: 3P0
  //s=1: 3P1
  //s=1: 3P2
  //note that for s=0 the p-waves are Pauli blocked, for s=1 these are the s and d waves
  if (pwaves) {
    Kitty.SetNumChannels(4);
    if (dwaves) Kitty.SetNumPW(0, 3);
    else Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 2);
    Kitty.SetNumPW(2, 2);
    Kitty.SetNumPW(3, 2);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetSpin(2, 1);
    Kitty.SetSpin(3, 1);
    Kitty.SetChannelWeight(0, 3. / 12.);
    Kitty.SetChannelWeight(1, 1. / 12.);
    Kitty.SetChannelWeight(2, 3. / 12.);
    Kitty.SetChannelWeight(3, 5. / 12.);
  } else {
    //important: even with the p-waves switched off, physics wise the spin 1 state still exists and
    //the p-waves are there, just in there asymptotic state (free wave). To include this in the computation,
    //CATS still needs a second channel, even if it is `empty`!
    Kitty.SetNumChannels(2);
    if (dwaves) Kitty.SetNumPW(0, 3);
    else Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 0);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetChannelWeight(0, 1.0 / 4.0);
    Kitty.SetChannelWeight(1, 3.0 / 4.0);
//    Kitty.SetChannelWeight(0, 1.0 );
//    Kitty.SetChannelWeight(1, 0.0 );
  }
  //to set up the strong interaction, one can use the predefined functions available in DLM_Potentials.h
  //the main idea is to always pass the function fDlmPot but with different input parameters, based on which the interaction is set up automatically
  //this works not only for pp, but many other systems are included.
  //The fDlmPot is only used as an interface to easily use different potentials, that are hard coded in DLM_Potentials.h
  //Feel free to expand the data base of this file if you think others will benefit from it. For further details contact Dimi
  //The input parameters are by default 9 and are defined as follows:
  //0: potential flag (defines which potential to use, see the enumerators in DLM_Potentials.h for more info)
  //1: a second flag, that can be used if needed (depending on the definition of the potential)
  //2: total isospin
  //3: 2 x isospin of particle 1
  //4: 2 x isospin of particle 2
  //5: total spin
  //6: l quantum number
  //7: j quantum number
  CATSparameters cPars_pp_1S0(CATSparameters::tPotential, 8, true);
  cPars_pp_1S0.SetParameter(0, NN_ReidV8); //choose the AV18
  cPars_pp_1S0.SetParameter(1, v18_Coupled3P2); //default option, which takes the 3P2 channel from a coupled-channel computation, but in CATS only the first diagonal potential elements is used
  cPars_pp_1S0.SetParameter(2, 1);
  cPars_pp_1S0.SetParameter(3, 1);
  cPars_pp_1S0.SetParameter(4, 1);
  cPars_pp_1S0.SetParameter(5, 0);
  cPars_pp_1S0.SetParameter(6, 0);
  cPars_pp_1S0.SetParameter(7, 0);
  //copy all settings from cPars_pp_1S0, and just change quantum numbers s,l,j
  CATSparameters cPars_pp_3P0(cPars_pp_1S0);
  cPars_pp_3P0.SetParameter(5, 1);
  cPars_pp_3P0.SetParameter(6, 1);
  cPars_pp_3P0.SetParameter(7, 0);
  CATSparameters cPars_pp_3P1(cPars_pp_1S0);
  cPars_pp_3P1.SetParameter(5, 1);
  cPars_pp_3P1.SetParameter(6, 1);
  cPars_pp_3P1.SetParameter(7, 1);
  CATSparameters cPars_pp_3P2(cPars_pp_1S0);
  cPars_pp_3P2.SetParameter(5, 1);
  cPars_pp_3P2.SetParameter(6, 1);
  cPars_pp_3P2.SetParameter(7, 2);
  CATSparameters cPars_pp_1D2(cPars_pp_1S0);
  cPars_pp_1D2.SetParameter(5, 0);
  cPars_pp_1D2.SetParameter(6, 2);
  cPars_pp_1D2.SetParameter(7, 2);
  //plug in the strong potential for each channel and partial wave
  //the arguments are: #WhichChannel,#WhichPartialWave,#PotentialFunction,#PotentialParameters
  Kitty.SetShortRangePotential(0, 0, fDlmPot, cPars_pp_1S0);
  if (pwaves) {
    Kitty.SetShortRangePotential(1, 1, fDlmPot, cPars_pp_3P0);
    Kitty.SetShortRangePotential(2, 1, fDlmPot, cPars_pp_3P1);
    Kitty.SetShortRangePotential(3, 1, fDlmPot, cPars_pp_3P2);
  }
  if (dwaves) {
    Kitty.SetShortRangePotential(0, 2, fDlmPot, cPars_pp_1D2);
  }
  //if later on you would like to switch some contribution off, this can be done with:
  //Kitty.RemoveShortRangePotential(#WhichChannel,#WhichPartialWave);
}