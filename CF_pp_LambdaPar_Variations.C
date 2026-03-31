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

double GetRandomElement(const std::vector<double>& vec, TRandom3& rand) {
    int index = rand.Integer(vec.size());
    return vec[index];
}

void CF_pp_var_lambdaPar(){

    //******************* input constants ******************
    const double massProton = ParticleMass::MassProton;
    const double mu12=massProton/2.; // m1*m2/(m1+m2)
    const double SourceRadValue=1.249; 
    const double SourceRadValue_pL = 1.2;

    const double lambdaParGen = LambdaParameters_pp::LambdaPar_pp_Gen;
    const double lambdaParpp_lam = LambdaParameters_pp::LambdaPar_pp_lam;

    const int nBins=250;
    const double kMin=0.;
    const double kMax=1000.;

    //******************* input files *************************
    //plambda
    DLM_CommonAnaFunctions AnalysisObject_pL;
    AnalysisObject_pL.SetCatsFilesFolder("/home/feriorob/PhD_analysis/analysis/initial_test/CernBox");
    TFile *InputFileRes_pl = new TFile("inputfiles/ALICE_pp_13TeV_MEpL.root");
    TH2F* hResolution_pl = (TH2F*)InputFileRes_pl -> Get("h_RESO_pL_MeV"); //resolution matrix for pL

    //pp
    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder("/home/feriorob/PhD_analysis/analysis/initial_test/CernBox");
    TFile *InputFileME = new TFile("inputfiles/CFOutput_pp.root"); 
    TH1F *hPhaseSpace = (TH1F*)((TList*)((TList*)InputFileME->Get("PairDist"))
                  ->FindObject("Pair"))
                  ->FindObject("MEDist_Particle0_Particle0_clone");
    if (!hPhaseSpace) {
        std::cerr << "Error: Histogram 'PairDist/Pair/MEDist_Particle0_Particle0_clone' not found in file." << std::endl;
    }
    hPhaseSpace->GetXaxis()->SetLimits(hPhaseSpace->GetXaxis()->GetXmin()*1000.,hPhaseSpace->GetXaxis()->GetXmax()*1000.);

    TFile *InputFileRes = new TFile("inputfiles/ALICE_pp_13TeV_MEpp.root");
    TFile *InputFileResidual_pL = new TFile("inputfiles/Decay_matrices_2020.root"); //residual matrix for feed-down from pL to pp
    TH2F *hResidual_pp_pL = (TH2F*)InputFileResidual_pL -> Get("hRes_pp_pL");
    TH2F* hResolution_pp = (TH2F*)InputFileRes -> Get("h_RESO_pp_MeV"); //resolution matrix for pp

    TH1F *hData = (TH1F*)InputFileME -> Get("hCk_ReweightedMeV_0"); //data histogram for pp
    //***************************** random variations *******************************
    std::vector<double> LambdaParsCoeff = {0.8, 0.84, 0.88, 0.92, 0.96, 1.00, 1.04, 1.08, 1.12, 1.16, 1.20};

    //***************************** CF plambda *****************************************
    CATS CF_plambda;
    CF_plambda.SetMomBins(nBins,kMin,kMax);
    CF_plambda.SetQ1Q2(0);
    CF_plambda.SetQuantumStatistics(false);
    AnalysisObject_pL.SetUpCats_pL(CF_plambda,"Chiral_Coupled_SPD","Gauss",0,0);

    double CuspWeight = 0.33;//0.54 tiene in conto del picchettino 
    CF_plambda.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s) // 
    CF_plambda.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
    CF_plambda.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
    CF_plambda.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
    CF_plambda.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)
    CF_plambda.SetAnaSource(0, SourceRadValue_pL);
    CF_plambda.SetNotifications(CATS::nError);

    CF_plambda.KillTheCat();
    DLM_Ck Ck_pL(1,0, CF_plambda);
    Ck_pL.Update();

    //***************************** genuine pp CF ********************************
    CATS CF_pp_gen;
    CF_pp_gen.SetMomBins(nBins, kMin, kMax);
    CATSparameters* cParSource = new CATSparameters(CATSparameters::tSource, 1, true);
    //cParSource->SetParameter(0, SourceRadValue);
    //KittyAV18.SetAnaSource(GaussSource, *cParSource);
    CF_pp_gen.SetQ1Q2(1); // 1 same charge; 0 neutral charge; -1
    CF_pp_gen.SetQuantumStatistics(true); // set it to true if you want quantum statistics
    CF_pp_gen.SetRedMass(mu12);
    //CF_pp_gen.SetQ1Q2(1);
    //CF_pp_gen.SetPdgId(2212, 2212);
    AnalysisObject.SetUpCats_pp(CF_pp_gen,"AV18","Gauss",0,0);
    CF_pp_gen.SetAnaSource(0, SourceRadValue);
    CF_pp_gen.SetNotifications(CATS::nError);

    //CATS_pp_AV181_SwaveOnly(CF_pp_gen, true, true); // true true means p- and d-wabes included
    CF_pp_gen.KillTheCat();

    DLM_Ck Ck_pp(1,8, CF_pp_gen);
    Ck_pp.Update();


    //**************************+ loop on variations *******************************
    int nRepetition = 500;
    TRandom3 rand(0); // random seed, can be changed for different variations

    TF1* pol3_fit_nosignal = new TF1("pol3_fit_nosignal","[0]+[1]*x+[2]*x*x+[3]*x*x*x", 300, kMax); //side band region for fitting the baseline    
    hData->Fit(pol3_fit_nosignal,"Q, S, N, R, M");
    double bkg_baseline = pol3_fit_nosignal->GetParameter(0);
    double bkg_slope = pol3_fit_nosignal->GetParameter(1);
    double bkg_quad = pol3_fit_nosignal->GetParameter(2);
    double bkg_cubic = pol3_fit_nosignal->GetParameter(3);

    double kfitMin = 0;
    double kfitMax = 375;

    std::vector<std::vector<double>> final_fit_values(nRepetition, std::vector<double>(nBins));
    TH1F* SystErr_point5 = new TH1F("SystErr_point5", "SystErr_point5", 2500,0.,5.);
    TH1F* SystErr_point100= new TH1F("SystErr_point100", "SystErr_point100", 2500,0.,5.);
    TH1F* SystErr_point25= new TH1F("SystErr_point25", "SystErr_point25", 2500,0.,5.);
    TH1F* SystErr_point50= new TH1F("SystErr_point50", "SystErr_point50", 2500,0.,5.);

    for (int index_var=0; index_var<nRepetition; index_var++){

        double v_lGEN = GetRandomElement(LambdaParsCoeff, rand);
        double v_lpL = GetRandomElement(LambdaParsCoeff, rand);
        double v_lpSp = GetRandomElement(LambdaParsCoeff, rand);
        double v_lTOT = GetRandomElement(LambdaParsCoeff, rand);
     
        double lam_flat = v_lTOT*(15.7+19.0+17.6)/(100-0.6);
        if(lam_flat >1) lam_flat = 1.;
    
        mattermost CkDec_pL("pL",1,Ck_pL,NULL);
        CkDec_pL.AddContribution(0,lam_flat,DLM_CkDecomposition::cFeedDown);
        CkDec_pL.Update();

        double lam_pp[] = {v_lGEN*0.67,v_lpL*0.203,v_lpSp*0.116,1.-v_lGEN*0.67-v_lpL*0.203-v_lpSp*0.116};

        if (1.-v_lGEN*0.67-v_lpL*0.203-v_lpSp*0.116<0) {
            lam_pp[3] = 0.;
        }


        DLM_CkDecomposition CkDec_pp("pp",3,Ck_pp,hResolution_pp);

        CkDec_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);//this is how you add the pL feed-down
        CkDec_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
        CkDec_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);

        CkDec_pp.Update();
        CkDec_pp.AddPhaseSpace(hPhaseSpace);

        //fit
        DLM_CkDecomposition* GraphForFit;
        GraphForFit = &CkDec_pp;
        TF1* fit_pol3 = new TF1(
            "pol3",
            [GraphForFit](double* x, double* par) {
                double Mom = x[0];
                return (par[0]+par[1]*Mom+par[2]*Mom*Mom+par[3]*Mom*Mom*Mom)*GraphForFit->EvalCk(Mom);
                //return (par[0]+par[1]*Mom*Mom+par[2]*Mom*Mom*Mom)*GraphForFit->EvalCk(Mom);
                },
            kfitMin, kfitMax, 4
          );
        fit_pol3->SetParameter(0,bkg_baseline);
        fit_pol3->SetParameter(1,bkg_slope);
        fit_pol3->SetParameter(2,bkg_quad);
        fit_pol3->SetParameter(3,bkg_cubic);
        hData->Fit(fit_pol3,"Q, S, N, R, M");
        double chi2 = fit_pol3->GetChisquare();
        double ndf = fit_pol3->GetNDF();
        cout<<"Variation "<<index_var<<": chi2/ndf = "<<chi2<<"/"<<ndf<<endl;

        TGraph* fitGraph = FitToGraph(fit_pol3, kfitMin, kfitMax, nBins);
        for (int i=0; i<nBins; i++){
            final_fit_values[index_var][i] = fitGraph->GetY()[i];
        }
        delete fit_pol3;
        delete fitGraph;
    }
    TGraph* gSystErr = new TGraph();
    for (int point_index=0; point_index<nBins; point_index++){
        TH1F* hPointValues = new TH1F(Form("hPointError_%d", point_index), Form("hPointError_%d", point_index), 2500,0.,5.);
            for (int rep_index=0; rep_index<nRepetition; rep_index++){
                hPointValues->Fill(final_fit_values[rep_index][point_index]);
                if(point_index==5) SystErr_point5->Fill(final_fit_values[rep_index][point_index]);
                if(point_index==25) SystErr_point25->Fill(final_fit_values[rep_index][point_index]);
                if(point_index==50) SystErr_point50->Fill(final_fit_values[rep_index][point_index]);
                if(point_index==100) SystErr_point100->Fill(final_fit_values[rep_index][point_index]);
                }
            double sysErr = hPointValues->GetStdDev();
            gSystErr->SetPoint(point_index, kMin + point_index*(kMax-kMin)/nBins, sysErr);
            delete hPointValues;
    
  }
    TFile* outputFile = new TFile("/home/feriorob/PhD_analysis/analysis/output_files/pp_CF_SystErrors.root", "RECREATE");
    outputFile->cd();
    gSystErr->Write("gSystErr");
    SystErr_point5->GetXaxis()->SetTitle("C(k*) at k*=20 MeV/c");
    SystErr_point5->Write("SystErr_point5");
    SystErr_point25->GetXaxis()->SetTitle("C(k*) at k*=100MeV/c");
    SystErr_point25->Write("SystErr_point25");
    SystErr_point50->GetXaxis()->SetTitle("C(k*) at k*=200MeV/c");
    SystErr_point50->Write("SystErr_point50");
    SystErr_point100->GetXaxis()->SetTitle("C(k*) at k*=400 MeV/c");
    SystErr_point100->Write("SystErr_point100");
    outputFile->Close();
}