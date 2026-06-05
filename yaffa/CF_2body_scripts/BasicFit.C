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

void CATS_pp_AV181_SwaveOnly(CATS & Kitty, const bool & pwaves, const bool & dwaves);


void BasicFit()
{
    const double massProton = ParticleMass::MassProton;
    const double mu=massProton/2.; // m1*m2/(m1+m2)
    const double SourceRadValue=1.25;
    const double SourceRadValue_pL = 1.2;
    
    const double lambdaParGen = LambdaParameters_pp::LambdaPar_pp_Gen;
    const double lambdaParpp_lam = LambdaParameters_pp::LambdaPar_pp_lam;

    const int nBins=250;
    const double kMin=0.;
    const double kMax=1000.;

    //********************************CF plambda ************************************
    CATS CF_plambda;
    CF_plambda.SetMomBins(nBins,kMin,kMax);
    CF_plambda.SetQ1Q2(0);
    CF_plambda.SetQuantumStatistics(false); 

    DLM_CommonAnaFunctions AnalysisObject_pL;
    AnalysisObject_pL.SetCatsFilesFolder("/home/feriorob/PhD_analysis/analysis/initial_test/CernBox");
    TFile *InputFileRes_pl = new TFile("inputfiles/ALICE_pp_13TeV_MEpL.root");
    AnalysisObject_pL.SetUpCats_pL(CF_plambda,"Chiral_Coupled_SPD","Gauss",0,0);//NLO 13 cutoff 600 spd waves
//    AnalysisObject_pL.SetUpCats_pL(CF_proton_proton,"Chiral_Coupled_SPD","Gauss",-11600,0);//LO 13 cutoff 600
    double CuspWeight = 0.33;//0.54 tiene in conto del picchettino 
    CF_plambda.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s) // 
    CF_plambda.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
    CF_plambda.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
    CF_plambda.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
    CF_plambda.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)
    CF_plambda.SetAnaSource(0, SourceRadValue_pL);
    CF_plambda.SetNotifications(CATS::nError);

    CF_plambda.KillTheCat();

    TGraph *gCF_plambda = new TGraph(nBins);
    for (size_t k_index=0; k_index < nBins; k_index++){
        double CFvalue = CF_plambda.GetCorrFun(k_index);
        double kstar_value = CF_plambda.GetMomentum(k_index);
        gCF_plambda->SetPoint(k_index, kstar_value, CFvalue);
    }

    DLM_Ck Ck_pL(1,0, CF_plambda);
    Ck_pL.Update();

    double lam_flat = (15.7+19.0+17.6)/(100-0.6);
    DLM_CkDecomposition CkDec_pL("pL",1,Ck_pL,NULL);
    CkDec_pL.AddContribution(0,lam_flat,DLM_CkDecomposition::cFeedDown);
    CkDec_pL.Update();

    TGraph *genModelGraph_pL = new TGraph(nBins);
    for (size_t k_index=0; k_index < nBins; k_index++){
        genModelGraph_pL->SetPoint(k_index,CkDec_pL.GetCk()->GetBinCenter(0,k_index),CkDec_pL.EvalCk(CkDec_pL.GetCk()->GetBinCenter(0,k_index)));
    }

    //***************************** CF pp *****************************************
    CATS KittyAV18;
    KittyAV18.SetMomBins(nBins,kMin,kMax);
    CATSparameters* cParSource = new CATSparameters(CATSparameters::tSource, 1, true);
    cParSource->SetParameter(0, SourceRadValue); 
    KittyAV18.SetAnaSource(GaussSource, *cParSource);

    KittyAV18.SetQ1Q2(1);
    KittyAV18.SetQuantumStatistics(true);
    KittyAV18.SetRedMass(mu);
    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder("/home/feriorob/PhD_analysis/analysis/initial_test/CernBox");

    //CATS_pp_AV181_SwaveOnly(KittyAV18, true, true);
    AnalysisObject.SetUpCats_pp(KittyAV18,"AV18","Gauss",0,0);
    KittyAV18.SetAnaSource(0, SourceRadValue);
    //KittyAV18.SetNotifications(CATS::nError);
    KittyAV18.KillTheCat();

    //build graphs without decomposition
    TGraph* gCF_Gen = new TGraph(nBins);
    //TGraph* gCF_Gen_lambdaPar = new TGraph(nBins);
    for (size_t k_index=0; k_index < nBins; k_index++){
        double CFvalue = KittyAV18.GetCorrFun(k_index);
        double kstar_value = KittyAV18.GetMomentum(k_index);
        gCF_Gen->SetPoint(k_index, kstar_value, CFvalue);
    }

     //decomposition
    DLM_Ck Ck_pp(1,8, KittyAV18);
    Ck_pp.Update();

    double lambdaParFlat = 1.-lambdaParGen-lambdaParpp_lam- 0.011;

    TFile *InputFileME = new TFile("inputfiles/CFOutput_pp.root"); 

    TH1F *hPhaseSpace = (TH1F*)((TList*)((TList*)InputFileME->Get("PairDist"))
                  ->FindObject("Pair"))
                  ->FindObject("MEDist_Particle0_Particle0_clone");
    if (!hPhaseSpace) {
        std::cerr << "Error: Histogram 'PairDist/Pair/MEDist_Particle0_Particle0_clone' not found in file." << std::endl;
    }
    hPhaseSpace->GetXaxis()->SetLimits(hPhaseSpace->GetXaxis()->GetXmin()*1000.,hPhaseSpace->GetXaxis()->GetXmax()*1000.);

    TFile *InputFileRes = new TFile("inputfiles/ALICE_pp_13TeV_MEpp.root");
    TFile *InputFileResidual_pL = new TFile("inputfiles/Decay_matrices_2020.root");
    TH2F *hResidual_pp_pL = (TH2F*)InputFileResidual_pL -> Get("hRes_pp_pL");
    TH2F* hResolution_pp = (TH2F*)InputFileRes -> Get("h_RESO_pp_MeV");
    //hResolution_pp->RebinX(4);
    //hResolution_pp->RebinY(4);
    //hResidual_pp_pL->RebinX(4);
    //hResidual_pp_pL->RebinY(4);

    DLM_CkDecomposition CkDec_pp_Tot("pp",3,Ck_pp,hResolution_pp);


    CkDec_pp_Tot.AddContribution(0,lambdaParFlat,DLM_CkDecomposition::cFeedDown);
    CkDec_pp_Tot.AddContribution(1,lambdaParpp_lam,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
    CkDec_pp_Tot.AddContribution(2,0.011,DLM_CkDecomposition::cFake);
    CkDec_pp_Tot.Update();
    CkDec_pp_Tot.AddPhaseSpace(0,hPhaseSpace);
    CkDec_pp_Tot.Update();

    //building the graph with decomposition
    TGraph* genModelGraph = new TGraph(nBins);

    for (size_t k_index=0; k_index < nBins; k_index++){
        genModelGraph->SetPoint(k_index,CkDec_pp_Tot.GetCk()->GetBinCenter(0,k_index),CkDec_pp_Tot.EvalCk(CkDec_pp_Tot.GetCk()->GetBinCenter(0,k_index)));
    }

    genModelGraph->SetLineColor(kGreen+2);
    genModelGraph->SetLineWidth(2);

    double kfitMin = 4;
    double kfitMax = 375.;
    TH1F *hData = (TH1F*)InputFileME -> Get("hCk_ReweightedMeV_0");
    TGraphErrors *gData = HistToGraph(hData, kMin, kMax, nBins);
    TF1 * pol3_nolinear = new TF1("pol3_nolinear", "[0]+[1]*x*x/1000000.+[2]*x*x*x/1000000000.+[3]*x/1000.", 0., 600.);  // baseline fit
    //TF1 * pol3_nolinear = new TF1("pol3_nolinear", "[0]+[1]*x*x/1000000.+[2]*x*x*x/1000000000.", 300., 800.);  // baseline fit
    pol3_nolinear->SetParameter(0,1.);
    //pol3_nolinear->SetParameter(3,0.);
    hData->Fit(pol3_nolinear," S, R, M,N");
    double bkg_baseline = pol3_nolinear->GetParameter(0);
    double bkg_quad = pol3_nolinear->GetParameter(1);
    double bkg_cubic = pol3_nolinear->GetParameter(2);
    double linear = pol3_nolinear->GetParameter(3);
    std::cout << "chi2/ndf baseline fit: " << pol3_nolinear->GetChisquare() << "/" << pol3_nolinear->GetNDF() << std::endl;
    std::cout << "Baseline fit parameters: " << std::endl;
    std::cout << "Baseline: " << bkg_baseline << std::endl;
    std::cout << "Quadratic: " << bkg_quad << std::endl;
    std::cout << "Cubic: " << bkg_cubic << std::endl;
    std::cout << "Linear: " << linear << std::endl;
    pol3_nolinear->SetLineColor(kOrange-3);
    pol3_nolinear->SetLineWidth(2);
    pol3_nolinear->SetLineStyle(7);

   
    TF1* fit_pp = new TF1("fit_pp",[&](double* x, double* par) {
            double MOM = x[0];
            return (par[0] + par[1]*MOM*MOM/1000000.+ par[2]*MOM*MOM*MOM/1000000000.+par[3]*MOM/1000.) *CkDec_pp_Tot.EvalCk(MOM);
            //return (par[0] + par[1]*MOM*MOM/1000000.+ par[2]*MOM*MOM*MOM/1000000000.) *CkDec_pp_Tot.EvalCk(MOM);
            //return par[0]*(1+par[1]*MOM*MOM+par[2]*MOM*MOM*MOM) *CkDec_pp_Tot.EvalCk(MOM);
        },
        kfitMin, kfitMax, 4);

    fit_pp->SetNpx(1024);

    /*if(fit_pp->GetParameter(0) > 0.) {
        fit_pp->SetParLimits(0, 0.1*bkg_baseline, 10*bkg_baseline);
        //fit_pp->SetParameter(0, bkg_baseline);
    } else {
        fit_pp->SetParLimits(0, 10*bkg_baseline, 0.1*bkg_baseline);
        //fit_pp->SetParameter(0, bkg_baseline);
    }
    if(fit_pp->GetParameter(1) > 0.) {
        fit_pp->SetParLimits(1, 0.1*bkg_quad, 10*bkg_quad);
        //fit_pp->SetParameter(1, bkg_quad);
    } else {
        fit_pp->SetParLimits(1, 10*bkg_quad, 0.1*bkg_quad);
        //fit_pp->SetParameter(1, bkg_quad);
    }
    if(fit_pp->GetParameter(2) > 0.) {
        fit_pp->SetParLimits(2, 0.1*bkg_cubic, 10*bkg_cubic);
        fit_pp->SetParameter(2, bkg_cubic);
    } else {
        fit_pp->SetParLimits(2, 10*bkg_cubic, 0.1*bkg_cubic);
        fit_pp->SetParameter(2, bkg_cubic);
    }
    if(fit_pp->GetParameter(3) > 0.) {
        fit_pp->SetParLimits(3, 0.1*linear, 10*linear);
        fit_pp->SetParameter(3, linear);
    } else {
        fit_pp->SetParLimits(3, 10*linear, 0.1*linear);
        fit_pp->SetParameter(3, linear);
    }*/
    
    hData->Fit(fit_pp,"MR+", "", kfitMin, kfitMax);
    std::cout << "Fit parameters: " << std::endl;
    std::cout << "Baseline: " << fit_pp->GetParameter(0) << std::endl;
    std::cout << "Quadratic: " << fit_pp->GetParameter(1) << std::endl;
    std::cout << "Cubic: " << fit_pp->GetParameter(2) << std::endl;
    std::cout << "Linear: " << fit_pp->GetParameter(3) << std::endl;

    std::cout << "chi2/ndf pp fit: " << fit_pp->GetChisquare() << "/" << fit_pp->GetNDF() << std::endl;
    std::cout<<"chi2 hist fit: "<<fit_pp->GetChisquare()<<"/"<<fit_pp->GetNDF()<<std::endl;
    for (int i=0; i<gData->GetN(); i++){
        gData->SetPointError(i, 0,gData->GetErrorY(i));
    }
    //gData->Fit(fit_pp,"QMRN+", "", kfitMin, kfitMax);
    //std::cout<<"chi2/ndf with x errors =0: "<<fit_pp->GetChisquare()<<"/"<<fit_pp->GetNDF()<<std::endl;
    TGraph* fitGraph = FitToGraph(fit_pp, kMin, kMax, nBins);
    fitGraph->SetLineColor(kBlue);
    fitGraph->SetLineWidth(2);

    TF1* pol3_updated = new TF1("pol3_updated", "[0]+[1]*x*x/1000000.+[2]*x*x*x/1000000000.+[3]*x/1000.", 4., 600.);
    //TF1* pol3_updated = new TF1("pol3_updated", "[0]+[1]*x*x/1000000.+[2]*x*x*x/1000000000.", 4., 800.);
    //TF1* pol3_updated = new TF1("pol3_updated", "[0]*(1+[1]*x*x+[2]*x*x*x)", 4., 800.);

    pol3_updated->SetParameter(0, fit_pp->GetParameter(0));
    pol3_updated->SetParameter(1, fit_pp->GetParameter(1));
    pol3_updated->SetParameter(2, fit_pp->GetParameter(2));
    pol3_updated->SetParameter(3, fit_pp->GetParameter(3));

    pol3_updated->SetLineColor(kRed);
    pol3_updated->SetLineWidth(2);

    TFile *outputFile = new TFile("/home/feriorob/PhD_analysis/analysis/output_files/CF_pp_outputs/CF_pp_mid_steps/BasicFitLacyComparison1.root", "RECREATE");
    outputFile->cd();
    gCF_Gen->Write("gCF_Gen");
    genModelGraph->Write("genModelGraph");
    gCF_plambda->Write("gCF_plambda");
    genModelGraph_pL->Write("genModelGraph_pL");
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->cd();
    gData->SetMarkerStyle(20);
    gData->SetMarkerSize(1);
    gData->SetMarkerColor(kBlack);
    gData->SetLineColor(kBlack);
    gData->Draw("E");
    genModelGraph->Draw("Same");
    fitGraph->Draw("Same");
    pol3_updated->Draw("Same");
    pol3_nolinear->Draw("Same");
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(gData, "Data", "lep");
    legend->AddEntry(pol3_nolinear, "Baseline Fit", "l");
    legend->AddEntry(genModelGraph, "Generated CF", "l");
    legend->AddEntry(fitGraph, "Fit with CF", "l");
    legend->AddEntry(pol3_updated, "Updated Baseline Fit", "l");
    legend->Draw();
    c1->Write();

    outputFile->Close();
}

void CATS_pp_AV181_SwaveOnly(CATS & Kitty, const bool & pwaves, const bool & dwaves) {
  //the 4 channels for pp are:
  //s=0: 1S0 + 1D2
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