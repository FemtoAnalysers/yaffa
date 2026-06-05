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
#include "TGraphErrors.h"



/*void SetParValues(TRandom &rand, std::vector<double> &lambda_variations, double &lGenCoeff, double &lppLamCoeff, double &lFlatCoeff, double &lFakeCoeff, double &r0_val, bool indip_var=true){
    int size = lambda_variations.size();
    if(indip_var){
        lGenCoeff=lambda_variations[rand.Integer(size)];
        lppLamCoeff=lambda_variations[rand.Integer(size)];
        lFlatCoeff=lambda_variations[rand.Integer(size)];
        lFakeCoeff=lambda_variations[rand.Integer(size)];
        r0_val = r0_variations[rand.Integer(43)];
    }
    else{
        lGenCoeff=lambda_variations[rand.Integer(size)];
        lppLamCoeff=lGenCoeff;
        lFlatCoeff=lGenCoeff;
        lFakeCoeff=lGenCoeff;
        r0_val = r0_variations[rand.Integer(43)];
    }
}*/



void CFpp_fitR0() // "main" function
{
    const double massProton = ParticleMass::MassProton;
    const double mu12=massProton/2.; // m1*m2/(m1+m2)
    const double SourceRadValue_pL = 1.2;

    const double lParGen = LambdaParameters_pp::LambdaPar_pp_Gen;
    const double lParpp_lam = LambdaParameters_pp::LambdaPar_pp_lam;
    const double lParFlat = LambdaParameters_pp::LambdaPar_pp_Flat;
    const double lParFake = LambdaParameters_pp::LambdaPar_pp_Fake;

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
    hData->SetName("hData");
    TGraphErrors *gData = HistToGraph(hData, kMin, kMax, nBins);
    gData->SetName("gData");

    //***************************** CF plambda *****************************************
    CATS CF_plambda;
    CF_plambda.SetMomBins(nBins,kMin,kMax);
    CF_plambda.SetQ1Q2(0);
    CF_plambda.SetQuantumStatistics(false);
    AnalysisObject_pL.SetUpCats_pL(CF_plambda,"Chiral_Coupled_SPD","Gauss",0,0);

    double CuspWeight = 0.33;//0.54 to take account of peak
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

    double lam_flat = (15.7+19.0+17.6)/(100-0.6);
    DLM_CkDecomposition CkDec_pL("pL",1,Ck_pL,NULL);
    CkDec_pL.AddContribution(0,lam_flat,DLM_CkDecomposition::cFeedDown);
    CkDec_pL.Update();

    //***************************** CF pp with variations on lPars and r0 *****************************************
    //std::vector<double> lambda_variations = {0.9, 1.0, 1.1}; //10% variations on lambda parameters
    const int n_r0var=34; // for start evaluation
    //const int n_r0var=16; // after posterior results
    //double r0_variations[n_r0var] = {1.19, 1.20, 1.21, 1.215, 1.22, 1.225, 1.23,1.232, 1.234, 1.236, 1.238, 1.24, 1.242, 1.244, 1.246, 1.248, 1.25, 1.252, 1.254, 1.256, 1.258, 1.26, 1.262, 1.264, 1.266, 1.268, 1.27, 1.275, 1.28, 1.285, 1.29, 1.295, 1.3, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36,1.37, 1.38,1.39, 1.40}; //
    //double r0_variations[n_r0var] ={1.23,1.232, 1.234, 1.236, １.
    //const int n_r0var=1;
    //double r0_variations[n_r0var] ={1.23};  
    double r0_variations[n_r0var] ={1.20,1.205,1.21,1.215,1.22,1.225,1.23,1.235,1.24,1.245,1.25,1.255,1.26,1.265,1.27,1.275,1.28,1.285,1.29,1.295,1.3,1.305,1.31,1.315,1.32,1.325,1.33,1.335,1.34,1.345,1.35, 1.355, 1.36, 1.365}; // variations of r0 to test in the fit for pol 2, pol3
    //double r0_variations[n_r0var] ={1.15,1.155,1.16,1.165,1.17,1.175,1.18,1.185,1.19,1.195,1.20,1.205,1.21,1.215,1.22,1.225,1.23,1.235,1.24,1.245,1.25,1.255,1.26,1.265,1.27,1.275,1.28,1.285,1.29,1.295,1.3,1.305,1.31,1.315};
    std::vector<double> lambda_variations = {0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.0,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.2}; //10% variations on lambda parameters
    int nRep=lambda_variations.size();
    TH2F* hChi2Map = new TH2F("hChi2Map","Chi2 map for variations of r0 and lambda parameters; r0 (fm); lambda variation", n_r0var, 1.1975, 1.3675, nRep, 0.795, 1.205); // to store the chi2 values for each variation of r0 and lambda parameters, for the posterior evaluation
    //TH2F* hChi2Map = new TH2F("hChi2Map","Chi2 map for variations of r0 and lambda parameters; r0 (fm); lambda variation", n_r0var, 1.1475, 1.3175, nRep, 0.795, 1.205); // to store the chi2 values for each variation of r0 and lambda parameters, for the posterior evaluation
    
    //TRandom3 rand(2); //fixed seed for reproducibility
     //3;// or -+1 or 0 for lambdaPar
     // to store the chi2 values for each variation of r0 and lambda parameters


    double r0_val=0;
    /*double lGenCoeff=0;
    double lppLamCoeff=0;
    double lFlatCoeff=0;
    double lFakeCoeff=0;*/
    double lCoeff=0.;
    
    double LGenVal=0;
    double LppLVal=0;
    double LFlatVal=0;
    double LFakeVal=0;

    double kfitMin = 0;
    double kfitMax = 375.;
    TFile *outputFile = new TFile("/home/feriorob/PhD_analysis/analysis/output_files/CF_pp_outputs/CF_pp_mid_steps/Chi2MapPol4Wolt_20percent.root", "RECREATE");
    outputFile->cd();
    
    //NO PREFIT!!!!
    TGraph *gPosterior_r0 = new TGraph(); // to store the posterior distribution of r0 values after the variation

    // loop over variations of r0 and lambda parameters, compute the CF for each variation,

    for(int index_r=0; index_r<n_r0var; index_r++){
        r0_val = r0_variations[index_r];
        double posterior_sum = 0;
        CATS CF_pp_gen;
        CF_pp_gen.SetMomBins(nBins, kMin, kMax);
        CF_pp_gen.SetQ1Q2(1);
        CF_pp_gen.SetQuantumStatistics(true);
        DLM_CommonAnaFunctions AnalysisObject;
        AnalysisObject.SetCatsFilesFolder("/home/feriorob/PhD_analysis/analysis/initial_test/CernBox");
        Use_P_D_waves(CF_pp_gen, true, true);// true true means p- and d-wabes included
        AnalysisObject.SetUpCats_pp(CF_pp_gen,"AV18","Gauss",0,0);
        CF_pp_gen.SetAnaSource(0, r0_val);
        CF_pp_gen.KillTheCat();
        DLM_Ck Ck_pp(1,8, CF_pp_gen);
        Ck_pp.Update();
        for(int index_rep=0; index_rep<nRep; index_rep++){
            std::cout<<"Starting repetition "<<index_rep+1<<" out of "<<nRep<<"for r0 = "<<r0_val<<std::endl;
            //SetParValues(rand, lambda_variations, lGenCoeff, lppLamCoeff, lFlatCoeff, lFakeCoeff,r0_val, indip_var);
            lCoeff = lambda_variations[index_rep];
            LppLVal=lCoeff*lParpp_lam;
            LFlatVal=lCoeff*lParFlat;
            LFakeVal=lCoeff*lParFake;
            if(LFakeVal+LppLVal+LFlatVal>1) { // no physics case
                LFakeVal = lParFake;
                LppLVal = lParpp_lam;
                LFlatVal = lParFlat;
            }
            LGenVal = 1. - LFakeVal - LppLVal - LFlatVal; // to ensure that the sum of the lambda parameters is 1, we fix the fake lambda to be the complement of the others, to high variation of 10% on lGen
           
          

            DLM_CkDecomposition CkDec_pp("pp",3,Ck_pp,hResolution_pp);

            CkDec_pp.AddContribution(0,LppLVal,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);//this is how you add the pL feed-down
            CkDec_pp.AddContribution(1,LFlatVal,DLM_CkDecomposition::cFeedDown);
            CkDec_pp.AddContribution(2,LFakeVal,DLM_CkDecomposition::cFake);

            CkDec_pp.AddPhaseSpace(0,hPhaseSpace);
            CkDec_pp.Update();
            //DLM_CkDecomposition CkDec_pp = ComputeTot_ppCF(hResolution_pp,hResidual_pp_pL,hPhaseSpace, CkDec_PL, LppLVal, LFlatVal, LFakeVal, r0_val, nBins, kMin, kMax);

            TF1* fit_pp = new TF1("fit_pp", // fit with a neutered pol3, no linear term
                [&](double* x, double* par) {
                    double MOM = x[0];
                    //
                    return (par[0]+par[1]*MOM*MOM+par[2]*MOM*MOM*MOM+par[3]*MOM*MOM*MOM*MOM) *CkDec_pp.EvalCk(MOM);
                }, 
                kfitMin, kfitMax, 4);



            hData->Fit(fit_pp,"MRNQ");
            std::cout<<"parameters of fit: "<<std::endl;
            /*std::cout<<"Parameter 0: "<<fit_pp->GetParameter(0)<<std::endl;
            std::cout<<"Parameter 1: "<<fit_pp->GetParameter(1)<<std::endl;
            std::cout<<"Parameter 2: "<<fit_pp->GetParameter(2)<<std::endl;
            //std::cout<<"Parameter 3: "<<fit_pp->GetParameter(3)<<std::endl;

            std::cout<<"par fit errors: "<<std::endl;
            std::cout<<"Parameter 0 error: "<<fit_pp->GetParError(0)<<std::endl;
            std::cout<<"Parameter 1 error: "<<fit_pp->GetParError(1)<<std::endl;
            std::cout<<"Parameter 2 error: "<<fit_pp->GetParError(2)<<std::endl;
           \ std::cout<<"Parameter 3 error: "<<fit_pp->GetParError(3)<<std::endl;*/


            std::cout<<"Fit results for repetition "<<index_rep+1<<": "<<std::endl;
            double chi2 = fit_pp->GetChisquare();
            cout<<"Chi2 = "<<chi2<<std::endl;
            double ndf = fit_pp->GetNDF();
            hChi2Map->Fill(r0_val, lCoeff, chi2);
            posterior_sum +=TMath::Exp(-0.5*chi2);
            std::cout<<"posterior value for r0 = "<<r0_val<<" is "<<posterior_sum<<std::endl;
        }

        gPosterior_r0->SetPoint(index_r, r0_val, posterior_sum/nRep);
    }

    gPosterior_r0->SetName("gPosterior_r0");
    gPosterior_r0->Write();   
    hChi2Map->Write();
    


    /*for(int i=0; i<fitGraphs.size(); i++){
        TCanvas *cFit = new TCanvas(Form("cFit_%d", i), Form("cFit_%d", i), 800, 600);
        fitGraphs[i]->SetTitle(Form("Fit for repetition %d", i+1));
        fitGraphs[i]->GetXaxis()->SetTitle("k* (MeV/c)");
        fitGraphs[i]->GetYaxis()->SetTitle("C(k*)");
        hData->SetMarkerStyle(20);
        hData->SetMarkerColor(kBlack);
        hData->SetLineColor(kBlack);
        hData->Draw("E");
        fitGraphs[i]->Draw("L");
        modelGraphs[i]->Draw("L same");
        fitFunctions[i]->Draw("L same");
        cFit->Write();
    }*/
    outputFile->Close();
}