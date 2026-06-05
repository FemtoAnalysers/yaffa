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

void FinalResultsOnlyOneBaseline()
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

    //***************************** CF proton-proton max and min values *****************************************
    std::vector<double> lambda_variations = {0.8,1.0,1.2}; //lambdapar
    //std::vector<double> lambda_variations = {1.0}; //lambdapar
    const int nRep=3;
    double lCOeff=0;
    std::vector<double> r0Values = {1.2350,1.2535};
    const int n_r0var = r0Values.size();
    double r0_val=0;
    //not parameters limits
    std::vector<double> fit_range_max = {375., 325.,425.}; //max fit range for the two r0 values
    const int nIntervals = fit_range_max.size();


    TFile *outfile = new TFile("/home/feriorob/PhD_analysis/analysis/output_files/CF_pp_outputs/CF_pp_mid_steps/PreFinalResult_20Perc_Pol2.root", "RECREATE");
    outfile->cd();

    double LGenVal=0;
    double LppLVal=0;
    double LFlatVal=0;
    double LFakeVal=0;

    double kfitMin = 0;
    double kfitMax = 375.;//std::min(fit_range_max[0], fit_range_max[1]); // to have the same fit range for the two r0 values, we take the minimum of the two max fit ranges. This is to ensure that the fit results are comparable, as the fit range can have an impact on the fit results. For the final result, we can then choose to use different fit ranges for the two r0 values, based on the stability of the fit results in different ranges.
    //double kfitMax = 375.;

    std::vector<TGraph*> genModelGraphs;
    std::vector<TGraph*> fitGraphs;
    std::vector<TF1*> fitBaselineFunc;
    TGraph* genGraphMin = new TGraph(nBins);
    TGraph* genGraphMax = new TGraph(nBins);
    TGraph* fitGraphMin = new TGraph(nBins);
    TGraph* fitGraphMax = new TGraph(nBins);
    TGraph* baselineMin = new TGraph(nBins);
    TGraph* baselineMax = new TGraph(nBins);
    TH2F* hCf_kstar = new TH2F("hCf_kstar","hCf_kstar", nBins, kMin, kMax, 3000, 0., 4.); // histogram to store the fit results for all the variations of the parameters, to be used for the evaluation of the systematic uncertainty on the final result  
    double probs[2] = {0.16, 0.84}; // quantiles for the evaluation of the systematic uncertainty on the final result, corresponding to 68% confidence interval
    double quantiles[2];   

    for(int index_r=0; index_r<r0Values.size(); index_r++){
         r0_val = r0Values[index_r];
         for(int index_rep=0; index_rep<nRep;index_rep++){
            double lCoeff = lambda_variations[index_rep];
            LppLVal=lCoeff*lParpp_lam;
            LFlatVal=lCoeff*lParFlat;
            LFakeVal=lCoeff*lParFake;
            if(LFakeVal+LppLVal+LFlatVal>1) { // no physics case
                LFakeVal = lParFake;
                LppLVal = lParpp_lam;
                LFlatVal = lParFlat;
            }
            LGenVal = 1. - LFakeVal - LppLVal - LFlatVal;

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

            DLM_CkDecomposition CkDec_pp("pp",3,Ck_pp,hResolution_pp);

            CkDec_pp.AddContribution(0,LppLVal,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);//this is how you add the pL feed-down
            CkDec_pp.AddContribution(1,LFlatVal,DLM_CkDecomposition::cFeedDown);
            CkDec_pp.AddContribution(2,LFakeVal,DLM_CkDecomposition::cFake);

            CkDec_pp.AddPhaseSpace(0,hPhaseSpace);
            CkDec_pp.Update();

            TGraph *genModelGraph = new TGraph(nBins);
            for(int index_point=0; index_point<nBins; index_point++){
                genModelGraph->SetPoint(index_point, CkDec_pp.GetCk()->GetBinCenter(0,index_point),CkDec_pp.EvalCk(CkDec_pp.GetCk()->GetBinCenter(0,index_point)));
            }
            genModelGraphs.push_back(genModelGraph);
            for(int index_interval=0; index_interval<nIntervals; index_interval++){
                kfitMax = fit_range_max[index_interval];
                    TF1* fit_pp = nullptr;
                    fit_pp = new TF1("fit_pp", // fit with a neutered pol3 complete
                        [&](double* x, double* par) {
                            double MOM = x[0];
                        return (par[0] + par[1]*MOM*MOM/1000000. +  par[2]*MOM/1000.) *CkDec_pp.EvalCk(MOM);
                        //return (par[0] + par[1]*MOM*MOM/1000000.) *CkDec_pp.EvalCk(MOM);
                        //return (par[0] + par[1]*MOM*MOM/1000000.) *CkDec_pp.EvalCk(MOM); //pol2
                    }, 
                    kfitMin, kfitMax, 2);
             
                    hData->Fit(fit_pp,"Q, S, N, R, M");
                    fitBaselineFunc.push_back(fit_pp); // save the baseline function for later use in the final plot of the baseline variation
                     // to keep the same size of the two vectors for the baseline functions, we push a nullptr in the vector of the fit functions without linear term when the fit function has the linear ter
                    //std::cout<<"Fit results for r0: "<<r0_val<<" and lambda variation: "<<lCoeff<<std::endl;
                    //std::cout<<"Baseline: "<<fit_pp->GetParameter(0)<<" Quadratic: "<<fit_pp->GetParameter(1)<<" Cubic: "<<fit_pp->GetParameter(2)<<" Linear: "<<fit_pp->GetParameter(3)<<std::endl;
                
                    TGraph* fitGraph = FitToGraph(fit_pp, 0., 1000., nBins);
                    fitGraphs.push_back(fitGraph);

                    //delete fit_pp;
                    //delete genModelGraph;
                    //delete fitGraph;
                    std::cout<<"finished fit for r0: "<<r0_val<<" and lambda variation: "<<lCoeff<<" and fit range max: "<<kfitMax<<"and fit type: "<<std::endl;
                    
            }
        }
    }
    std::cout<<"numbero of graph"<<fitGraphs.size()<<std::endl;
    std::cout<<"numero di modelli"<<genModelGraphs.size()<<std::endl;
    std::cout<<"finished loop on variations"<<std::endl;
    // take minimun and maximun of each point of CF and model
    for(int index_point=0; index_point<nBins; index_point++){
        double genMin = 999;
        double genMax = -999;
        double baselineMinVal = 999;
        double baselineMaxVal = -999;
        for(int index_graph=0; index_graph<fitGraphs.size(); index_graph++){
            if(index_graph < genModelGraphs.size()){
                double genVal = genModelGraphs[index_graph]->Eval(genModelGraphs[index_graph]->GetX()[index_point]);
                if(genVal<genMin) genMin=genVal;
                if(genVal>genMax) genMax=genVal;
            }
            double fitVal = fitGraphs[index_graph]->Eval(fitGraphs[index_graph]->GetX()[index_point]);
            double x=fitGraphs[index_graph]->GetX()[index_point];
            
            hCf_kstar->Fill(x, fitVal); // fill the histogram with the fit results for all the variations of the parameters, to be used for the evaluation of the systematic uncertainty on the final resu
            // if the fit type is the one with the linear term
            TF1* fitFuncBaseline= new TF1("fitFuncBaseline", "[0]+[1]*x*x/1000000. + [2]*x*x*x/1000000000. ", kfitMin, kfitMax);
            //TF1* fitFuncBaseline= new TF1("fitFuncBaseline", "[0]+[1]*x*x/1000000.  + [2]*x/1000.", kfitMin, kfitMax);
            fitFuncBaseline->SetParameter(0, fitBaselineFunc[index_graph]->GetParameter(0));
            fitFuncBaseline->SetParameter(1, fitBaselineFunc[index_graph]->GetParameter(1));
            fitFuncBaseline->SetParameter(2, fitBaselineFunc[index_graph]->GetParameter(2));
            //fitFuncBaseline->SetParameter(3, fitBaselineFunc[index_graph]->GetParameter(3));
            double baselineVal = fitFuncBaseline->Eval(x);
            if(baselineVal<baselineMinVal) baselineMinVal=baselineVal;
            if(baselineVal>baselineMaxVal) baselineMaxVal=baselineVal;
            delete fitFuncBaseline; // delete the fit function to avoid memory leak

        }
        genGraphMin->SetPoint(index_point, genModelGraphs[0]->GetX()[index_point], genMin);
        genGraphMax->SetPoint(index_point, genModelGraphs[0]->GetX()[index_point], genMax);
        baselineMin->SetPoint(index_point, fitGraphs[0]->GetX()[index_point], baselineMinVal);
        baselineMax->SetPoint(index_point, fitGraphs[0]->GetX()[index_point], baselineMaxVal);
        TH1D* yProjCF= hCf_kstar->ProjectionY(Form("yProjCF_%d", index_point), index_point+1, index_point+1);
        yProjCF->GetQuantiles(2, quantiles, probs);
        fitGraphMin->SetPoint(index_point, fitGraphs[0]->GetX()[index_point], quantiles[0]);
        fitGraphMax->SetPoint(index_point, fitGraphs[0]->GetX()[index_point], quantiles[1]);
        delete yProjCF; // delete the projection histogram to avoid
    }
    genGraphMin->SetName("genGraphMin");
    genGraphMax->SetName("genGraphMax");
    fitGraphMin->SetName("fitGraphMin");
    fitGraphMax->SetName("fitGraphMax");
    genGraphMin->Write("genGraphMin");
    genGraphMax->Write("genGraphMax");
    fitGraphMin->Write("fitGraphMin");
    fitGraphMax->Write("fitGraphMax");
    baselineMin->SetName("baselineMin");
    baselineMax->SetName("baselineMax");
    baselineMin->Write("baselineMin");
    baselineMax->Write("baselineMax");
    hCf_kstar->Write("hCf_kstar");
    outfile->Close();
    std::cout<<"finished saving graphs"<<std::endl;
}