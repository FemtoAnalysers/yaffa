//R__LOAD_LIBRARY($HOME/DLM/install/lib/libCATS.so) // good for sunrise setting. To change in case of different installation path

#include "/home/feriorob/PhD_analysis/analysis/constants.h" // contains the constants used in the analysis, such as the mass of the proton, lambdapar etc.
#include "/home/feriorob/PhD_analysis/analysis/utils.cpp"
//#include "/home/feriorob/PhD_analysis/analysis/utils_physics.cpp"

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

TH2F* hWF_SC = nullptr;
TH2F* hWF_MC = nullptr;


double GaussSource(double r, double r0)
{
    return 4.*TMath::Pi()*r*r*
           pow(4.*TMath::Pi()*r0*r0,-1.5)*
           exp(-(r*r)/(4.*r0*r0));
}

double GetCorrFunc(double k, double r0, TH2F* hWaveFun)
{
    double source_int=0.;
    double sum = 0.;
    int ik = hWaveFun->GetYaxis()->FindBin(k);
    for(int ir=1; ir<=hWaveFun->GetNbinsX(); ++ir)
    {
        double r  = hWaveFun->GetXaxis()->GetBinCenter(ir);
        double dr = hWaveFun->GetXaxis()->GetBinWidth(ir);
        double WF_val = hWaveFun->GetBinContent(ir,ik); // data are in GeV, WF in MeV
        double S = GaussSource(r,r0);
        sum += WF_val*S*dr;
        source_int += S*dr;
    }
    if(source_int>1 + 1e-4){
        std::cerr << "Warning: source integral is greater than 1: " << source_int << std::endl;
        sum /= source_int;
    }
    return sum;
}

void CF_pipi_SingleFit() // take a TH2 with WF value, integrate the source with 3 Gauss and fit them
{
    double LamParGen= 0.8;
    double fracResCut= 0.122;
    LamParGen = LamParGen*(1-fracResCut)*(1-fracResCut); // lambda flat value
    //double lParFlat=1. -LamParGen; // lambda flat value
    double lParFlat_singleGauss = 0.4;
    double lParFlat = 0.228; // for systematic evaluation
    double p_val=0.3679;
    double p_val_min=0.3551;
    double p_val_max=0.3828;

    TFile *InputFileWF_SC = new TFile("/home/feriorob/PhD_analysis/analysis/analysis_output_files/CoulomStudies/WFpipi_SC.root","read");
    TH2F *h_WF_sc = (TH2F*)InputFileWF_SC -> Get("hWFpipi"); //WF histogram for pp
    hWF_SC = h_WF_sc; // assign to global variable for use in the fit function

    TFile *InputFileWF_MC = new TFile("/home/feriorob/PhD_analysis/analysis/analysis_output_files/CoulomStudies/WFpipi_MC.root","read");
    TH2F *h_WF_mc = (TH2F*)InputFileWF_MC -> Get("hWFpipi"); 
    hWF_MC = h_WF_mc; // assign to global variable for use in the fit function

    TFile *InputDataFile = new TFile("/home/feriorob/data/CoulombStudies/pipiRun1data.root","read");

    TH1F *hData_qinv_SC = (TH1F*)InputDataFile -> Get("Table 1/Hist1D_y1");  //table 1= same charge
    TH1F *hSystErrData_SC = (TH1F*)InputDataFile -> Get("Table 1/Hist1D_y1_e2");
    TH1F *hStatErrData_SC = (TH1F*)InputDataFile -> Get("Table 1/Hist1D_y1_e1");

    TH1F *hData_qinv_MC = (TH1F*)InputDataFile -> Get("Table 2/Hist1D_y1");  //table 1= same charge
    TH1F *hSystErrData_MC = (TH1F*)InputDataFile -> Get("Table 2/Hist1D_y1_e2");
    TH1F *hStatErrData_MC = (TH1F*)InputDataFile -> Get("Table 2/Hist1D_y1_e1");

        //put errors on bin contents
    for (int idexBin=1; idexBin<=hData_qinv_SC->GetNbinsX(); idexBin++){ 
        double systErr = hSystErrData_SC->GetBinContent(idexBin);
        double statErr = hStatErrData_SC->GetBinContent(idexBin);
        double totalErr = sqrt(systErr*systErr + statErr*statErr);
        hData_qinv_SC->SetBinError(idexBin, totalErr);
    }

    for (int idexBin=1; idexBin<=hData_qinv_MC->GetNbinsX(); idexBin++){ 
        double systErr = hSystErrData_MC->GetBinContent(idexBin);
        double statErr = hStatErrData_MC->GetBinContent(idexBin);
        double totalErr = sqrt(systErr*systErr + statErr*statErr);
        hData_qinv_MC->SetBinError(idexBin, totalErr);
    }

    TH1F *hData_SC = new TH1F("hData","hData",hData_qinv_SC->GetNbinsX(),hData_qinv_SC->GetXaxis()->GetXmin()*1000./2.,hData_qinv_SC->GetXaxis()->GetXmax()*1000./2.);
    for (int idexBin=1; idexBin<=hData_qinv_SC->GetNbinsX(); idexBin++){
       hData_SC->SetBinContent(idexBin,hData_qinv_SC->GetBinContent(idexBin));
       hData_SC->SetBinError(idexBin,hData_qinv_SC->GetBinError(idexBin));
    }

    TH1F *hData_MC = new TH1F("hData","hData",hData_qinv_MC->GetNbinsX(),hData_qinv_MC->GetXaxis()->GetXmin()*1000./2.,hData_qinv_MC->GetXaxis()->GetXmax()*1000./2.);
    for (int idexBin=1; idexBin<=hData_qinv_MC->GetNbinsX(); idexBin++){
       hData_MC->SetBinContent(idexBin,hData_qinv_MC->GetBinContent(idexBin));
       hData_MC->SetBinError(idexBin,hData_qinv_MC->GetBinError(idexBin));
    }

    TFile *outputFile = new TFile("/home/feriorob/PhD_analysis/analysis/analysis_output_files/CoulomStudies/systStudyOutput/CF_pipi_SingleFit_SingleGauss_25.root", "RECREATE");

    TF1* fit_pipi_SC = new TF1(
        "fit_pipi_SC",
        [](double* x, double* par)
        {
            double k = x[0]; 
            double baseline = par[0];
            double p_val    = par[1];
            double r_p      = par[2];
            double r_r      = par[3];
            double lambdaFlat = par[4];
            double r_pr = sqrt((r_p*r_p + r_r*r_r)/2.); // effective radius for the mixed pairs

            double CF_p  = GetCorrFunc(k,r_p,hWF_SC);
            double CF_r  = GetCorrFunc(k,r_r,hWF_SC);
            double CF_pr = GetCorrFunc(k,r_pr,hWF_SC);
            //return baseline * ((1. - LFlatVal) * ( p_val*p_val*CF_p + 2.*p_val*(1-p_val)*CF_pr + (1-p_val)*(1-p_val)*CF_r ) + LFlatVal); // source is sum of 3 gauss, possible to separate the integral in 3 different CF
            return baseline * ((1. - lambdaFlat) * ( p_val*p_val*CF_p + 2.*p_val*(1-p_val)*CF_pr + (1-p_val)*(1-p_val)*CF_r ) + lambdaFlat); // source is sum of 3 gauss, possible to separate the integral in 3 different CF
        },
        0.,
        hData_SC->GetXaxis()->GetXmax(),
        5
    );

    std::cout<<"Lambda flat value: "<<lParFlat<<std::endl;

    fit_pipi_SC->SetParNames("Norm", "p", "r_p", "r_r", "lambdaFlat");
    fit_pipi_SC->SetParameters(
        1.0,    // Norm
        p_val,   // p
        6.501,     // r_p (fm)
        14.74,    // r_r (fm)
        lParFlat    // lambdaFlat
    );

    fit_pipi_SC->SetParLimits(0, 0.99, 1.01);  // Norm
    fit_pipi_SC->SetParLimits(1, p_val_min, p_val_max);  // p
    //fit_pipi_SC->FixParameter(1, p_val);  // p
    //fit_pipi_SC->FixParameter(0, 1.);  // p
    //fit_pipi_SC->FixParameter(2, 6.);      // r_p
    fit_pipi_SC->SetParLimits(2, 3., 9.);      // r_p
    //fit_pipi_SC->SetParLimits(2, 6, 7.);      // r_p
    fit_pipi_SC->SetParLimits(3, 11., 17.);    // r_r
    //fit_pipi_SC->FixParameter(3, 14.74);    // r_r
    fit_pipi_SC->FixParameter(4, lParFlat);  // lambdaFlat

    // -------------------- the actual fit --------------------
    TFitResultPtr fitResult = hData_SC->Fit(fit_pipi_SC, "MRNS", "", 0., hData_SC->GetXaxis()->GetXmax());

    std::cout << "Fit status: " << fitResult->Status() << " (0 = OK)" << std::endl;

    double chi2 = fit_pipi_SC->GetChisquare();
    int    ndf  = fit_pipi_SC->GetNDF();
    std::cout << "Chi2 = " << chi2 << std::endl;
    std::cout << "NDF  = " << ndf  << std::endl;
    std::cout << "Chi2/NDF = " << (ndf > 0 ? chi2/ndf : -1.) << std::endl;

    double norm_SC   = fit_pipi_SC->GetParameter(0);
    double p_SC      = fit_pipi_SC->GetParameter(1);
    double rp_SC     = fit_pipi_SC->GetParameter(2);
    double rr_SC     = fit_pipi_SC->GetParameter(3);
    double lambda_SC = fit_pipi_SC->GetParameter(4);
    double rpr_SC    = sqrt((rp_SC*rp_SC + rr_SC*rr_SC)/2.);

    std::cout << "Norm       = " << norm_SC   << " +/- " << fit_pipi_SC->GetParError(0) << std::endl;
    std::cout << "p_val      = " << p_SC      << " +/- " << fit_pipi_SC->GetParError(1) << std::endl;
    std::cout << "r_p        = " << rp_SC     << " +/- " << fit_pipi_SC->GetParError(2) << std::endl;
    std::cout << "r_r        = " << rr_SC     << " +/- " << fit_pipi_SC->GetParError(3) << std::endl;
    std::cout << "r_pr       = " << rpr_SC    << std::endl;
    std::cout << "lambdaFlat = " << lambda_SC << " +/- " << fit_pipi_SC->GetParError(4) << std::endl;

    // correlation matrix (useful to spot degeneracies, e.g. p_val vs r_p/r_r)
    TMatrixDSym cov = fitResult->GetCovarianceMatrix();
    for (int i = 0; i < 5; ++i) {
        for (int j = i+1; j < 5; ++j) {
            double corr = cov(i,j) / sqrt(cov(i,i)*cov(j,j));
            std::cout << fit_pipi_SC->GetParName(i) << " vs " << fit_pipi_SC->GetParName(j)
                       << " : " << corr << std::endl;
        }
    }

        TF1* fit_pipi_MC = new TF1(
        "fit_pipi_MC",
        [](double* x, double* par)
        {
            double k = x[0]; 
            double baseline = par[0];
            double p_val    = par[1];
            double r_p      = par[2];
            double r_r      = par[3];
            double lambdaFlat = par[4];
            double r_pr = sqrt((r_p*r_p + r_r*r_r)/2.); // effective radius for the mixed pairs

            double CF_p  = GetCorrFunc(k,r_p,hWF_SC);
            double CF_r  = GetCorrFunc(k,r_r,hWF_SC);
            double CF_pr = GetCorrFunc(k,r_pr,hWF_SC);
            //return baseline * ((1. - LFlatVal) * ( p_val*p_val*CF_p + 2.*p_val*(1-p_val)*CF_pr + (1-p_val)*(1-p_val)*CF_r ) + LFlatVal); // source is sum of 3 gauss, possible to separate the integral in 3 different CF
            return baseline * ((1. - lambdaFlat) * ( p_val*p_val*CF_p + 2.*p_val*(1-p_val)*CF_pr + (1-p_val)*(1-p_val)*CF_r ) + lambdaFlat); // source is sum of 3 gauss, possible to separate the integral in 3 different CF
        },
        0.,
        hData_MC->GetXaxis()->GetXmax(),
        5
    );

    fit_pipi_MC->SetParLimits(0, 0.99, 1.01);  // Norm
    fit_pipi_MC->SetParLimits(1, p_val_min, p_val_max);  // p
    //fit_pipi_MC->FixParameter(1, p_val);  // p
    //fit_pipi_MC->FixParameter(0, 1.);  // p
    //fit_pipi_MC->FixParameter(2, 6.);      // r_p
    fit_pipi_MC->SetParLimits(2, 3., 9.);      // r_p
    //fit_pipi_MC->SetParLimits(2, 6, 7.);      // r_p
    fit_pipi_MC->SetParLimits(3, 11., 17.);    // r_r
    //fit_pipi_MC->FixParameter(3, 14.74);    // r_r
    fit_pipi_MC->FixParameter(4, lParFlat);  // lambdaFlat

    // -------------------- the actual fit --------------------
    TFitResultPtr fitResultMC = hData_MC->Fit(fit_pipi_MC, "MRNS", "", 0., hData_MC->GetXaxis()->GetXmax());

    std::cout << "Fit status: " << fitResultMC->Status() << " (0 = OK)" << std::endl;

    double chi2MC = fit_pipi_MC->GetChisquare();
    int    ndfMC  = fit_pipi_MC->GetNDF();
    std::cout << "Chi2 = " << chi2MC << std::endl;
    std::cout << "NDF  = " << ndfMC  << std::endl;
    std::cout << "Chi2/NDF = " << (ndfMC > 0 ? chi2MC/ndfMC : -1.) << std::endl;

    double norm_MC   = fit_pipi_MC->GetParameter(0);
    double p_MC      = fit_pipi_MC->GetParameter(1);
    double rp_MC     = fit_pipi_MC->GetParameter(2);
    double rr_MC     = fit_pipi_MC->GetParameter(3);
    double lambda_MC = fit_pipi_MC->GetParameter(4);
    double rpr_MC    = sqrt((rp_MC*rp_MC + rr_MC*rr_MC)/2.);

    std::cout << "Norm       = " << norm_MC   << " +/- " << fit_pipi_MC->GetParError(0) << std::endl;
    std::cout << "p_val      = " << p_MC      << " +/- " << fit_pipi_MC->GetParError(1) << std::endl;
    std::cout << "r_p        = " << rp_MC     << " +/- " << fit_pipi_MC->GetParError(2) << std::endl;
    std::cout << "r_r        = " << rr_MC     << " +/- " << fit_pipi_MC->GetParError(3) << std::endl;
    std::cout << "r_pr       = " << rpr_MC    << std::endl;
    std::cout << "lambdaFlat = " << lambda_MC << " +/- " << fit_pipi_MC->GetParError(4) << std::endl;

    // --------------------single gaussian source function for comparison--------------------
    TF1* fSingleGauss_SC = new TF1("SC_fSingleGauss",
        [=](double* x, double* par){
            double k = x[0];
            double baseline = par[0];
            double r_single = par[1];
            double lambdaFlat = par[2];

            double CF_single = GetCorrFunc(k, r_single, hWF_SC);
            return baseline * ((1. - lambdaFlat) * CF_single + lambdaFlat);
        }, 0., hData_SC->GetXaxis()->GetXmax(), 3);

        fSingleGauss_SC->SetParNames("Baseline", "r_single", "lambdaFlat");
        fSingleGauss_SC->SetParameters(
            1.0,    // Baseline
            8.,  // r_single (fm)
            lParFlat_singleGauss // lambdaFlat
        );

        fSingleGauss_SC->SetParLimits(0, 0.99, 1.01);  // Baseline
        fSingleGauss_SC->SetParLimits(1, 3., 16.);     // r_single
        fSingleGauss_SC->FixParameter(2, lParFlat_singleGauss);  // lambdaFlat

        TFitResultPtr fitResultSingle = hData_SC->Fit(fSingleGauss_SC, "MRNS", "", 0., hData_SC->GetXaxis()->GetXmax());
        std::cout << "Fit status (single gauss): " << fitResultSingle->Status() << " (0 = OK)" << std::endl;

        std::cout << "Single Gaussian Fit Results:" << std::endl;
        std::cout << "Baseline = " << fSingleGauss_SC->GetParameter(0) << " +/- " << fSingleGauss_SC->GetParError(0) << std::endl;
        std::cout << "r_single = " << fSingleGauss_SC->GetParameter(1) << " +/- " << fSingleGauss_SC->GetParError(1) << std::endl;
        std::cout << "lambdaFlat = " << fSingleGauss_SC->GetParameter(2) << " +/- " << fSingleGauss_SC->GetParError(2) << std::endl;
        std::cout << "Chi2/NDF (single gauss) = " << fSingleGauss_SC->GetChisquare() / fSingleGauss_SC->GetNDF() << std::endl;

        double r_single= fSingleGauss_SC->GetParameter(1);
        double flat_single= fSingleGauss_SC->GetParameter(2);

    // -------------------- decomposed components (weighted) --------------------
    TF1* fPrim_SC = new TF1("SC_fPrim",
        [=](double* x, double*) {
            return norm_SC*(1.-lambda_SC) * p_SC*p_SC * GetCorrFunc(x[0], rp_SC, hWF_SC);
        }, 0., hData_SC->GetXaxis()->GetXmax(), 0);

    TF1* fMix_SC = new TF1("SC_fMix",
        [=](double* x, double*) {
            return norm_SC*(1.-lambda_SC) * 2.*p_SC*(1.-p_SC) * GetCorrFunc(x[0], rpr_SC, hWF_SC);
        }, 0., hData_SC->GetXaxis()->GetXmax(), 0);

    TF1* fRes_SC = new TF1("SC_fRes",
        [=](double* x, double*) {
            return norm_SC*(1.-lambda_SC) * (1.-p_SC)*(1.-p_SC) * GetCorrFunc(x[0], rr_SC, hWF_SC);
        }, 0., hData_SC->GetXaxis()->GetXmax(), 0);

    TF1* fFlat_SC = new TF1("SC_fFlat", Form("%f", norm_SC*lambda_SC), 0., hData_SC->GetXaxis()->GetXmax());

    // -------------------- decomposed components (no weight, for shape comparison) --------------------
    TF1* fPrim_SC_noWeight = new TF1("SC_fPrim_noWeight",
        [=](double* x, double*) { return norm_SC*(1.-lambda_SC) * GetCorrFunc(x[0], rp_SC, hWF_SC); },
        0., hData_SC->GetXaxis()->GetXmax(), 0);

    TF1* fMix_SC_noWeight = new TF1("SC_fMix_noWeight",
        [=](double* x, double*) { return norm_SC*(1.-lambda_SC) * GetCorrFunc(x[0], rpr_SC, hWF_SC); },
        0., hData_SC->GetXaxis()->GetXmax(), 0);

    TF1* fRes_SC_noWeight = new TF1("SC_fRes_noWeight",
        [=](double* x, double*) { return norm_SC*(1.-lambda_SC) * GetCorrFunc(x[0], rr_SC, hWF_SC); },
        0., hData_SC->GetXaxis()->GetXmax(), 0);

    // -------------------- source functions S(r) --------------------
    TF1* fSp_SC = new TF1("SC_fSp",
        [=](double* x, double*) { return p_SC*p_SC*GaussSource(x[0], rp_SC); }, 0., 100., 0);

    TF1* fSpr_SC = new TF1("SC_fSpr",
        [=](double* x, double*) { return 2.*p_SC*(1.-p_SC)*GaussSource(x[0], rpr_SC); }, 0., 100., 0);

    TF1* fSr_SC = new TF1("SC_fSr",
        [=](double* x, double*) { return (1.-p_SC)*(1.-p_SC)*GaussSource(x[0], rr_SC); }, 0., 100., 0);

    TF1* fSp_SC_noWeight = new TF1("SC_fSp_noWeight",
        [=](double* x, double*) { return GaussSource(x[0], rp_SC); }, 0., 100., 0);

    TF1* fSpr_SC_noWeight = new TF1("SC_fSpr_noWeight",
        [=](double* x, double*) { return GaussSource(x[0], rpr_SC); }, 0., 100., 0);

    TF1* fSr_SC_noWeight = new TF1("SC_fSr_noWeight",
        [=](double* x, double*) { return GaussSource(x[0], rr_SC); }, 0., 100., 0);

    TF1* TotalSource_SC = new TF1("SC_TotalSource",
        [=](double* x, double*) {
            return p_SC*p_SC*GaussSource(x[0], rp_SC)
                 + 2.*p_SC*(1.-p_SC)*GaussSource(x[0], rpr_SC)
                 + (1.-p_SC)*(1.-p_SC)*GaussSource(x[0], rr_SC);
        }, 0., 100., 0);

    TF1* SingleSource_SC = new TF1("SC_SingleSource",
        [=](double* x, double*) { return GaussSource(x[0], r_single); }, 0., 100., 0);

    // -------------------- canvas: fit + data (weighted components) --------------------
    TCanvas *cFitData = new TCanvas("cFitData_SC","cFitData_SC",800,600);
    cFitData->cd();
    gPad->DrawFrame(0., 0., hData_SC->GetXaxis()->GetXmax(), 5., "Fit of CF for #pi^{+}#pi^{+} #oplus #pi^{-}#pi^{-}; k* (MeV/c); C(k*)");

    SetHistInfo(hData_SC, kBlack, 20, 1.0, "Data #pi^{+}#pi^{+} #oplus #pi^{-}#pi^{-}", "k* (MeV/c)", "C(k*)", 1);
    hData_SC->Draw("E same");

    SetTF1info(fit_pipi_SC, kRed, 20, 1.0, "Fit with CF", "r (fm)", "S(r)", 1);
    fit_pipi_SC->Draw("L same");
    SetTF1info(fPrim_SC, kViolet-1, 20, 1.0, "Primary-Primary", "r (fm)", "S(r)", 2);
    fPrim_SC->Draw("L same");
    SetTF1info(fMix_SC, kViolet+9, 20, 1.0, "Primary-Resonance", "r (fm)", "S(r)", 4);
    fMix_SC->Draw("L same");
    SetTF1info(fRes_SC, kViolet-8, 20, 1.0, "Resonance-Resonance", "r (fm)", "S(r)", 9);
    fRes_SC->Draw("L same");
    SetTF1info(fFlat_SC, kMagenta, 20, 1.0, "Flat", "r (fm)", "S(r)", 3);
    fFlat_SC->Draw("L same");

    TLine *line_SC = new TLine(0., norm_SC, hData_SC->GetXaxis()->GetXmax(), norm_SC);
    line_SC->SetLineColor(kGray);
    line_SC->SetLineStyle(2);
    line_SC->Draw("Same");

    TLegend *legend_SC = new TLegend(0.6, 0.65, 0.9, 0.9);
    legend_SC->AddEntry(hData_SC, "Data #pi^{+}#pi^{+} #oplus #pi^{-}#pi^{-}", "lep");
    legend_SC->AddEntry(fit_pipi_SC, "Fit with CF", "l");
    legend_SC->AddEntry(fPrim_SC, "Primary-Primary", "l");
    legend_SC->AddEntry(fMix_SC, "Primary-Resonance", "l");
    legend_SC->AddEntry(fRes_SC, "Resonance-Resonance", "l");
    legend_SC->AddEntry(fFlat_SC, "Flat", "l");
    legend_SC->AddEntry(line_SC, "Baseline", "l");
    legend_SC->SetBorderSize(0);
    legend_SC->Draw();
    cFitData->Write();


    TCanvas *cFitData_MC = new TCanvas("cFitData_MC","cFitData_MC",800,600);
    cFitData_MC->cd();
    gPad->DrawFrame(0., 0., hData_MC->GetXaxis()->GetXmax(), 5., "Fit of CF for #pi^{-}#pi^{+} #oplus #pi^{+}#pi^{-} (MC); k* (MeV/c); C(k*)");

    SetHistInfo(hData_MC, kBlack, 20, 1.0, "Data #pi^{-}#pi^{+} #oplus #pi^{+}#pi^{-}", "k* (MeV/c)", "C(k*)", 1);
    hData_MC->Draw("E same");

    SetTF1info(fit_pipi_MC, kOrange, 20, 1.0, "Fit with CF", "r (fm)", "S(r)", 1);
    fit_pipi_MC->Draw("L same");

    TLegend *legend_MC = new TLegend(0.6, 0.65, 0.9, 0.9);
    legend_MC->AddEntry(hData_MC, "Data #pi^{-}#pi^{+} #oplus #pi^{+}#pi^{-}", "lep");
    legend_MC->AddEntry(fit_pipi_MC, "Fit with CF", "l");
    legend_MC->SetBorderSize(0);
    legend_MC->Draw();
    cFitData_MC->Write();


    // -------------------- canvas: fit + data (no-weight components) --------------------
    TCanvas *cFitData_noWeight = new TCanvas("cFitData_SC_noWeight","cFitData_SC_noWeight",800,600);
    cFitData_noWeight->cd();
    gPad->DrawFrame(0., 0., hData_SC->GetXaxis()->GetXmax(), 5., "Fit of CF (no weight) for #pi^{+}#pi^{+} #oplus #pi^{-}#pi^{-}; k* (MeV/c); C(k*)");

    SetHistInfo(hData_SC, kBlack, 20, 1.0, "Data #pi^{+}#pi^{+} #oplus #pi^{-}#pi^{-}", "k* (MeV/c)", "C(k*)", 1);
    hData_SC->Draw("E same");
    SetTF1info(fit_pipi_SC, kRed, 20, 1.0, "Fit with CF", "r (fm)", "S(r)", 1);
    fit_pipi_SC->Draw("L same");
    SetTF1info(fPrim_SC_noWeight, kViolet-1, 20, 1.0, "Primary-Primary", "r (fm)", "S(r)", 2);
    fPrim_SC_noWeight->Draw("L same");
    SetTF1info(fMix_SC_noWeight, kViolet+9, 20, 1.0, "Primary-Resonance", "r (fm)", "S(r)", 4);
    fMix_SC_noWeight->Draw("L same");
    SetTF1info(fRes_SC_noWeight, kViolet-8, 20, 1.0, "Resonance-Resonance", "r (fm)", "S(r)", 9);
    fRes_SC_noWeight->Draw("L same");
    SetTF1info(fFlat_SC, kMagenta, 20, 1.0, "Flat", "r (fm)", "S(r)", 3);
    fFlat_SC->Draw("L same");

    TLegend *legend_noWeight_SC = new TLegend(0.6, 0.65, 0.9, 0.9);
    legend_noWeight_SC->AddEntry(hData_SC, "Data #pi^{+}#pi^{+} #oplus #pi^{-}#pi^{-}", "lep");
    legend_noWeight_SC->AddEntry(fit_pipi_SC, "Fit with CF", "l");
    legend_noWeight_SC->AddEntry(fPrim_SC_noWeight, "Primary-Primary", "l");
    legend_noWeight_SC->AddEntry(fMix_SC_noWeight, "Primary-Resonance", "l");
    legend_noWeight_SC->AddEntry(fRes_SC_noWeight, "Resonance-Resonance", "l");
    legend_noWeight_SC->AddEntry(fFlat_SC, "Flat", "l");
    legend_noWeight_SC->SetBorderSize(0);
    legend_noWeight_SC->Draw();
    cFitData_noWeight->Write();

    // -------------------- canvas: source decomposition --------------------
    TCanvas *cCanvasSource_SC = new TCanvas("cCanvasSource_SC","cCanvasSource_SC",800,600);
    cCanvasSource_SC->cd();
    cCanvasSource_SC->DrawFrame(0., 0., 100., 0.25, "Source functions for #pi^{+}#pi^{+} #oplus #pi^{-}#pi^{-}; r (fm); S(r)");
    SetTF1info(fSp_SC, kViolet-1, 20, 1.0, "Primary-Primary", "r (fm)", "S(r)", 2);
    fSp_SC->Draw("L same");
    SetTF1info(fSpr_SC, kViolet+9, 20, 1.0, "Primary-Resonance", "r (fm)", "S(r)", 4);
    fSpr_SC->Draw("L same");
    SetTF1info(fSr_SC, kViolet-8, 20, 1.0, "Resonance-Resonance", "r (fm)", "S(r)", 9);
    fSr_SC->Draw("L same");
    SetTF1info(TotalSource_SC, kOrange+7, 20, 1.0, "Total Source", "r (fm)", "S(r)", 1);
    TotalSource_SC->Draw("L same");
    SetTF1info(SingleSource_SC, kGreen+2, 20, 1.0, "Single Gaussian Source", "r (fm)", "S(r)", 1);
    SingleSource_SC->Draw("L same");

    TLegend *legendSource = new TLegend(0.6, 0.7, 0.9, 0.9);
    legendSource->AddEntry(TotalSource_SC, "Total Source", "l");
    legendSource->AddEntry(fSp_SC, "Primary-Primary", "l");
    legendSource->AddEntry(fSpr_SC, "Primary-Resonance", "l");
    legendSource->AddEntry(fSr_SC, "Resonance-Resonance", "l");
    legendSource->AddEntry(SingleSource_SC, "Single Gaussian Source", "l");
    legendSource->SetBorderSize(0);
    legendSource->Draw();
    cCanvasSource_SC->Write();

    TCanvas *cSingleGaussSourceFit = new TCanvas("cSingleGaussSourceFit_SC","cSingleGaussSourceFit_SC",800,600);
    cSingleGaussSourceFit->cd();
    gPad->DrawFrame(0., 0., 100., 0.25, "Single Gaussian Source Fit for #pi^{+}#pi^{+} #oplus #pi^{-}#pi^{-}; r (fm); S(r)");
    SetTF1info(fSingleGauss_SC, kRed, 20, 1.0, "Single Gaussian Fit", "r (fm)", "S(r)", 1);
    hData_SC->Draw("E same");
    fSingleGauss_SC->Draw("L same");
    TLine *line_single = new TLine(0., flat_single, 35., flat_single);
    line_single->SetLineColor(kMagenta);
    line_single->SetLineStyle(2);
    line_single->Draw("Same");
    
    TLine *line_single_baseline = new TLine(0., fSingleGauss_SC->GetParameter(0), 35., fSingleGauss_SC->GetParameter(0));
    line_single_baseline->SetLineColor(kGray);
    line_single_baseline->SetLineStyle(2);
    line_single_baseline->Draw("Same");
    TLegend *legendSingleGauss = new TLegend(0.6, 0.7, 0.9, 0.9);
    legendSingleGauss->AddEntry(fSingleGauss_SC, "Single Gaussian Fit", "l");
    legendSingleGauss->AddEntry(line_single, "Flat", "l");
    legendSingleGauss->AddEntry(line_single_baseline, "Baseline", "l");
    legendSingleGauss->SetBorderSize(0);
    legendSingleGauss->Draw();
    cSingleGaussSourceFit->Write();


    std::cout << "Done. Saving output file..." << std::endl;
    outputFile->Close();
}