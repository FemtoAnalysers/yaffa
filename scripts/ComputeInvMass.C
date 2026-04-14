#include <math.h>

// #include "Fitting.h"
#include "TF1.h"
#include "TROOT.h"
#include "TSystem.h"
// #include "src/common.h"

const double mL = 1115.683;
const double mK = 493.677;
const double RedMass = (mL * mK) / (mL + mK);
const double EnThrLK = mL + mK;
const double MassL = 1115.683;   // Mass of Daughter 1 in MeV
const double MassKmin = 493.68;  // Mass of Daughter 2 in MeV
const double MassXi = 1321.71;
const double MassPi = 139.57;

void SetStyle(bool graypalette, bool title) {
    const int NCont = 255;
    gStyle->Reset("Plain");
    gStyle->SetNumberContours(NCont);
    gStyle->SetOptTitle(title);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    gStyle->SetFuncColor(kGreen);
    gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.045, "x");
    gStyle->SetLabelSize(0.045, "y");
    gStyle->SetLabelOffset(0.01, "y");
    gStyle->SetLabelOffset(0.01, "x");
    gStyle->SetLabelColor(kBlack, "xyz");
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetTitleOffset(1.25, "y");
    gStyle->SetTitleOffset(1.2, "x");
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendBorderSize(0);
}
//------------------------------------------
double fit_FemtoSillWeightedLK(double *x, double *par) {
    double t = *x;
    // par: using Eq(77)  s-1 < s < s_2  from F. Giacosa Eur. Phys. J. A (2021) s->E
    //      [0] width 1 decay into πΞ
    //      [1] width 2 decay into ΛK
    //      [2] mass of Ξ(1620)
    //      [3] weight

    double RedMass = MassL * MassKmin / (MassL + MassKmin);
    double enThr1 = MassXi + MassPi;
    double enThr2 = MassL + MassKmin;

    double arg1 = 0.;
    double numerator = 0.;
    double arg3 = 0.;
    double denominator = 0.;
    double denominator2 = 0.;
    double Sill1 = 0.;
    double Sill2 = 0.;
    double Sill = 0.;

    if (t >= enThr2) {
        arg1 = (2. * t) / TMath::Pi();
        numerator = par[1] * TMath::Sqrt(t * t - enThr2 * enThr2);  // decay into ΛK-
        denominator2 =
            pow((t * t - par[2] * par[2]), 2) + pow(numerator + par[0] * TMath::Sqrt(t * t - enThr1 * enThr1), 2);
        Sill1 = arg1 * numerator / denominator2;
    }

    Sill = par[3] * (Sill1 + Sill2);

    return Sill;
}

double fit_WeightedBreitWigner(double *x, double *par) {
    double t = *x;  // energy

    // par:
    //      [0] mass
    //      [1] width
    //      [2] weight
    return (par[2] * TMath::BreitWigner(t, par[0], par[1]));
}

double fit_WeightedGauss_Exclude(double *x, double *par) {
    double t = *x;  // energy
    double xmin = par[3];
    double xmax = par[4];
    if (t < xmin && t > xmax) {
        TF1::RejectPoint();
        return (double)0;
    }

    // par:
    //      [0] mass
    //      [1] width
    //      [2] weight
    return (par[2] * TMath::Gaus(t, par[0], par[1], kTRUE));
}

double fit_WeightedGauss(double *x, double *par) {
    double t = *x;  // energy
    // par:
    //      [0] mass
    //      [1] width
    //      [2] weight
    return (par[2] * TMath::Gaus(t, par[0], par[1], kTRUE));
}

double fit_axis(double *x) {
    double t = *x;  // energy

    return sqrt(t * t + mL * mL) + sqrt(t * t + mK * mK);
}

//------------------------------------------
void ComputeInvMass(std::string inDataFileName, std::string inMCFileName, std::string oFileName, int REBINMC = 1) {

    SetStyle(false, false);

    TFile *OutputFile = new TFile(oFileName.data(), "recreate");

    /// Accessing resonances properties
    /// Reading the resonance properties from the fit
    TFile *InputTree;

    TString fName_Output;

    TString fName_OutputFigs;
    TString fName_Output_AncestorsFigs;

    fName_Output = homeFolder + OutputFolderLednickySillGaussERE;
    fName_OutputFigs = homeFolder + OutputFolderLednickySillGaussERE;

    std::cout << "Output folder is:" << fName_Output << std::endl;

    int TypeOfFit = 4;
    bool WithOmega = false;
    bool IsEffRangeZero = false;
    int UsePbPb = 0;
    InputTree = TFile::Open(
        fName_Output +
        TString::Format("OutputTree_rebinMC2rebinDATA1_Fit%i_Back0_Omega%i_EffRangeZero%i_UsePbPb%i_Resonance.root",
                        TypeOfFit, int(WithOmega), int(IsEffRangeZero), UsePbPb));

    // TH1F *hMassOmegafFullFemto = (TH1F *)(InputTree->Get("hMassOmegafFullFemto"));
    // TH1F *hMassXi1fFullFemto = (TH1F *)(InputTree->Get("hMassXi1fFullFemto"));
    // TH1F *hMassXi2fFullFemto = (TH1F *)(InputTree->Get("hMassXi2fFullFemto"));
    // TH1F *hGammaOmegafFullFemto = (TH1F *)(InputTree->Get("hGammaOmegafFullFemto"));
    // TH1F *hGammaXi1fFullFemto = (TH1F *)(InputTree->Get("hGammaXi1fFullFemto"));
    // TH1F *hGammaXi2fFullFemto = (TH1F *)(InputTree->Get("hGammaXi2fFullFemto"));
    // TH1F *hHeightOmegafFullFemto = (TH1F *)(InputTree->Get("hHeightOmegafFullFemto"));
    // TH1F *hHeightXi1fFullFemto = (TH1F *)(InputTree->Get("hHeightXi1fFullFemto"));
    // TH1F *hHeightXi2fFullFemto = (TH1F *)(InputTree->Get("hHeightXi2fFullFemto"));

    // TH1F *hMassOmegafFullFemto_stat = (TH1F *)(InputTree->Get("hMassOmegafFullFemto_stat"));
    // TH1F *hMassXi1fFullFemto_stat = (TH1F *)(InputTree->Get("hMassXi1fFullFemto_stat"));
    // TH1F *hMassXi2fFullFemto_stat = (TH1F *)(InputTree->Get("hMassXi2fFullFemto_stat"));
    // TH1F *hGammaOmegafFullFemto_stat = (TH1F *)(InputTree->Get("hGammaOmegafFullFemto_stat"));
    // TH1F *hGammaXi1fFullFemto_stat = (TH1F *)(InputTree->Get("hGammaXi1fFullFemto_stat"));
    // TH1F *hGammaXi2fFullFemto_stat = (TH1F *)(InputTree->Get("hGammaXi2fFullFemto_stat"));
    // TH1F *hHeightOmegafFullFemto_stat = (TH1F *)(InputTree->Get("hHeightOmegafFullFemto_stat"));
    // TH1F *hHeightXi1fFullFemto_stat = (TH1F *)(InputTree->Get("hHeightXi1fFullFemto_stat"));
    // TH1F *hHeightXi2fFullFemto_stat = (TH1F *)(InputTree->Get("hHeightXi2fFullFemto_stat"));

    // TH1F *hMassXi1620fFullFemto = (TH1F *)(InputTree->Get("hMassXi1620fFullFemto"));
    // TH1F *hGammaXi1620XiPifFullFemto = (TH1F *)(InputTree->Get("hGammaXi1620XiPifFullFemto"));
    // TH1F *hGammaXi1620LKfFullFemto = (TH1F *)(InputTree->Get("hGammaXi1620LKfFullFemto"));
    // TH1F *hHeightXi1620fFullFemto = (TH1F *)(InputTree->Get("hHeightXi1620fFullFemto"));

    // TH1F *hMassXi1620fFullFemto_stat = (TH1F *)(InputTree->Get("hMassXi1620fFullFemto_stat"));
    // TH1F *hGammaXi1620XiPifFullFemto_stat = (TH1F *)(InputTree->Get("hGammaXi1620XiPifFullFemto_stat"));
    // TH1F *hGammaXi1620LKfFullFemto_stat = (TH1F *)(InputTree->Get("hGammaXi1620LKfFullFemto_stat"));
    // TH1F *hHeightXi1620fFullFemto_stat = (TH1F *)(InputTree->Get("hHeightXi1620fFullFemto_stat"));
    // TH1F *hWeightfFullFemto_stat = (TH1F *)(InputTree->Get("hWeightfFullFemto_stat"));

    // double MassOmega_mean = hMassOmegafFullFemto_stat->GetMean();
    // double MassOmegaEnergy_mean =
    //     sqrt(MassOmega_mean * MassOmega_mean + mL * mL) + sqrt(MassOmega_mean * MassOmega_mean + mK * mK);

    // double MassOmega_sig_stat = hMassOmegafFullFemto_stat->GetStdDev();
    // double MassOmega_sig_tot = hMassOmegafFullFemto->GetStdDev();
    // double MassOmega_sig_syst =
    //     TMath::Sqrt(MassOmega_sig_tot * MassOmega_sig_tot - MassOmega_sig_stat * MassOmega_sig_stat);

    // double GammaOmega_mean = hGammaOmegafFullFemto_stat->GetMean();
    // double GammaOmega_sig_stat = hGammaOmegafFullFemto_stat->GetStdDev();
    // double GammaOmega_sig_tot = hGammaOmegafFullFemto->GetStdDev();
    // double GammaOmega_sig_syst =
    //     TMath::Sqrt(GammaOmega_sig_tot * GammaOmega_sig_tot - GammaOmega_sig_stat * GammaOmega_sig_stat);

    // double HeightOmega_mean = hHeightOmegafFullFemto_stat->GetMean();
    // double HeightOmega_sig_stat = hHeightOmegafFullFemto_stat->GetStdDev();
    // double HeightOmega_sig_tot = hHeightOmegafFullFemto->GetStdDev();
    // double HeightOmega_sig_syst =
    //     TMath::Sqrt(HeightOmega_sig_tot * HeightOmega_sig_tot - HeightOmega_sig_stat * HeightOmega_sig_stat);

    // /// Ξ(1690)
    // double MassXi1_mean = hMassXi1fFullFemto_stat->GetMean();
    // double MassXi1Energy_mean =
    //     sqrt(MassXi1_mean * MassXi1_mean + mL * mL) + sqrt(MassXi1_mean * MassXi1_mean + mK * mK);
    // double MassXi1_sig_stat = hMassXi1fFullFemto_stat->GetStdDev();
    // double MassXi1_sig_tot = hMassXi1fFullFemto->GetStdDev();
    // double MassXi1_sig_syst = TMath::Sqrt(MassXi1_sig_tot * MassXi1_sig_tot - MassXi1_sig_stat * MassXi1_sig_stat);

    // double GammaXi1_mean = hGammaXi1fFullFemto_stat->GetMean();
    // double GammaXi1_sig_stat = hGammaXi1fFullFemto_stat->GetStdDev();
    // double GammaXi1_sig_tot = hGammaXi1fFullFemto->GetStdDev();
    // double GammaXi1_sig_syst = TMath::Sqrt(GammaXi1_sig_tot * GammaXi1_sig_tot - GammaXi1_sig_stat * GammaXi1_sig_stat);

    // double HeightXi1_mean = hHeightXi1fFullFemto_stat->GetMean();
    // double HeightXi1_sig_stat = hHeightXi1fFullFemto_stat->GetStdDev();
    // double HeightXi1_sig_tot = hHeightXi1fFullFemto->GetStdDev();
    // double HeightXi1_sig_syst =
    //     TMath::Sqrt(HeightXi1_sig_tot * HeightXi1_sig_tot - HeightXi1_sig_stat * HeightXi1_sig_stat);

    // /// Ξ(1820)
    // double MassXi2_mean = hMassXi2fFullFemto_stat->GetMean();
    // double MassXi2Energy_mean =
    //     sqrt(MassXi2_mean * MassXi2_mean + mL * mL) + sqrt(MassXi2_mean * MassXi2_mean + mK * mK);
    // double MassXi2_sig_stat = hMassXi2fFullFemto_stat->GetStdDev();
    // double MassXi2_sig_tot = hMassXi2fFullFemto->GetStdDev();
    // double MassXi2_sig_syst = TMath::Sqrt(MassXi2_sig_tot * MassXi2_sig_tot - MassXi2_sig_stat * MassXi2_sig_stat);

    // double GammaXi2_mean = hGammaXi2fFullFemto_stat->GetMean();
    // double GammaXi2_sig_stat = hGammaXi2fFullFemto_stat->GetStdDev();
    // double GammaXi2_sig_tot = hGammaXi2fFullFemto->GetStdDev();
    // double GammaXi2_sig_syst = TMath::Sqrt(GammaXi2_sig_tot * GammaXi2_sig_tot - GammaXi2_sig_stat * GammaXi2_sig_stat);

    // double HeightXi2_mean = hHeightXi2fFullFemto_stat->GetMean();
    // double HeightXi2_sig_stat = hHeightXi2fFullFemto_stat->GetStdDev();
    // double HeightXi2_sig_tot = hHeightXi2fFullFemto->GetStdDev();
    // double HeightXi2_sig_syst =
    //     TMath::Sqrt(HeightXi2_sig_tot * HeightXi2_sig_tot - HeightXi2_sig_stat * HeightXi2_sig_stat);

    // /// Ξ(1820)
    // double MassXi1620_mean = hMassXi1620fFullFemto_stat->GetMean();

    // /// Data
    // TString fName_Data = homeFolder + foldernameData + TString::Format("CFOutput_ALK_0.root");
    // TFile *FileData;
    // FileData = new TFile(fName_Data, "read");
    // TList *List1 = (TList *)FileData->Get("PairDist");
    // TList *List1A = (TList *)FileData->Get("AntiPairDist");
    // TList *List2 = (TList *)List1->FindObject("Pair");
    // TList *List2A = (TList *)List1A->FindObject("Pair");
    // auto hSEPair_GeV = (TH1F *)List2->FindObject("SEDist_Particle1_Particle2_clone");
    // auto hSEPair = (TH1F *)TransformToMeV1D(hSEPair_GeV);
    // auto hSEAntiPair_GeV = (TH1F *)List2A->FindObject("SEDist_Particle0_Particle3_clone");
    // auto hSEAntiPair = (TH1F *)TransformToMeV1D(hSEAntiPair_GeV);
    // auto hSE = (TH1F *)hSEPair->Clone("hSEData");
    // hSE->Add(hSEAntiPair);

    // auto hMEPair_GeV = (TH1F *)List2->FindObject("MEDist_Particle1_Particle2_clone");
    // auto hMEPair = (TH1F *)TransformToMeV1D(hMEPair_GeV);
    // auto hMEAntiPair_GeV = (TH1F *)List2A->FindObject("MEDist_Particle0_Particle3_clone");
    // auto hMEAntiPair = (TH1F *)TransformToMeV1D(hMEAntiPair_GeV);
    // auto hME = (TH1F *)hMEPair->Clone("hMEData");
    // hME->Add(hMEAntiPair);

    // /// MC
    // TString fName_MC = homeFolder + foldernameMC + TString::Format("CFOutput_ALK_0.root");
    // TFile *FileMC;
    // FileMC = new TFile(fName_MC, "read");
    // TList *List1_MC = (TList *)FileMC->Get("PairDist");
    // TList *List1A_MC = (TList *)FileMC->Get("AntiPairDist");
    // TList *List2_MC = (TList *)List1_MC->FindObject("Pair");
    // TList *List2A_MC = (TList *)List1A_MC->FindObject("Pair");
    // auto hSEPairMC_GeV = (TH1F *)List2_MC->FindObject("SEDist_Particle1_Particle2_clone");
    // auto hSEPairMC = (TH1F *)TransformToMeV1D(hSEPairMC_GeV);
    // auto hSEAntiPairMC_GeV = (TH1F *)List2A_MC->FindObject("SEDist_Particle0_Particle3_clone");
    // auto hSEAntiPairMC = (TH1F *)TransformToMeV1D(hSEAntiPairMC_GeV);
    // auto hSEMC = (TH1F *)hSEPairMC->Clone("hSEMC");
    // hSEMC->Add(hSEAntiPairMC);

    // auto hMEPairMC_GeV = (TH1F *)List2_MC->FindObject("MEDist_Particle1_Particle2_clone");
    // auto hMEPairMC = (TH1F *)TransformToMeV1D(hMEPairMC_GeV);
    // auto hMEAntiPairMC_GeV = (TH1F *)List2A_MC->FindObject("MEDist_Particle0_Particle3_clone");
    // auto hMEAntiPairMC = (TH1F *)TransformToMeV1D(hMEAntiPairMC_GeV);
    // auto hMEMC = (TH1F *)hMEPairMC->Clone("hMEMC");
    // hMEMC->Add(hMEAntiPairMC);

    // double norm1 = 500.;
    // double norm2 = 800.;

    // /// Same Event manipulation
    // /// Data
    // auto hSEDataNorm = (TH1F *)hSE->Clone("hSEDataNorm");

    // auto hMEDataNorm = (TH1F *)hME->Clone("hMEDataNorm");
    // hMEDataNorm->Scale(hSEDataNorm->Integral(hSEDataNorm->FindBin(norm1), hSEDataNorm->FindBin(norm2)) /
    //                    hMEDataNorm->Integral(hMEDataNorm->FindBin(norm1), hMEDataNorm->FindBin(norm2)));

    // auto hSEDataMinusMEDataNorm = (TH1F *)hSEDataNorm->Clone("hSEDataMinusMEDataNorm");
    // hSEDataMinusMEDataNorm->Sumw2();
    // hSEDataMinusMEDataNorm->Add(hMEDataNorm, -1.);

    // /// MonteCarlo
    // auto hSEMCNorm = (TH1F *)hSEMC->Clone("hSEMCNorm");
    // hSEMCNorm->Scale(hSEDataNorm->Integral(hSEDataNorm->FindBin(norm1), hSEDataNorm->FindBin(norm2)) /
    //                  hSEMCNorm->Integral(hSEMCNorm->FindBin(norm1), hSEMCNorm->FindBin(norm2)));

    // auto hMEMCNorm = (TH1F *)hMEMC->Clone("hMEMCNorm");
    // // hMEMCNorm->Scale(1. / hMEMCNorm->Integral(hMEMCNorm->FindBin(norm1), hMEMCNorm->FindBin(norm2)));
    // hMEMCNorm->Scale(hSEDataNorm->Integral(hSEDataNorm->FindBin(norm1), hSEDataNorm->FindBin(norm2)) /
    //                  hMEMCNorm->Integral(hMEMCNorm->FindBin(norm1), hMEMCNorm->FindBin(norm2)));

    // auto hSEMCMinusMEMCNorm = (TH1F *)hSEMCNorm->Clone("hSEMCMinusMEMCNorm");
    // hSEMCMinusMEMCNorm->Sumw2();
    // hSEMCMinusMEMCNorm->Add(hMEMCNorm, -1.);

    // /// Fitting the residual background with Pol3
    // TF1 *fPol3_MC = nullptr;
    // TF1 *fFullPol3_MC = nullptr;
    // double CutOm_min = 200.;
    // double CutOm_max = 220.;
    // double upperFitRange = 500.;
    // TH1F *h_MCBoot;
    // TH1F *h_IMMinusMC;
    // TRandom3 rangen(5);

    // TTree *outTree = new TTree(TString::Format("Tree"), TString::Format("Tree"));
    // outTree->Branch("fFullPol3_MC", "TF1", &fFullPol3_MC, sizeof(TF1));

    // OutputFile->cd();
    // unsigned uSyst = 0;
    // float Value;
    // float Error;
    // float RanVal;
    // h_MCBoot = (TH1F *)hSEMCMinusMEMCNorm->Clone(TString::Format("h_MCBoot"));
    // h_MCBoot->Rebin(REBINMC);
    // h_MCBoot->Scale(1. / double(REBINMC));
    // // if (uSyst == 0)
    // // {
    // h_MCBoot->Write();
    // // }
    // // for (unsigned uBin = 0; uBin < h_MCBoot->GetNbinsX(); uBin++)
    // // {
    // //   Value = h_MCBoot->GetBinContent(uBin + 1);
    // //   Error = h_MCBoot->GetBinError(uBin + 1);
    // //   // RanVal = rangen.Gaus(Value, Error);
    // //   RanVal = Value;
    // //   h_MCBoot->SetBinContent(uBin + 1, RanVal);
    // //   h_MCBoot->SetBinError(uBin + 1, Error);
    // // }
    // fPol3_MC = new TF1(TString::Format("fPol3_MC_var%i", uSyst), fit_Pol4MC_Exclude, 0., upperFitRange, 7);
    // fPol3_MC->SetLineColor(kGreen + 2);
    // fPol3_MC->FixParameter(5, CutOm_min);
    // fPol3_MC->FixParameter(6, CutOm_max);
    // fPol3_MC->SetParameter(0, 0.5);
    // // fPol3_MC->SetParLimits(0, 0.1, 5.);
    // fPol3_MC->SetParameter(1, -0.001);
    // fPol3_MC->SetParameter(2, 0.0001);
    // fPol3_MC->SetParameter(3, -0.0001);
    // fPol3_MC->SetParameter(4, +0.0001);
    // h_MCBoot->Fit(fPol3_MC, "Q,S, N, R, M");
    // fPol3_MC->Write();

    // fFullPol3_MC = new TF1(TString::Format("fFullPol3_MC_var%i", uSyst), fit_Pol4MC, 0., upperFitRange, 7);
    // fFullPol3_MC->SetLineColor(kYellow + 3);
    // for (unsigned uPar = 0; uPar < 7; uPar++) {
    //     fFullPol3_MC->FixParameter(uPar, fPol3_MC->GetParameter(uPar));
    // }
    // fFullPol3_MC->Write();

    // h_IMMinusMC = (TH1F *)hSEDataMinusMEDataNorm->Clone(TString::Format("h_IMMinusMC_var%i", uSyst));
    // h_IMMinusMC->Rebin(REBINMC);
    // h_IMMinusMC->Scale(1. / double(REBINMC));
    // float ValueData;
    // float ValueMCFit;
    // for (unsigned uBin = 0; uBin < h_IMMinusMC->GetNbinsX(); uBin++) {
    //     double MOM = h_IMMinusMC->GetBinCenter(uBin + 1);
    //     if (MOM < 0. || MOM > upperFitRange) {
    //         h_IMMinusMC->SetBinContent(uBin + 1, 0.);
    //         h_IMMinusMC->SetBinError(uBin + 1, 0.);
    //         continue;
    //     }

    //     ValueData = h_IMMinusMC->GetBinContent(uBin + 1);
    //     ValueMCFit = fFullPol3_MC->Eval(MOM);
    //     Error = h_IMMinusMC->GetBinError(uBin + 1);
    //     h_IMMinusMC->SetBinContent(uBin + 1, ValueData - ValueMCFit);
    //     h_IMMinusMC->SetBinError(uBin + 1, Error);
    // }
    // h_IMMinusMC->Write();

    // /// Writing as a function of the energy

    // double LastBink = h_IMMinusMC->GetXaxis()->GetBinUpEdge(h_IMMinusMC->GetNbinsX());
    // double MaxEn = sqrt(LastBink * LastBink + mL * mL) + sqrt(LastBink * LastBink + mK * mK);

    // const int nBins = h_IMMinusMC->GetNbinsX();
    // Float_t *BinBoundaries = new Float_t[nBins + 1];

    // BinBoundaries[0] = EnThrLK;
    // for (unsigned uBin = 1; uBin < h_IMMinusMC->GetNbinsX() + 1; uBin++) {
    //     // BinBoundaries[uBin] = pow(h_IMMinusMC->GetXaxis()->GetBinUpEdge(uBin), 2.) / (2. * RedMass) + EnThrLK;
    //     BinBoundaries[uBin] =
    //         sqrt(h_IMMinusMC->GetXaxis()->GetBinUpEdge(uBin) * h_IMMinusMC->GetXaxis()->GetBinUpEdge(uBin) + mL * mL) +
    //         sqrt(h_IMMinusMC->GetXaxis()->GetBinUpEdge(uBin) * h_IMMinusMC->GetXaxis()->GetBinUpEdge(uBin) + mK * mK);
    // }

    // TH1F *h_IMMinusMCEn =
    //     new TH1F(TString::Format("h_IMMinusMCEn_var%i", uSyst), TString::Format("h_IMMinusMCEn_var%i", uSyst),
    //              h_IMMinusMC->GetNbinsX(), BinBoundaries);

    // for (unsigned uBin = 0; uBin < h_IMMinusMC->GetNbinsX(); uBin++) {
    //     double MOM = h_IMMinusMCEn->GetBinCenter(uBin + 1);
    //     double En = sqrt(MOM * MOM + mL * mL) + sqrt(MOM * MOM + mK * mK);

    //     ValueData = h_IMMinusMC->GetBinContent(uBin + 1);
    //     Error = h_IMMinusMC->GetBinError(uBin + 1);
    //     h_IMMinusMCEn->SetBinContent(uBin + 1, ValueData);
    //     h_IMMinusMCEn->SetBinError(uBin + 1, Error);
    // }

    // h_IMMinusMCEn->Write();

    // h_IMMinusMCEn->SetFillColor(kBlack);
    // h_IMMinusMCEn->SetMarkerColor(kBlack);
    // h_IMMinusMCEn->SetLineColor(kBlack);
    // h_IMMinusMCEn->SetLineWidth(2);
    // h_IMMinusMCEn->SetMarkerStyle(kOpenCircle);
    // h_IMMinusMCEn->SetMarkerSize(0.5);

    // h_IMMinusMC->SetFillColor(kBlack);
    // h_IMMinusMC->SetMarkerColor(kBlack);
    // h_IMMinusMC->SetLineColor(kBlack);
    // h_IMMinusMC->SetLineWidth(2);
    // h_IMMinusMC->SetMarkerStyle(kOpenCircle);
    // h_IMMinusMC->SetMarkerSize(0.5);
    // h_IMMinusMC->SetFillColor(kGray + 1);

    // double xmin = 0.5;
    // double xmax = 0.74;
    // double ymax = 0.665;
    // double ymin = 0.665;
    // unsigned NumRows = 1;

    // double MOM = 500.;
    // double EnMax = sqrt(MOM * MOM + mL * mL) + sqrt(MOM * MOM + mK * mK);
    // double EnMin = sqrt(0. * 0. + mL * mL) + sqrt(0. * 0. + mK * mK);

    // auto canvas = new TCanvas("cIm", "cIm", 800, 500);
    // gStyle->SetOptStat(0);
    // gStyle->SetOptTitle(0);
    // h_IMMinusMCEn->GetXaxis()->SetNdivisions(515);
    // h_IMMinusMCEn->GetXaxis()->SetRangeUser(BinBoundaries[0], EnMax);
    // h_IMMinusMCEn->GetXaxis()->SetTitle("M_{#Lambda#minus K^{-}#oplus#bar{#Lambda}#minus K^{+}}(MeV/#it{c}^{2})");
    // h_IMMinusMCEn->GetYaxis()->SetTitle("a. u.");
    // h_IMMinusMCEn->Draw();
    // canvas->Update();

    // TLine *line = new TLine(EnThrLK, 0, BinBoundaries[125], 0);
    // line->SetLineStyle(2);
    // line->Draw("same");

    // TLatex BeamText;
    // BeamText.SetTextSize(gStyle->GetTextSize() * 0.90);
    // BeamText.SetNDC(kTRUE);
    // BeamText.DrawLatex(0.50, 0.765, "ALICE Preliminary");
    // BeamText.DrawLatex(0.50, 0.715, "pp #sqrt{#it{s}} = 13 TeV");
    // BeamText.DrawLatex(0.50, 0.665, "High Mult. (0#minus0.17% INEL > 0)");
    // canvas->SaveAs(TString::Format("InvMassSpectra_Resonances.pdf"));
    // delete canvas;
    // delete line;

    // canvas = new TCanvas("cZoom", "cZoom", 800, 500);
    // gStyle->SetOptStat(0);
    // gStyle->SetOptTitle(0);
    // h_IMMinusMCEn->GetXaxis()->SetNdivisions(515);
    // h_IMMinusMCEn->GetXaxis()->SetRangeUser(EnThrLK, 1690.);
    // h_IMMinusMCEn->GetXaxis()->SetTitle("M_{#Lambda#minus K^{-}#oplus#bar{#Lambda}#minus K^{+}}(MeV/#it{c}^{2})");
    // h_IMMinusMCEn->Draw();

    // line = new TLine(EnThrLK, 0, 1690., 0);
    // line->Draw("same");

    // BeamText.SetTextSize(gStyle->GetTextSize() * 0.90);
    // BeamText.SetNDC(kTRUE);
    // BeamText.DrawLatex(0.30, 0.765, "ALICE pp #sqrt{#it{s}} = 13 TeV");
    // BeamText.DrawLatex(0.30, 0.715, "High Mult. (0#minus0.17% INEL > 0)");
    // canvas->SaveAs(TString::Format("InvMassSpectra_ResonancesZoom.pdf"));
    // delete canvas;
    // delete line;

    // // #########################################################
    // /// Plot IM as function of k*
    // double xmin_beam = 0.18;

    // canvas = new TCanvas("cIm", "cIm", 800, 500);
    // gStyle->SetOptStat(0);
    // gStyle->SetOptTitle(0);
    // canvas->SetTopMargin(0.15);
    // h_IMMinusMC->GetXaxis()->SetNdivisions(515);
    // h_IMMinusMC->GetXaxis()->SetRangeUser(0, 500.);
    // h_IMMinusMC->GetYaxis()->SetRangeUser(-6200, 35000.);

    // h_IMMinusMC->GetXaxis()->SetTitle("#it{k*} (MeV/#it{c})");
    // h_IMMinusMC->GetYaxis()->SetTitle("a. u.");
    // h_IMMinusMC->Draw();
    // canvas->Update();

    // // Omega
    // double MassOmega_meanEn =
    //     sqrt(MassOmega_mean * MassOmega_mean + mL * mL) + sqrt(MassOmega_mean * MassOmega_mean + mK * mK);
    // TLine *lineOm = new TLine(MassOmega_mean, -6200, MassOmega_mean, 35000);
    // lineOm->SetLineStyle(7);
    // lineOm->SetLineWidth(2);
    // lineOm->SetLineColor(kGreen + 2);
    // lineOm->SetLineColorAlpha(kGreen + 2, 0.5);
    // lineOm->Draw("same");
    // TLatex BeamTextOm;
    // BeamTextOm.SetTextSize(gStyle->GetTextSize() * 0.70);
    // BeamTextOm.SetNDC(kFALSE);
    // BeamTextOm.SetTextColor(kGreen + 2);
    // BeamTextOm.DrawLatex(149, 31000, TString::Format("#Omega: %.1f", MassOmega_meanEn));

    // // Xi 1690
    // double MassXi1_meanEn = sqrt(MassXi1_mean * MassXi1_mean + mL * mL) + sqrt(MassXi1_mean * MassXi1_mean + mK * mK);
    // TLine *lineXi1 = new TLine(MassXi1_mean, -6200, MassXi1_mean, 35000);
    // lineXi1->SetLineStyle(7);
    // lineXi1->SetLineWidth(2);
    // lineXi1->SetLineColor(kAzure + 1);
    // lineXi1->SetLineColorAlpha(kAzure + 1, 0.5);
    // lineXi1->Draw("same");

    // TLatex BeamTextXi1;
    // BeamTextXi1.SetTextSize(gStyle->GetTextSize() * 0.70);
    // BeamTextXi1.SetNDC(kFALSE);
    // BeamTextXi1.SetTextColor(kAzure + 1);
    // BeamTextXi1.DrawLatex(248, 31000, TString::Format("#Xi(1690): %.1f", MassXi1_meanEn));

    // // Xi 1820
    // double MassXi2_meanEn = sqrt(MassXi2_mean * MassXi2_mean + mL * mL) + sqrt(MassXi2_mean * MassXi2_mean + mK * mK);
    // TLine *lineXi2 = new TLine(MassXi2_mean, -6200, MassXi2_mean, 35000);
    // lineXi2->SetLineStyle(7);
    // lineXi2->SetLineWidth(2);
    // lineXi2->SetLineColor(kOrange + 1);
    // lineXi2->SetLineColorAlpha(kOrange + 1, 0.5);
    // lineXi2->Draw("same");

    // TLatex BeamTextXi2;
    // BeamTextXi2.SetTextSize(gStyle->GetTextSize() * 0.70);
    // BeamTextXi2.SetNDC(kFALSE);
    // BeamTextXi2.SetTextColor(kOrange + 1);
    // BeamTextXi2.DrawLatex(402, 31000, TString::Format("#Xi(1820): %.1f", MassXi2_meanEn));

    // // Xi 1620
    // double kstar_Xi1620 = sqrt(pow(MassXi1620_mean, 4) - 2 * MassXi1620_mean * MassXi1620_mean * mL * mL + pow(mL, 4) -
    //                            2 * MassXi1620_mean * MassXi1620_mean * mK * mK - 2 * mL * mL * mK * mK + pow(mK, 4)) /
    //                       (2 * MassXi1620_mean);
    // TLine *lineXiex = new TLine(kstar_Xi1620, -6200, kstar_Xi1620, 12000);
    // lineXiex->SetLineStyle(7);
    // lineXiex->SetLineWidth(2);
    // lineXiex->SetLineColor(kPink + 2);
    // lineXiex->SetLineColorAlpha(kPink + 2, 0.5);
    // lineXiex->Draw("same");

    // TLine *lineXiex2 = new TLine(kstar_Xi1620, 24000, kstar_Xi1620, 35000);
    // lineXiex2->SetLineStyle(7);
    // lineXiex2->SetLineWidth(2);
    // lineXiex2->SetLineColor(kPink + 2);
    // lineXiex2->SetLineColorAlpha(kPink + 2, 0.5);
    // lineXiex2->Draw("same");

    // TLatex BeamTextXi1620;
    // BeamTextXi1620.SetTextSize(gStyle->GetTextSize() * 0.70);
    // BeamTextXi1620.SetNDC(kFALSE);
    // BeamTextXi1620.SetTextColor(kPink + 2);
    // BeamTextXi1620.DrawLatex(81, 28000, TString::Format("#Xi(1620): %.1f", MassXi1620_mean));

    // TLegend *legend = new TLegend(xmin_beam, 0.515 - 0.05, 0.30, 0.515);  // lbrt
    // legend->SetBorderSize(0);
    // legend->SetTextFont(42);
    // legend->SetTextSize(gStyle->GetTextSize() * 0.7);
    // legend->SetFillColor(kWhite);
    // TH1F *hCk_Fake = (TH1F *)h_IMMinusMC->Clone("hCk_Fake");
    // hCk_Fake->SetName("hCk_Fake");
    // hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());
    // legend->AddEntry(hCk_Fake, "#Lambda#minus K^{-} #oplus #bar{#Lambda}#minus K^{+}", "fp");
    // legend->Draw("same");

    // line = new TLine(0, 0, 500., 0);
    // line->SetLineStyle(2);
    // line->Draw("same");

    // TF1 *f1 = nullptr;
    // // f1 = new TF1("fit_axis", "sqrt(x*x+[0]*[0])+sqrt(x*x+[1]*[1])", 0., 500.);
    // f1 = new TF1("fit_axis", "sqrt(x*x*x*x-2*x*x*[0]*[0]+pow([0],4)-2*x*x*[1]*[1]-2*[0]*[0]*[1]*[1]+pow([1],4))/(2*x)",
    //              EnMin, EnMax);
    // f1->SetParameter(0, mL);
    // f1->SetParameter(1, mK);
    // TGaxis *axis =
    //     new TGaxis(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax(), "fit_axis", 507, "-");
    // axis->SetLabelSize(0.045);
    // axis->SetTextFont(42);
    // axis->SetLabelFont(42);
    // axis->SetLabelOffset(0.01);
    // axis->SetLabelColor(kBlack);
    // axis->SetTitleSize(0.05);
    // axis->SetTitleOffset(1.2);  // 1.2
    // axis->SetTextSizePixels(26);
    // axis->SetTitle("M_{#Lambda#minus K^{-}#oplus#bar{#Lambda}#minus K^{+}}(MeV/#it{c}^{2})");
    // axis->Draw();

    // BeamText.SetTextSize(gStyle->GetTextSize() * 0.70);
    // BeamText.SetNDC(kTRUE);
    // BeamText.DrawLatex(xmin_beam, 0.615, "ALICE pp #sqrt{#it{s}} = 13 TeV");
    // BeamText.DrawLatex(xmin_beam, 0.565, "High Mult. (0#minus0.17% INEL > 0)");
    // canvas->SaveAs(TString::Format("InvMassSpectraKSTAR_Resonances.pdf"));
    // delete canvas;
    // delete line;

    // // }/// end bootstrap

    // /// Reading the resonance properties from the fit
    // /*  TFile *InputTree;

    //   TString fName_Output;

    //   TString fName_OutputFigs;
    //   TString fName_Output_AncestorsFigs;

    //   fName_Output = homeFolder + OutputFolderLednickySillGaussERE;
    //   fName_OutputFigs = homeFolder + OutputFolderLednickySillGaussERE;

    //   std::cout << "Output folder is:" << fName_Output << std::endl;

    //   int TypeOfFit = 4;
    //   bool WithOmega = false;
    //   bool IsEffRangeZero = false;
    //   int UsePbPb = 0;
    //   InputTree = TFile::Open(fName_Output +
    //   TString::Format("OutputTree_rebinMC2rebinDATA1_Fit%i_Back0_Omega%i_EffRangeZero%i_UsePbPb%i_Resonance.root",
    //   TypeOfFit, int(WithOmega), int(IsEffRangeZero), UsePbPb));

    //   TH1F *hMassOmegafFullFemto = (TH1F *)(InputTree->Get("hMassOmegafFullFemto"));
    //   TH1F *hMassXi1fFullFemto = (TH1F *)(InputTree->Get("hMassXi1fFullFemto"));
    //   TH1F *hMassXi2fFullFemto = (TH1F *)(InputTree->Get("hMassXi2fFullFemto"));
    //   TH1F *hGammaOmegafFullFemto = (TH1F *)(InputTree->Get("hGammaOmegafFullFemto"));
    //   TH1F *hGammaXi1fFullFemto = (TH1F *)(InputTree->Get("hGammaXi1fFullFemto"));
    //   TH1F *hGammaXi2fFullFemto = (TH1F *)(InputTree->Get("hGammaXi2fFullFemto"));
    //   TH1F *hHeightOmegafFullFemto = (TH1F *)(InputTree->Get("hHeightOmegafFullFemto"));
    //   TH1F *hHeightXi1fFullFemto = (TH1F *)(InputTree->Get("hHeightXi1fFullFemto"));
    //   TH1F *hHeightXi2fFullFemto = (TH1F *)(InputTree->Get("hHeightXi2fFullFemto"));

    //   TH1F *hMassOmegafFullFemto_stat = (TH1F *)(InputTree->Get("hMassOmegafFullFemto_stat"));
    //   TH1F *hMassXi1fFullFemto_stat = (TH1F *)(InputTree->Get("hMassXi1fFullFemto_stat"));
    //   TH1F *hMassXi2fFullFemto_stat = (TH1F *)(InputTree->Get("hMassXi2fFullFemto_stat"));
    //   TH1F *hGammaOmegafFullFemto_stat = (TH1F *)(InputTree->Get("hGammaOmegafFullFemto_stat"));
    //   TH1F *hGammaXi1fFullFemto_stat = (TH1F *)(InputTree->Get("hGammaXi1fFullFemto_stat"));
    //   TH1F *hGammaXi2fFullFemto_stat = (TH1F *)(InputTree->Get("hGammaXi2fFullFemto_stat"));
    //   TH1F *hHeightOmegafFullFemto_stat = (TH1F *)(InputTree->Get("hHeightOmegafFullFemto_stat"));
    //   TH1F *hHeightXi1fFullFemto_stat = (TH1F *)(InputTree->Get("hHeightXi1fFullFemto_stat"));
    //   TH1F *hHeightXi2fFullFemto_stat = (TH1F *)(InputTree->Get("hHeightXi2fFullFemto_stat"));

    //   TH1F *hMassXi1620fFullFemto = (TH1F *)(InputTree->Get("hMassXi1620fFullFemto"));
    //   TH1F *hGammaXi1620XiPifFullFemto = (TH1F *)(InputTree->Get("hGammaXi1620XiPifFullFemto"));
    //   TH1F *hGammaXi1620LKfFullFemto = (TH1F *)(InputTree->Get("hGammaXi1620LKfFullFemto"));
    //   TH1F *hHeightXi1620fFullFemto = (TH1F *)(InputTree->Get("hHeightXi1620fFullFemto"));

    //   TH1F *hMassXi1620fFullFemto_stat = (TH1F *)(InputTree->Get("hMassXi1620fFullFemto_stat"));
    //   TH1F *hGammaXi1620XiPifFullFemto_stat = (TH1F *)(InputTree->Get("hGammaXi1620XiPifFullFemto_stat"));
    //   TH1F *hGammaXi1620LKfFullFemto_stat = (TH1F *)(InputTree->Get("hGammaXi1620LKfFullFemto_stat"));
    //   TH1F *hHeightXi1620fFullFemto_stat = (TH1F *)(InputTree->Get("hHeightXi1620fFullFemto_stat"));
    //   TH1F *hWeightfFullFemto_stat = (TH1F *)(InputTree->Get("hWeightfFullFemto_stat"));

    //   double MassOmega_mean = hMassOmegafFullFemto_stat->GetMean();
    //   double MassOmegaEnergy_mean = sqrt(MassOmega_mean * MassOmega_mean + mL * mL) + sqrt(MassOmega_mean *
    //   MassOmega_mean + mK * mK);

    //   double MassOmega_sig_stat = hMassOmegafFullFemto_stat->GetStdDev();
    //   double MassOmega_sig_tot = hMassOmegafFullFemto->GetStdDev();
    //   double MassOmega_sig_syst = TMath::Sqrt(MassOmega_sig_tot * MassOmega_sig_tot - MassOmega_sig_stat *
    //   MassOmega_sig_stat);

    //   double GammaOmega_mean = hGammaOmegafFullFemto_stat->GetMean();
    //   double GammaOmega_sig_stat = hGammaOmegafFullFemto_stat->GetStdDev();
    //   double GammaOmega_sig_tot = hGammaOmegafFullFemto->GetStdDev();
    //   double GammaOmega_sig_syst = TMath::Sqrt(GammaOmega_sig_tot * GammaOmega_sig_tot - GammaOmega_sig_stat *
    //   GammaOmega_sig_stat);

    //   double HeightOmega_mean = hHeightOmegafFullFemto_stat->GetMean();
    //   double HeightOmega_sig_stat = hHeightOmegafFullFemto_stat->GetStdDev();
    //   double HeightOmega_sig_tot = hHeightOmegafFullFemto->GetStdDev();
    //   double HeightOmega_sig_syst = TMath::Sqrt(HeightOmega_sig_tot * HeightOmega_sig_tot - HeightOmega_sig_stat *
    //   HeightOmega_sig_stat);

    //   /// Ξ(1690)
    //   double MassXi1_mean = hMassXi1fFullFemto_stat->GetMean();
    //   double MassXi1Energy_mean = sqrt(MassXi1_mean * MassXi1_mean + mL * mL) + sqrt(MassXi1_mean * MassXi1_mean + mK *
    //   mK); double MassXi1_sig_stat = hMassXi1fFullFemto_stat->GetStdDev(); double MassXi1_sig_tot =
    //   hMassXi1fFullFemto->GetStdDev(); double MassXi1_sig_syst = TMath::Sqrt(MassXi1_sig_tot * MassXi1_sig_tot -
    //   MassXi1_sig_stat * MassXi1_sig_stat);

    //   double GammaXi1_mean = hGammaXi1fFullFemto_stat->GetMean();
    //   double GammaXi1_sig_stat = hGammaXi1fFullFemto_stat->GetStdDev();
    //   double GammaXi1_sig_tot = hGammaXi1fFullFemto->GetStdDev();
    //   double GammaXi1_sig_syst = TMath::Sqrt(GammaXi1_sig_tot * GammaXi1_sig_tot - GammaXi1_sig_stat *
    //   GammaXi1_sig_stat);

    //   double HeightXi1_mean = hHeightXi1fFullFemto_stat->GetMean();
    //   double HeightXi1_sig_stat = hHeightXi1fFullFemto_stat->GetStdDev();
    //   double HeightXi1_sig_tot = hHeightXi1fFullFemto->GetStdDev();
    //   double HeightXi1_sig_syst = TMath::Sqrt(HeightXi1_sig_tot * HeightXi1_sig_tot - HeightXi1_sig_stat *
    //   HeightXi1_sig_stat);

    //   /// Ξ(1820)
    //   double MassXi2_mean = hMassXi2fFullFemto_stat->GetMean();
    //   double MassXi2Energy_mean = sqrt(MassXi2_mean * MassXi2_mean + mL * mL) + sqrt(MassXi2_mean * MassXi2_mean + mK *
    //   mK); double MassXi2_sig_stat = hMassXi2fFullFemto_stat->GetStdDev(); double MassXi2_sig_tot =
    //   hMassXi2fFullFemto->GetStdDev(); double MassXi2_sig_syst = TMath::Sqrt(MassXi2_sig_tot * MassXi2_sig_tot -
    //   MassXi2_sig_stat * MassXi2_sig_stat);

    //   double GammaXi2_mean = hGammaXi2fFullFemto_stat->GetMean();
    //   double GammaXi2_sig_stat = hGammaXi2fFullFemto_stat->GetStdDev();
    //   double GammaXi2_sig_tot = hGammaXi2fFullFemto->GetStdDev();
    //   double GammaXi2_sig_syst = TMath::Sqrt(GammaXi2_sig_tot * GammaXi2_sig_tot - GammaXi2_sig_stat *
    //   GammaXi2_sig_stat);

    //   double HeightXi2_mean = hHeightXi2fFullFemto_stat->GetMean();
    //   double HeightXi2_sig_stat = hHeightXi2fFullFemto_stat->GetStdDev();
    //   double HeightXi2_sig_tot = hHeightXi2fFullFemto->GetStdDev();
    //   double HeightXi2_sig_syst = TMath::Sqrt(HeightXi2_sig_tot * HeightXi2_sig_tot - HeightXi2_sig_stat *
    //   HeightXi2_sig_stat);

    //   TF1 *fGaussOmega = nullptr;
    //   double EMinOm = 0.;
    //   double EMaxOm = hME->GetXaxis()->GetBinUpEdge(hME->GetNbinsX());
    //   std::cout << "EMaxOm = " << EMaxOm << std::endl;
    //   fGaussOmega = new TF1(TString::Format("fGaussOmega"), fit_WeightedGauss, EMinOm, EMaxOm, 3);
    //   fGaussOmega->SetNpx(1000000);
    //   fGaussOmega->SetLineColor(kGreen + 2);
    //   fGaussOmega->FixParameter(0, MassOmega_mean);
    //   fGaussOmega->FixParameter(1, GammaOmega_mean);
    //   fGaussOmega->FixParameter(2, 1.);
    //   std::cout << "Int.fGaussOmega = " << fGaussOmega->Integral(180, 240) << std::endl;

    //   auto hMEDataNormFullRange = (TH1F *)hME->Clone("hMEDataNormFullRange");
    //   hMEDataNormFullRange->Scale(1. / hMEDataNormFullRange->Integral());

    //   TH1F *hOmegaMENormFullRange = new TH1F(TString::Format("hOmegaMENormFullRange"),
    //   TString::Format("hOmegaMENormFullRange"), hMEDataNormFullRange->GetNbinsX(), 0.,
    //   hMEDataNormFullRange->GetXaxis()->GetBinUpEdge(hMEDataNormFullRange->GetNbinsX()));
    //   hOmegaMENormFullRange->FillRandom("fGaussOmega", 1000000);
    //   hOmegaMENormFullRange->Scale(1. / hOmegaMENormFullRange->Integral());
    //   hOmegaMENormFullRange->Divide(hMEDataNormFullRange);

    //   double TheIntegral_Omega = hOmegaMENormFullRange->Integral(40., 5000.);
    //   std::cout << "k Omega=" << TheIntegral_Omega << std::endl;
    //   // h_IMMinusMCEn->Fit(fGaussOmega, "Q,S, N, R, M");
    //   hOmegaMENormFullRange->Scale(1. / hOmegaMENormFullRange->Integral());
    //   std::cout << "Int.hOmegaMENormFullRange = " << hOmegaMENormFullRange->Integral() << std::endl;

    //   TF1 *fBWXi1690 = nullptr;
    //   fBWXi1690 = new TF1(TString::Format("fBWXi1690"), fit_WeightedBreitWigner, EMinOm, EMaxOm, 3);
    //   fBWXi1690->SetNpx(1000000);
    //   fBWXi1690->SetLineColor(kRed + 2);
    //   fBWXi1690->FixParameter(0, MassXi1_mean);
    //   fBWXi1690->FixParameter(1, GammaXi1_mean);
    //   fBWXi1690->FixParameter(2, 1.);
    //   std::cout << "Int.fBWXi1690 = " << fBWXi1690->Integral(EMinOm, EMaxOm) << std::endl;

    //   TF1 *fBWXi1820 = nullptr;
    //   fBWXi1820 = new TF1(TString::Format("fBWXi1820"), fit_WeightedBreitWigner, EMinOm, EMaxOm, 3);
    //   fBWXi1820->SetNpx(1000000);
    //   fBWXi1820->SetLineColor(kRed + 2);
    //   fBWXi1820->FixParameter(0, MassXi2_mean);
    //   fBWXi1820->FixParameter(1, GammaXi2_mean);
    //   fBWXi1820->FixParameter(2, 1.);
    //   std::cout << "Int.fBWXi1820 = " << fBWXi1820->Integral(EMinOm, EMaxOm) << std::endl;

    //   TH1F *hXi1MENormFullRange = new TH1F(TString::Format("hXi1MENormFullRange"),
    //   TString::Format("hXi1MENormFullRange"), hMEDataNormFullRange->GetNbinsX(), 0.,
    //   hMEDataNormFullRange->GetXaxis()->GetBinUpEdge(hMEDataNormFullRange->GetNbinsX()));
    //   hXi1MENormFullRange->FillRandom("fBWXi1690", 1000000);
    //   hXi1MENormFullRange->Scale(1. / hXi1MENormFullRange->Integral());
    //   hXi1MENormFullRange->Divide(hMEDataNormFullRange);

    //   double TheIntegral_Xi1 = hXi1MENormFullRange->Integral(40., 5000.);
    //   std::cout << "k Xi 1690=" << TheIntegral_Xi1 << std::endl;

    //   TH1F *hXi2MENormFullRange = new TH1F(TString::Format("hXi2MENormFullRange"),
    //   TString::Format("hXi2MENormFullRange"), hMEDataNormFullRange->GetNbinsX(), 0.,
    //   hMEDataNormFullRange->GetXaxis()->GetBinUpEdge(hMEDataNormFullRange->GetNbinsX()));
    //   hXi2MENormFullRange->FillRandom("fBWXi1820", 1000000);
    //   hXi2MENormFullRange->Scale(1. / hXi2MENormFullRange->Integral());
    //   hXi2MENormFullRange->Divide(hMEDataNormFullRange);

    //   double TheIntegral_Xi2 = hXi2MENormFullRange->Integral(40., 5000.);
    //   std::cout << "k Xi 1820=" << TheIntegral_Xi2 << std::endl;
    // */
    // OutputFile->cd();
    // // fGaussOmega->Write();
    // // hMEDataNormFullRange->Write();
    // // hOmegaMENormFullRange->Write();
    // // fBWXi1690->Write();
    // // hXi1MENormFullRange->Write();
    // // fBWXi1820->Write();
    // // hXi2MENormFullRange->Write();
    // //-----------------------------------------------------

    // hSE->Write();
    // hME->Write();

    // hSEMC->Write();
    // hMEMC->Write();

    // hSEDataNorm->Write();
    // hMEDataNorm->Write();
    // hSEMCNorm->Write();
    // hMEMCNorm->Write();
    // hSEDataMinusMEDataNorm->Write();
    // hSEMCMinusMEMCNorm->Write();

    // delete OutputFile;
    // delete fPol3_MC;
}