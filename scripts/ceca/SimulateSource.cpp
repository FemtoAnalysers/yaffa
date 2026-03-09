/*
 * Script to simulate the source function
 */

// C++ headers
#include <omp.h>
#include <unistd.h>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

// Third party headers
#include <boost/algorithm/string.hpp>

// ROOT headers
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGenPhaseSpace.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"

// DLM headers
#include "CATS.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "CECA.h"
#include "CommonAnaFunctions.h"
#include "DLM_CppTools.h"
#include "DLM_Histo.h"
#include "DLM_HistoAnalysis.h"
#include "DLM_MathFunctions.h"
#include "DLM_MultiFit.h"
#include "DLM_OmpTools.h"
#include "DLM_Potentials.h"
#include "DLM_Random.h"
#include "DLM_RootWrapper.h"
#include "DLM_Source.h"
#include "TREPNI.h"

double GaussianSource(double* x, double* par) {
    // Variables
    double rstar = x[0];

    // Parameters
    double rzero = par[0];

    // Calculation
    double jacobian = 4 * Pi * rstar * rstar;
    double norm = 1. / pow(4. * Pi * rzero * rzero, 1.5);
    double expo = exp(-(rstar * rstar) / (4. * rzero * rzero));

    return jacobian * norm * expo;
}

double ScaledGauss(double* x, double* par) { return par[1] * GaussSourceTF1(x, par); }

double GaussSourceMean(double* x, double* Pars) {
    double& size = Pars[0];
    TF1 f_gauss("f_gauss", "[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]", 0,
                256);
    f_gauss.FixParameter(0, 1);
    f_gauss.FixParameter(1, size);
    return f_gauss.Mean(0, 256);
}

double GaussFromMean(const double mean) {
    TH1F hHist("hHist", "hHist", 1, 0, 1);
    hHist.SetBinContent(1, mean);
    hHist.SetBinError(1, mean * 0.01);
    TF1 fHist("fHist", GaussSourceMean, 0, 1, 1);
    fHist.SetParameter(0, mean);
    fHist.SetParLimits(0, mean * 0.1, mean * 2.);
    hHist.Fit(&fHist, "Q, S, N, R, M");
    return fHist.GetParameter(0);
}

// for Lambda, more like pT in 0.4 --> inf
DLM_Histo<float>* GetPtEta_13TeV(TString FileNameIn, TString GraphNameIn, const double pT_min, const double pT_max,
                                 const double EtaCut) {
    TGraphAsymmErrors* gSpectrum;
    TFile file_in(FileNameIn, "read");
    gSpectrum = (TGraphAsymmErrors*)file_in.Get(GraphNameIn);
    if (!gSpectrum) printf("ISSUE with gSpectrum\n");
    // gROOT->cd();
    double* BinRange = new double[gSpectrum->GetN() + 1];
    double* BinCenter = new double[gSpectrum->GetN()];
    double* BinContent = new double[gSpectrum->GetN()];
    // printf("Iter over %u\n",gSpectrum->GetN());
    for (unsigned uBin = 0; uBin < gSpectrum->GetN(); uBin++) {
        // printf(" -- %u\n",uBin);
        double pT, Yield;
        gSpectrum->GetPoint(uBin, pT, Yield);
        pT *= 1000;

        double pT_low = pT - gSpectrum->GetErrorXlow(uBin) * 1000.;
        double pT_high = pT + gSpectrum->GetErrorXhigh(uBin) * 1000.;

        BinCenter[uBin] = 0.5 * (pT_high + pT_low);
        if (BinCenter[uBin] < pT_min || BinCenter[uBin] > pT_max)
            BinContent[uBin] = 0;
        else
            BinContent[uBin] = Yield;

        // if(uBin) BinRange[uBin] = BinRange[uBin-1];
        // else BinRange[uBin] = pT_low;
        BinRange[uBin] = pT_low;

        if (uBin == gSpectrum->GetN() - 1) {
            BinRange[uBin + 1] = pT_high;
        }
    }

    // for(unsigned uBin=0; uBin<=gSpectrum->GetN(); uBin++){
    // printf(" BinRange[%u] = %.f\n",uBin,BinRange[uBin]);
    // }

    DLM_Histo<float>* dlm_pT_eta = new DLM_Histo<float>();
    dlm_pT_eta->SetUp(2);
    dlm_pT_eta->SetUp(0, gSpectrum->GetN(), BinRange);
    dlm_pT_eta->SetUp(1, 1, -EtaCut, EtaCut);
    dlm_pT_eta->Initialize();
    for (unsigned uBin = 0; uBin < gSpectrum->GetN(); uBin++) {
        dlm_pT_eta->SetBinContent(uBin, BinContent[uBin]);
    }

    file_in.Close();
    delete[] BinRange;
    delete[] BinCenter;
    delete[] BinContent;
    return dlm_pT_eta;
}

int main(int argc, const char** argv) {
    // Parse arguments form command line
    std::string cfg_file;
    if (argc == 1) {
        cfg_file = "config.yaml";
    } else {
        cfg_file = argv[1];
    }

    // Load configuration from file
    YAML::Node cfg = YAML::LoadFile(cfg_file);

    unsigned NUM_CPU = cfg["ncpu"].as<unsigned>() ? cfg["ncpu"].as<unsigned>() : omp_get_max_threads();
    std::string oFileName = cfg["ofile"].as<std::string>();
    std::string system = cfg["system"].as<std::string>();
    const unsigned globalTimeout = cfg["glob_timeout"].as<unsigned>();
    const unsigned threadTimeout = cfg["thread_timeout"].as<unsigned>();
    double HadronSize = cfg["hadron_size"].as<double>();
    double HadronSlope = cfg["hadron_slope"].as<double>();
    double EtaCut = cfg["eta_cut"].as<double>();
    const bool PROTON_RESO = cfg["enable_resonances"].as<bool>();
    const double frac_protons = cfg["frac_prim"]["p"].as<double>();;
    const bool EQUALIZE_TAU = cfg["equalize_tau"].as<bool>();
    const unsigned Multiplicity = cfg["mult"].as<unsigned>();
    const double femto_region = cfg["femto_region"].as<double>();;
    const unsigned target_yield = cfg["target_yield"].as<unsigned>();  // originally 4M
    const int EffFix = 9001;
    // Effectively remove the lorentz boost effect by setting the particle's masses to 1 TeV
    const bool removeBoost = cfg["remove_boost"].as<bool>();


    // Default parameters:
    // double HadronSize = 0;   // 0.75
    // double HadronSlope = 0;  // 0.2
    // const double EtaCut = 0.8;
    // const bool PROTON_RESO = true;
    // const double frac_protons = 35.78;
    // const bool EQUALIZE_TAU = true;
    // const unsigned Multiplicity = 2;
    // const double femto_region = 100;
    // const unsigned target_yield = 512 * 1000 / 64.;  // originally 4M
    // const int EffFix = -1402;
    // const bool removeBoost = false; // effectively remove the lorentz boost effect by setting the particle's masses
    // to 1 TeV

    // Reproduce simple gaussian-like source in CECA paper (only rd = 0.85 fm):
    // double HadronSize = 0;   // 0.75
    // double HadronSlope = 0;  // 0.2
    // const double EtaCut = 0.8;
    // const bool PROTON_RESO = false;
    // const double frac_protons = 35.78;
    // const bool EQUALIZE_TAU = true;
    // const unsigned Multiplicity = 2;
    // const double femto_region = 100;
    // const unsigned target_yield = 512 * 1000 / 64.;  // originally 4M
    // const int EffFix = 9000;
    // const bool removeBoost = false; // effectively remove the lorentz boost effect by setting the particle's masses
    // to 1 TeV

    // Reproduce simple 1fm source with and remove lorentz boost


    // Reproduce simple 1fm source - verify non gaussianity of source due to lorentz boost
    // double HadronSize = 0;   // 0.75
    // double HadronSlope = 0;  // 0.2
    // const double EtaCut = 0.8;
    // const bool PROTON_RESO = false;
    // const double frac_protons = 35.78;
    // const bool EQUALIZE_TAU = true;
    // const unsigned Multiplicity = 2;
    // const double femto_region = 100;
    // const unsigned target_yield = 512 * 1000 / 64.;  // originally 4M
    // const int EffFix = 9002;
    // const bool removeBoost = false; // effectively remove the lorentz boost effect by setting the particle's masses
    // to 1 TeV

    // we run to either reproduce the core of 0.97,
    // or the upper limit of reff = 1.06+0.04
    // this leads to a 10% difference in the SP core source
    double rSP_core = 0;
    double rSP_dispZ = 0;
    double rSP_hadr = 0;
    double rSP_hadrZ = 0;
    double rSP_hflc = 0;
    double rSP_tau = 0;
    bool tau_prp = true;
    double rSP_tflc = 0;
    double rSP_ThK = 0;

    // default is true
    bool rSP_FixedHadr = true;

    // default is 0
    float rSP_FragBeta = 0;

    if (system == "pp") {
        // jaime issue, to get his stuff, rSP is 0, rSP_hflc = 0.15
        rSP_core = EffFix ? 0.915 * 1.100 : 0.915;
        rSP_core = 1.12;
        rSP_dispZ = rSP_core;
        if (EffFix == -1) {
            rSP_core = 0.2;
            rSP_dispZ = 0.2;
            // rSP_hadr = 2.7;
            // rSP_tau = 1.7;
            // rSP_hadr = 3.55;
            // rSP_tau = 2.3;
            rSP_hadr = 2.3;
            rSP_tau = 3.55;
            rSP_hflc = 0.00;
        } else if (EffFix == -2) {
            rSP_core = 0.25;
            rSP_dispZ = 0.25;
            rSP_hadr = 3.55;
            rSP_tau = 2.3;
            rSP_hflc = 0.1;
        } else if (EffFix == -3) {
            // this works well
            // rSP_core = 0.35;
            // rSP_hadr = 4.50;
            // rSP_tau = 1.0;
            // rSP_hflc = 0.0;
            // HadronSize = 0.8;
            // HadronSlope = 0.2;
            // what works is all zero, but this is 0.25 and the tau is 3.10 exponential
            // what kind of works (but not perfect) is to have tau of c.a. 4.5, 15% fluctuations
            // the former leads to exp source (also the core), while for the latter it is a Gauss
            rSP_core = 0.0;
            rSP_dispZ = 0.0;
            rSP_hadr = 0.0;
            rSP_tau = 4.5;
            rSP_tflc = 0.15;
            tau_prp = false;
            rSP_hflc = 0.0;
            HadronSize = 0.0;
            HadronSlope = 0.0;
            rSP_ThK = 0.0;
            // small increase of scaling at large mT values
            rSP_FixedHadr = true;
            // can also increase the radius.
            // Big effect on Tau, where it creates strong scaling (at low mT values)
            rSP_FragBeta = 1.0;

            // rSP_core = 0.25;
            // rSP_hadr = 3.55;
            // rSP_tau = 1.0;
        }
        // some random tests
        else if (EffFix == -900) {
            rSP_core = 0.32;
            rSP_dispZ = 0.32;
            rSP_hadr = 2.3 / 2.;
            rSP_hadrZ = 2.3 * 4.;
            rSP_tau = 3.55 * 0;
            rSP_hflc = 0.00;
        }
        // for FEMTUM
        else if (EffFix == -1001) {
            // rSP_core = 0.80;
            // rSP_dispZ = 0.80;
            // these values are for the ceca paper, work well with 1 fm effective
            rSP_core = 0.85;
            rSP_dispZ = 0.85;
            rSP_hadr = 0.00;
            rSP_tau = 0.00;
        } else if (EffFix == -1002) {
            rSP_core = 1.20;
            rSP_dispZ = 1.20;
            rSP_hadr = 0.00;
            rSP_tau = 0.00;
        } else if (EffFix == -1010) {
            rSP_core = 0.00;
            rSP_dispZ = 0.00;
            // rSP_hadr = 4.00;
            // these values are for the ceca paper, work well with 1 fm effective
            rSP_hadr = 4.2;
            rSP_tau = 0.00;
        } else if (EffFix == -1020) {
            rSP_core = 0.00;
            rSP_dispZ = 0.00;
            rSP_hadr = 6.00;
            rSP_tau = 0.00;
        } else if (EffFix == -1100) {
            rSP_core = 0.00;
            rSP_dispZ = 0.00;
            rSP_hadr = 0.00;
            // rSP_tau = 4.00;
            // these values are for the ceca paper, work well with 1 fm effective
            rSP_tau = 4.2;
        } else if (EffFix == -1200) {
            rSP_core = 0.00;
            rSP_dispZ = 0.00;
            rSP_hadr = 0.00;
            rSP_tau = 6.00;
        } else if (EffFix == -1300) {
            rSP_core = 0.25;
            rSP_dispZ = 0.00;
            rSP_hadr = 7.00;
            rSP_tau = 0.00;
            // rSP_FragBeta = 1.0;
            // tau_prp = false;
            rSP_ThK = 100;
        }
        // examples for the animation for hadron 2023
        else if (EffFix == -1400) {
            rSP_core = 0.176;
            rSP_dispZ = 0.176;
            rSP_hadr = 0.00;
            rSP_tau = 0.00;
        } else if (EffFix == -1401) {
            rSP_core = 0.176;
            rSP_dispZ = 0.176;
            rSP_hadr = 2.68;
            rSP_tau = 0.00;
        } else if (EffFix == -1402) {
            // this so happens to be exactly as the usmani fit in the ceca paper
            rSP_core = 0.176;
            rSP_dispZ = 0.176;
            rSP_hadr = 2.68;
            rSP_tau = 3.76;
        } else if (EffFix == -1403) {
            rSP_core = 0.176;
            rSP_dispZ = 0.176;
            rSP_hadr = 2.70;
            rSP_tau = 3.78;
        } else if (EffFix == -1404) {
            rSP_core = 0.176;
            rSP_dispZ = 0.176;
            rSP_hadr = 2.24;
            rSP_tau = 3.13;
        } else if (EffFix == -1405) {
            rSP_core = 0.176;
            rSP_dispZ = 0.176;
            rSP_hadr = 1.96;
            rSP_tau = 2.74;
        } else if (EffFix == -1406) {
            rSP_core = 0.176;
            rSP_dispZ = 0.176;
            rSP_hadr = 1.44;
            rSP_tau = 2.02;
        } else if (EffFix == -1407) {
            rSP_core = 0.176;
            rSP_dispZ = 0.176;
            rSP_hadr = 2.00;
            rSP_tau = 3.78;
        }
        // these values are for the ceca paper -> with NLO
        else if (EffFix == -1408) {
            rSP_core = 0.288;
            rSP_dispZ = 0.288;
            rSP_hadr = 3.23;
            rSP_tau = 3.26;
        } else if (EffFix == 9000) {
            rSP_core = 0.85;
            rSP_dispZ = rSP_core;
            rSP_hadr = 0;
            rSP_tau = 0;
        } else if (EffFix == 9001 || EffFix == 9002) {
            rSP_core = 1.;
            rSP_dispZ = rSP_core;
            rSP_hadr = 0;
            rSP_tau = 0;
        }
    }

    TREPNI Database(0);
    Database.SetSeed(11);
    std::vector<TreParticle*> ParticleList;
    ParticleList.push_back(Database.NewParticle("Proton"));
    ParticleList.push_back(Database.NewParticle("PrimProton"));
    ParticleList.push_back(Database.NewParticle("Lambda"));
    ParticleList.push_back(Database.NewParticle("Kaon"));
    ParticleList.push_back(Database.NewParticle("Deuteron"));
    ParticleList.push_back(Database.NewParticle("Pion"));
    ParticleList.push_back(Database.NewParticle("PionFSI"));
    ParticleList.push_back(Database.NewParticle("Photon"));
    ParticleList.push_back(Database.NewParticle("Dch"));
    ParticleList.push_back(Database.NewParticle("Dch_star"));
    ParticleList.push_back(Database.NewParticle("PionReso"));
    ParticleList.push_back(Database.NewParticle("ProtonReso"));
    ParticleList.push_back(Database.NewParticle("LambdaReso"));
    ParticleList.push_back(Database.NewParticle("KaonReso"));
    ParticleList.push_back(Database.NewParticle("DeuteronReso"));

    TFile fOutput(oFileName.data(), "recreate");

    DLM_Histo<float>* dlm_pT_p = NULL;
    DLM_Histo<float>* dlm_pT_d = NULL;
    DLM_Histo<float> dlm_pT_eta_p;
    DLM_Histo<float> dlm_pT_eta_d;
    TH1F* h_pT_p_all = NULL;
    TH1F* h_pT_d_all = NULL;
    TH1F* h_pT_p = NULL;
    TH1F* h_pT_ap = NULL;
    TH1F* h_pT_d = NULL;
    TH1F* h_pT_ad = NULL;

    TString FilePath;
    if (system == "pp") {
        FilePath = TString::Format("~/ph/sw/yaffa/input/ptshapes/", "");
        TString FileNameP1, FileNameAP1, FileNameP2, FileNameAP2;
        FileNameP1 = "p_HMpp13TeV_reco.root";
        FileNameAP1 = "pbar_HMpp13TeV_reco.root";
        FileNameP2 = "p_HMpp13TeV_reco.root";
        FileNameAP2 = "pbar_HMpp13TeV_reco.root";

        TFile file_p(FilePath + FileNameP1, "read");
        h_pT_p = (TH1F*)file_p.Get("pTDist_after");
        if (!h_pT_p) printf("ISSUE with h_pT_p\n");
        fOutput.cd();
        h_pT_p_all = (TH1F*)h_pT_p->Clone("h_pT_p_all");

        TFile file_ap(FilePath + FileNameAP1, "read");
        h_pT_ap = (TH1F*)file_ap.Get("pTDist_after");
        if (!h_pT_ap) printf("ISSUE with h_pT_ap\n");
        h_pT_p_all->Add(h_pT_ap);

        TFile file_d(FilePath + FileNameP2, "read");
        h_pT_d = (TH1F*)file_d.Get("pTDist_after");
        if (!h_pT_d) printf("ISSUE with h_pT_d\n");
        fOutput.cd();
        h_pT_d_all = (TH1F*)h_pT_d->Clone("h_pT_d_all");

        TFile file_ad(FilePath + FileNameAP2, "read");
        h_pT_ad = (TH1F*)file_ad.Get("pTDist_after");
        if (!h_pT_ad) printf("ISSUE with h_pT_ad\n");
        h_pT_d_all->Add(h_pT_ad);
    }
    // N.B. here we set the Kaons to the proton histo

    if (h_pT_p_all) {
        dlm_pT_p = Convert_TH1F_DlmHisto(h_pT_p_all);
        dlm_pT_p->RescaleAxis(0, 1000, false);
    }
    if (h_pT_d_all) {
        dlm_pT_d = Convert_TH1F_DlmHisto(h_pT_d_all);
        dlm_pT_d->RescaleAxis(0, 1000, false);
    }

    double* BinRange = NULL;
    double axis[2];

    if (dlm_pT_p) {
        dlm_pT_eta_p.SetUp(2);
        BinRange = dlm_pT_p->GetBinRange(0);
        dlm_pT_eta_p.SetUp(0, dlm_pT_p->GetNbins(), BinRange);
        delete[] BinRange;
        dlm_pT_eta_p.SetUp(1, 1, -EtaCut, EtaCut);
        dlm_pT_eta_p.Initialize();
        for (unsigned uBin = 0; uBin < dlm_pT_p->GetNbins(); uBin++) {
            dlm_pT_eta_p.SetBinContent(uBin, 0, dlm_pT_p->GetBinContent(uBin));
            // printf("b%u %.3e\n",uBin,dlm_pT_p->GetBinContent(uBin));
        }
    }

    if (dlm_pT_d) {
        dlm_pT_eta_d.SetUp(2);
        BinRange = dlm_pT_d->GetBinRange(0);
        dlm_pT_eta_d.SetUp(0, dlm_pT_d->GetNbins(), BinRange);
        delete[] BinRange;
        dlm_pT_eta_d.SetUp(1, 1, -EtaCut, EtaCut);
        dlm_pT_eta_d.Initialize();
        for (unsigned uBin = 0; uBin < dlm_pT_d->GetNbins(); uBin++) {
            dlm_pT_eta_d.SetBinContent(uBin, 0, dlm_pT_d->GetBinContent(uBin));
        }
    }

    TH2F* hSampleQA_p = new TH2F("hSampleQA_p", "hSampleQA_p", 64, 0, 5000, 64, -1, 1);
    DLM_Random* dlmRandom = new DLM_Random(1);
    if (dlm_pT_p) {
        for (unsigned uIter = 0; uIter < 100 * 1000; uIter++) {
            double axisValues[2];
            dlm_pT_eta_p.SampleYield(axisValues, false, dlmRandom);
            hSampleQA_p->Fill(axisValues[0], axisValues[1]);
        }
    }

    fOutput.cd();
    hSampleQA_p->Write();
    if (h_pT_p_all) h_pT_p_all->Write();
    if (h_pT_d_all) h_pT_d_all->Write();

    for (TreParticle* prt : ParticleList) {
        if (prt->GetName() == "Proton" || prt->GetName() == "PrimProton") {
            // if(prt->GetName()=="PrimProton")prt->SetMass(2.*Mass_p);
            // else prt->SetMass(Mass_p);
            if (removeBoost) {
                prt->SetMass(1000000);
            } else {
                prt->SetMass(Mass_p);
            }
            if (system == "pp")
                prt->SetAbundance(frac_protons + (100. - frac_protons) * (!PROTON_RESO));
            else if (system == "pP")
                prt->SetAbundance(100.);
            else
                prt->SetAbundance(0);
            prt->SetRadius(HadronSize);
            prt->SetRadiusSlope(HadronSlope);
            if (dlm_pT_p)
                prt->SetPtEtaPhi(dlm_pT_eta_p);
            else
                prt->SetPtPz(0.85 * prt->GetMass(), 0.85 * prt->GetMass());
        } else if (prt->GetName() == "ProtonReso") {
            if (removeBoost) {
                prt->SetMass(1000000 * 1362. / 0.938);
            } else {
                prt->SetMass(1362);
            }
            if (system == "pp" || system == "pP")
                prt->SetAbundance((100. - frac_protons) * PROTON_RESO);
            else
                prt->SetAbundance(0);
            prt->SetWidth(hbarc / 1.65);
            prt->SetRadius(HadronSize);
            prt->SetRadiusSlope(HadronSlope);
            if (dlm_pT_p)
                prt->SetPtEtaPhi(dlm_pT_eta_p);
            else
                prt->SetPtPz(0.85 * prt->GetMass(), 0.85 * prt->GetMass());

            prt->NewDecay();
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
            prt->GetDecay(0)->SetBranching(100);
        }
    }

    std::vector<std::string> ListOfParticles;
    if (system == "pp") {
        ListOfParticles.push_back("Proton");
        ListOfParticles.push_back("Proton");
    }
    if (system == "pP") {
        ListOfParticles.push_back("Proton");
        ListOfParticles.push_back("PrimProton");
    }

    CECA ceca(Database, ListOfParticles);
    for (int iThread = 0; iThread < NUM_CPU; iThread++) {
        ceca.SetSeed(iThread, iThread + 1);
    }

    printf("d %.3f %.3f\n", rSP_core, rSP_dispZ);
    printf("h %.3f %.3f\n", rSP_hadr, rSP_hadrZ);
    printf("t %.3f %.3f\n", rSP_tau, float(tau_prp));
    printf("tf %.3f\n", rSP_tflc);
    printf("tk %.3f\n", rSP_ThK);
    printf("fh %.3f\n", float(rSP_FixedHadr));
    printf("fb %.3f\n", rSP_FragBeta);

    ceca.SetDisplacementZ(rSP_dispZ);
    ceca.SetDisplacementT(rSP_core);
    ceca.SetHadronizationZ(rSP_hadrZ);
    ceca.SetHadronizationT(rSP_hadr);
    ceca.SetHadrFluctuation(rSP_hflc);
    ceca.SetTau(rSP_tau, tau_prp);
    ceca.SetTauFluct(rSP_tflc);
    ceca.SetThermalKick(rSP_ThK);
    ceca.SetFixedHadr(rSP_FixedHadr);
    ceca.SetFragmentBeta(rSP_FragBeta);

    // ceca.SetDisplacement(3.55*0.25);
    // ceca.SetHadronizationT(3.55);
    // ceca.SetHadrFluctuation(0.0);
    // ceca.SetTau(0);

    // ceca.SetDisplacement(1);
    // ceca.SetHadronizationT(0);
    // ceca.SetTau(0);

    ceca.SetTargetStatistics(target_yield);
    ceca.SetEventMult(Multiplicity);
    ceca.SetSourceDim(2);
    ceca.SetDebugMode(true);
    ceca.SetThreadTimeout(threadTimeout);
    ceca.SetGlobalTimeout(globalTimeout);
    ceca.EqualizeFsiTime(EQUALIZE_TAU);
    ceca.SetFemtoRegion(femto_region);
    ceca.GHETTO_EVENT = true;
    if (system == "pp") {
        // ceca paper
        ceca.Ghetto_NumMtBins = 10;
        ceca.Ghetto_MtBins = new double[ceca.Ghetto_NumMtBins + 1];
        ceca.Ghetto_MtBins[0] = 930;   // avg  983 ( 985)
        ceca.Ghetto_MtBins[1] = 1020;  // avg 1054 (1055)
        ceca.Ghetto_MtBins[2] = 1080;  // avg 1110 (1110)
        ceca.Ghetto_MtBins[3] = 1140;  // avg 1168 (1170)
        ceca.Ghetto_MtBins[4] = 1200;  // avg 1228 (1230)
        ceca.Ghetto_MtBins[5] = 1260;  // avg 1315 (1315)
        ceca.Ghetto_MtBins[6] = 1380;  // avg 1463 (1460)
        ceca.Ghetto_MtBins[7] = 1570;  // avg 1681 (1680)
        ceca.Ghetto_MtBins[8] = 1840;  // avg 1923 (1920)
        ceca.Ghetto_MtBins[9] = 2030;  // avg 2303 (2300)
        ceca.Ghetto_MtBins[10] = 4500;

        ceca.Ghetto_NumMomBins = 150;
        ceca.Ghetto_MomMin = 0;
        ceca.Ghetto_MomMax = 600;

        ceca.Ghetto_NumRadBins = 2048;
        ceca.Ghetto_RadMin = 0;
        ceca.Ghetto_RadMax = 64;
    }

    ceca.GoBabyGo(NUM_CPU);

    // ceca.Ghetto_kstar_rstar_mT->QuickWrite(BaseFileName + ".Ghetto_kstar_rstar_mT", true);

    // return;
    double TotPairs = ceca.GhettoPrimReso[0] + ceca.GhettoPrimReso[1] + ceca.GhettoPrimReso[2] + ceca.GhettoPrimReso[3];
    double TotPP = double(ceca.GhettoPrimReso[0]) / TotPairs;
    double TotPR = double(ceca.GhettoPrimReso[1]) / TotPairs;
    double TotRP = double(ceca.GhettoPrimReso[2]) / TotPairs;
    double TotRR = double(ceca.GhettoPrimReso[3]) / TotPairs;

    double FemtoPairs = ceca.GhettoFemtoPrimReso[0] + ceca.GhettoFemtoPrimReso[1] + ceca.GhettoFemtoPrimReso[2] +
                        ceca.GhettoFemtoPrimReso[3];
    double FemtoPP = double(ceca.GhettoFemtoPrimReso[0]) / FemtoPairs;
    double FemtoPR = double(ceca.GhettoFemtoPrimReso[1]) / FemtoPairs;
    double FemtoRP = double(ceca.GhettoFemtoPrimReso[2]) / FemtoPairs;
    double FemtoRR = double(ceca.GhettoFemtoPrimReso[3]) / FemtoPairs;

    printf("     Total  Femto\n");
    printf("PP%6.2f%% %6.2f\n", TotPP * 100., FemtoPP * 100.);
    printf("PR%6.2f%% %6.2f\n", TotPR * 100., FemtoPR * 100.);
    printf("RP%6.2f%% %6.2f\n", TotRP * 100., FemtoRP * 100.);
    printf("RR%6.2f%% %6.2f\n", TotRR * 100., FemtoRR * 100.);

    ceca.GhettoFemto_rstar->ComputeError();
    ceca.GhettoFemto_rstar->ScaleToIntegral();
    ceca.GhettoFemto_rstar->ScaleToBinSize();

    CATS cat;
    cat.SetMomBins(80, 0, 320);
    cat.SetNotifications(CATS::nWarning);
    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles", "").Data());
    AnalysisObject.SetUpCats_pp(cat, "AV18", "NULL", 0, 0);
    DLM_HistoSource ppHistoSource(*ceca.GhettoFemto_rstar);
    cat.SetAnaSource(CatsSourceForwarder, &ppHistoSource, 0);
    cat.SetUseAnalyticSource(true);
    cat.SetAutoNormSource(false);
    cat.SetNormalizedSource(true);
    cat.KillTheCat();
    TGraph Ck_pp;
    Ck_pp.SetName("Ck_pp");
    Ck_pp.SetLineColor(kBlue);
    Ck_pp.SetLineWidth(5);
    for (unsigned uBin = 0; uBin < cat.GetNumMomBins(); uBin++) {
        Ck_pp.SetPoint(uBin, cat.GetMomentum(uBin), cat.GetCorrFun(uBin));
    }
    fOutput.cd();
    Ck_pp.Write();

    TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(ceca.GhettoFemto_rstar, "GhettoFemto_rstar");
    fOutput.cd();
    h_GhettoFemto_rstar->SetLineWidth(3);
    h_GhettoFemto_rstar->SetLineColor(kAzure);
    h_GhettoFemto_rstar->Write();
    TF1* fit_rstar = new TF1(
        "fit_rstar", "[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]", 0, 16);
    fit_rstar->SetLineColor(kAzure);
    fit_rstar->SetLineWidth(3);
    // plot the fit to the total dist as an example
    // below we find the same for each mT, there we will just get out the reff_Ceca
    double reff_Ceca = Get_reff(h_GhettoFemto_rstar);
    printf("reff_Ceca = %f\n", reff_Ceca);
    fOutput.cd();
    fit_rstar->FixParameter(0, 1);
    fit_rstar->FixParameter(1, reff_Ceca);
    fit_rstar->Write();

    ceca.GhettoFemto_rcore->ComputeError();
    ceca.GhettoFemto_rcore->ScaleToIntegral();
    ceca.GhettoFemto_rcore->ScaleToBinSize();
    TH1F* h_GhettoFemto_rcore = Convert_DlmHisto_TH1F(ceca.GhettoFemto_rcore, "GhettoFemto_rcore");
    fOutput.cd();
    h_GhettoFemto_rcore->SetLineWidth(3);
    h_GhettoFemto_rcore->SetLineColor(kBlack);
    h_GhettoFemto_rcore->Write();
    TF1* fit_rcore = new TF1(
        "fit_rcore", "[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]", 0, 16);
    fit_rcore->SetLineColor(kBlack);
    fit_rcore->SetLineWidth(3);
    double rcore_Ceca = Get_reff(h_GhettoFemto_rcore);
    printf("rcore_Ceca = %f\n", rcore_Ceca);
    fOutput.cd();
    fit_rcore->FixParameter(0, 1);
    fit_rcore->FixParameter(1, rcore_Ceca);
    fit_rcore->Write();

    ceca.Ghetto_kstar->ComputeError();
    ceca.Ghetto_kstar->ScaleToIntegral();
    ceca.Ghetto_kstar->ScaleToBinSize();
    TH1F* h_Ghetto_kstar = Convert_DlmHisto_TH1F(ceca.Ghetto_kstar, "Ghetto_kstar");
    ceca.Ghetto_kstar_rstar->ComputeError();
    TH2F* h_Ghetto_kstar_rstar = Convert_DlmHisto_TH2F(ceca.Ghetto_kstar_rstar, "Ghetto_kstar_rstar");

    ceca.Ghetto_kstar_rstar_PP->ComputeError();
    TH2F* h_Ghetto_kstar_rstar_PP = Convert_DlmHisto_TH2F(ceca.Ghetto_kstar_rstar_PP, "Ghetto_kstar_rstar_PP");
    ceca.Ghetto_kstar_rstar_PR->ComputeError();
    TH2F* h_Ghetto_kstar_rstar_PR = Convert_DlmHisto_TH2F(ceca.Ghetto_kstar_rstar_PR, "Ghetto_kstar_rstar_PR");
    ceca.Ghetto_kstar_rstar_RP->ComputeError();
    TH2F* h_Ghetto_kstar_rstar_RP = Convert_DlmHisto_TH2F(ceca.Ghetto_kstar_rstar_RP, "Ghetto_kstar_rstar_RP");
    ceca.Ghetto_kstar_rstar_RR->ComputeError();
    TH2F* h_Ghetto_kstar_rstar_RR = Convert_DlmHisto_TH2F(ceca.Ghetto_kstar_rstar_RR, "Ghetto_kstar_rstar_RR");

    TH2F* h_Ghetto_mT_rstar = Convert_DlmHisto_TH2F(ceca.Ghetto_mT_rstar, "Ghetto_mT_rstar");
    ceca.GhettoFemto_mT_rstar->ComputeError();
    TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(ceca.GhettoFemto_mT_rstar, "GhettoFemto_mT_rstar");

    // fOutput.cd();
    // h_GhettoFemto_mT_rstar->Write();
    ceca.GhettoFemto_mT_rcore->ComputeError();
    TH2F* h_GhettoFemto_mT_rcore = Convert_DlmHisto_TH2F(ceca.GhettoFemto_mT_rcore, "GhettoFemto_mT_rcore");
    TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(ceca.GhettoFemto_mT_kstar, "GhettoFemto_mT_kstar");
    TH2F* h_Ghetto_mT_costh = Convert_DlmHisto_TH2F(ceca.Ghetto_mT_costh, "Ghetto_mT_costh");
    TH2F* h_GhettoSP_pT_th = Convert_DlmHisto_TH2F(ceca.GhettoSP_pT_th, "GhettoSP_pT_th");
    TH1F* h_GhettoSP_pT_1 = Convert_DlmHisto_TH1F(ceca.GhettoSP_pT_1, "GhettoSP_pT_1");
    TH1F* h_GhettoSP_pT_2 = Convert_DlmHisto_TH1F(ceca.GhettoSP_pT_2, "GhettoSP_pT_2");
    // ceca.GhettoFemto_rstar->ComputeError();
    // TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(ceca.GhettoFemto_rstar,"GhettoFemto_rstar");
    ceca.Ghetto_PP_AngleRcP1->ComputeError();
    ceca.Ghetto_PP_AngleRcP2->ComputeError();
    ceca.Ghetto_PP_AngleP1P2->ComputeError();
    ceca.Ghetto_RP_AngleRcP1->ComputeError();
    ceca.Ghetto_PR_AngleRcP2->ComputeError();
    ceca.Ghetto_RR_AngleRcP1->ComputeError();
    ceca.Ghetto_RR_AngleRcP2->ComputeError();
    ceca.Ghetto_RR_AngleP1P2->ComputeError();
    TH1F* h_Ghetto_PP_AngleRcP1 = Convert_DlmHisto_TH1F(ceca.Ghetto_PP_AngleRcP1, "Ghetto_PP_AngleRcP1");
    TH1F* h_Ghetto_PP_AngleRcP2 = Convert_DlmHisto_TH1F(ceca.Ghetto_PP_AngleRcP2, "Ghetto_PP_AngleRcP2");
    TH1F* h_Ghetto_PP_AngleP1P2 = Convert_DlmHisto_TH1F(ceca.Ghetto_PP_AngleP1P2, "Ghetto_PP_AngleP1P2");
    TH1F* h_Ghetto_RP_AngleRcP1 = Convert_DlmHisto_TH1F(ceca.Ghetto_RP_AngleRcP1, "Ghetto_RP_AngleRcP1");
    TH1F* h_Ghetto_PR_AngleRcP2 = Convert_DlmHisto_TH1F(ceca.Ghetto_PR_AngleRcP2, "Ghetto_PR_AngleRcP2");
    TH1F* h_Ghetto_RR_AngleRcP1 = Convert_DlmHisto_TH1F(ceca.Ghetto_RR_AngleRcP1, "Ghetto_RR_AngleRcP1");
    TH1F* h_Ghetto_RR_AngleRcP2 = Convert_DlmHisto_TH1F(ceca.Ghetto_RR_AngleRcP2, "Ghetto_RR_AngleRcP2");
    TH1F* h_Ghetto_RR_AngleP1P2 = Convert_DlmHisto_TH1F(ceca.Ghetto_RR_AngleP1P2, "Ghetto_RR_AngleP1P2");

    ceca.GhettoSPr_X->ComputeError();
    ceca.GhettoSPr_X->ScaleToIntegral();
    ceca.GhettoSPr_X->ScaleToBinSize();
    TH1F* h_GhettoSPr_X = Convert_DlmHisto_TH1F(ceca.GhettoSPr_X, "GhettoSPr_X");

    ceca.GhettoSPr_Y->ComputeError();
    ceca.GhettoSPr_Y->ScaleToIntegral();
    ceca.GhettoSPr_Y->ScaleToBinSize();
    TH1F* h_GhettoSPr_Y = Convert_DlmHisto_TH1F(ceca.GhettoSPr_Y, "GhettoSPr_Y");

    ceca.GhettoSPr_Z->ComputeError();
    ceca.GhettoSPr_Z->ScaleToIntegral();
    ceca.GhettoSPr_Z->ScaleToBinSize();
    TH1F* h_GhettoSPr_Z = Convert_DlmHisto_TH1F(ceca.GhettoSPr_Z, "GhettoSPr_Z");

    ceca.GhettoSPr_Rho->ComputeError();
    ceca.GhettoSPr_Rho->ScaleToIntegral();
    ceca.GhettoSPr_Rho->ScaleToBinSize();
    TH1F* h_GhettoSPr_Rho = Convert_DlmHisto_TH1F(ceca.GhettoSPr_Rho, "GhettoSPr_Rho");

    ceca.GhettoSPr_R->ComputeError();
    ceca.GhettoSPr_R->ScaleToIntegral();
    ceca.GhettoSPr_R->ScaleToBinSize();
    TH1F* h_GhettoSPr_R = Convert_DlmHisto_TH1F(ceca.GhettoSPr_R, "GhettoSPr_R");

    ceca.GhettoSP_X->ComputeError();
    ceca.GhettoSP_X->ScaleToIntegral();
    ceca.GhettoSP_X->ScaleToBinSize();
    TH1F* h_GhettoSP_X = Convert_DlmHisto_TH1F(ceca.GhettoSP_X, "GhettoSP_X");

    ceca.GhettoSP_Y->ComputeError();
    ceca.GhettoSP_Y->ScaleToIntegral();
    ceca.GhettoSP_Y->ScaleToBinSize();
    TH1F* h_GhettoSP_Y = Convert_DlmHisto_TH1F(ceca.GhettoSP_Y, "GhettoSP_Y");

    ceca.GhettoSP_Z->ComputeError();
    ceca.GhettoSP_Z->ScaleToIntegral();
    ceca.GhettoSP_Z->ScaleToBinSize();
    TH1F* h_GhettoSP_Z = Convert_DlmHisto_TH1F(ceca.GhettoSP_Z, "GhettoSP_Z");

    ceca.GhettoSP_Rho->ComputeError();
    ceca.GhettoSP_Rho->ScaleToIntegral();
    ceca.GhettoSP_Rho->ScaleToBinSize();
    TH1F* h_GhettoSP_Rho = Convert_DlmHisto_TH1F(ceca.GhettoSP_Rho, "GhettoSP_Rho");

    ceca.GhettoSP_R->ComputeError();
    ceca.GhettoSP_R->ScaleToIntegral();
    ceca.GhettoSP_R->ScaleToBinSize();
    TH1F* h_GhettoSP_R = Convert_DlmHisto_TH1F(ceca.GhettoSP_R, "GhettoSP_R");

    ceca.GhettoFemto_mT_mTwrong->ComputeError();
    TH2F* h_GhettoFemto_mT_mTwrong = Convert_DlmHisto_TH2F(ceca.GhettoFemto_mT_mTwrong, "GhettoFemto_mT_mTwrong");

    ceca.Ghetto_mT_mTwrong->ComputeError();
    TH2F* h_Ghetto_mT_mTwrong = Convert_DlmHisto_TH2F(ceca.Ghetto_mT_mTwrong, "Ghetto_mT_mTwrong");

    ceca.GhettoFemto_pT1_pT2->ComputeError();
    TH2F* h_GhettoFemto_pT1_pT2 = Convert_DlmHisto_TH2F(ceca.GhettoFemto_pT1_pT2, "GhettoFemto_pT1_pT2");

    ceca.GhettoFemto_pT1_div_pT->ComputeError();
    TH1F* h_GhettoFemto_pT1_div_pT = Convert_DlmHisto_TH1F(ceca.GhettoFemto_pT1_div_pT, "GhettoFemto_pT1_div_pT");

    TGraphErrors g_GhettoFemto_mT_rstar;
    g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
    g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
    g_GhettoFemto_mT_rstar.SetMarkerSize(1);
    g_GhettoFemto_mT_rstar.SetLineWidth(3);

    TGraphErrors g_GhettoFemto_mT_rstar_G;
    g_GhettoFemto_mT_rstar_G.SetName("g_GhettoFemto_mT_rstar_G");
    g_GhettoFemto_mT_rstar_G.SetMarkerStyle(20);
    g_GhettoFemto_mT_rstar_G.SetMarkerSize(1);
    g_GhettoFemto_mT_rstar_G.SetLineWidth(6);
    g_GhettoFemto_mT_rstar_G.SetLineColor(kAzure);

    TGraphErrors g_GhettoFemto_mT_rcore;
    g_GhettoFemto_mT_rcore.SetName("g_GhettoFemto_mT_rcore");
    g_GhettoFemto_mT_rcore.SetMarkerStyle(20);
    g_GhettoFemto_mT_rcore.SetMarkerSize(1);
    g_GhettoFemto_mT_rcore.SetLineWidth(3);

    TGraphErrors g_GhettoFemto_mT_rcore_G;
    g_GhettoFemto_mT_rcore_G.SetName("g_GhettoFemto_mT_rcore_G");
    g_GhettoFemto_mT_rcore_G.SetMarkerStyle(20);
    g_GhettoFemto_mT_rcore_G.SetMarkerSize(1);
    g_GhettoFemto_mT_rcore_G.SetLineWidth(6);
    g_GhettoFemto_mT_rcore_G.SetLineColor(kBlack);

    unsigned uPointRS = 0;
    unsigned uPointRC = 0;
    for (unsigned uBin = 0; uBin < h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++) {
        TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"), uBin + 1, uBin + 1);
        double Mean = hProj->GetMean();
        double Err = hProj->GetStdDev();
        double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin + 1);
        if (Mean && Err && hProj->GetEntries() > 128) {
            hProj->Scale(1. / hProj->Integral(), "width");
            g_GhettoFemto_mT_rstar.SetPoint(uPointRS, mT * 0.001, Mean);
            // g_GhettoFemto_mT_rstar.SetPointError(uPointRS,0,Err);
            g_GhettoFemto_mT_rstar.SetPointError(uPointRS, 0, 0);

            reff_Ceca = Get_reff(hProj);
            // printf("mT%u reff_Ceca=%.3f\n",uBin,reff_Ceca);
            g_GhettoFemto_mT_rstar_G.SetPoint(uPointRS, mT * 0.001, reff_Ceca);
            g_GhettoFemto_mT_rstar_G.SetPointError(uPointRS, 0, 0);

            uPointRS++;
        }
        delete hProj;

        hProj = (TH1F*)h_GhettoFemto_mT_rcore->ProjectionY(TString::Format("hProj"), uBin + 1, uBin + 1);
        Mean = hProj->GetMean();
        Err = hProj->GetStdDev();
        mT = h_GhettoFemto_mT_rcore->GetXaxis()->GetBinCenter(uBin + 1);
        if (Mean && Err && hProj->GetEntries() > 256) {
            hProj->Scale(1. / hProj->Integral(), "width");
            g_GhettoFemto_mT_rcore.SetPoint(uPointRC, mT * 0.001, Mean);
            g_GhettoFemto_mT_rcore.SetPointError(uPointRC, 0, 0);

            rcore_Ceca = Get_reff(hProj);
            // printf("    rcore_Ceca=%.3f\n",rcore_Ceca);
            g_GhettoFemto_mT_rcore_G.SetPoint(uPointRC, mT * 0.001, rcore_Ceca);
            g_GhettoFemto_mT_rcore_G.SetPointError(uPointRC, 0, 0);

            uPointRC++;
        }
        delete hProj;
    }

    const double MeanToGauss = 1. / 2.3;
    double mT_first = g_GhettoFemto_mT_rstar.GetPointY(0);
    printf("r(%.2f GeV) = %.2f\n", g_GhettoFemto_mT_rstar.GetPointX(0) * 0.001, mT_first * MeanToGauss);
    double mT_2GeV = g_GhettoFemto_mT_rstar.Eval(2000);
    printf("r(%.2f GeV) = %.2f\n", 2., mT_2GeV * MeanToGauss);

    h_GhettoFemto_rstar->Scale(1. / h_GhettoFemto_rstar->Integral(), "width");

    double lowerlimit;
    double upperlimit;
    double alpha = cfg["alpha"].as<double>();
    GetCentralInterval(*h_GhettoFemto_rstar, alpha, lowerlimit, upperlimit, true);

    unsigned lowerbin = h_GhettoFemto_rstar->FindBin(lowerlimit);
    unsigned upperbin = h_GhettoFemto_rstar->FindBin(upperlimit);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n", lowerlimit, upperlimit);
    TF1* fSource = new TF1("fSource", ScaledGauss, 0, 32, 2);
    fSource->SetParameter(0, h_GhettoFemto_rstar->GetMean() / 2.3);
    fSource->SetParLimits(0, 0.5, 5.0);
    fSource->FixParameter(1, 1.0);

    h_GhettoFemto_rstar->Fit(fSource, "Q, S, R, M +", "", lowerlimit, upperlimit);
    fOutput.cd();
    h_GhettoFemto_rstar->Write("hRStar");
    printf("The effective Gaussian size is %.3f +/- %.3f fm, λ = %.3f +/- %.3f\n", fSource->GetParameter(0),
           fSource->GetParError(0), fSource->GetParameter(1), fSource->GetParError(1));

    GetCentralInterval(*h_GhettoFemto_rstar, 0.9, lowerlimit, upperlimit, true);
    TF1* fitDG_rstar = new TF1("fitDG_rstar", NormDoubleGaussSourceTF1, 0.4, 16.0, 4);
    fitDG_rstar->FixParameter(3, 1);
    fitDG_rstar->FixParameter(0, fit_rcore->GetParameter(1));  // 3.30932e-01

    fitDG_rstar->SetParameter(1, 1.30);      // 1.23015e+00
    fitDG_rstar->SetParLimits(1, 0.9, 5.0);  // second source
    fitDG_rstar->FixParameter(2, FemtoPP);

    h_GhettoFemto_rstar->Fit(fitDG_rstar, "S, N, R, M");
    printf("DG = %.3f x %.3ffm + %.3f x %.3ffm\n", fitDG_rstar->GetParameter(2), fitDG_rstar->GetParameter(0),
           1. - fitDG_rstar->GetParameter(2), fitDG_rstar->GetParameter(1));

    TGraph gRadKstar;
    gRadKstar.SetName("gRadKstar");
    gRadKstar.SetLineColor(kRed + 1);
    gRadKstar.SetLineWidth(3);

    TGraph gMeanRadKstar;
    gMeanRadKstar.SetName("gMeanRadKstar");
    gMeanRadKstar.SetLineColor(kRed + 1);
    gMeanRadKstar.SetLineWidth(3);

    // evaluate the radius based on the mean of the source distribution,
    // setting the actual r value to that of a Gauss source with the same mean
    TGraph gGhettoRadKstar;
    gGhettoRadKstar.SetName("gGhettoRadKstar");
    gGhettoRadKstar.SetLineColor(kRed + 1);
    gGhettoRadKstar.SetLineWidth(3);

    TH1D** hkstar_rstar;
    hkstar_rstar = new TH1D*[h_Ghetto_kstar_rstar->GetXaxis()->GetNbins()];

    for (unsigned uMom = 0; uMom < h_Ghetto_kstar_rstar->GetXaxis()->GetNbins(); uMom++) {
        double kstar = h_Ghetto_kstar_rstar->GetXaxis()->GetBinCenter(uMom + 1);
        hkstar_rstar[uMom] = NULL;
        if (kstar > 780) continue;
        hkstar_rstar[uMom] =
            h_Ghetto_kstar_rstar->ProjectionY(TString::Format("hkstar_rstar_%.0f", kstar), uMom + 1, uMom + 1);
        hkstar_rstar[uMom]->Scale(1. / hkstar_rstar[uMom]->Integral(), "width");
        GetCentralInterval(*hkstar_rstar[uMom], 0.9, lowerlimit, upperlimit, true);
        fSource->SetParameter(0, hkstar_rstar[uMom]->GetMean() / 2.3);
        fSource->SetParLimits(0, hkstar_rstar[uMom]->GetMean() / 4., hkstar_rstar[uMom]->GetMean());
        hkstar_rstar[uMom]->Fit(fSource, "Q, S, N, R, M", "", lowerlimit, upperlimit);
        gRadKstar.SetPoint(uMom, kstar, fSource->GetParameter(0));
        gMeanRadKstar.SetPoint(uMom, kstar, hkstar_rstar[uMom]->GetMean());
        double gfm = GaussFromMean(hkstar_rstar[uMom]->GetMean());
        gGhettoRadKstar.SetPoint(uMom, kstar, gfm);
    }

    h_Ghetto_PP_AngleRcP1->Scale(1. / h_Ghetto_PP_AngleRcP1->Integral(), "width");
    h_Ghetto_PP_AngleRcP2->Scale(1. / h_Ghetto_PP_AngleRcP2->Integral(), "width");
    h_Ghetto_PP_AngleP1P2->Scale(1. / h_Ghetto_PP_AngleP1P2->Integral(), "width");
    h_Ghetto_RP_AngleRcP1->Scale(1. / h_Ghetto_RP_AngleRcP1->Integral(), "width");
    h_Ghetto_PR_AngleRcP2->Scale(1. / h_Ghetto_PR_AngleRcP2->Integral(), "width");
    h_Ghetto_RR_AngleRcP1->Scale(1. / h_Ghetto_RR_AngleRcP1->Integral(), "width");
    h_Ghetto_RR_AngleRcP2->Scale(1. / h_Ghetto_RR_AngleRcP2->Integral(), "width");
    h_Ghetto_RR_AngleP1P2->Scale(1. / h_Ghetto_RR_AngleP1P2->Integral(), "width");

    fOutput.cd();
    // h_Ghetto_rstar->Write();
    // fit_rstar->Write();
    // h_Ghetto_rcore->Write();
    // fit_rcore->Write();

    TH1F* hAxisSource = new TH1F("hAxisSource", "hAxisSource", 128, 0, 8);
    hAxisSource->SetStats(false);
    hAxisSource->SetTitle("");
    hAxisSource->GetXaxis()->SetTitle("r (fm)");
    hAxisSource->GetXaxis()->SetTitleSize(0.06);
    hAxisSource->GetXaxis()->SetLabelSize(0.06);
    hAxisSource->GetXaxis()->SetTitleOffset(1.3);
    hAxisSource->GetXaxis()->SetLabelOffset(0.02);
    hAxisSource->GetYaxis()->SetTitle("4#pir*^{2}S(r) (1/fm)");
    hAxisSource->GetYaxis()->SetTitleSize(0.06);
    hAxisSource->GetYaxis()->SetLabelSize(0.06);
    hAxisSource->GetYaxis()->SetTitleOffset(1.00);
    hAxisSource->GetYaxis()->SetRangeUser(0.001, 1.50);
    hAxisSource->Write();

    TH1F* hAxisMt = new TH1F("hAxisMt", "hAxisMt", 128, 1.05, 2.25);
    hAxisMt->SetStats(false);
    hAxisMt->SetTitle("");
    hAxisMt->GetXaxis()->SetTitle("m_{T} (GeV)");
    hAxisMt->GetXaxis()->SetTitleSize(0.06);
    hAxisMt->GetXaxis()->SetLabelSize(0.06);
    hAxisMt->GetXaxis()->SetTitleOffset(1.3);
    hAxisMt->GetXaxis()->SetLabelOffset(0.02);
    hAxisMt->GetYaxis()->SetTitle("r_{Gauss} (fm)");
    hAxisMt->GetYaxis()->SetTitleSize(0.06);
    hAxisMt->GetYaxis()->SetLabelSize(0.06);
    hAxisMt->GetYaxis()->SetTitleOffset(1.00);
    hAxisMt->GetYaxis()->SetRangeUser(0.4, 1.85);
    hAxisMt->Write();

    h_Ghetto_mT_mTwrong->Write();
    h_GhettoFemto_mT_mTwrong->Write();
    h_GhettoFemto_pT1_pT2->Write();
    h_GhettoFemto_pT1_div_pT->Write();

    h_Ghetto_kstar_rstar->GetXaxis()->SetRangeUser(0, 1200);
    h_Ghetto_kstar_rstar_PP->GetXaxis()->SetRangeUser(0, 1200);
    h_Ghetto_kstar_rstar_PR->GetXaxis()->SetRangeUser(0, 1200);
    h_Ghetto_kstar_rstar_RP->GetXaxis()->SetRangeUser(0, 1200);
    h_Ghetto_kstar_rstar_RR->GetXaxis()->SetRangeUser(0, 1200);

    double UpRU = 100;
    h_Ghetto_kstar_rstar->GetYaxis()->SetRangeUser(0, UpRU);
    h_Ghetto_kstar_rstar_PP->GetYaxis()->SetRangeUser(0, UpRU);
    h_Ghetto_kstar_rstar_PR->GetYaxis()->SetRangeUser(0, UpRU);
    h_Ghetto_kstar_rstar_RP->GetYaxis()->SetRangeUser(0, UpRU);
    h_Ghetto_kstar_rstar_RR->GetYaxis()->SetRangeUser(0, UpRU);

    h_Ghetto_kstar_rstar->GetXaxis()->SetTitle("k* (MeV)");
    h_Ghetto_kstar_rstar_PP->GetXaxis()->SetTitle("k* (MeV)");
    h_Ghetto_kstar_rstar_PR->GetXaxis()->SetTitle("k* (MeV)");
    h_Ghetto_kstar_rstar_RP->GetXaxis()->SetTitle("k* (MeV)");
    h_Ghetto_kstar_rstar_RR->GetXaxis()->SetTitle("k* (MeV)");

    h_Ghetto_kstar_rstar->GetYaxis()->SetTitle("r* (fm)");
    h_Ghetto_kstar_rstar_PP->GetYaxis()->SetTitle("r* (fm)");
    h_Ghetto_kstar_rstar_PR->GetYaxis()->SetTitle("r* (fm)");
    h_Ghetto_kstar_rstar_RP->GetYaxis()->SetTitle("r* (fm)");
    h_Ghetto_kstar_rstar_RR->GetYaxis()->SetTitle("r* (fm)");

    TCanvas* cSource = new TCanvas("cSource", "cSource", 1);
    cSource->cd(0);
    cSource->SetCanvasSize(1280, 720);
    cSource->SetMargin(0.15, 0.05, 0.2, 0.05);  // lrbt
    cSource->SetLogy(true);
    hAxisSource->Draw("axis");
    h_GhettoFemto_rcore->Draw("same");
    h_GhettoFemto_rstar->Draw("same");
    fit_rcore->Draw("same");
    fit_rstar->Draw("same");

    TCanvas* cMt = new TCanvas("cMt", "cMt", 1);
    cMt->cd(0);
    cMt->SetCanvasSize(1280, 720);
    cMt->SetMargin(0.15, 0.05, 0.2, 0.05);  // lrbt
    hAxisMt->Draw("axis");
    g_GhettoFemto_mT_rcore_G.Draw("same");
    g_GhettoFemto_mT_rstar_G.Draw("same");

    fOutput.cd();
    cSource->Write();
    cMt->Write();
    // h_GhettoFemto_rstar->Write();
    fSource->Write();
    fitDG_rstar->Write();
    h_Ghetto_kstar->Write();
    h_Ghetto_kstar_rstar->Write();
    h_Ghetto_kstar_rstar_PP->Write();
    h_Ghetto_kstar_rstar_PR->Write();
    h_Ghetto_kstar_rstar_RP->Write();
    h_Ghetto_kstar_rstar_RR->Write();
    h_Ghetto_mT_rstar->Write();
    h_GhettoFemto_mT_rstar->Write();
    g_GhettoFemto_mT_rstar.Write();
    g_GhettoFemto_mT_rstar_G.Write();
    g_GhettoFemto_mT_rcore.Write();
    g_GhettoFemto_mT_rcore_G.Write();
    h_GhettoFemto_mT_kstar->Write();
    h_Ghetto_mT_costh->Write();
    h_GhettoSP_pT_th->Write();
    h_GhettoSP_pT_1->Write();
    h_GhettoSP_pT_2->Write();
    h_Ghetto_PP_AngleRcP1->Write();
    h_Ghetto_PP_AngleRcP2->Write();
    h_Ghetto_PP_AngleP1P2->Write();
    h_Ghetto_RP_AngleRcP1->Write();
    h_Ghetto_PR_AngleRcP2->Write();
    h_Ghetto_RR_AngleRcP1->Write();
    h_Ghetto_RR_AngleRcP2->Write();
    h_Ghetto_RR_AngleP1P2->Write();
    gRadKstar.Write();
    gMeanRadKstar.Write();
    gGhettoRadKstar.Write();
    h_GhettoSPr_X->Write();
    h_GhettoSPr_Y->Write();
    h_GhettoSPr_Z->Write();
    h_GhettoSPr_Rho->Write();
    h_GhettoSPr_R->Write();
    h_GhettoSP_X->Write();
    h_GhettoSP_Y->Write();
    h_GhettoSP_Z->Write();
    h_GhettoSP_Rho->Write();
    h_GhettoSP_R->Write();
    for (unsigned uMom = 0; uMom < h_Ghetto_kstar_rstar->GetXaxis()->GetNbins(); uMom++) {
        if (hkstar_rstar[uMom]) {
            hkstar_rstar[uMom]->Write();
            delete hkstar_rstar[uMom];
        }
    }

    return 0;
}