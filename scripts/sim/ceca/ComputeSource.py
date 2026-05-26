#!/bin/env python3
from pathlib import Path
import argparse
import sys
from dotenv import load_dotenv
import numpy as np
import os

from ROOT import RDataFrame, TChain, TFile, gROOT, TGraphAsymmErrors, EnableImplicitMT, RDF

env_path = Path(__file__).resolve().parent.parent.parent.parent / ".env"
if not load_dotenv(dotenv_path=env_path):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
sys.path.append(f'{YAFFA_PATH}/src/python')

from yaffa import logger as log
from yaffa import utils

BINNING_MT = (30, 1000, 2500) # um = MeV
BINNING_SOURCE = (200, 0, 20) # um = fm
BINNING_MOMENTUM = (2000, 0, 2000) # um = MeV/c

def DefineVariables(df, m1, m2, m3, arbitraryMass):
    '''
    Define variables for single particles, pairs, and triplets.

    The origin (according to CECA implementation) is specified as a binary number. For single particles
      1: primary
      0: secondary

    Origin for pairs and triplets is stored as an integer that can be accessed via bits: e.g. to check if the second
    particle in a triplet is primary, do
      origin == 0b010
      
    '''
    for idx in [1, 2, 3]:
        df = df \
            .Define(f"t{idx}", f"x{idx}.T()") \
            .Define(f"x{idx}_x", f"x{idx}.X()") \
            .Define(f"x{idx}_y", f"x{idx}.Y()") \
            .Define(f"x{idx}_z", f"x{idx}.Z()") \
            .Define(f'origin{idx}', f'std::isnan(p{idx}_mother.Px())') \

    # Define variables for pairs
    for idx1, idx2 in [(1, 2), (1, 3), (2, 3)]:
        df = df \
            .Define(f"beta{idx1}{idx2}", f"(p{idx1}+p{idx2}).BoostVector()") \
            .Define(f"origin{idx1}{idx2}", f"(origin{idx1} << 1) + origin{idx2}") \
            .Define(f"p{idx1}_com{idx1}{idx2}", f"TLorentzVector tmp = p{idx1}; tmp.Boost(-beta{idx1}{idx2}); return tmp;") \
            .Define(f"p{idx2}_com{idx1}{idx2}", f"TLorentzVector tmp = p{idx2}; tmp.Boost(-beta{idx1}{idx2}); return tmp;") \
            .Define(f"x{idx1}_com{idx1}{idx2}", f"TLorentzVector tmp = x{idx1}; tmp.Boost(-beta{idx1}{idx2}); return tmp;") \
            .Define(f"x{idx2}_com{idx1}{idx2}", f"TLorentzVector tmp = x{idx2}; tmp.Boost(-beta{idx1}{idx2}); return tmp;") \
            .Define(f'beta{idx1}x_com{idx1}{idx2}', f'p{idx1}_com{idx1}{idx2}.Px()/p{idx1}_com{idx1}{idx2}.E()') \
            .Define(f'beta{idx1}y_com{idx1}{idx2}', f'p{idx1}_com{idx1}{idx2}.Py()/p{idx1}_com{idx1}{idx2}.E()') \
            .Define(f'beta{idx1}z_com{idx1}{idx2}', f'p{idx1}_com{idx1}{idx2}.Pz()/p{idx1}_com{idx1}{idx2}.E()') \
            .Define(f'beta{idx2}x_com{idx1}{idx2}', f'p{idx2}_com{idx1}{idx2}.Px()/p{idx2}_com{idx1}{idx2}.E()') \
            .Define(f'beta{idx2}y_com{idx1}{idx2}', f'p{idx2}_com{idx1}{idx2}.Py()/p{idx2}_com{idx1}{idx2}.E()') \
            .Define(f'beta{idx2}z_com{idx1}{idx2}', f'p{idx2}_com{idx1}{idx2}.Pz()/p{idx2}_com{idx1}{idx2}.E()') \
            .Define(f"tmax{idx1}{idx2}", f"std::max(x{idx1}_com{idx1}{idx2}.T(), x{idx2}_com{idx1}{idx2}.T())") \
            .Define(f"x{idx1}_com{idx1}{idx2}_prop", f"TLorentzVector tmp(x{idx1}_com{idx1}{idx2}.X() + beta{idx1}x_com{idx1}{idx2} * (tmax{idx1}{idx2} - x{idx1}_com{idx1}{idx2}.T()), x{idx1}_com{idx1}{idx2}.Y() + beta{idx1}y_com{idx1}{idx2} * (tmax{idx1}{idx2} - x{idx1}_com{idx1}{idx2}.T()), x{idx1}_com{idx1}{idx2}.Z() + beta{idx1}z_com{idx1}{idx2} * (tmax{idx1}{idx2} - x{idx1}_com{idx1}{idx2}.T()), tmax{idx1}{idx2}); return tmp;") \
            .Define(f"x{idx2}_com{idx1}{idx2}_prop", f"TLorentzVector tmp(x{idx2}_com{idx1}{idx2}.X() + beta{idx2}x_com{idx1}{idx2} * (tmax{idx1}{idx2} - x{idx2}_com{idx1}{idx2}.T()), x{idx2}_com{idx1}{idx2}.Y() + beta{idx2}y_com{idx1}{idx2} * (tmax{idx1}{idx2} - x{idx2}_com{idx1}{idx2}.T()), x{idx2}_com{idx1}{idx2}.Z() + beta{idx2}z_com{idx1}{idx2} * (tmax{idx1}{idx2} - x{idx2}_com{idx1}{idx2}.T()), tmax{idx1}{idx2}); return tmp;") \
            .Define(f'ptot_com{idx1}{idx2}', f'(p{idx2}_com{idx1}{idx2} + p{idx1}_com{idx1}{idx2}).P()') \
            .Define(f'kstar{idx1}{idx2}', f'0.5 * (p{idx2}_com{idx1}{idx2} - p{idx1}_com{idx1}{idx2}).P()') \
            .Define(f'rstar{idx1}{idx2}', f'(x{idx2}_com{idx1}{idx2}_prop - x{idx1}_com{idx1}{idx2}_prop).P()') \
            .Define(f"mT{idx1}{idx2}", f"(p{idx1}+p{idx2}).Mt()/2") \

    mtot = m1 + m2 + m3;
    mu12 = m1 * m2 / (m1 + m2); # Reduced mass of particles 1 and 2
    mu3_12 = m3 * (m1 + m2) / mtot; # Reduced mass of particle 3 wrt 1 and 2

    # Based on Eq. 1.80 of 3B notes
    alpha = 4 * m3 * m3 / pow(m1 + m3, 2. ) + 4 * m3 * m3 / pow(m2 + m3, 2.) + 4
    beta = 4 * m3 * mtot / (m1 + m2) * (m2 / pow(m2 + m3, 2) - m1 / pow(m1 + m3, 2))
    gamma = 4 * mtot * mtot / pow(m1 + m2, 2) *  (m1 * m1 / pow(m1 + m3, 2) + m2 * m2 / pow(m2 + m3, 2))

    df = df \
        .Define('beta', '(p1+p2+p3).BoostVector()') \
        .Define('origin', '(origin1 << 2) + (origin2 << 1) + origin3') \
        .Define('p1_com', 'TLorentzVector tmp = p1; tmp.Boost(-beta); return tmp;') \
        .Define('p2_com', 'TLorentzVector tmp = p2; tmp.Boost(-beta); return tmp;') \
        .Define('p3_com', 'TLorentzVector tmp = p3; tmp.Boost(-beta); return tmp;') \
        .Define('x1_com', 'TLorentzVector tmp = x1; tmp.Boost(-beta); return tmp;') \
        .Define('x2_com', 'TLorentzVector tmp = x2; tmp.Boost(-beta); return tmp;') \
        .Define('x3_com', 'TLorentzVector tmp = x3; tmp.Boost(-beta); return tmp;') \
        .Define('beta1x', 'p1_com.Px()/p1_com.E()') \
        .Define('beta1y', 'p1_com.Py()/p1_com.E()') \
        .Define('beta1z', 'p1_com.Pz()/p1_com.E()') \
        .Define('beta2x', 'p2_com.Px()/p2_com.E()') \
        .Define('beta2y', 'p2_com.Py()/p2_com.E()') \
        .Define('beta2z', 'p2_com.Pz()/p2_com.E()') \
        .Define('beta3x', 'p3_com.Px()/p3_com.E()') \
        .Define('beta3y', 'p3_com.Py()/p3_com.E()') \
        .Define('beta3z', 'p3_com.Pz()/p3_com.E()') \
        .Define('tmax', 'std::max(std::max(x1_com.T(), x2_com.T()), x3_com.T())') \
        .Define('x1_com_prop', 'TLorentzVector tmp(x1_com.X() + beta1x * (tmax - x1_com.T()), x1_com.Y() + beta1y * (tmax - x1_com.T()), x1_com.Z() + beta1z * (tmax - x1_com.T()), tmax); return tmp;') \
        .Define('x2_com_prop', 'TLorentzVector tmp(x2_com.X() + beta2x * (tmax - x2_com.T()), x2_com.Y() + beta2y * (tmax - x2_com.T()), x2_com.Z() + beta2z * (tmax - x2_com.T()), tmax); return tmp;') \
        .Define('x3_com_prop', 'TLorentzVector tmp(x3_com.X() + beta3x * (tmax - x3_com.T()), x3_com.Y() + beta3y * (tmax - x3_com.T()), x3_com.Z() + beta3z * (tmax - x3_com.T()), tmax); return tmp;') \
        .Define('mT', '(p1 + p2 + p3).Mt() / 3') \
        .Define('k12', f'{m2} / ({m1 + m2}) * p1_com.Vect() - {m1} / ({m1 + m2}) * p2_com.Vect()') \
        .Define('k3_12', f' {(m1 + m2) / mtot}  * p3_com.Vect() - {m3 / mtot} * (p1_com.Vect() + p2_com.Vect())') \
        .Define('Q3', f'std::sqrt({alpha} * (k12 * k12) + 2 * {beta} * (k12 * k3_12) + {gamma} * (k3_12 * k3_12))') \
        .Define('Q', f'std::sqrt({arbitraryMass} / {mu12} * k12 * k12 + {arbitraryMass} / {mu3_12} * k3_12 * k3_12)') \
        .Define('r12', 'x1_com_prop.Vect() - x2_com_prop.Vect()') \
        .Define('r3_12', f'-{m1 / (m1 + m2)} * x1_com_prop.Vect() - {m2 / (m1 + m2)} * x2_com_prop.Vect() + x3_com_prop.Vect()') \
        .Define('hyp_rad', f'std::sqrt(({mu12} * r12 * r12 + {mu3_12} * r3_12 * r3_12) / {arbitraryMass})') \
        .Define('hyp_angle', f'std::atan2(std::sqrt({mu3_12} * r3_12 * r3_12), std::sqrt({mu12} * r12 * r12))') \

    return df

def BookParticleHistograms(df, idx):
    return {
        f'hT{idx}' : df.Histo1D((f'hT{idx}', f";t_{{idx}} (fm/c);Counts", *BINNING_SOURCE), f't{idx}'),
        f'hX{idx}' : df.Histo1D((f'hX{idx}', f";x_{{idx}} (fm);Counts", *BINNING_SOURCE), f'x{idx}_x'),
        f'hY{idx}' : df.Histo1D((f'hY{idx}', f";y_{{idx}} (fm);Counts", *BINNING_SOURCE), f'x{idx}_y'),
        f'hZ{idx}' : df.Histo1D((f'hZ{idx}', f";z_{{idx}} (fm);Counts", *BINNING_SOURCE), f'x{idx}_z'),
    }
    
def BookPairHistograms(df, idx1, idx2, max_kstar):
    # If the column contains nan it was not supposed to be used. Use then an empty dataframe and ensure empty but properly defined histograms
    df = df.Filter(f"!std::isnan(x{idx1}.X())").Filter(f"!std::isnan(x{idx2}.X())")
    df_femto = df.Filter(f'kstar{idx1}{idx2} < {max_kstar}')

    return {
        f"hKStar{idx1}{idx2}" : df.Histo1D((f"hKStar{idx1}{idx2}", f";k*_{{({idx1}, {idx2})}};Counts", *BINNING_MOMENTUM), f'kstar{idx1}{idx2}'),
        f"hTMax{idx1}{idx2}" : df_femto.Histo1D((f"hTMax{idx1}{idx2}", f";t_{{max}}^{{({idx1}, {idx2})}};Counts", *BINNING_SOURCE), f'tmax{idx1}{idx2}'),
        f"hRStar{idx1}{idx2}" : df_femto.Histo1D((f"hRStar{idx1}{idx2}", f";r*_{{({idx1}, {idx2})}} (fm);Counts", *BINNING_SOURCE), f'rstar{idx1}{idx2}'),
        f"hRStar{idx1}{idx2}_pp" : df_femto.Filter(f'origin{idx1}{idx2} == 0b11').Histo1D((f"hRStar{idx1}{idx2}pp", f";r*_{{({idx1}, {idx2})}}^{{pp}} (fm);Counts", *BINNING_SOURCE), f'rstar{idx1}{idx2}'),
        f"hRStar{idx1}{idx2}_ps" : df_femto.Filter(f'origin{idx1}{idx2} == 0b10').Histo1D((f"hRStar{idx1}{idx2}ps", f";r*_{{({idx1}, {idx2})}}^{{ps}} (fm);Counts", *BINNING_SOURCE), f'rstar{idx1}{idx2}'),
        f"hRStar{idx1}{idx2}_sp" : df_femto.Filter(f'origin{idx1}{idx2} == 0b01').Histo1D((f"hRStar{idx1}{idx2}sp", f";r*_{{({idx1}, {idx2})}}^{{sp}} (fm);Counts", *BINNING_SOURCE), f'rstar{idx1}{idx2}'),
        f"hRStar{idx1}{idx2}_ss" : df_femto.Filter(f'origin{idx1}{idx2} == 0b00').Histo1D((f"hRStar{idx1}{idx2}ss", f";r*_{{({idx1}, {idx2})}}^{{ss}} (fm);Counts", *BINNING_SOURCE), f'rstar{idx1}{idx2}'),
        f"hRStarVsMt{idx1}{idx2}" : df_femto.Histo2D((f"hRStarVsMt{idx1}{idx2}", f";m_{{T}}^{{({idx1}, {idx2})}} (GeV/#it{{c}});r*_{{({idx1}, {idx2})}} (fm);Counts",*BINNING_MT, *BINNING_SOURCE), f'mT{idx1}{idx2}', f'rstar{idx1}{idx2}'),
    }

def BookTripletHistograms(df, max_Q3):
    # If the column contains nan it was not supposed to be used. Use then an empty dataframe and ensure empty but properly defined histograms
    df = df.Filter("!std::isnan(x2.X())").Filter("!std::isnan(x3.X())")

    df_femto = df.Filter(f'Q3 < {max_Q3}')

    return {
        'hQ' : df.Histo1D((f'hQ', ';Q (GeV/#it{c});Counts', *BINNING_MOMENTUM), 'Q'),
        'hQ3' : df.Histo1D((f'hQ3', ';Q_{3} (GeV/#it{c});Counts', *BINNING_MOMENTUM), 'Q3'),
        'hHypRad' : df_femto.Histo1D((f'hHypRad', ';#rho (fm);Counts', *BINNING_SOURCE), 'hyp_rad'),
        'hHypRad_ppp' : df_femto.Filter(f'origin == 0b111').Histo1D((f'hHypRad_ppp', ';#rho_{{ppp}} (fm);Counts', *BINNING_SOURCE), 'hyp_rad'),
        'hHypRad_pps' : df_femto.Filter(f'origin == 0b110').Histo1D((f'hHypRad_pps', ';#rho_{{pps}} (fm);Counts', *BINNING_SOURCE), 'hyp_rad'),
        'hHypRad_psp' : df_femto.Filter(f'origin == 0b101').Histo1D((f'hHypRad_psp', ';#rho_{{psp}} (fm);Counts', *BINNING_SOURCE), 'hyp_rad'),
        'hHypRad_pss' : df_femto.Filter(f'origin == 0b100').Histo1D((f'hHypRad_pss', ';#rho_{{pss}} (fm);Counts', *BINNING_SOURCE), 'hyp_rad'),
        'hHypRad_spp' : df_femto.Filter(f'origin == 0b011').Histo1D((f'hHypRad_spp', ';#rho_{{spp}} (fm);Counts', *BINNING_SOURCE), 'hyp_rad'),
        'hHypRad_sps' : df_femto.Filter(f'origin == 0b010').Histo1D((f'hHypRad_sps', ';#rho_{{sps}} (fm);Counts', *BINNING_SOURCE), 'hyp_rad'),
        'hHypRad_ssp' : df_femto.Filter(f'origin == 0b001').Histo1D((f'hHypRad_ssp', ';#rho_{{ssp}} (fm);Counts', *BINNING_SOURCE), 'hyp_rad'),
        'hHypRad_sss' : df_femto.Filter(f'origin == 0b000').Histo1D((f'hHypRad_sss', ';#rho_{{sss}} (fm);Counts', *BINNING_SOURCE), 'hyp_rad'),
        'hHypAngle' : df_femto.Histo1D((f'hHypAngle', ';#varphi (rad);Counts', 200, 0, np.pi / 2), 'hyp_angle'),
        'hHypRadVsMt' : df_femto.Histo2D(('hHypRadVsMt', ';m_{T}^{3B} (GeV/#it{c});#rho (fm)', *BINNING_MT, *BINNING_SOURCE), 'mT', 'hyp_rad'),
        'hRStar12VsMt' : df_femto.Histo2D((f'hRStar12VsMt', ';m_{T}^{3B} (GeV/#it{c});r*_{(1,2)} (fm)', *BINNING_MT, *BINNING_SOURCE), 'mT', 'rstar12'),
        'hRStar13VsMt' : df_femto.Histo2D((f'hRStar13VsMt', ';m_{T}^{3B} (GeV/#it{c});r*_{(1,3)} (fm)', *BINNING_MT, *BINNING_SOURCE), 'mT', 'rstar13'),
        'hRStar23VsMt' : df_femto.Histo2D((f'hRStar23VsMt', ';m_{T}^{3B} (GeV/#it{c});r*_{(2,3)} (fm)', *BINNING_MT, *BINNING_SOURCE), 'mT', 'rstar23'),
        'hPair12MtVsTripletMt' : df_femto.Histo2D((f'hPair12MtVsTripletMt', ';m_{T}^{3B} (GeV/#it{c});m_{T}^{(1,2)} (GeV/#it{c});Counts', *BINNING_MT, *BINNING_MT), 'mT', 'mT12'),
        'hPair13MtVsTripletMt' : df_femto.Histo2D((f'hPair13MtVsTripletMt', ';m_{T}^{3B} (GeV/#it{c});m_{T}^{(1,3)} (GeV/#it{c});Counts', *BINNING_MT, *BINNING_MT), 'mT', 'mT13'),
        'hPair23MtVsTripletMt' : df_femto.Histo2D((f'hPair23MtVsTripletMt', ';m_{T}^{3B} (GeV/#it{c});m_{T}^{(2,3)} (GeV/#it{c});Counts', *BINNING_MT, *BINNING_MT), 'mT', 'mT23'),
    }

def ProcessParticle(hists, idx, dir):
    log.info(f'Processing particle {idx}...')

    subdir = dir.mkdir(f'part{idx}')
    subdir.cd()
    
    for hist in hists[f'part{idx}'].values():
        hist.Write()

    dir.cd()

def ProcessPair(hists, idx1, idx2, dir):
    log.info(f'Processing pair ({idx1}, {idx2})...')

    subdir = dir.mkdir(f'pair{idx1}{idx2}')
    subdir.cd()

    for hist in hists[f'pair{idx1}{idx2}'].values():
        hist.Write()

    pRStarVsMt = hists[f'pair{idx1}{idx2}'][f"hRStarVsMt{idx1}{idx2}"].ProfileX()
    pRStarVsMt.SetTitle(';m_{T} (GeV/#it{c});r* (fm)')
    pRStarVsMt.Write(f'pRStarVsMt{idx1}{idx2}')

    dir.cd()

def ProcessTriplet(hists, dir):
    log.info(f'Processing triplet...')

    subdir = dir.mkdir('triplet')
    subdir.cd()

    for key, hist in hists['triplet'].items():
        if key in ['hPair12MtVsTripletMt', 'hPair13MtVsTripletMt','hPair23MtVsTripletMt', 'hRStar12VsMt', 'hRStar13VsMt', 'hRStar23VsMt']:
            continue
        hist.Write()

    pHypRadVsMt = hists['triplet']['hHypRadVsMt'].ProfileX('pHypRadVsMt')
    pHypRadVsMt.SetTitle(';m_{T}^{3B} (GeV/#it{c});#LT#rho#GT (fm)')
    pHypRadVsMt.Write(f'gHypRadVsMt')

    subsubdir = dir.mkdir('triplet/slices_mT')
    subsubdir.cd()
    hHypRadVsMt = utils.analysis.SliceVertically(hists['triplet']['hHypRadVsMt'], name='hHypRad_mT')
    for hist in hHypRadVsMt:
        hist.Write()
    subdir.cd()

    hRStarVsMt = hists['triplet']['hRStar12VsMt'].GetValue().Clone("hRStarVsMt")
    hRStarVsMt.Add(hists['triplet']['hRStar13VsMt'].GetValue())
    hRStarVsMt.Add(hists['triplet']['hRStar23VsMt'].GetValue())
    hRStarVsMt.SetTitle(';m_{T}^{3B} (GeV/#it{c});r* (fm)')
    hRStarVsMt.Write()

    hPairMtVsTripletMt = hists['triplet']['hPair12MtVsTripletMt'].GetValue().Clone("hPairMtVsTripletMt")
    hPairMtVsTripletMt.Add(hists['triplet']['hPair13MtVsTripletMt'].GetValue())
    hPairMtVsTripletMt.Add(hists['triplet']['hPair23MtVsTripletMt'].GetValue())
    hPairMtVsTripletMt.SetTitle(';m_{T}^{3B} (GeV/#it{c});m_{T}^{2B} (GeV/#it{c});Counts')
    hPairMtVsTripletMt.Write()

    pPairMtVsTripletMt = hPairMtVsTripletMt.ProfileX('pPairMtVsTripletMt')
    pPairMtVsTripletMt.SetTitle(';m_{T}^{3B} (GeV/#it{c});m_{T}^{2B} (GeV/#it{c})')
    pPairMtVsTripletMt.Write('pPairMtVsTripletMt')

    pRStarVsMt = hRStarVsMt.ProfileX('pRStarVsMt')
    pRStarVsMt.SetTitle(';m_{T}^{3B} (GeV/#it{c});#LT r*#GT (fm)')
    pRStarVsMt.Write(f'gRStarVsMt123')

    pExpectedHypRadVsMt = pRStarVsMt.Clone('pExpectedHypRadVsMt')
    pExpectedHypRadVsMt.SetTitle(';m_{T}^{3B} (GeV/#it{c});#LT#rho#GT_{exp} (fm)')
    pExpectedHypRadVsMt.Scale(15 * np.pi / 32)
    pExpectedHypRadVsMt.Write()

    # Compute expected mT scaling of 3B via mT sclaing of 2B 
    hRStarVsPairMt = hists['pair12']['hRStarVsMt12'].GetValue().Clone("hRStarVsPairMt")
    hRStarVsPairMt.Add(hists['pair13']['hRStarVsMt13'].GetValue())
    hRStarVsPairMt.Add(hists['pair23']['hRStarVsMt23'].GetValue())
    pRStarVsPairMt = hRStarVsPairMt.ProfileX('pRStarVsPairMt')
    pRStarVsPairMt.Write()
    avgPairMts = []
    for iBin in range(pPairMtVsTripletMt.GetNbinsX()):
        entries = pPairMtVsTripletMt.GetBinEntries(iBin + 1) > 0
        if entries > 0:
            avgPairMts.append(pPairMtVsTripletMt.GetBinContent(iBin + 1))
            if entries < 30:
                log.warning('Less than 30 entries found: unstable! Consider rebinning or running more statistics')
        else:
            log.warning('Not enough statistics for profiling')
            avgPairMts.append(None)
    rStars = [pRStarVsPairMt.GetBinContent(pRStarVsPairMt.FindBin(mT)) if mT is not None else None for mT in avgPairMts]
    rStarUncs = [pRStarVsPairMt.GetBinError(pRStarVsPairMt.FindBin(mT)) if mT is not None else None for mT in avgPairMts]
    hExpHypRadFrom2BVsMt = pHypRadVsMt.ProjectionX().Clone('hExpHypRadFrom2BVsMt')
    hExpHypRadFrom2BVsMt.Reset()
    for iBin, (rStar, rStarUnc) in enumerate(zip(rStars, rStarUncs)):
        if rStar is None or rStarUnc is None:
            continue
        hExpHypRadFrom2BVsMt.SetBinContent(iBin + 1, rStar)
        hExpHypRadFrom2BVsMt.SetBinError(iBin + 1, rStarUnc)
    hExpHypRadFrom2BVsMt.Scale(15 * np.pi / 32)
    hExpHypRadFrom2BVsMt.Write()

    dir.cd()

if __name__ == '__main__':
    gROOT.SetBatch(True)
    EnableImplicitMT()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('ofile')
    parser.add_argument('-f', help=' Fraction of the dataset to analyze', default=1.0, type=float)
    parser.add_argument('--max-kstar', help=' Max k* (MeV)', default=100, type=float)
    parser.add_argument('--max-Q3', help=' Max Q3 (MeV)', default=600, type=float)
    args = parser.parse_args()

    tree = TChain('tEvents')
    tree.Add(args.infile)
    nEntries = tree.GetEntries()

    df = RDataFrame(tree)
    RDF.Experimental.AddProgressBar(df)
    df = df.Filter(f"rdfentry_ < {int(args.f * nEntries)}")

    df = DefineVariables(df, 938.272, 938.272, 938.272, 938.272 / 2)

    oFile = TFile(args.ofile, 'recreate')
    oFile.cd()

    hists = {}
    hists['part1'] = BookParticleHistograms(df, 1)
    hists['part2'] = BookParticleHistograms(df, 2)
    hists['part3'] = BookParticleHistograms(df, 3)
    hists['pair12'] = BookPairHistograms(df, 1, 2, args.max_kstar)
    hists['pair13'] = BookPairHistograms(df, 1, 3, args.max_kstar)
    hists['pair23'] = BookPairHistograms(df, 2, 3, args.max_kstar)
    hists['triplet'] = BookTripletHistograms(df, args.max_Q3)

    log.info('The loop over the dataset has begun. This might take some time...')
    RDF.RunGraphs([g for group in hists.values() for g in group.values()])

    ProcessParticle(hists, 1, oFile)
    ProcessParticle(hists, 2, oFile)
    ProcessParticle(hists, 3, oFile)

    ProcessPair(hists, 1, 2, oFile)
    ProcessPair(hists, 1, 3, oFile)
    ProcessPair(hists, 2, 3, oFile)

    ProcessTriplet(hists, oFile)

    oFile.Close()
