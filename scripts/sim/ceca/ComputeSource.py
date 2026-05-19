#!/bin/env python3
from pathlib import Path
import argparse
import sys
from dotenv import load_dotenv
import numpy as np
import os

from ROOT import RDataFrame, TChain, TFile, gROOT, TGraphAsymmErrors

env_path = Path(__file__).resolve().parent.parent.parent.parent / ".env"
if not load_dotenv(dotenv_path=env_path):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
sys.path.append(f'{YAFFA_PATH}/src/python')

from yaffa import logger as log
from yaffa import utils

femtoLimit = 100

def GetMtScaling(gRadiusVsMt):
    mTbins = range(0, 4001, 100)
    slices = utils.analysis.SliceVertically(gRadiusVsMt, mTbins, name='hRStar_mT')
    means = [hist.GetMean() for hist in slices]
    errors = [hist.GetMeanError() for hist in slices]
    
    x = np.array([(mTbins[i] + mTbins[i + 1]) / 2 for i in range(len(mTbins) - 1)], dtype='double')
    y = np.array(means, dtype='double')
    yUnc = np.array(errors, dtype='double')

    # TODO: Now using standard mean error (assume gaussianity). Think of a better way (Bootstrap?)
    gMtScaling = TGraphAsymmErrors(len(mTbins) - 1, x, y, 0, 0, yUnc, yUnc)

    return gMtScaling

def ProcessParticle(df, idx, dir):
    log.info(f'Processing particle {idx}...')
    df_tmp = df \
        .Define(f"t{idx}", f"x{idx}.T()") \
        .Define(f"x{idx}_x", f"x{idx}.X()") \
        .Define(f"x{idx}_y", f"x{idx}.Y()") \
        .Define(f"x{idx}_z", f"x{idx}.Z()") \

    dir.mkdir(f'part{idx}').cd()
    
    df_tmp.Histo1D((f'hT{idx}', f";t_{idx} (fm/c);Counts", 200, 0, 20), f't{idx}').Write()
    df_tmp.Histo1D((f'hX{idx}', f";x_{idx} (fm/c);Counts", 200, 0, 20), f'x{idx}_x').Write()
    df_tmp.Histo1D((f'hY{idx}', f";y_{idx} (fm/c);Counts", 200, 0, 20), f'x{idx}_y').Write()
    df_tmp.Histo1D((f'hZ{idx}', f";z_{idx} (fm/c);Counts", 200, 0, 20), f'x{idx}_z').Write()

    dir.cd()

def ProcessPair(df, idx1, idx2, dir):
    # If the column contains nan it was not supposed to be used. Use then an empty dataframe and ensure empty but properly defined histograms
    df = df.Filter(f"!std::isnan(x{idx1}.X())").Filter(f"!std::isnan(x{idx2}.X())")

    log.info(f'Processing pair ({idx1}, {idx2})...')

    dir.mkdir(f'pair{idx1}{idx2}').cd()

    df_tmp = df \
        .Define(f'is_p{idx1}_primary', f'std::isnan(p{idx1}_mother.Px())') \
        .Define(f'is_p{idx2}_primary', f'std::isnan(p{idx2}_mother.Px())') \
        .Define("beta", f"(p{idx1}+p{idx2}).BoostVector()") \
        .Define(f"p{idx1}_com", f"TLorentzVector tmp = p{idx1}; tmp.Boost(-beta); return tmp;") \
        .Define(f"p{idx2}_com", f"TLorentzVector tmp = p{idx2}; tmp.Boost(-beta); return tmp;") \
        .Define(f"x{idx1}_com", f"TLorentzVector tmp = x{idx1}; tmp.Boost(-beta); return tmp;") \
        .Define(f"x{idx2}_com", f"TLorentzVector tmp = x{idx2}; tmp.Boost(-beta); return tmp;") \
        .Define(f'beta{idx1}x', f'p{idx1}_com.Px()/p{idx1}_com.E()') \
        .Define(f'beta{idx1}y', f'p{idx1}_com.Py()/p{idx1}_com.E()') \
        .Define(f'beta{idx1}z', f'p{idx1}_com.Pz()/p{idx1}_com.E()') \
        .Define(f'beta{idx2}x', f'p{idx2}_com.Px()/p{idx2}_com.E()') \
        .Define(f'beta{idx2}y', f'p{idx2}_com.Py()/p{idx2}_com.E()') \
        .Define(f'beta{idx2}z', f'p{idx2}_com.Pz()/p{idx2}_com.E()') \
        .Define("tmax", f"std::max(x{idx1}_com.T(), x{idx2}_com.T())") \
        .Define(f"x{idx1}_com_prop", f"TLorentzVector tmp(x{idx1}_com.X() + beta{idx1}x * (tmax - x{idx1}_com.T()), x{idx1}_com.Y() + beta{idx1}y * (tmax - x{idx1}_com.T()), x{idx1}_com.Z() + beta{idx1}z * (tmax - x{idx1}_com.T()), tmax); return tmp;") \
        .Define(f"x{idx2}_com_prop", f"TLorentzVector tmp(x{idx2}_com.X() + beta{idx2}x * (tmax - x{idx2}_com.T()), x{idx2}_com.Y() + beta{idx2}y * (tmax - x{idx2}_com.T()), x{idx2}_com.Z() + beta{idx2}z * (tmax - x{idx2}_com.T()), tmax); return tmp;") \
        .Define('ptot_com', f'(p{idx2}_com + p{idx1}_com).P()') \
        .Define('kstar', f'0.5 * (p{idx2}_com - p{idx1}_com).P()') \
        .Define('rstar', f'(x{idx2}_com_prop - x{idx1}_com_prop).P()') \
        .Define("mT", f"(p{idx1}+p{idx2}).Mt()/2") \

    df_tmp.Histo1D((f"hTMax{idx1}{idx2}", f";t_{max}^{({idx1}, {idx2})}", 200, 0, 20), 'tmax').Write()
    df_tmp.Filter(f'kstar < {femtoLimit}').Histo1D((f"hRStar{idx1}{idx2}", f";r*_{({idx1}, {idx2})}", 200, 0, 20), 'rstar').Write()
    df_tmp.Histo1D((f"hKStar{idx1}{idx2}", f";k*_{({idx1}, {idx2})}", 200, 0, 2000), 'kstar').Write()
    df_tmp.Histo1D((f"hPTotCOM{idx1}{idx2}", f";k*_{({idx1}, {idx2})}", 200, 0, 0.1), 'ptot_com').Write()

    hRStarVsMt = df_tmp.Filter(f'kstar < {femtoLimit}').Histo2D((f"hRStarVsMt{idx1}{idx2}", f";m_{{T}}^{({idx1}, {idx2})}", 82, 950, 5050, 200, 0, 20), 'mT', 'rstar')
    hRStarVsMt.Write()
    gMt = GetMtScaling(hRStarVsMt)
    gMt.Write(f'gMt{idx1}{idx2}')

    dir.cd()

if __name__ == '__main__':
    gROOT.SetBatch(True)

    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('ofile')
    parser.add_argument('-f', help=' Fraction of the dataset to analyze', default=1.0, type=float)
    args = parser.parse_args()

    tree = TChain("tEvents")
    tree.Add(args.infile)
    nEntries = tree.GetEntries()

    df = RDataFrame(tree)
    df = df.Range(int(args.f * nEntries))

    oFile = TFile(args.ofile, 'recreate')
    oFile.cd()

    ProcessParticle(df, 1, oFile)
    ProcessParticle(df, 2, oFile)
    ProcessParticle(df, 3, oFile)

    ProcessPair(df, 1, 2, oFile)
    ProcessPair(df, 1, 3, oFile)
    ProcessPair(df, 2, 3, oFile)

    oFile.Close()
