#!/bin/env python3
from pathlib import Path
import argparse
import sys
from dotenv import load_dotenv
import numpy as np
import os

from ROOT import RDataFrame, TCanvas, TFile, gDirectory, gInterpreter, TF1, gROOT, TGraphErrors, TGraphAsymmErrors

env_path = Path(__file__).resolve().parent.parent.parent.parent / ".env"
if not load_dotenv(dotenv_path=env_path):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
sys.path.append(f'{YAFFA_PATH}/src/python')

from yaffa import logger as log
from yaffa import utils

gROOT.SetBatch(True)


parser = argparse.ArgumentParser()
parser.add_argument('infile')
parser.add_argument('ofile')
args = parser.parse_args()


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

def ProcessPair(df, idx1, idx2):
    df_tmp = df \
        .Define("tmax", f"std::max(p{idx1}.T(), p{idx2}.T())") \
        .Define(f'is_p{idx1}_primary', f'std::isnan(p{idx1}_mother.Px())') \
        .Define(f'is_p{idx2}_primary', f'std::isnan(p{idx2}_mother.Px())') \
        .Define("beta", f"(p{idx1}+p{idx2}).BoostVector()") \
        .Define(f"p{idx1}_com", f"TLorentzVector tmp = p{idx1}; tmp.Boost(-beta); return tmp;") \
        .Define(f"p{idx2}_com", f"TLorentzVector tmp = p{idx2}; tmp.Boost(-beta); return tmp;") \
        .Define(f"x{idx1}_com", f"TLorentzVector tmp = x{idx1}; tmp.Boost(-beta); return tmp;") \
        .Define(f"x{idx2}_com", f"TLorentzVector tmp = x{idx2}; tmp.Boost(-beta); return tmp;") \
        .Define(f"x{idx1}_com_prop", f"TLorentzVector tmp(x{idx1}_com.X() + beta.X() * (tmax - x{idx1}_com.T()), x{idx1}_com.Y() + beta.Y() * (tmax - x{idx1}_com.T()), x{idx1}_com.Z() + beta.Z() * (tmax - x{idx1}_com.T()), tmax); return tmp;") \
        .Define('kstar', f'0.5 * (p{idx2}_com - p{idx1}_com).P()') \
        .Define('rstar', f'(x{idx2}_com - x{idx1}_com).P()') \
        .Define("mT", f"(p{idx1}+p{idx2}).Mt()/2") \

    c = TCanvas('c', '', 600, 600)

    hRStarVsMt = df_tmp.Filter('kstar < 100').Histo2D((f"hRStarVsMt{idx1}{idx2}", f";m_{{T}}^{({idx1}, {idx2})}", 82, 950, 5050, 200, 0, 20), 'mT', 'rstar')
    hRStarVsMt.Draw()
    c.SaveAs('chRStarVsMt.pdf')

    gMt = GetMtScaling(hRStarVsMt)
    gMt.Write(f'gMt{idx1}{idx2}')

df = RDataFrame('tEvents', args.infile)

oFile = TFile(args.ofile, 'recreate')

ProcessPair(df, 1, 2)
ProcessPair(df, 1, 3)
ProcessPair(df, 2, 3)

oFile.Close()
