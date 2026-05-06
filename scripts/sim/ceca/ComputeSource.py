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

from yaffa import logger as log

sys.path.append(f'{YAFFA_PATH}/src/python')

from yaffa import utils


parser = argparse.ArgumentParser()
parser.add_argument('infile')
parser.add_argument('ofile')
args = parser.parse_args()


def GetMtScaling(gRadiusVsMt):
    mTbins = range(0, 4001, 100)
    slices = utils.analysis.SliceVertically(gRadiusVsMt, mTbins, name='hRStar_mT')
    means = [hist.GetMean() for hist in slices]
    errors = [hist.GetMeanError() for hist in slices]
    
    x = np.array([(mTbins[i] + mTbins[i + 1]) / 2 for i in range(len(mTbins) -1)], dtype='double')
    y = np.array(means, dtype='double')
    yUnc = np.array(errors, dtype='double')

    # TODO: Now using standard mean error (assume gaussianity). Think of a better way (Bootstrap?)
    gMtScaling = TGraphAsymmErrors(len(mTbins) - 1, x, y, 0, 0, yUnc, yUnc)

    return gMtScaling

def ProcessPair(df):
    df_tmp = df \
        .Define('x1x', 'x1.Px()') \
        .Define('x1y', 'x1.Py()') \
        .Define('x1z', 'x1.Pz()') \
        .Define('x2t', 'x2.T()') \
        .Define('x2x', 'x2.Px()') \
        .Define('x2y', 'x2.Py()') \
        .Define('x2z', 'x2.Pz()') \
        .Define('x3t', 'x3.T()') \
        .Define('x3x', 'x3.Px()') \
        .Define('x3y', 'x3.Py()') \
        .Define('x3z', 'x3.Pz()') \
        .Define('p1e', 'p1.E()') \
        .Define('p1x', 'p1.Px()') \
        .Define('p1y', 'p1.Py()') \
        .Define('p1z', 'p1.Pz()') \
        .Define('p2e', 'p2.E()') \
        .Define('p2x', 'p2.Px()') \
        .Define('p2y', 'p2.Py()') \
        .Define('p2z', 'p2.Pz()') \
        .Define('p3e', 'p3.E()') \
        .Define('p3x', 'p3.Px()') \
        .Define('p3y', 'p3.Py()') \
        .Define('p3z', 'p3.Pz()') \
        .Define("tmax", "std::max(p1.T(), p2.T())") \
        .Define('is_p1_primary', 'std::isnan(p1_mother.Px())') \
        .Define('is_p2_primary', 'std::isnan(p2_mother.Px())') \
        .Define('is_p3_primary', 'std::isnan(p3_mother.Px())') \
        .Define("beta12", "(p1+p2).BoostVector()") \
        .Define("p1_com12", "TLorentzVector tmp = p1; tmp.Boost(-beta12); return tmp;") \
        .Define("p2_com12", "TLorentzVector tmp = p2; tmp.Boost(-beta12); return tmp;") \
        .Define("x1_com12", "TLorentzVector tmp = x1; tmp.Boost(-beta12); return tmp;") \
        .Define("x2_com12", "TLorentzVector tmp = x2; tmp.Boost(-beta12); return tmp;") \
        .Define('p1_comx', 'p1_com12.Px()') \
        .Define("x1_com12_prop", "TLorentzVector tmp(x1_com12.X() + beta12.X() * (tmax - x1_com12.T()), x1_com12.Y() + beta12.Y() * (tmax - x1_com12.T()), x1_com12.Z() + beta12.Z() * (tmax - x1_com12.T()), tmax); return tmp;") \
        .Define('p1_com12_propx', 'x1_com12_prop.Px()') \
        .Define('kstar12', '0.5 * (p2_com12 - p1_com12).P()') \
        .Define('rstar12', '(x2_com12 - x1_com12).P()') \
        .Define("mT12", "(p1+p2).Mt()/2") \

    c = TCanvas('c', '', 600, 600)

    hRStarVsMt12 = df_tmp.Filter('kstar12 < 100').Histo2D(("hRStarVsMt12", "", 82, 950, 5050, 200, 0, 20), 'mT12', 'rstar12')
    hRStarVsMt12.Draw()
    c.SaveAs('chRStarVsMt12.pdf')

    gMt = GetMtScaling(hRStarVsMt12)
    gMt.Write('gMt')

df = RDataFrame('tEvents', args.infile)

oFile = TFile(args.ofile, 'recreate')

ProcessPair(df)

oFile.Close()
