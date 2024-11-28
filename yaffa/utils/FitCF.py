'''
Script to fit femtoscopic correlation functions.
'''

import os
import argparse
import yaml

from ROOT import TFile, TCanvas, gInterpreter, gROOT
gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/Observable.h"')
gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/SuperFitter.h"')
from ROOT import Observable, SuperFitter

from yaffa import utils

def FitCF(cfg): # pylint disable:missing-function-docstring
    '''Fit the correlation functions.

    Args:
        cfg (dict): configuration of the fit
    '''
    inFile = TFile(cfg['fits'][0]['infile'])
    hObs = utils.io.Load(inFile, cfg['fits'][0]['name'])
    hObs.SetDirectory(0)
    oObs = Observable(hObs)
    inFile.Close()

    fitter = SuperFitter(oObs, 0, 0.5)
    fitter.Add('pol1', [
        ("a0", 1.05, 1, -1),
        ("a1", -0.1, 1, -1),
    ])

    # Add template
    templFile = TFile(cfg['fits'][0]['terms'][4]['file'])
    hTemplate = utils.io.Load(templFile, cfg['fits'][0]['terms'][4]['path'])

    # input()
    templFile1 = TFile(cfg['fits'][0]['terms'][5]['file'])
    hTemplate1 = utils.io.Load(templFile1, cfg['fits'][0]['terms'][5]['path'])

    # print(hTemplate)
    fitter.Add('Sigma1385_dir', hTemplate, [
        ("norm_Sigma1385_dir", 0.1, 0, -2.5),
    ])
    fitter.Add('Sigma1385_indir', hTemplate1, [
        ("norm_Sigma1385_indir", 5, 0, -2.5),
    ])

    cFit = TCanvas('cFit', '', 600, 600)
    hTemplate1.Draw()

    fitter.Fit("poll1 + Sigma1385_dir + Sigma1385_indir", 'MR+')

    cFit.DrawFrame(0, 0.99, 0.5, 1.12)
    fitter.Draw('same')

    cFit.SaveAs('cFit.pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', default='cfg_fit.yml', nargs='?')
    args = parser.parse_args()

    with open(args.cfg) as file:
        config = yaml.safe_load(file)

    gROOT.SetBatch(True)

    FitCF(config)
