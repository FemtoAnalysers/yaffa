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
    pars = [
        ("a0", 1.1, -10, 10),
        ("a1", 0.1, -10, 20),
        ("a2", 0.1, -10, 20),
    ]
    fitter.Add('pol2', pars)

    cFit = TCanvas('cFit', '', 600, 600)

    fitter.Fit("pol2", 'MR+')

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
