'''
Script to fit femtoscopic correlation functions.
'''

import argparse
import yaml
import os

from ROOT import TFile, TCanvas, gInterpreter
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

    fitter = SuperFitter(oObs)

    cFit = TCanvas('cFit', '', 600, 600)
    fitter.Fit()
    cFit.DrawFrame(0, 0.99, 0.5, 1.12)
    fitter.Draw('same')

    cFit.SaveAs('cFit.pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', default='cfg_fit.yml', nargs='?')
    args = parser.parse_args()

    with open(args.cfg) as file:
        config = yaml.safe_load(file)

    FitCF(config)
