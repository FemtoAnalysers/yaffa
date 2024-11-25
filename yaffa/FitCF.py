'''
Script to fit femtoscopic correlation functions.
'''

import argparse
import yaml


from ROOT import TFile, TCanvas, gInterpreter
gInterpreter.ProcessLine('#include "utils/Observable.h"')
from ROOT import Observable

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

    cFit = TCanvas('cFit', '', 600, 600)
    cFit.DrawFrame(0, 0.99, 0.5, 1.12)
    oObs.Draw('same')
    cFit.SaveAs('cFit.pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', default='cfg_fit.yml', nargs='?')
    args = parser.parse_args()

    with open(args.cfg) as file:
        config = yaml.safe_load(file)

    FitCF(config)
