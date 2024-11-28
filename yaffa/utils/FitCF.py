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

    for fitCfg in cfg['fits']:
        inFile = TFile(fitCfg['infile'])
        hObs = utils.io.Load(inFile, fitCfg['path'])
        hObs.SetDirectory(0)
        oObs = Observable(hObs)
        inFile.Close()

        fitter = SuperFitter(oObs, fitCfg['fitrange'][0][0], fitCfg['fitrange'][0][1])

        # Add template to the fitter
        for term in fitCfg['terms']:
            templFile = TFile(term['file'])
            hTemplate = utils.io.Load(templFile, term['path'])
            hTemplate.SetDirectory(0)
            templFile.Close()

            fitter.Add(term['name'], hTemplate, term['params'])

        fitter.Fit(fitCfg['model'], 'MR+')

        cFit = TCanvas('cFit', '', 600, 600)
        cFit.DrawFrame(0, 0.99, 0.6, 1.12)
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
