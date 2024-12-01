'''
Script to fit femtoscopic correlation functions.
'''

import sys
import os
import argparse
import yaml

from ROOT import TFile, TCanvas, gInterpreter, gROOT, TH1, TF1
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
            if templFileName := term.get('file'):
                templFile = TFile(templFileName)
                hTemplate = utils.io.Load(templFile, term['path'])
                print("Adding ", hTemplate.GetName())
                if isinstance(hTemplate, TH1):
                    print("Is TH1")
                    hTemplate.SetDirectory(0)
                    hTemplate = utils.analysis.ChangeUnits(hTemplate, term.get('units', 1))
                    fitter.Add(term['name'], hTemplate, term['params'])

                elif isinstance(hTemplate, TF1):
                    unitMult = term.get('unit_mult', 1)
                    print('mmimi', hTemplate.Eval(200))
                    fitter.Add(term['name'], hTemplate.Clone(), term['params'], unitMult)

                else:
                    print('Type not implemented. Exit!')
                    sys.exit()
                templFile.Close()
            else:
                fitter.Add(term['func'], term['params'])

        fitter.Fit(fitCfg['model'], 'MR+')

        cFit = TCanvas('cFit', '', 600, 600)
        # cFit.DrawFrame(0, 0.8, 0.5, 1.2)

        cFit.DrawFrame(0, 0.9, 0.5, 1.12)
        fitter.Draw(fitCfg['draw_recipes'])
        cFit.SaveAs('cFit.pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', default='cfg_fit.yml', nargs='?')
    args = parser.parse_args()

    with open(args.cfg) as file:
        config = yaml.safe_load(file)

    gROOT.SetBatch(True)

    FitCF(config)
