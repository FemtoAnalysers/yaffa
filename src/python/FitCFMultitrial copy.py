# pylint: skip-file
'''
Script to fit femtoscopic correlation functions.
'''

import sys
import os
import argparse
import yaml
import tabulate


from yaffa import utils

def FitCF(cfgFits, cfgMultitrials): # pylint disable:missing-function-docstring
    '''Fit the correlation functions.

    Args:
        cfg (dict): configuration of the fit
    '''

    oFile = TFile(f"{cfgFits['ofile']}.root", 'recreate')

    for cfgFit, cfgMultitrial in zip(cfgFits['fits'], cfgMultitrials['fits']):
        multifitter = SuperFitterMultitrial()
        for inFileName in cfgMultitrial['infiles']:
            inFile = TFile(inFileName)
            hObs = utils.io.Load(inFile, cfgFit['path'])
            hObs.SetDirectory(0)
            oObs = Observable(hObs)
            inFile.Close()

            multifitter.AddCF(oObs)

            # break
        
        # for template

        # oObs, cfgFit['fitrange'][0][0], cfgFit['fitrange'][0][1])
        multifitter.SetDrawRange(*cfgFit['drawrange'])

        print("receee", cfgFit['draw_recipes'])
        multifitter.SetDrawRecipes(cfgFit['draw_recipes'])

        # Add template to the multifitter
        for iTerm, term in enumerate(cfgFit['terms']):
            if templFileName := term.get('file'):
                templFile = TFile(templFileName)
                template = utils.io.Load(templFile, term['path'])
                if isinstance(template, TH1):
                    hTemplate = utils.analysis.ChangeUnits(template, term.get('unit_mult', 1))
                    hTemplate.SetDirectory(0)
                    multifitter.Add(term['name'], hTemplate)
                else:
                    print('Type not implemented. Exit!')
                    sys.exit()
                templFile.Close()
            else:
                multifitter.Add(term['name'], term['func'])

            multifitter.SetPars(term['name'], term['params'])
        multifitter.FitMultitrials(cfgFit['model'], 'MR+')
        # multifitter.Fit(cfgFit['model'], 'MR+')

        # cFit = TCanvas('cFit', '', 600, 600)
        # cFit.DrawFrame(*cfgFit['frame'], ';#it{k}* (GeV/c);#it{C}(#it{k}*)')
        # cFit.SaveAs(f'{cfg["ofile"]}.pdf')
        # cFit.SaveAs(f'{cfg["ofile"]}.root')

    oFile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfgFit')
    parser.add_argument('cfgMultitrial')
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('-x', default=False, action='store_true', help='plot the canvas')
    args = parser.parse_args()

    utils.style.SetStyle()

    from ROOT import TFile, TCanvas, gInterpreter, gROOT, TH1
    gInterpreter.ProcessLine(f'#define DO_DEBUG {1 if args.debug else 0}')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/Observable.h"')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/SuperFitterMultitrial.h"')
    from ROOT import Observable, SuperFitterMultitrial # plint: disable=ungrouped-imports

    with open(args.cfgFit) as file:
        configFit = yaml.safe_load(file)

    with open(args.cfgMultitrial) as file:
        configMultitrial = yaml.safe_load(file)

    gROOT.SetBatch(not args.x)

    FitCF(configFit, configMultitrial)
