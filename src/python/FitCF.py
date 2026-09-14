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
from yaffa import logger as log

from dotenv import load_dotenv
from pathlib import Path

env_path = Path(__file__).resolve().parents[2] / ".env"
print(f'Loading env from {env_path}')
if not load_dotenv(dotenv_path=env_path, verbose=True, override=True):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
if not YAFFA_PATH:
    print("\033[33mWARNING: Path to yaffa is empty, something might break!\033[0m")




def FitCF(cfg): # pylint disable:missing-function-docstring
    '''Fit the correlation functions.

    Args:
        cfg (dict): configuration of the fit
    '''

    terms = {}
    hObsList = []
    # oFile = TFile('ancestors_LPiplus.root', 'create')
    fitter = SuperFitter()
    fitter.SetFitRange(cfg['fits'][0]['fitrange'])
    fitter.SetDrawRange(*cfg['fits'][0]['drawrange'])

    for iFit, fitCfg in enumerate(cfg['fits']):
        inFile = TFile(fitCfg['infile'])
        obs = utils.io.Load(inFile, fitCfg['path'])
        if isinstance(obs, TGraphMultiErrors):
            # Only the statistical uncertainties enter the fit. The x errors are dropped, otherwise ROOT uses the effective chi2
            hObs = TGraphErrors(obs.GetN())
            hObs.SetName(obs.GetName())
            for iPoint in range(obs.GetN()):
                hObs.SetPoint(iPoint, obs.GetPointX(iPoint) * fitCfg.get('unit_mult', 1), obs.GetPointY(iPoint))
                hObs.SetPointError(iPoint, 0, obs.GetErrorY(iPoint, 0))
        else:
            hObs = utils.analysis.ChangeUnits(obs, fitCfg.get('unit_mult', 1))
            hObs.SetDirectory(0)
        hObsList.append(hObs)
        oObs = Observable(hObs)
        inFile.Close()

        fitter.AddObservable(oObs)

        # Add template to the fitter
        for iTerm, term in enumerate(fitCfg['terms']):
            if templFileName := term.get('file'):
                templFile = TFile(templFileName)
                template = utils.io.Load(templFile, term['path'])
                if isinstance(template, TH1):
                    hTemplate = utils.analysis.ChangeUnits(template, term.get('unit_mult', 1))
                    hTemplate.SetDirectory(0)
                    fitter.Add(iFit, term['name'], hTemplate, term['params'])
                elif isinstance(template, (TGraph, TGraphErrors)):
                    fitter.Add(iFit, term['name'], template, term['params'], term.get('unit_mult', 1))
                elif isinstance(template, TF1):
                    fitter.Add(iFit, term['name'], template, term['params'], 1)

                # elif isinstance(template, TF1):
                    # Ccnvert to hist
                    # hTemplate = TH1D(f"h{iTerm}", "", 500, 0, 2)
                    # for iBin in range(500):
                    #     bc = template.Eval(hTemplate.GetBinCenter(iBin + 1) * term.get('unit_mult', 1))
                    #     print("bc", bc)
                    #     hTemplate.SetBinContent(iBin + 1, bc)
                    # hTemplate.SetDirectory(0)

                    # c = TCanvas()
                    # hTemplate.Draw()

                    # hTemplate.Write()
                    # # oFile.cd()
                    # hTemplate.Write(f'hCF_{iTerm}')
                    # c.SaveAs(f'test{iTerm}.png')


                    # unitMult = term.get('unit_mult', 1)
                    # SetOwnership(template, False)
                    # template.SetName(f'f{iTerm}')

                    # print('zezez', template)
                    # fitter.Add(term['name'], template, term['params'], term.get('unit_mult', 1))
                    # fitter.Add(term['name'], hTemplate, term['params'])
                else:
                    raise ValueError("Type not implemented")
                templFile.Close()
            else:
                fitter.Add(iFit, term['name'], term['func'], term['params'])

        fitter.SetModel(iFit, fitCfg['model'])
    fitter.Fit('MR+')

    oFileName = cfg["ofile"]
    if suffix := cfg['suffix']:
        oFileName = f'{oFileName}_{suffix}'

    oFile = TFile(f'{oFileName}.root', 'recreate')
    panels = utils.style.GetNPanels(len(cfg['fits']))
    cFit = TCanvas('cFit', '', 600 * panels[0], 600 * panels[1])
    cFit.Divide(*panels)
    for iFit in range(len(cfg['fits'])):
        cFit.cd(iFit + 1)
        cFit.DrawFrame(*cfg['fits'][iFit]['frame'], ';#it{k}* (GeV/c);#it{C}(#it{k}*)')
        fitter.Draw(iFit, cfg['fits'][iFit]['draw_recipes'], cfg['fits'][iFit]['datalabel'], cfg['fits'][iFit]['header'])
        terms[iFit] = fitter.GetTerms()
        fitter.GetFitFunction(iFit).Write()

    cFit.SaveAs(f'{oFileName}.pdf')

    # Save fit parameters to file
    fFit = fitter.GetFitFunction()
    colNames = ['Parameter', 'Name', 'Value', 'Error', 'Step size', 'Derivative']
    pars = []
    result = {
        'chi2ndf': 999 # fFit.GetChisquare() / fFit.GetNDF(),
    }
    for iPar in range(fFit.GetNpar()):
        result[fFit.GetParName(iPar)] = fFit.GetParameter(iPar)
        pars.append([
            iPar,
            fFit.GetParName(iPar),
            f'{fFit.GetParameter(iPar):+10.5e}',
            f'{fFit.GetParError(iPar):+10.5e}',
            '-',  # Placeholder for "Step size"
            '-'   # Placeholder for "Derivative"
        ])

    a0_re_idx = fFit.GetParNumber('a0_re')
    a0_im_idx = fFit.GetParNumber('a0_im')

    gScatLen = TGraphErrors(1)
    gScatLen.SetName('gScatLen')
    gScatLen.SetPoint(0, fFit.GetParameter(a0_re_idx), fFit.GetParameter(a0_im_idx))
    gScatLen.SetPointError(0, fFit.GetParError(a0_re_idx), fFit.GetParError(a0_im_idx))
    gScatLen.Write()

    table = tabulate.tabulate(pars, headers=colNames, tablefmt='pipe', floatfmt=".5e")  # "grid" is one of many styles
    with open(f'{oFileName}_parameters.txt', "w") as file:
        file.write(table)

    for hObs in hObsList:
        hObs.Write()
    for idx, _ in enumerate(cfg['fits']):
        gGenCF = fitter.GetGenuineCF(idx, cfg['fits'][idx]['gencf']) # explicit cast to int for some reason
        gGenCF.SetName(f'gGenCF{idx}')
        gGenCF.Write()
    cFit.Write()

    for iTerm, termlist in terms.items():
        oFile.mkdir(f'term{iTerm}')
        oFile.cd(f'term{iTerm}')
        for t in termlist:
            t.Write()
    oFile.Close()

    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', default='cfg_fit.yml', nargs='?')
    parser.add_argument('--debug', type=int, default=10)
    parser.add_argument('-x', default=False, action='store_true', help='plot the canvas')
    args = parser.parse_args()

    from ROOT import TF1, TFile, TCanvas, gInterpreter, gROOT, TH1, TGraph, TGraphErrors, TGraphMultiErrors
    gInterpreter.ProcessLine(f'#undef DEBUG_LEVEL')
    gInterpreter.ProcessLine(f'#define DEBUG_LEVEL {args.debug}')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/src/cpp/Observable.h"')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/src/cpp/SuperFitter.h"')
    from ROOT import Observable, SuperFitter # plint: disable=ungrouped-imports

    

    utils.style.SetStyle()

    with open(args.cfg) as file:
        config = yaml.safe_load(file)

    gROOT.SetBatch(not args.x)

    FitCF(config)
