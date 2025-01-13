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

def FitCF(cfg): # pylint disable:missing-function-docstring
    '''Fit the correlation functions.

    Args:
        cfg (dict): configuration of the fit
    '''

    # oFile = TFile('ancestors_LPiplus.root', 'create')


    for fitCfg in cfg['fits']:
        inFile = TFile(fitCfg['infile'])
        hObs = utils.io.Load(inFile, fitCfg['path'])
        hObs.SetDirectory(0)
        oObs = Observable(hObs)
        inFile.Close()

        fitter = SuperFitter(oObs, fitCfg['fitrange'][0][0], fitCfg['fitrange'][0][1])

        # Add template to the fitter
        for iTerm, term in enumerate(fitCfg['terms']):
            if templFileName := term.get('file'):
                templFile = TFile(templFileName)
                template = utils.io.Load(templFile, term['path'])
                if isinstance(template, TH1):
                    hTemplate = utils.analysis.ChangeUnits(template, term.get('unit_mult', 1))
                    hTemplate.SetDirectory(0)
                    fitter.Add(term['name'], hTemplate, term['params'])

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
                    print('Type not implemented. Exit!')
                    sys.exit()
                templFile.Close()
            else:
                fitter.Add(term['name'], term['func'], term['params'])

        fitter.Fit(fitCfg['model'], 'MR+')

        cFit = TCanvas('cFit', '', 600, 600)
        cFit.DrawFrame(0, 0.98, fitCfg['drawrange'][1] * 0.9999, 1.11999, ';#it{k}* (GeV/c);#it{C}(#it{k}*)')
        fitter.Draw(fitCfg['draw_recipes'])
        cFit.SaveAs(f'{cfg["ofile"]}.pdf')

        # Save fit parameters to file
        fFit = fitter.GetFitFunction()
        colNames = ['Parameter', 'Name', 'Value', 'Error', 'Step size', 'Derivative']
        pars = []
        for iPar in range(fFit.GetNpar()):
            pars.append([
                iPar,
                fFit.GetParName(iPar),
                f'{fFit.GetParameter(iPar):+10.5e}',
                f'{fFit.GetParError(iPar):+10.5e}',
                '-',  # Placeholder for "Step size"
                '-'   # Placeholder for "Derivative"
            ])        
        table = tabulate.tabulate(pars, headers=colNames, tablefmt='pipe', floatfmt=".5e")  # "grid" is one of many styles
        with open(f'{cfg["ofile"]}_parameters.txt', "w") as file:
            file.write(table)
        # oFile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', default='cfg_fit.yml', nargs='?')
    parser.add_argument('--debug', default=False, action='store_true')
    args = parser.parse_args()

    utils.style.SetStyle()

    from ROOT import TFile, TCanvas, gInterpreter, gROOT, TH1, TH1D, TF1
    gInterpreter.ProcessLine(f'#define DO_DEBUG {1 if args.debug else 0}')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/Observable.h"')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/SuperFitter.h"')
    from ROOT import Observable, SuperFitter # plint: disable=ungrouped-imports

    with open(args.cfg) as file:
        config = yaml.safe_load(file)

    gROOT.SetBatch(True)

    FitCF(config)
