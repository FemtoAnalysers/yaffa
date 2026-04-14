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

    oFile = TFile(f"{cfg['ofile']}.root", 'recreate')

    sfmt = SuperFitterMultitrial()
    for fitCfg in cfg['fits']:
        sfmt.SetDrawRange(*fitCfg['drawrange'])

        for inFileName in fitCfg['infiles']:
            inFile = TFile(inFileName)
            hObs = utils.io.Load(inFile, fitCfg['path'])
            hObs.SetDirectory(0)
            oObs = Observable(hObs)
            inFile.Close()

            sfmt.AddCF(oObs)
        
        # = SuperFitter(oObs, fitCfg['fitrange'][0][0], fitCfg['fitrange'][0][1])

    #     # Add template to the sfmt
        for iTerm, term in enumerate(fitCfg['terms']):
            if templFileName := term.get('file'):
                templFile = TFile(templFileName)
                template = utils.io.Load(templFile, term['path'])
                if isinstance(template, TH1):
                    hTemplate = utils.analysis.ChangeUnits(template, term.get('unit_mult', 1))
                    hTemplate.SetDirectory(0)
                    sfmt.Add(term['name'], hTemplate, term['params'])
                else:
                    print('Type not implemented. Exit!')
                    sys.exit()
                templFile.Close()
            else:
                sfmt.Add(term['name'], term['func'], term['params'])

        # sfmt.AddSourceParameters([[1, 2], [3, 4]])
        # sfmt.AddSourceParameters([[1, 2], [3, 4]])
        # sfmt.AddFitRange([[1, 2], [3, 4]])
        sfmt.FitMultitrials(fitCfg['model'], 'MR+')

    #     cFit = TCanvas('cFit', '', 600, 600)
    #     cFit.DrawFrame(*fitCfg['frame'], ';#it{k}* (GeV/c);#it{C}(#it{k}*)')
    #     sfmt.Draw(fitCfg['draw_recipes'])
    #     cFit.SaveAs(f'{cfg["ofile"]}.pdf')
    #     cFit.SaveAs(f'{cfg["ofile"]}.root')

    #     # Save fit parameters to file
    #     fFit = sfmt.GetFitFunction()
    #     colNames = ['Parameter', 'Name', 'Value', 'Error', 'Step size', 'Derivative']
    #     pars = []
    #     for iPar in range(fFit.GetNpar()):
    #         pars.append([
    #             iPar,
    #             fFit.GetParName(iPar),
    #             f'{fFit.GetParameter(iPar):+10.5e}',
    #             f'{fFit.GetParError(iPar):+10.5e}',
    #             '-',  # Placeholder for "Step size"
    #             '-'   # Placeholder for "Derivative"
    #         ])        
    #     table = tabulate.tabulate(pars, headers=colNames, tablefmt='pipe', floatfmt=".5e")  # "grid" is one of many styles
    #     with open(f'{cfg["ofile"]}_parameters.txt', "w") as file:
    #         file.write(table)
    oFile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', default='cfg_fit.yml', nargs='?')
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('-x', default=False, action='store_true', help='plot the canvas')
    args = parser.parse_args()

    # utils.style.SetStyle()

    from ROOT import TFile, TCanvas, gInterpreter, gROOT, TH1
    gInterpreter.ProcessLine(f'#define DO_DEBUG {1 if args.debug else 0}')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/Observable.h"')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/SuperFitterMultitrial.h"')
    from ROOT import Observable, SuperFitterMultitrial # plint: disable=ungrouped-imports

    with open(args.cfg) as file:
        config = yaml.safe_load(file)

    # gROOT.SetBatch(not args.x)

    FitCF(config)
