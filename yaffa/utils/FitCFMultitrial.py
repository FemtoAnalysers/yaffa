# pylint: skip-file
'''
Script to fit femtoscopic correlation functions.
'''

import sys
import os
import argparse
import yaml
import random
from rich import print
from ROOT import TFile, TNtuple

from yaffa import utils
from yaffa.utils.FitCF import FitCF

random.seed(1)

def UpdateConfig(cfg, key, value):
    keys = key.split('.')
    cfgTmp = cfg

    for key in keys[:-1]:
        if 'fits-' in key:
            cfgTmp = cfgTmp['fits'][int(key[5:])]
        elif 'terms-' in key:
            cfgTmp = next((term for term in cfgTmp['terms'] if term['name'] == key[6:]))
        else:
            cfgTmp = cfgTmp[key]

    cfgTmp[keys[-1]] = value


def FitCFMultitrial(config, cfgVariations):
    '''Fit the correlation functions.

    Args:
        cfg (dict): configuration of the fit
    '''
    utils.style.SetStyle()

    oFile = cfgVariations.pop('ofile')
    suffix = cfgVariations.pop('suffix')
    if suffix:
        oFile = f'{oFile}_{suffix}'


    # Build variation of the default configuration
    NIter = 20
    results = []
    for iVar in range(NIter):
        print(f'Iteration {iVar}/{NIter}')
        # Update configuration
        for key, valueList, in cfgVariations.items():
            config.update({'ofile': oFile, 'suffix': f'var{iVar}'})
            UpdateConfig(config, key, random.choice(valueList))
        
        # Fit using varied config
        results.append(FitCF(config))
        results[-1].pop('LPi_a0_im')
        results[-1].pop('LPi_d0')
        results[-1].pop('LPi_r1')
        results[-1].pop('LPi_r2')
        results[-1].pop('LPi_lambda')
        results[-1].pop('LPiplus_p0')
        results[-1].pop('LPiplus_p1')
        results[-1].pop('LPiplus_p2')
        results[-1].pop('chi2ndf')
        
    file = TFile(f"{oFile}.root", "recreate")

    branches = results[0].keys()

    # Create TNtuple
    ntuple = TNtuple(
        "tResults",                    # name
        "Multitrial results",            # title
        ":".join(branches)
    )

    # Fill with data
    for result in results:
        ntuple.Fill(*[result[branch] for branch in branches if branch not in ['LPi_a0_im']])

    # Write to file
    ntuple.Write()
    file.Close()

    print(f'FOutput saved to {oFile}.root')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfgFit')
    parser.add_argument('cfgMultitrial')
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('-x', default=False, action='store_true', help='plot the canvas')
    args = parser.parse_args()
    
    from ROOT import gInterpreter, gROOT
    gInterpreter.ProcessLine(f'#define DEBUG_LEVEL {1 if args.debug else 0}')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/Observable.h"')
    gInterpreter.ProcessLine(f'#include "{os.environ.get("YAFFA")}/yaffa/utils/SuperFitter.h"')

    with open(args.cfgFit) as file:
        configFit = yaml.safe_load(file)

    with open(args.cfgMultitrial) as file:
        configMultitrial = yaml.safe_load(file)

    gROOT.SetBatch(not args.x)

    cfgVariations = configMultitrial.copy()
    del cfgVariations['ofile']
    del cfgVariations['suffix']
    FitCFMultitrial(configFit, configMultitrial)
