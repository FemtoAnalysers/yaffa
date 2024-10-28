'''
Script to compbine different templates based on branching ratios, acceptance and thermal weights

Usage:
python3 BuildCocktail.py cfg.yml
'''

import argparse
import itertools
import math
import yaml
import numpy as np

from ROOT import TFile, TCanvas  # pylint: disable=no-name-in-module

from yaffa import logger as log

def ApplyUncertainty(*args): # pylint: disable=inconsistent-return-statements
    '''
    Applies the uncertainty to a value

    Parameters
    ----------
    args : tuple
        Variable length argument list. Each argument should be of type `float`. The first value is considered to be the
        "central" one, the second is lower uncertainty and the third is the upper uncertainty. If the third value is
        missing, the second value is assumed to be a symmetric uncertainty. If the second value is missing, `value` is
        interpreted to be without uncertainty.

    Returns
    -------
    list
        The value and its variations within uncertainties

    Returns
    -------
    Examples:
        no uncertainty:         [10]       -> [10]
        symmetric uncertainty:  [10, 1]    -> [10, 9, 11]
        asymmetric uncertainty: [10, 1, 2] -> [10, 9, 12]
    '''

    if len(args) == 1:
        return args
    if len(args) == 2:
        return args[0], args[0] - abs(args[1]), args[0] + abs(args[1])
    if len(args) == 3:
        return args[0], args[0] - abs(args[1]), args[0] + abs(args[2])

    log.critical('Invalid number of parameters')

def BuildCocktail():
    '''
    Function to build a cocktail using different tempalates and weighting factors.
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', default='')
    parser.add_argument('--central', default=False, action='store_true')
    args = parser.parse_args()

    # Load yaml file
    with open(args.cfg, 'r') as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError:
            log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

    templates = []
    # Load histograms
    for iEntry, entry in enumerate(cfg['cocktail']):
        inFile = TFile(entry['infile'])
        templates.append({})
        templates[-1]['template'] = inFile.Get(entry['name'])
        templates[-1]['template'].SetDirectory(0)
        templates[-1]['template'].Rebin(cfg['rebin'])
        templates[-1]['template'].Scale(1./templates[-1]['template'].GetEntries())
        templates[-1]['weight'] = entry['weight']
        templates[-1]['acceptance'] = entry['acceptance']

    # Convert BR from percentage to absolute values
    for iEntry, entry in enumerate(cfg['cocktail']):
        for iChn, chn in enumerate(entry['bratio']):
            if args.central:
                cfg['cocktail'][iEntry]['bratio'][iChn] = [cfg['cocktail'][iEntry]['bratio'][iChn][0]]
            cfg['cocktail'][iEntry]['bratio'][iChn] = [value / 100 for value in cfg['cocktail'][iEntry]['bratio'][iChn]]

    # Convert BR uncertainty to upper and lower values
    for iEntry, entry in enumerate(cfg['cocktail']):
        for iChn, chn in enumerate(entry['bratio']):
            uncs = list(ApplyUncertainty(*chn)) # Substitute [value, unc] with [value, value-unc, value+unc]
            cfg['cocktail'][iEntry]['bratio'][iChn] = uncs

        # Compute all the possible combinations of branching ratios
        brs = itertools.product(*cfg['cocktail'][iEntry]['bratio'])

        # Compute the total BRs
        cfg['cocktail'][iEntry]['bratio'] = [math.prod(br) for br in brs]

    oFile = TFile(cfg['ofile'], 'recreate')

    cSummary = TCanvas("cSummary", "", 600, 600)
    variations = []
    # Make cocktails
    for iBR, br in enumerate(itertools.product(*(comp['bratio'] for comp in cfg['cocktail']))):
        print('Making cocktail with BRs:', br)

        hCocktail = templates[0]['template'].Clone(f'hCocktail{iBR}')
        hCocktail.Reset()
        for template, bratio in zip(templates, br):
            hCocktail.Add(template['template'], bratio * template['weight'] * template['acceptance'])

        variations.append([hCocktail.GetBinContent(iBin + 1) for iBin in range(hCocktail.GetNbinsX())])
        hCocktail.Write()

    hSummary = hCocktail.Clone('hSummary')
    hSummary.Reset()
    hSummary.SetFillColor(4)
    hSummary.SetLineColor(4)
    hSummary.SetMarkerColor(4)

    averages = np.mean(np.array(variations), axis=0)
    errors = np.std(np.array(variations), axis=0)
    for iPoint, (avg, error) in enumerate(zip(averages, errors)):
        hSummary.SetBinContent(iPoint + 1, avg)
        hSummary.SetBinError(iPoint + 1, error)
    hSummary.Write()
    hSummary.Draw('e3')
    cSummary.Write()
    oFile.Close()

    print(f'Output saved in {cfg["ofile"]}')

if __name__ == '__main__':
    BuildCocktail()
