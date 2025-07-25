# pylint: skip-file

'''
Script to compute the raw correlation function.
The output file is RawCF_suffix.root

Usage:
python3 ComputeRawCF.py cfg.yml

'''
import os
import re
import argparse
import yaml
from rich import print # pylint: disable=redefined-builtin

from ROOT import TFile, TH1D, TH2D

from yaffa import logger as log
from yaffa.utils.io import Load, GetKeyNames


def LoadMultVsKstarAncestorFemtoDream(inFile, **kwargs):
    '''Load Mult vs k* ancestors from FD'''
    suffix = kwargs['suffix']
    regions = kwargs['regions']
    combs = kwargs['combs']

    if regions != ['sgn']:
        log.fatal('Only signal region is implemented.')

    if combs != ['p02', 'p03', 'p12', 'p13']:
        log.fatal('Combinations not implemented')

    hDistr = {}
    for comb in combs:
        hDistr[comb] = {}
        for region in regions:
            log.info('Loading %s %s', comb, region)
            fdcomb = f'Particle{comb[1]}_Particle{comb[2]}'
            folder = f'HMResults{suffix}/HMResults{suffix}/{fdcomb}/'

            for ancestor in ['Common', 'NonCommon']:
                hAncestor = Load(inFile, f'{folder}/SEMultDist{ancestor}_{fdcomb}')
                if hAncestor == None: # pylint: disable=singleton-comparison
                    log.info('Ancestors were not found')
                    return {}

                log.info('Loading %s ancestor', ancestor)
                hDistr[comb][f'{region}/{ancestor}'] = hAncestor
                hDistr[comb][f'{region}/{ancestor}'].SetDirectory(0)

            iMt = 0
            while True:
                hSECommon = Load(inFile, f'{folder}/SEmTMultCommon_{iMt}_{fdcomb}')
                hSENonCommon = Load(inFile,f'{folder}/SEmTMultNonCommon_{iMt}_{fdcomb}')

                if hSECommon == None or hSENonCommon == None:  # pylint: disable=singleton-comparison
                    break

                log.info('Loading mT %d bins', iMt)

                hDistr[comb][f'{region}/mT{iMt}/Common'] = hSECommon
                hDistr[comb][f'{region}/mT{iMt}/Common'].SetDirectory(0)
                hDistr[comb][f'{region}/mT{iMt}/NonCommon'] = hSENonCommon
                hDistr[comb][f'{region}/mT{iMt}/NonCommon'].SetDirectory(0)

                iMt += 1

    return hDistr




def LoadMultVsKstarFemtoDream(inFile, **kwargs):
    '''Load Mult vs k* from FD'''

    suffix = kwargs['suffix']
    regions = kwargs['regions']
    combs = kwargs['combs']

    if regions != ['sgn']:
        log.fatal('Only signal region is implemented.')

    if combs != ['p02', 'p03', 'p12', 'p13']:
        log.fatal('Combinations not implemented')

    hSE = {}
    hME = {}
    for comb in combs:
        hSE[comb] = {}
        hME[comb] = {}
        for region in regions:
            log.info('Loading %s %s', comb, region)

            fdcomb = f'Particle{comb[1]}_Particle{comb[2]}'
            folder = f'HMResults{suffix}/HMResults{suffix}/{fdcomb}/'

            hSE[comb][region] = Load(inFile, f'{folder}/SEMultDist_{fdcomb}')
            hME[comb][region] = Load(inFile, f'{folder}/MEMultDist_{fdcomb}')
            hSE[comb][region].SetDirectory(0)
            hME[comb][region].SetDirectory(0)


            iMt = 0
            while True:
                hSEmT = Load(inFile, f'{folder}/SEmTMult_{iMt}_{fdcomb}')
                hMEmT = Load(inFile, f'{folder}/MEmTMult_{iMt}_{fdcomb}')

                if hSEmT == None or hMEmT == None: # pylint: disable=singleton-comparison
                    break

                log.info('Loading mT %d bins', iMt)

                hSE[comb][f'{region}/mT{iMt}'] = hSEmT
                hSE[comb][f'{region}/mT{iMt}'].SetDirectory(0)
                hME[comb][f'{region}/mT{iMt}'] = hMEmT
                hME[comb][f'{region}/mT{iMt}'].SetDirectory(0)

                iMt += 1
    return hSE, hME


def LoadMultVsKstarPythia(inFile, **kwargs):
    '''Load Mult vs k* from Pythia simulations'''
    
    suffix = kwargs['suffix']
    regions = kwargs['regions']
    combs = kwargs['combs']

    if regions != ['sgn']:
        log.fatal('Only signal region is implemented.')

    hSE = {}
    hME = {}
    for comb in combs:
        hSE[comb] = {}
        hME[comb] = {}
        for region in regions:
            log.info('Loading %s %s', comb, region)

            hSE[comb][region] = Load(inFile, f'{comb}/hSE')
            hME[comb][region] = Load(inFile, f'{comb}/hME')
            hSE[comb][region].SetDirectory(0)
            hME[comb][region].SetDirectory(0)

    return hSE, hME


def LoadMultVsKstar(inFile, **kwargs):
    ''' Load Mult vs k*'''
    suffix = kwargs['suffix']
    hSE = {}
    hME = {}

    print(kwargs)
    print(GetKeyNames(inFile))
    if f'HMResults{suffix}' in GetKeyNames(inFile): # Make correlation functions from FemtoDream
        hSE, hME = LoadMultVsKstarFemtoDream(inFile, **kwargs)
        hAnc = LoadMultVsKstarAncestorFemtoDream(inFile, **kwargs)

        if hAnc:
            for comb in hSE:
                hSE[comb].update(hAnc[comb])
    else:
        hSE, hME = LoadMultVsKstarPythia(inFile, **kwargs)

    if not hSE:
        raise ValueError('Same Event could not be loaded properly')

    if not hME:
        raise ValueError('Mixed Event could not be loaded properly')
    return hSE, hME

def Reweight(hSE, hME, normRange = None, multRange = None, name=None):
    '''reweight'''

    suffix = f"_{name}" if name else ''

    hSEMultSlices = []
    hMEMultSlices = []
    hCFMultSlices = []

    nBins = hME.GetXaxis().GetNbins()
    xMin = hME.GetXaxis().GetXmin()
    xMax = hME.GetXaxis().GetXmax()
    hMERew = TH1D(f'hMERew{suffix}', ';#it{k}* (GeV/#it{c});Counts', nBins, xMin, xMax)

    nBins = hME.GetYaxis().GetNbins() + 2
    xMin = hME.GetYaxis().GetXmin()
    xMax = hME.GetYaxis().GetXmax()
    hWeights = TH1D(f'hWeights{suffix}', ';Mult bin (a.u.); Weight', nBins, xMin, xMax)

    startBinMultRew, endBinMultRew = multRange if multRange is not None else [0, nBins]
    for iBin in range(startBinMultRew, endBinMultRew): # Loop over underflow, all bins, and overflow
        hSEbinmult = hSE.ProjectionX(f'hSEdistr_{iBin}', iBin, iBin)
        hMEbinmult = hME.ProjectionX(f'hMEdistr_{iBin}', iBin, iBin)
        hSEMultSlices.append(hSEbinmult)
        hMEMultSlices.append(hMEbinmult)

        firstBin = hSEbinmult.FindBin(normRange[0] * 1.0001)
        lastBin = hSEbinmult.FindBin(normRange[1] * 0.9999)
        if hSEbinmult.Integral(firstBin, lastBin) > 0 and hMEbinmult.Integral() > 0:
            hMERew.Add(hMEbinmult, hSEbinmult.Integral()/hMEbinmult.Integral())
            hWeights.SetBinContent(iBin, hSEbinmult.Integral()/hMEbinmult.Integral())

            # Compute the CFs for each multiplicity bin
            norm = hMEbinmult.Integral(firstBin, lastBin) / hSEbinmult.Integral(firstBin, lastBin)
            hCFbinmult = norm * hSEbinmult / hMEbinmult
            hCFbinmult.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
            hCFMultSlices.append(hCFbinmult)
        else:
            hCFMultSlices.append(hSEbinmult.Clone(f'hCF_multbin{iBin}'))
            hCFMultSlices[-1].Reset()

    return hMERew, hWeights, (hSEMultSlices, hMEMultSlices, hCFMultSlices)


def ProjectDistr(hDistrMult):
    '''Project distributions
    '''
    log.info('Projecting distributions...')
    hDistr = {}
    for comb in hDistrMult:
        hDistr[comb] = {}

        for region in hDistrMult[comb]:
            hDistr[comb][region] = hDistrMult[comb][region].ProjectionX(f'{comb}SEdistr', 1, hDistrMult[comb][region].GetNbinsX())
            hDistr[comb][region].SetDirectory(0)

    return hDistr


def Combine(hDistr, recipes):
    '''
    Sums histograms based on the given recipes.
    
    Args:
        hDistr (dict(dict(TH1))): histograms to be summed

        recipes dict(list(str)): instructions on how to combine the histograms
    '''

    print(recipes)
    print(hDistr)
    
    for target, (comb1, comb2) in recipes.items():
        hDistr[target] = {}
        for region in hDistr[comb1]:
            log.info("Summing %s = %s + %s for region '%s'", target, comb1, comb2, region)
            hDistr[target][region] = hDistr[comb1][region] +  hDistr[comb1][region]

    return hDistr
    
    
def main(cfg): # pylint: disable=too-many-statements
    '''
    main function

    Args:
        cfg (dict): configuration for the calculation of the CFs
    '''

    regions = ['sgn']
    sumRecipe = cfg['combine']

    inFile = TFile(cfg['infile'])
    combs = [comb for comb in GetKeyNames(inFile) if re.match('p[0-9][0-9]', comb)]
    if not combs:
        combs = ['p02', 'p03', 'p12', 'p13']

    # Define the output file
    oFileBaseName = 'RawCF'
    if cfg['suffix'] != '' and cfg['suffix'] is not None:
        oFileBaseName += f'_{cfg["suffix"]}'
    oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

    # Open the output file
    oFile = TFile.Open(oFileName, 'recreate')

    # 5 bins: sgn, sgn_Common, sgn_NonCommon, sbl, sbr. If some are not present they are left empty
    nBins = len(combs) + int(len(combs)/2)
    hFemtoPairs = TH2D("hFemtoPairs", "Pairs in femto region", 5, 0, 5, nBins, 0, nBins)
    hFemtoPairs.GetXaxis().SetBinLabel(1, "sgn")
    hFemtoPairs.GetXaxis().SetBinLabel(2, "sgn/Common")
    hFemtoPairs.GetXaxis().SetBinLabel(3, "sgn/NonCommon")
    hFemtoPairs.GetXaxis().SetBinLabel(4, "sbl")
    hFemtoPairs.GetXaxis().SetBinLabel(5, "sbr")
    hFemtoPairs.GetYaxis().SetLabelSize(50)
    hFemtoPairs.GetXaxis().SetLabelSize(50)

    # Load 2D Mult-vs-kstar same- and mixed-event distributions
    log.info("Loading Mult vs k* histograms")
    hSEmultk, hMEmultk = LoadMultVsKstar(inFile, suffix=cfg['runsuffix'], regions=regions, combs=combs)

    log.info("Computing mult-reweighted ME")
    hMErew = {}
    hWeightsRew = {}
    for comb in combs:
        oFile.mkdir(comb)
        oFile.cd(comb)

        hMErew[comb] = {}
        hWeightsRew[comb] = {}
        for region in hSEmultk[comb]:
            log.info("Reweight for %s %s", comb, region)
            regionME = region.replace('/Common', '').replace('/NonCommon', '')
            log.info('reweight %s %s', comb, region)
            hMERew, hWeights, slices = Reweight(hSEmultk[comb][region], hMEmultk[comb][regionME], cfg['norm'], name=f'{comb}_{region}')

            for iSlice, (hSE, hME, hCF) in enumerate(zip(*slices)):
                log.info("Slice: %d", iSlice)
                
                oFile.mkdir(f'{comb}/{region}/multbins/{iSlice}')
                oFile.cd(f'{comb}/{region}/multbins/{iSlice}')

                hCF.Write(f'hCF_multbin{iSlice}')
                hSE.Write(f'hSE_multbin{iSlice}')
                hME.Write(f'hME_multbin{iSlice}')

                oFile.cd(comb)

            hMErew[comb][region] = hMERew
            hWeightsRew[comb][region] = hWeights

            
    log.info("Projecting 2D histograms")
    hSE = ProjectDistr(hSEmultk)
    hME = ProjectDistr(hMEmultk)

    print(hSE)
    print(hME)

    hSE = Combine(hSE, sumRecipe)
    hME = Combine(hME, sumRecipe)
    hMErew = Combine(hMErew, sumRecipe)
    hWeightsRew = Combine(hWeightsRew, sumRecipe)

    print(hSE)
    print(hME)

    # Compute the CF and write to file
    for iComb, comb in enumerate(hSE):
        for region in hSE[comb]:
            rebin = round(float(cfg['binwidth']) / (hSE[comb][region].GetBinWidth(1) * 1000))
            hSE[comb][region].Rebin(rebin)
            # remove '/Common' and '/NonCommon' since the ME is the same as the inclusive ancestor
            regionME = region.replace('/Common', '').replace('/NonCommon', '')
            hME[comb][regionME].Rebin(rebin)
            hMErew[comb][regionME].Rebin(rebin)

            # Count number of pairs in the femto region
            regionBin = hFemtoPairs.GetXaxis().FindBin(region)
            firstBin = hSE[comb][region].FindBin(0.0001)
            lastBin = hSE[comb][region].FindBin(0.2*0.9999)
            hFemtoPairs.GetYaxis().SetBinLabel(iComb + 1, comb)
            hFemtoPairs.SetBinContent(regionBin + 1, iComb + 1, hSE[comb][region].Integral(firstBin, lastBin))

            if cfg['norm'] is None:  # if not specified, normalize to the yields
                firstBinNorm = 1
                lastBinNorm = hMErew[comb][region].GetNbinsX()
            elif isinstance(cfg['norm'], list):
                firstBinNorm = hSE[comb][region].FindBin(cfg['norm'][0]*1.0001)
                lastBinNorm = hSE[comb][region].FindBin(cfg['norm'][1]*0.9999)
            else:
                log.critical('Normalization method not implemented')

            YieldSE = hSE[comb][region].Integral(firstBinNorm, lastBinNorm)
            norm = hME[comb][regionME].Integral(firstBinNorm, lastBinNorm) / YieldSE
            normRew = hMErew[comb][regionME].Integral(firstBinNorm, lastBinNorm) / YieldSE

            hSE[comb][region].Sumw2()
            hME[comb][regionME].Sumw2() # Just to trigger the same draw option as for hSE

            hCF = norm * hSE[comb][region] / hME[comb][regionME]
            hCFrew = normRew * hSE[comb][region] / hMErew[comb][regionME]

            oFile.mkdir(f'{comb}/{region}')
            oFile.cd(f'{comb}/{region}')

            hSE[comb][region].SetName('hSE')
            hSE[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')
            hSE[comb][region].Write()

            hME[comb][regionME].SetName('hME')
            hME[comb][regionME].SetTitle(';#it{k}* (GeV/#it{c});Counts')
            hME[comb][regionME].Write()

            hMErew[comb][region].SetName('hMErew')
            hMErew[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')
            hMErew[comb][region].Write('hMErew')

            hWeightsRew[comb][regionME].SetName('hWeightsRew')
            hWeightsRew[comb][regionME].SetTitle(';Mult bin;Weight')
            hWeightsRew[comb][regionME].Write()

            hCF.SetName('hCF')
            hCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
            hCF.Write()

            hCFrew.SetName('hCFrew')
            hCFrew.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
            hCFrew.Write()

    oFile.cd()
    hFemtoPairs.Write()
    oFile.Close()
    print(f'output saved in {oFileName}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg', default='')
    parser.add_argument('--debug', default=False, action='store_true')
    args = parser.parse_args()

    if args.debug:
        log.setLevel(1)

    # Load yaml file
    with open(args.cfg, "r") as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

    main(config)
