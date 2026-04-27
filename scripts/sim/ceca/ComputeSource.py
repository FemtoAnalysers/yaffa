import os
import sys
from dotenv import load_dotenv
from pathlib import Path
import numpy as np

env_path = Path(__file__).resolve().parent.parent.parent.parent / ".env"
if not load_dotenv(dotenv_path=env_path):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")

from ROOT import TFile, gDirectory, gInterpreter, TF1, gROOT, TGraphErrors, TGraphAsymmErrors
gInterpreter.Declare(f'#include "{YAFFA_PATH}/src/cpp/RootFunctions.hxx"')
from ROOT import SourceAAA, SourceCountsAAA, SourceCountsGauss

from yaffa import logger as log

sys.path.append(f'{YAFFA_PATH}/src/python')

from yaffa import utils, Chain, FitResult

def Fit(obj, func, range, pars, options=''):
    '''
    Create fit function
    '''
    
    fFit = TF1(f'fFit', func, *range, len(pars))
    fFit.SetNpx(10000)
    for iPar, par in enumerate(pars):
        fFit.SetParName(iPar, par[0])
        fFit.SetParameter(iPar, (par[1] + par[2]) / 2)
        fFit.SetParLimits(iPar, par[1], par[2])

    result = FitResult(obj.GetName(), obj.Fit(fFit, 'QLS'))
    if 'Q' not in options:
        print(result)

    return result


def Test2B3BConsistency(hRStarInTriplets, hRho):
    # r* of pairs in triplets
    norm = hRStarInTriplets.GetEntries() * hRStarInTriplets.GetBinWidth(1)
    resultRStar = Fit(hRStarInTriplets, SourceCountsGauss, [0, 8], [['norm', norm / 2, norm * 2] ,['r0', 0.1, 5]])
    avgRStar = hRStarInTriplets.GetMean()
    r0 = resultRStar.pars[1][1]
    expectedAvgRStar = r0 * 4 / np.sqrt(np.pi)
    delta = (expectedAvgRStar - avgRStar) / expectedAvgRStar
    print(f"r0 (fit)   = {r0:.2f} fm")
    print(f"<r*> (exp) = {expectedAvgRStar:.2f} fm")
    print(f"<r*>       = {avgRStar:.2f} fm")
    print(f"delta      = {delta * 100:.2f} %")

    # Hyper-radius of triplets
    norm = hRho.GetEntries() * hRho.GetBinWidth(1)
    resultRho = Fit(hRho, SourceCountsAAA, [0, 8], [['norm', norm / 2, norm * 2] ,['rho0', 0.1, 5]])
    avgRho = hRho.GetMean()
    rho0 = resultRho.pars[1][1]
    expectedAvgRho = rho0 * 15 * np.sqrt(np.pi) / 16
    delta = (expectedAvgRho - avgRho) / expectedAvgRho
    print(f"rho0 (fit)  = {rho0:.2f} fm")
    print(f"<rho> (exp) = {expectedAvgRho:.2f} fm")
    print(f"<rho>       = {avgRho:.2f} fm")
    print(f"delta       = {delta * 100:.2f} %")

    # Compatibility between 2B and 3B
    expectedAvgRho = avgRStar * 15 * np.pi / 32
    deviation = (avgRho - expectedAvgRho) / expectedAvgRho
    print(f"\n<rho> (exp) = {expectedAvgRho:.2f} fm")
    print(f"<rho>       = {avgRho:.2f} fm")
    print(f"delta       = {deviation * 100:.2f} %")


def GetMtScaling(hist, dirName):
    file = gDirectory.GetFile()
    if not file:
        log.error('No file is currently opened: output of mT scaling processing will not be saved!')

    if not file.IsWritable():
        log.error('Current file is not writable: output of mT scaling processing will not be saved!')


    nPoints = len(cfg['mt_bins']) - 1

    file.mkdir(dirName)
    file.cd(dirName)

    chRStar = Chain(utils.analysis.SliceVertically(hist, cfg['mt_bins'], name='hRStar_mT'))
    chRStar.Scale(1. / chRStar.Integral('width'))

    chRStar.Write()
    file.cd('../')

    chRStar.collect(lambda h : h.GetMean(), id='avgRStar')
    chRStar.collect(lambda h : h.GetMeanError(), id='avgRStarUnc')

    # TODO: Now using bin center for mT -> use center of gravity instead
    x = np.array([(cfg['mt_bins'][i] + cfg['mt_bins'][i + 1]) / 2 for i in range(nPoints)], dtype='double')
    y = np.array(chRStar.results['avgRStar'], dtype='double')
    yUnc = np.array(chRStar.results['avgRStarUnc'], dtype='double')

    # TODO: Now using standard mean error (assume gaussianity). Think of a better way (Bootstrap?)
    gMtScaling = TGraphAsymmErrors(nPoints, x, y, 0, 0, yUnc, yUnc)

    return gMtScaling

def main(cfg:dict):
    '''
    Compute the source size as a function of mT
    '''
    # Open input file
    inFile = TFile(cfg['infile'])
    oFile = TFile(cfg['ofile'], 'recreate')

    # Load histograms from CECA simulation
    hRStarVsMt = inFile.Get('hRStarVsMt')
    hRho = inFile.Get('hPhiVsRho').ProjectionX()
    hRStarInTriplets = inFile.Get('hRStarInTriplets')
    hFemtoRho = inFile.Get('hFemtoPhiVsRho').ProjectionX()
    hFemtoRStarFemtoPairsInTripletsVsMt = inFile.Get('hFemtoRStarFemtoPairsInTripletsVsMt')

    # hRho.SetName('hRho')
    # hFemtoRho.SetName('hFemtoRho')
    if hRStarVsMt and hRStarVsMt.GetEntries():
        gMtScaling = GetMtScaling(hRStarVsMt, 'rStar_slices')
        gMtScaling.SetName('gMtScaling')
        gMtScaling.SetTitle(';m_{T} (MeV); <r*> (fm)')
        utils.style.SetObjectStyle(gMtScaling)
        gMtScaling.Write()

        gMtScaling3BExpectedFrom2B = utils.analysis.ScaleGraph(gMtScaling, 15 * np.pi / 16, name='gMtScaling3BExpectedFrom2B')
        gMtScaling3BExpectedFrom2B.Write()

    if hFemtoRStarFemtoPairsInTripletsVsMt and hFemtoRStarFemtoPairsInTripletsVsMt.GetEntries():
        gMtScalingInTriplets = GetMtScaling(hFemtoRStarFemtoPairsInTripletsVsMt, 'rStarInTriplets_slices')
        gMtScalingInTriplets.SetName('gMtScalingInTriplets')
        gMtScalingInTriplets.SetTitle(';m_{T} (MeV); <r*> (fm)')
        utils.style.SetObjectStyle(gMtScalingInTriplets)
        gMtScalingInTriplets.Write()


    # Test2B3BConsistency(hRStarInTriplets, hRho)
    # Test2B3BConsistency(hFemtoRStarInTriplets, hFemtoRho)

    # mT scaling for 3B
    # hRhoVsMt = inFile.Get('hRhoVsMt')
    # hRhoVsMt.Write()

    # oFile.mkdir('rho_vs_mt')
    # oFile.cd('rho_vs_mt')
    # hists = Chain(utils.analysis.SliceVertically(hRhoVsMt, cfg['mt_bins'], name='hRho_mT'))
    # hists.apply(Fit, SourceCountsAAA, [0, 8], [['norm', 1, 100000], ['rho0', 0.1, 5]], id='results')

    # gMtScaling, hists = GetMtScaling(hRhoVsMt, SourceCountsAAA)
    # gMtScaling.Write('gRhoVsMt')
    # hists.apply(lambda h : h.Write())
    # oFile.cd('/')

    # hFemtoRhoVsMt = inFile.Get('hFemtoRhoVsMt')
    # hFemtoRhoVsMt.Write()

    # oFile.mkdir('femto_rho_vs_mt')
    # oFile.cd('femto_rho_vs_mt')
    # gMtScaling, hists = GetMtScaling(hFemtoRhoVsMt, SourceCountsAAA)
    # gMtScaling.Write('gFemtoRhoVsMt')
    # hists.apply(lambda h : h.Write())
    # oFile.cd('/')

    # hRho.Write()
    # hRStarInTriplets.Write()
    # hFemtoRho.Write()
    # hFemtoRStarInTriplets.Write()

    oFile.Close()

    print(f'Output saved to {cfg["ofile"]}')
    
if __name__ == '__main__':
    import argparse
    import yaml
    
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg')
    args = parser.parse_args()

    with open(args.cfg) as file:
        cfg = yaml.safe_load(file)

    main(cfg)
