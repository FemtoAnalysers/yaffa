import sys

from ROOT import TFile, gInterpreter, TF1, gROOT, TGraphErrors
gInterpreter.Declare('#include "../../src/RootFunctions.hxx"')
from ROOT import SourceAAA, SourceGauss

sys.path.append('../../src')

from utils import SliceVertically, FitHistList

def Write(objects):
    for object in objects:
        object.Write()

def GetMtScaling(hSources, mTMins, mTMaxs, func, name):
    gRStar0VsMt = TGraphErrors()
    gRStar0VsMt.SetName(name)

    if func == 'average':
        gRStar0VsMt.SetTitle(';m_{T} (GeV);#LT r_{0} #GT (fm)')
        for iBin, hRStar in enumerate(hSources):
            rStarAvg = hRStar.GetMean()
            rStarAvgUnc = hRStar.GetMeanError()
            gRStar0VsMt.SetPoint(iBin, (mTMaxs[iBin] + mTMins[iBin]) / 2, rStarAvg)
            gRStar0VsMt.SetPointError(iBin, (mTMaxs[iBin] - mTMins[iBin]) / 2, rStarAvgUnc)
        return gRStar0VsMt
    
    for iBin, hRStar in enumerate(hSources):
        fSource = TF1(f'fSource_{mTMins[iBin]:.0f}_{mTMins[iBin]:.0f}', func, 0, 20, 1)
        fSource.SetNpx(10000)
        fSource.SetParameter(0, 1)
        fSource.SetParLimits(0, 0.1, 5)

        hRStar.Scale(1./hRStar.GetEntries() / hRStar.GetBinWidth(1))
        hRStar.Fit(fSource, 'MR+L')
        hRStar.SetTitle(';r* (fm);Counts')
        fSource.Write()

        # Fit results
        rho0 = fSource.GetParameter(0)
        rho0unc = fSource.GetParError(0)

        gRStar0VsMt.SetPoint(iBin, (mTMaxs[iBin] + mTMins[iBin]) / 2, rho0)
        gRStar0VsMt.SetPointError(iBin, (mTMaxs[iBin] - mTMins[iBin]) / 2, rho0unc)
    return gRStar0VsMt

def main(cfg:dict):
    '''
    Compute the source size as a function of mT
    '''

    mTMins = cfg['mt_bins'][:-1]
    mTMaxs = cfg['mt_bins'][1:]

    # Open input file
    inFile = TFile(cfg['infile'])

    oFile = TFile(cfg['ofile'], 'recreate')

    # Load hitograms from file
    hRhoVsMt = inFile.Get('hRhoVsMt')
    # if hRhoVsMt.Integral() > 0:
    #     hMt3B = hRhoVsMt.ProjectionX('hMt3B')
    #     hMt3B.Write()

    #     hRhos = SliceVertically(hRhoVsMt, cfg['mt_bins'], name='hRho')

    #     gRho0VsMt = TGraphErrors()
    #     gRho0VsMt.SetName('gRho0VsMt')
    #     for iBin, hRho in enumerate(hRhos):
    #         fSource = TF1(f'fSource_{mTMins[iBin]:.0f}_{mTMins[iBin]:.0f}', SourceAAA, 0, 20, 1)
    #         fSource.SetNpx(10000)
    #         fSource.SetParameter(0, 1)
            
    #         hRho.Scale(1./hRho.GetEntries() / hRho.GetBinWidth(1))
    #         hRho.Fit(fSource, 'MR+L')
    #         hRho.SetTitle(';#rho* (fm);Counts')
    #         hRho.Write()
    #         fSource.Write()

    #         # Fit results
    #         rho0 = fSource.GetParameter(0)
    #         rho0unc = fSource.GetParError(0)

    #         gRho0VsMt.SetPoint(iBin, (mTMaxs[iBin] + mTMins[iBin]) / 2, rho0)
    #         gRho0VsMt.SetPointError(iBin, (mTMaxs[iBin] - mTMins[iBin]) / 2, rho0unc)
    #     gRho0VsMt.SetTitle(';m_{T} (GeV);#rho_{0} (fm)')
    #     gRho0VsMt.Write()
    
    # mT scaling for 2B
    hRStarVsMt = inFile.Get('hRStarVsMt')
    hMt = hRStarVsMt.ProjectionX('hMt')
    hMt.Write()
    
    hRStars = SliceVertically(hRStarVsMt, cfg['mt_bins'], name='hRStar')
    gRStar0VsMt = GetMtScaling(hRStars, mTMins, mTMaxs, SourceGauss, 'gRStarVsMt')
    Write(hRStars)
    gRStar0VsMt.Write()

    # Method = average
    gAvgRStar0VsMt = GetMtScaling(hRStars, mTMins, mTMaxs, 'average', 'gAvgRStarVsMt')
    gAvgRStar0VsMt.Write()

    oFile.Close()

if __name__ == '__main__':
    import argparse
    import yaml
    
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg')
    args = parser.parse_args()

    with open(args.cfg) as file:
        cfg = yaml.safe_load(file)

    main(cfg)
