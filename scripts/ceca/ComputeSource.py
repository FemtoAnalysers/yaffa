import sys

from ROOT import TFile, gInterpreter, TF1, gROOT, TGraphErrors
gInterpreter.Declare('#include "../../src/RootFunctions.hxx"')
from ROOT import SourceAAA

sys.path.append('../../src')

from utils import SliceVertically

def main(cfg:dict):
    mTMins = cfg['mt_bins'][:-1]
    mTMaxs = cfg['mt_bins'][1:]

    # Open input file
    inFile = TFile(cfg['infile'])

    oFile = TFile(cfg['ofile'], 'recreate')

    # Load hitograms from file
    hRhoVsMt = inFile.Get('hRhoVsMt')

    fSource = TF1('fSource', SourceAAA, 0, 20, 1)
    fSource.SetNpx(10000)
    fSource.SetParameter(0, 1)
    print("int: ", integ)

    hMt = hRhoVsMt.ProjectionX('hMt')
    hRhos = SliceVertically(hRhoVsMt, cfg['mt_bins'], name='hRho')

    gRho0VsMt = TGraphErrors()
    gRho0VsMt.SetName('gRho0VsMt')
    for hRho in hRhos:
        hRho.SetTitle(';#rho* (fm);Counts')
        hRho.Scale(1./hRho.GetEntries() / hRho.GetBinWidth(1))
        hRho.Fit(fSource, 'MR+')
        hRho.Write()
        rho0 = fSource.GetParameter(0)
        rho0unc = fSource.GetParError(0)

        integ = fSource.Integral(0, 20)
    
        iBin = gRho0VsMt.GetN()
        gRho0VsMt.SetPoint(iBin, (mTMaxs[iBin] + mTMins[iBin]) / 2, rho0)
        gRho0VsMt.SetPointError(iBin, (mTMaxs[iBin] - mTMins[iBin]) / 2, rho0unc)
    fSource.Write()
    hMt.Write()
    gRho0VsMt.Write()

if __name__ == '__main__':
    import argparse
    import yaml
    
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg')
    args = parser.parse_args()

    with open(args.cfg) as file:
        cfg = yaml.safe_load(file)

    main(cfg)
