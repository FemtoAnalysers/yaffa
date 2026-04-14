'''
Script to compute the acceptance of a decay chain

Usage:
python3 ComputeAcceptance.py Distr.root
'''

import argparse
import numpy as np
from rich import print  # pylint: disable=redefined-builtin

from ROOT import TFile, TF1, TCanvas, TLegend  # pylint: disable=no-name-in-module

from yaffa.utils.io import GetSubdirNames
from yaffa.utils.analysis import WeightedAverage
from yaffa.utils.style import SetStyle

def ComputeAcceptance(args): # pylint: disable=too-many-statements
    '''
    Compute the acceptance and efficiency
    '''

    inFile = TFile(args.file)
    nGen = inFile.Get('hEvt').GetEntries()
    print(f'Number of events: {nGen}')

    hPt0 = inFile.Get("qa0/hPt0")
    hPt1 = inFile.Get("qa1/hPt1")

    # Load efficiency 0
    if '.root' in args.eff0:
        inFileEff0 = TFile(args.eff0)
        aEff0 = inFileEff0.Get('hEff')
        aEff0.SetDirectory(0)
        edges = np.array([0, ] * aEff0.GetNbinsX(), dtype='d')
        aEff0.GetXaxis().GetLowEdge(edges)

        # Ensure correct binning
        # edges = np.concatenate((edges, [hPt0.GetXaxis().GetXmax()]))
        if edges[0] > 0.001: # Assuing the Eff is given vs pT (GeV/c)
            edges = np.concatenate(([0], edges))

        # Rebin
        hPt0 = hPt0.Rebin(len(edges) -1, "hPt0_rebin", edges)
        print(edges)
        hPt0.GetXaxis().GetLowEdge(edges)
        print(edges)

    else:
        aEff0 = TF1('fEff0', args.eff0, 0, 5)

    # Load efficiency 1
    if '.root' in args.eff1:
        inFileEff1 = TFile(args.eff1)
        aEff1 = inFileEff1.Get('hEff')
        aEff1.SetDirectory(0)
        edges = np.array([0, ] * aEff1.GetNbinsX(), dtype='d')
        aEff1.GetXaxis().GetLowEdge(edges)

        # Ensure correct binning
        edges = np.concatenate((edges, [hPt0.GetXaxis().GetXmax()]))
        if edges[0] > 0.001: # Assuing the Eff is given vs pT (GeV/c)
            edges = np.concatenate(([0], edges))

        # Rebin
        hPt0 = hPt0.Rebin(len(edges) -1, "hPt0_rebin", edges)
    else:
        aEff1 = TF1('fEff1', args.eff1, 0, 5)

    eff0 = WeightedAverage(aEff0, hPt0)
    eff1 = WeightedAverage(aEff1, hPt1)

    SetStyle()

    oFile = TFile('EffAcc.root', 'recreate')

    # Plot Efficiency
    cEff = TCanvas('cEff', '', 1200, 600)
    cEff.Divide(2, 1)
    cEff.cd(1)
    aEff0.SetStats(0)
    aEff0.Draw()
    aEff0.Write()

    cEff.cd(2)
    aEff0.SetStats(0)
    aEff1.Draw()
    aEff1.Write()

    cEff.SaveAs('cEff.pdf')

    # Plot pt distributions
    cPt0 = TCanvas('cPt0', '', 1200, 600)
    cPt0.Divide(2, 1)
    cPt0.cd(1)
    hPt0.GetXaxis().SetRangeUser(0, 4)
    hPt0.SetMarkerColor(4)
    hPt0.SetLineColor(4)
    hPt0.Draw()
    hPt0.Write()

    hPt0TimesEff = hPt0.Clone('hPt0TimesEff')
    hPt0TimesEff.Multiply(aEff0 if aEff0 else eff0)
    hPt0TimesEff.SetMarkerColor(2)
    hPt0TimesEff.SetLineColor(2)
    hPt0TimesEff.Draw('same')
    hPt0TimesEff.Write()

    leg0 = TLegend(0.5, 0.7, 0.9, 0.9)
    leg0.AddEntry(hPt0, 'Before')
    leg0.AddEntry(hPt0TimesEff, 'After')
    leg0.Draw()


    cPt0.cd(2)
    hPt1.GetXaxis().SetRangeUser(0, 4)
    hPt1.SetMarkerColor(4)
    hPt1.SetLineColor(4)
    hPt1.Draw()
    hPt1.Write()

    hPt1TimesEff = hPt1.Clone('hPt1TimesEff')
    hPt1TimesEff.Multiply(aEff1 if aEff1 else eff1)
    hPt1TimesEff.SetMarkerColor(2)
    hPt1TimesEff.SetLineColor(2)
    hPt1TimesEff.Draw('same')
    hPt1TimesEff.Write()

    leg1 = TLegend(0.5, 0.7, 0.9, 0.9)
    leg1.AddEntry(hPt1, 'Before')
    leg1.AddEntry(hPt1TimesEff, 'After')
    leg1.Draw()

    cPt0.SaveAs('cPt.pdf')

    oFile.Close()

    print(f'<eff> part0 = {eff0 * 100:.4f} %')
    print(f'<eff> part1 = {eff1 * 100:.4f} %')
    print(f'product = {eff0 * eff1 * 100:.4f} %')

    # Acceptance
    for pair in [name for name in GetSubdirNames(inFile) if name[0] == 'p']:
        nReco = inFile.Get(f'{pair}/hSE').GetEntries()
        print(f'pair: {pair}  Acc: {nReco/nGen:.5f}')


def main():
    '''
    Main function
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('file', default='Distr.root')
    parser.add_argument('--eff0', default='1')
    parser.add_argument('--eff1', default='1')
    args = parser.parse_args()

    ComputeAcceptance(args)

if __name__ == '__main__':
    main()
