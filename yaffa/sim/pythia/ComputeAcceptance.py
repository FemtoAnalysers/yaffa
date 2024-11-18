'''
Script to compute the acceptance of a decay chain

Usage:
python3 ComputeAcceptance.py Distr.root
'''

import argparse
from rich import print  # pylint: disable=redefined-builtin

from ROOT import TFile, TF1, TCanvas  # pylint: disable=no-name-in-module

from yaffa.utils.io import GetSubdirNames
from yaffa.utils.analysis import WeightedAverage
from yaffa.utils.style import SetStyle

def ComputeAcceptance(args):
    '''
    Compute the acceptance and efficiency
    '''

    inFile = TFile(args.file)
    nGen = inFile.Get('hEvt').GetEntries()
    print(f'Number of events: {nGen}')

    # Efficiency
    hPt0 = inFile.Get("qa0/hPt0")
    hPt1 = inFile.Get("qa1/hPt1")

    fEff = TF1('fEff', args.eff, 0, 5)
    eff0 = WeightedAverage(fEff, hPt0)
    eff1 = WeightedAverage(fEff, hPt1)

    SetStyle()
    cPt0 = TCanvas('cPt0', '', 1200, 600)
    cPt0.Divide(2, 1)
    pad1 = cPt0.cd(1)
    hPt0.SetLineColor(4)
    hPt0.GetXaxis().SetRangeUser(0, 4)
    hPt0.Draw()

    hPt0TimesEff = hPt0.Clone('hPt0TimesEff')
    hPt0TimesEff.Multiply(fEff)
    hPt0TimesEff.SetLineColor(2)
    hPt0TimesEff.Draw('same')
    pad1.BuildLegend()

    pad2 = cPt0.cd(2)
    hPt1.SetLineColor(4)
    hPt1.GetXaxis().SetRangeUser(0, 4)
    hPt1.Draw()

    hPt1TimesEff = hPt1.Clone('hPt1TimesEff')
    hPt1TimesEff.Multiply(fEff)
    hPt1TimesEff.SetLineColor(2)
    hPt1TimesEff.Draw('same')
    pad2.BuildLegend()

    cPt0.SaveAs('cPt0.pdf')


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
    parser.add_argument('eff', default='1')
    args = parser.parse_args()

    ComputeAcceptance(args)

if __name__ == '__main__':
    main()
