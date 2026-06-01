#!/bin/env python

import argparse

from ROOT import TF2, TCanvas, TFile, TF1, gInterpreter, gROOT
gInterpreter.Declare('#include "../src/cpp/RootFunctions.hxx"')
from ROOT import SourceAAAJC, SourceGauss


from yaffa import utils

utils.style.SetStyle()

gROOT.SetBatch(True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('system', choices=('aa', 'aaa'))
    parser.add_argument('ofile', default='Source.root', type=str)
    parser.add_argument('--r0', default=0, type=float)
    parser.add_argument('--jc', action='store_true', default=False, help='Use source function in jacobi coordinates (3B)')
    args = parser.parse_args()

    cSource = TCanvas('cSource', '', 600, 600)
    cSource.SetTopMargin(0.02)
    cSource.SetRightMargin(0.2)

    # Plot source function
    if args.system == 'aa':
        fSource = TF1('fSource', SourceGauss, 0, 8, 1)
        fSource.SetParameter(0, args.r0)
        fSource.SetNpx(1000)
        fSource.SetTitle(';r* (fm); 4#pi r*^{2} S(r*)')
        fSource.Draw('colz')

    elif args.system == 'aaa' and args.jc:
        cSource.SetRightMargin(0.22)
        fSource = TF2('fSource', SourceAAAJC, 0, 8, 0, 8, 1)
        fSource.SetParameter(0, args.r0)
        fSource.SetNpx(100)
        fSource.SetNpy(100)
        fSource.SetTitle(';r_{12} (fm);r_{3,12} (fm);(4#pi)^{2} r_{12}^{2} r_{3,12}^{2} S(r_{12}, r_{3,12})')
        fSource.GetZaxis().SetTitleOffset(1.7)
        fSource.Draw('colz')
        
    else:
        raise NotImplementedError('Source for this system is not implemented')

    # Save to file
    if '.root' in args.ofile:
        oFile = TFile(args.ofile, 'recreate')
        fSource.Write()
        oFile.Close()
    else:
        cSource.SaveAs(args.ofile)
