#!/bin/env python3

'''
Script to plot a wave function (.wf) as a function of the radius for some values of the momentum, and to compute the
correlation function with a Gaussian source.

Usage:
    python3 PlotWaveFunction.py input.wf [-n 5] [-k [0] 200] [--r0 1] [-o WaveFunction]
'''

import argparse
from pathlib import Path

import numpy as np

from ROOT import TCanvas, TFile, TGraph, TLegend, gInterpreter  # pylint: disable=import-error,no-name-in-module
from yaffa.utils import style


def main():
    '''
    Main function.
    '''
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('wf', help='wave function in the .wf format')
    parser.add_argument('-n', type=int, default=5, help='number of momentum values to plot')
    parser.add_argument('-k', type=float, nargs='+', default=None,
                        help='momentum range to plot (MeV/c): [min] max, min defaults to 0, default: full range')
    parser.add_argument('--r0', type=float, default=1., help='radius of the Gaussian source (fm), hyper-radius for 3B')
    parser.add_argument('-o', '--ofilename', default='WaveFunction', help='name of the output files, without extension')
    args = parser.parse_args()

    cppDir = Path(__file__).resolve().parents[1] / 'src/cpp'
    gInterpreter.Declare(f'#include "{cppDir / "WaveFunction.cpp"}"')
    gInterpreter.Declare(f'#include "{cppDir / "Functions.hxx"}"')
    from ROOT import WaveFunction, _SourcePdfGauss, _SourcePdfAAAHypRad  # pylint: disable=import-error,no-name-in-module,import-outside-toplevel

    wf = WaveFunction(args.wf)
    mom = np.array(wf.Momentum())
    radius = np.array(wf.Radius())
    values = np.array(wf.Values()).reshape(len(mom), len(radius))

    if wf.GetNBody() == 2:
        momLabel = '#it{k}*'
        radLabel = '#it{r}*'
        source = np.array([_SourcePdfGauss(r, args.r0) for r in radius])
    else:
        momLabel = '#it{Q}_{3}'
        radLabel = '#rho'
        source = np.array([_SourcePdfAAAHypRad(r, args.r0) for r in radius])

    style.SetStyle()

    # Wave function vs radius for n momentum values
    kMin, kMax = (0, mom[-1]) if args.k is None else ([0] + args.k)[-2:]
    iMoms = np.abs(mom[:, None] - np.linspace(kMin, kMax, args.n)).argmin(axis=0)
    yMax = 1.1 * values[iMoms][:, radius > 1].max()  # skip r < 1 fm, where the tables can be numerically unstable
    cWF = TCanvas('cWF', '', 600, 600)
    cWF.DrawFrame(radius[0], 0, radius[-1], yMax, f';{radLabel} (fm);|#Psi|^{{2}}')
    leg = TLegend(0.55, 0.9 - 0.05 * (args.n + 1), 0.9, 0.9)
    leg.SetHeader(f'{wf.System()} {wf.Potential()}')
    graphs = []
    for iMom in iMoms:
        gWF = TGraph(len(radius), radius, values[iMom].copy())
        gWF.SetLineWidth(2)
        gWF.Draw('same l plc')
        leg.AddEntry(gWF, f'{momLabel} = {mom[iMom]:.1f} MeV/#it{{c}}', 'l')
        graphs.append(gWF)
    leg.Draw()
    cWF.SaveAs(f'{args.ofilename}.pdf')

    # Correlation function with a Gaussian source
    cf = values @ source / source.sum()
    gCF = TGraph(len(mom), mom, cf)
    gCF.SetName('gCF')
    gCF.SetTitle(f';{momLabel} (MeV/#it{{c}});#it{{C}}({momLabel})')

    oFile = TFile(f'{args.ofilename}.root', 'recreate')
    gCF.Write()
    oFile.Close()
    print(f'Output saved in {args.ofilename}.pdf and {args.ofilename}.root')


if __name__ == '__main__':
    main()
