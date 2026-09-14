'''
Script to smear a TGraph or a wave function with a matrix, e.g. a momentum resolution matrix or a phase-space decay
matrix.

The matrix must have the input variable on the x axis (e.g. k*_gen, or k* of the parent pair) and the output variable on
the y axis (e.g. k*_reco, or k* of the daughter pair). All quantities are assumed to be in MeV.

Wave functions (.wf) are smeared in |psi|^2 at fixed radius: the smeared |psi|^2 at each output momentum is the average
of |psi|^2 over the input momenta, weighted with the corresponding row of the matrix. Folding it with a k*-independent
source gives the same correlation function as smearing the correlation function computed with the original wave function.

Usage:
    python3 Smear.py input.root:gCF output.root:gCF_smeared --matrix matrix.root:hPhaseSpace
    python3 Smear.py input.wf output.wf --matrix matrix.root:hResolutionMatrixME
'''

import argparse
from pathlib import Path

import numpy as np

from ROOT import TFile, TGraph, TH2, gInterpreter  # pylint: disable=import-error,no-name-in-module
from yaffa.utils.io import Load
from yaffa.utils.analysis import SmearGraph
from yaffa import logger as log


def Split(path):
    '''
    Split a "file.root:path/to/object" string into the file name and the path of the object inside the file.
    '''
    if ':' not in path:
        log.critical('Invalid path "%s", expected <file.root>:<path/to/object>', path)
    return path.rsplit(':', 1)


def LoadObject(path):
    '''
    Load an object from a "file.root:path/to/object" string.
    '''
    fileName, objPath = Split(path)
    inFile = TFile(fileName)
    if inFile.IsZombie():
        log.critical('Cannot open %s', fileName)

    obj = Load(inFile, objPath)
    if obj == None:  # pylint: disable=singleton-comparison
        log.critical('Object %s not found in %s', objPath, fileName)
    if hasattr(obj, 'SetDirectory'):
        obj.SetDirectory(0)
    inFile.Close()
    return obj


def SmearWaveFunction(wf, matrix, description=''):
    '''
    Smear |psi|^2 along the momentum axis with a matrix (x: input momentum, y: output momentum).

    The output momentum axis is given by the centers of the y bins of the matrix that are populated and lie within the
    momentum range of the wave function. |psi|^2 is linearly interpolated in momentum and held at the edge value
    outside the tabulated range.

    Parameters
    ----------
    wf : WaveFunction
        The wave function to be smeared
    matrix : TH2
        The smearing matrix
    description : str, optional
        Free-form notes written in the header of the smeared wave function

    Returns
    -------
    WaveFunction
        The smeared wave function
    '''
    from ROOT import WaveFunction  # pylint: disable=import-error,no-name-in-module,import-outside-toplevel

    mom = np.array(wf.Momentum())
    radius = np.array(wf.Radius())
    values = np.array(wf.Values()).reshape(len(mom), len(radius))

    nX = matrix.GetNbinsX()
    nY = matrix.GetNbinsY()
    xAxis = matrix.GetXaxis()
    yAxis = matrix.GetYaxis()
    kIn = np.array([xAxis.GetBinCenter(i + 1) for i in range(nX)])
    kOut = np.array([yAxis.GetBinCenter(i + 1) for i in range(nY)])

    # Matrix content without under/overflow: counts[iY, iX]
    counts = np.array([[matrix.GetBinContent(iX + 1, iY + 1) for iX in range(nX)] for iY in range(nY)])
    norm = counts.sum(axis=1)
    keep = (norm >= 1) & (kOut <= mom[-1])
    weights = counts[keep] / norm[keep, None]

    # |psi|^2 at the input momenta
    iHi = np.clip(np.searchsorted(mom, kIn), 1, len(mom) - 1)
    iLo = iHi - 1
    t = np.clip((kIn - mom[iLo]) / (mom[iHi] - mom[iLo]), 0, 1)
    valuesIn = (1 - t)[:, None] * values[iLo] + t[:, None] * values[iHi]

    smeared = weights @ valuesIn
    return WaveFunction(kOut[keep].tolist(), radius.tolist(), smeared.ravel().tolist(), wf.GetNBody(), wf.System(),
                        description)


def main():
    '''
    Main function.
    '''
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='graph to be smeared, as <file.root>:<path/to/graph>, or wave function as <file.wf>')
    parser.add_argument('output', help='where to save the smeared graph, as <file.root>:<name>, or <file.wf>')
    parser.add_argument('--matrix', required=True,
                        help='smearing matrix (x: input, y: output), as <file.root>:<path/to/TH2>')
    args = parser.parse_args()

    matrix = LoadObject(args.matrix)
    if not isinstance(matrix, TH2):
        log.critical('The smearing matrix must be a TH2, got %s', type(matrix))

    if args.input.endswith('.wf'):
        if not args.output.endswith('.wf'):
            log.critical('A smeared wave function must be saved in the .wf format, got "%s"', args.output)

        gInterpreter.Declare(f'#include "{Path(__file__).resolve().parents[1] / "src/cpp/WaveFunction.cpp"}"')
        from ROOT import WaveFunction  # pylint: disable=import-error,no-name-in-module,import-outside-toplevel

        wf = WaveFunction(args.input)
        description = f'smeared from {args.input}\nwith {args.matrix} (x: input momentum, y: output momentum)'
        SmearWaveFunction(wf, matrix, description).Save(args.output)
        print(f'Output saved in {args.output}')
        return

    graph = LoadObject(args.input)
    if not isinstance(graph, TGraph):
        log.critical('Smearing for type %s is not implemented. Only TGraph and .wf are supported.', type(graph))

    oFileName, oName = Split(args.output)
    title = f';{graph.GetXaxis().GetTitle()};{graph.GetYaxis().GetTitle()}'
    gSmeared = SmearGraph(graph, matrix, name=oName, title=title)

    oFile = TFile(oFileName, 'recreate')
    gSmeared.Write()
    oFile.Close()
    print(f'Output saved in {oFileName}')


if __name__ == '__main__':
    main()
