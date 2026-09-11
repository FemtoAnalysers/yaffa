'''
Script to smear a TGraph with a matrix, e.g. a momentum resolution matrix or a phase-space decay matrix.

The matrix must have the input variable on the x axis (e.g. k*_gen, or k* of the parent pair) and the output variable on
the y axis (e.g. k*_reco, or k* of the daughter pair). All quantities are assumed to be in MeV.

Usage:
    python3 Smear.py input.root:gCF output.root:gCF_smeared --matrix matrix.root:hPhaseSpace
'''

import argparse

from ROOT import TFile, TGraph, TH2  # pylint: disable=import-error,no-name-in-module
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


def main():
    '''
    Main function.
    '''
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='graph to be smeared, as <file.root>:<path/to/graph>')
    parser.add_argument('output', help='where to save the smeared graph, as <file.root>:<name>')
    parser.add_argument('--matrix', required=True,
                        help='smearing matrix (x: input, y: output), as <file.root>:<path/to/TH2>')
    args = parser.parse_args()

    graph = LoadObject(args.input)
    if not isinstance(graph, TGraph):
        log.critical('Smearing for type %s is not implemented. Only TGraph is supported.', type(graph))

    matrix = LoadObject(args.matrix)
    if not isinstance(matrix, TH2):
        log.critical('The smearing matrix must be a TH2, got %s', type(matrix))

    oFileName, oName = Split(args.output)
    title = f';{graph.GetXaxis().GetTitle()};{graph.GetYaxis().GetTitle()}'
    gSmeared = SmearGraph(graph, matrix, name=oName, title=title)

    oFile = TFile(oFileName, 'recreate')
    gSmeared.Write()
    oFile.Close()
    print(f'Output saved in {oFileName}')


if __name__ == '__main__':
    main()
