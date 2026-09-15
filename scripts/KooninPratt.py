import os
import numpy as np
from pathlib import Path

from ROOT import TGraph, TF1

from yaffa import utils
from yaffa import logger as log

from dotenv import load_dotenv
from pathlib import Path

env_path = Path(__file__).resolve().parent.parent / ".env"
print(f'Loading env from {env_path}')
if not load_dotenv(dotenv_path=env_path, verbose=True, override=True):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
if not YAFFA_PATH:
    print("\033[33mWARNING: Path to yaffa is empty, something might break!\033[0m")

from ROOT import gInterpreter, TFile
gInterpreter.Declare(f'#include "{YAFFA_PATH}/src/cpp/RootFunctions.hxx"')
from ROOT import _SourcePdfAAAHypRad, _SourcePdfGauss

utils.style.SetStyle()
from yaffa.utils.analysis import Convert

def ComputeSource(source, radii):
    first, second = source.split(':')

    if first == 'gaussAAA':
        source = [_SourcePdfAAAHypRad(radius, float(second)) for radius in radii]
    elif first == 'gauss2b':
        source = [_SourcePdfGauss(radius, float(second)) for radius in radii]    
    elif '.root' in first:
        inFile = TFile(first)
        hSource = inFile.Get(second)

        if isinstance(hSource, TF1):
            source = [hSource.Eval(radius) / hSource.GetParameter(0) for radius in radii]
        else:
            hSource.SetDirectory(0)
            source = [hSource.GetBinContent(hSource.FindBin(radius)) for radius in radii]

        inFile.Close()
    else:
        raise ValueError("Invalid source")

    return source

def main(ofile, wf, source=None):
    if not Path(wf).exists():
        log.error(f'File "{wf}" does not exist.')
        return

    if wf.endswith('.wf'):
        # The sources above already include the Jacobian, so only the tabulated |psi|^2 is taken from the file
        gInterpreter.Declare(f'#include "{YAFFA_PATH}/src/cpp/WaveFunction.cpp"')
        from ROOT import WaveFunction

        wfTable = WaveFunction(wf)
        radii = list(wfTable.Radius())
        momenta = np.array(wfTable.Momentum())
        wf = np.array(wfTable.Values()).reshape(len(momenta), len(radii))
    elif '.root' in wf:
        inFile = TFile(wf)
        hWF = inFile.Get('hWF')
        hWF.SetDirectory(0)
        inFile.Close()
        radii = [hWF.GetXaxis().GetBinCenter(iBin + 1) for iBin in range(hWF.GetNbinsX())]

        wf = Convert(hWF, 'numpy_array').T
        momenta = np.array([hWF.GetYaxis().GetBinCenter(iBin + 1) for iBin in range(hWF.GetNbinsY())])
    elif '.dat' in wf:
        with open(wf) as f:
            for line in f:
                if line.startswith("# radius "):
                    tokens = line.split()
                    momenta = np.array([float(x) for x in tokens[2:]])
                    break
        data = np.loadtxt(wf)

        radii = data[:, 0]
        wf = data[:, 1:].T

    oFile = TFile(ofile, 'recreate')

    for label, src in utils.io.Expand(source):
        if not label:
            suffix = ''
        elif label[0] == '_':
            suffix = label
        else:
            suffix = f'_{label}'

        sourceValues = np.array(ComputeSource(src, radii), dtype='d')

        gCF = TGraph(1)
        gCF.SetName(f'gCF{suffix}')
        gCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')

        gSource = TGraph(len(sourceValues))
        gSource.SetName(f'gSource{suffix}')
        for iPoint, (r, s) in enumerate(zip(radii, sourceValues)):
            gSource.SetPoint(iPoint, r, s)

        for iMomentum, (momentum, wf2) in enumerate(zip(momenta, wf)):
            gCF.SetPoint(iMomentum, momentum, sourceValues @ wf2 / sum(sourceValues))

        oFile.cd()
        gCF.Write()
        gSource.Write()

    oFile.Close()

    print(f'Output saved in {ofile}')

if __name__ == '__main__':    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('ofile', nargs='?', default='CF.pdf')
    parser.add_argument('--wf', help='Wave function')
    parser.add_argument('--source')
    args = parser.parse_args()

    main(args.ofile, args.wf, args.source)
