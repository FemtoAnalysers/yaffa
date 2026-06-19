import os
import numpy as np
from pathlib import Path

from ROOT import TGraph

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
from ROOT import _SourceAAA, _SourceGauss

utils.style.SetStyle()
from yaffa.utils.analysis import Convert

def ComputeSource(source, radii):
    first, second = source.split(':')

    if first == 'gaussAAA':
        source = [_SourceAAA(radius, float(second)) for radius in radii]
    elif first == 'gauss2b':
        source = [_SourceGauss(radius, float(second)) for radius in radii]    
    elif source:
        inFile = TFile(first)
        hSource = inFile.Get(second)
        hSource.SetDirectory(0)
        inFile.Close()

        source = [hSource.GetBinContent(hSource.FindBin(radius)) for radius in radii]    

    return source

def main(ofile, wf, source=None, radius=2.6):
    if not Path(wf).exists():
        log.error(f'File "{wf}" does not exist.')
        return

    gCF = TGraph(1)
    gCF.SetName('gCF')
    gCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')

    if '.root' in wf:
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

    source = np.array(ComputeSource(source, radii), dtype='d')

    for iMomentum, (momentum, wf2) in enumerate(zip(momenta, wf)):
        gCF.SetPoint(iMomentum, momentum, source @ wf2 / sum(source))

    oFile = TFile(ofile, 'recreate')
    gCF.Write()
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
