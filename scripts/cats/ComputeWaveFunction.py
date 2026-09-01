import os
import numpy as np
from ROOT import gSystem, gInterpreter, TH2D, TDatabasePDG, TGraph, TFile

from dotenv import load_dotenv
from pathlib import Path

env_path = Path(__file__).resolve().parents[2]/ ".env"
print(f'Loading env from {env_path}')
if not load_dotenv(dotenv_path=env_path, verbose=True, override=True):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
if not YAFFA_PATH:
    print("\033[33mWARNING: Path to yaffa is empty, something might break!\033[0m")

from ROOT import gInterpreter, TFile
# Include the .cpp (not just the header) so cling JIT-compiles the member
# definitions in WaveFunction.cpp -- there is no compiled libyaffa to link.
gInterpreter.Declare(f'#include "{YAFFA_PATH}/src/cpp/WaveFunction.cpp"')
from ROOT import WaveFunction

CATS_PATH = os.environ['CATS']

gSystem.Load(f'{CATS_PATH}/install/lib/libCATS.so')
gInterpreter.AddIncludePath(f'{CATS_PATH}/install/include')

gInterpreter.Declare('''
#include "CATS.h"
#include "CommonAnaFunctions.h"
''')

from ROOT import CATS
from ROOT import DLM_CommonAnaFunctions

RADIUS_STEP = 0.01
RADIUS_MAX = 20
KSTAR_STEP = 1
KSTAR_MAX = 500

def matrix_to_th2d(matrix, title=''):
    h = TH2D(
        'hWF', title,
        round(RADIUS_MAX / RADIUS_STEP), 0, RADIUS_MAX,
        round(KSTAR_MAX / KSTAR_STEP), 0, KSTAR_MAX
    )

    ny, nx = matrix.shape

    for iy in range(ny):
        for ix in range(nx):
            h.SetBinContent(ix + 1, iy + 1, float(matrix[iy, ix]))

    return h 

def get_wave_function(cats, channel):
    nk = cats.GetNumMomBins()
    nr = round(RADIUS_MAX / RADIUS_STEP)

    wf = np.empty((nk, nr))

    for ik in range(nk):
        for ir in range(nr):
            radius = RADIUS_STEP/2 + ir * RADIUS_STEP
            wf[ik, ir] = cats.EvalWaveFun2(ik, radius, channel)

    return wf

def reduced_mass(m1, m2):
    return m1 * m2 / (m1 + m2)

def compute_wave_function(system, oFile):
    n_kstar_bins = round(KSTAR_MAX / KSTAR_STEP)
    pdg = TDatabasePDG.Instance()

    cats = CATS()

    if system == 'pp':
        m1 = 1000. * pdg.GetParticle(2212).Mass()
        m2 = 1000. * pdg.GetParticle(2212).Mass()

        header = 'Wave function of proton-proton with Argonne v18 potential computed with CATS\n'
        cats.SetMomBins(n_kstar_bins, 0, KSTAR_MAX)
        cats.SetQ1Q2(1)
        cats.SetQuantumStatistics(True)
        cats.SetRedMass(reduced_mass(m1, m2))

        can = DLM_CommonAnaFunctions()
        can.SetUpCats_pp(cats, 'AV18', 'Gauss', 0, 0)
    else:
        raise RuntimeError('System not implemented')

    cats.SetAnaSource(0, 1)
    cats.KillTheCat()

    gCF = TGraph()

    for i in range(n_kstar_bins):
        gCF.SetPoint(i, cats.GetMomentum(i), cats.GetCorrFun(i))

    nChn = cats.GetNumChannels()
    if nChn == 2:
        weights = [0.25, 0.75]
    elif nChn == 4:
        weights = [3/12, 1/12, 3/12, 5/12]
    else:
        raise RuntimeError(f'Unexpected number of channels: {nChn}')

    wf = sum([get_wave_function(cats, iChn) * weights[iChn] for iChn in range(nChn)])

    hWF = matrix_to_th2d(wf)
    hWF.SetTitle('pp, AV18, |#psi|^{2};r (fm);k* (MeV/c);|#psi|^{2}')

    fout = TFile('wf.root', 'RECREATE')
    gCF.Write('gCF')
    hWF.Write('hWF')
    fout.Close()

    radius = RADIUS_STEP/2 + np.arange(wf.shape[1]) * RADIUS_STEP
    kstar = KSTAR_STEP/2 + np.arange(wf.shape[0]) * KSTAR_STEP

    # --- new: save via the WaveFunction C++ class (.wf standardized format) ---
    # wf has shape (nk, nr) == (row = momentum, column = r*), which is exactly
    # the class's row-major convention: values[iMom * nRad + iRad].
    wf_path = str(Path(oFile).with_suffix('.wf'))
    wff = WaveFunction(
        kstar.tolist(),        # momentum axis: k* (MeV/c), one entry per row
        radius.tolist(),       # r* axis (fm), one entry per column
        wf.ravel().tolist(),   # row-major values, iMom slow / iRad fast
        2,                     # nBody
        system,                # free-form system label, e.g. "pp"
    )
    wff.Save(wf_path)
    print(f'WaveFunction (.wf) written to {wf_path}')

    # --- old: legacy plain-text table kept for the time being (.dat format) ---
    dat_path = str(Path(oFile).with_suffix('.dat'))
    wf = np.column_stack([radius, wf.T])
    header += (
        f'radius_step = {RADIUS_STEP} fm\n'
        f'kstar_step = {KSTAR_STEP} MeV/c\n'
        'radius ' + ' '.join(f'{k:.3f}' for k in kstar)
    )

    np.savetxt(
        dat_path,
        wf,
        fmt='%.17e',
        header=header
    )
    print(f'legacy wave function (.dat) written to {dat_path}')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('system', choices=('pp'))
    parser.add_argument('oFile')
    args = parser.parse_args()
    
    compute_wave_function(args.system, args.oFile)
