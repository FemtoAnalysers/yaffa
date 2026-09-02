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

# Folder with the CATS input files (external wave functions, source
# distributions, ...). Only needed by the systems whose potential is read from
# disk rather than computed analytically.
CATS_FILES_PATH = os.getenv("CATS_FILES")

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
        title = 'pp, AV18, |#psi|^{2};r (fm);k* (MeV/c);|#psi|^{2}'
        cats.SetMomBins(n_kstar_bins, 0, KSTAR_MAX)
        cats.SetQ1Q2(1)
        cats.SetQuantumStatistics(True)
        cats.SetRedMass(reduced_mass(m1, m2))

        can = DLM_CommonAnaFunctions()
        can.SetUpCats_pp(cats, 'AV18', 'Gauss', 0, 0)
    elif system == 'pL':
        m1 = 1000. * pdg.GetParticle(2212).Mass()
        m2 = 1000. * pdg.GetParticle(3122).Mass()

        header = ('Wave function of proton-lambda with the chiral EFT NLO13(600) potential '
                  'including the coupled S, P and D channels, computed with CATS\n')
        title = 'p#Lambda, #chiEFT NLO13(600), |#psi|^{2};r (fm);k* (MeV/c);|#psi|^{2}'
        cats.SetMomBins(n_kstar_bins, 0, KSTAR_MAX)
        cats.SetQ1Q2(0)  # the lambda is neutral: no Coulomb
        cats.SetQuantumStatistics(False)
        cats.SetRedMass(reduced_mass(m1, m2))

        if not CATS_FILES_PATH:
            raise RuntimeError(
                'CATS_FILES is not set. Add it to the .env file: it must point to the folder '
                'with the CATS input files (the one containing Interaction/Haidenbauer/).'
            )

        can = DLM_CommonAnaFunctions()
        can.SetCatsFilesFolder(CATS_FILES_PATH)
        can.SetUpCats_pL(cats, 'Chiral_Coupled_SPD', 'Gauss', 0, 0)
    else:
        raise RuntimeError('System not implemented')

    cats.SetAnaSource(0, 1)
    cats.KillTheCat()

    gCF = TGraph()

    for i in range(n_kstar_bins):
        gCF.SetPoint(i, cats.GetMomentum(i), cats.GetCorrFun(i))

    # CATS builds the correlation function as C(k*) = Int S(r) sum_ch w_ch |psi_ch|^2,
    # so the same weights turn the per-channel wave functions into the tabulated one.
    # Channels with a vanishing weight (e.g. the p-waves of the coupled-channel
    # potentials, which exist only for some cutoffs) are skipped.
    channels = [
        (iChn, cats.GetChannelWeight(iChn))
        for iChn in range(cats.GetNumChannels())
        if cats.GetChannelWeight(iChn) > 0
    ]
    print(f'Summing {len(channels)} channels with weights: '
          + ', '.join(f'{iChn}: {w:.4f}' for iChn, w in channels))

    wf = sum([get_wave_function(cats, iChn) * w for iChn, w in channels])

    hWF = matrix_to_th2d(wf)
    hWF.SetTitle(title)

    root_path = str(Path(oFile).with_suffix('.root'))
    fout = TFile(root_path, 'RECREATE')
    gCF.Write('gCF')
    hWF.Write('hWF')
    fout.Close()
    print(f'correlation function and |#psi|^2 (.root) written to {root_path}')

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
    parser.add_argument('system', choices=('pp', 'pL'))
    parser.add_argument('oFile')
    args = parser.parse_args()
    
    compute_wave_function(args.system, args.oFile)
