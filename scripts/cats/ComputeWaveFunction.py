import os
import math
import numpy as np
from ROOT import gSystem, gInterpreter, TH2D, TDatabasePDG, TGraph, TFile

cats_dir = os.environ["CATS"]

gSystem.Load(f"{cats_dir}/install/lib/libCATS.so")
gInterpreter.AddIncludePath(f"{cats_dir}/install/include")

gInterpreter.Declare("""
#include "CATS.h"
#include "CommonAnaFunctions.h"
""")

from ROOT import CATS
from ROOT import DLM_CommonAnaFunctions

RADIUS_STEP = 0.01
RADIUS_MAX = 6.25
KSTAR_STEP = 0.4
KSTAR_MAX = 400.

def table_to_th2d(matrix):
    h = TH2D(
        "hWF", "",
        round(RADIUS_MAX / RADIUS_STEP), 0, RADIUS_MAX,
        round(KSTAR_MAX / KSTAR_STEP), 0, KSTAR_MAX
    )

    ny, nx = matrix.shape

    for iy in range(ny):
        for ix in range(nx):
            h.SetBinContent(ix + 1, iy + 1, float(matrix[iy, ix]))

    return h

def to_file(filename, matrix):
    kstar = KSTAR_STEP/2 + np.arange(matrix.shape[0]) * KSTAR_STEP
    radius = RADIUS_STEP/2 + np.arange(matrix.shape[1]) * RADIUS_STEP

    with open(filename, "w") as f:
        f.write("radius")
        for k in kstar:
            f.write(f"\t{k:.6f}")
        f.write("\n")

        for ir, r in enumerate(radius):
            f.write(f"{r:.6f}")
            for ik in range(matrix.shape[0]):
                f.write(f" {matrix[ik, ir]:.6f}")
            f.write("\n")

def get_wave_function(cats, weight=1., channel=-1):
    if channel < 0:
        n_channels = cats.GetNumChannels()

        if n_channels == 2:
            weights = [0.25, 0.75]
        elif n_channels == 4:
            weights = [3/12, 1/12, 3/12, 5/12]
        else:
            raise RuntimeError(f"Unexpected number of channels: {n_channels}")

        wf = get_wave_function(cats, weights[0], 0)

        for ich in range(1, n_channels):
            wf += get_wave_function(cats, weights[ich], ich)

        return wf

    nk = cats.GetNumMomBins()
    nr = round(RADIUS_MAX / RADIUS_STEP)

    wf = np.empty((nk, nr))

    for ik in range(nk):
        for ir in range(nr):
            radius = RADIUS_STEP/2 + ir * RADIUS_STEP
            wf[ik, ir] = (weight * cats.EvalWaveFun2(ik, radius, channel))

    return wf

def compute_wave_function(pdg1=2212, pdg2=2212, r0=1.25):
    pdg = TDatabasePDG.Instance()

    m1 = 1000. * pdg.GetParticle(pdg1).Mass()
    m2 = 1000. * pdg.GetParticle(pdg2).Mass()

    mu = m1 * m2 / (m1 + m2)
    n_kstar_bins = round(KSTAR_MAX / KSTAR_STEP)

    cats = CATS()
    cats.SetMomBins(n_kstar_bins, 0, KSTAR_MAX)
    cats.SetQ1Q2(1)
    cats.SetQuantumStatistics(True)
    cats.SetRedMass(mu)

    can = DLM_CommonAnaFunctions()

    if pdg1 == 2212 and pdg2 == 2212:
        header = "Wave function of proton-proton with Argonne v18 potential computed with CATS\n"
        can.SetUpCats_pp(cats, "AV18", "Gauss", 0, 0)
    else:
        raise RuntimeError("System not implemented")

    cats.SetAnaSource(0, r0)
    cats.KillTheCat()

    gCF = TGraph()

    for i in range(n_kstar_bins):
        gCF.SetPoint(i, cats.GetMomentum(i), cats.GetCorrFun(i))

    wf = get_wave_function(cats)

    to_file("wf.dat", wf)

    hWF = table_to_th2d(wf)
    hWF.SetTitle("pp, AV18, |#psi|^{2};r (fm);k* (MeV/c);|#psi|^{2}")

    fout = TFile("ppCF_WF_Source.root", "RECREATE")
    gCF.Write("gCF")
    hWF.Write("hWF")
    fout.Close()

    radius = RADIUS_STEP/2 + np.arange(wf.shape[1]) * RADIUS_STEP
    kstar = KSTAR_STEP/2 + np.arange(wf.shape[0]) * KSTAR_STEP

    table = np.column_stack([radius, wf.T])
    header += (
        f"radius step = {RADIUS_STEP} fm\n"
        f"kstar step = {KSTAR_STEP} MeV/c\n"
        "radius " + " ".join(f"{k:.3f}" for k in kstar)
    )

    np.savetxt(
        "wf2.dat",
        table,
        fmt="%.17e",
        header=header
    )

if __name__ == "__main__":
    compute_wave_function()
