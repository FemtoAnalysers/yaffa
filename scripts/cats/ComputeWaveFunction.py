import os
import math
import ROOT

# ------------------------------------------------------------------
# Load CATS
# ------------------------------------------------------------------

cats_dir = os.environ["CATS"]
ROOT.gSystem.Load(f"{cats_dir}/install/lib/libCATS.so")

ROOT.gInterpreter.AddIncludePath(f"{cats_dir}/install/include")

ROOT.gInterpreter.Declare("""
#include "CATS.h"
#include "CATSconstants.h"
#include "CommonAnaFunctions.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomp.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_Fitters.h"
#include "DLM_Histo.h"
#include "DLM_Potentials.h"
#include "DLM_ResponseMatrix.h"
#include "DLM_RootFit.h"
#include "DLM_RootWrapper.h"
#include "DLM_Source.h"
""")

# ------------------------------------------------------------------

RADIUS_STEP = 0.01
RADIUS_MAX = 6.25

KSTAR_STEP = 0.4
KSTAR_MAX = 400.0

# ------------------------------------------------------------------

def matrix_sum(m1, m2):
    return [
        [m1[i][j] + m2[i][j] for j in range(len(m1[i]))]
        for i in range(len(m1))
    ]


def table_to_th2d(matrix):

    h = ROOT.TH2D(
        "hWF",
        "",
        round(RADIUS_MAX / RADIUS_STEP),
        0,
        RADIUS_MAX,
        round(KSTAR_MAX / KSTAR_STEP),
        0,
        KSTAR_MAX,
    )

    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            h.SetBinContent(j + 1, i + 1, matrix[i][j])

    return h


def to_file(filename, matrix):

    with open(filename, "w") as out:

        out.write("radius")

        for ik in range(len(matrix)):
            out.write(f"\t{KSTAR_STEP/2 + KSTAR_STEP*ik:.6f}")

        out.write("\n")

        for ir in range(len(matrix[0])):

            out.write(f"{RADIUS_STEP/2 + RADIUS_STEP*ir:.6f}")

            for ik in range(len(matrix)):
                out.write(f" {matrix[ik][ir]:.6f}")

            out.write("\n")


def get_wave_function(cats, weight=1.0, channel=-1):

    if channel < 0:

        n_channels = cats.GetNumChannels()

        if n_channels == 2:
            weights = [0.25, 0.75]

        elif n_channels == 4:
            weights = [3/12, 1/12, 3/12, 5/12]

        else:
            raise RuntimeError(
                f"Unexpected number of channels: {n_channels}"
            )

        wf = get_wave_function(cats, weights[0], 0)

        for ich in range(1, n_channels):
            wf = matrix_sum(
                wf,
                get_wave_function(cats, weights[ich], ich)
            )

        return wf

    wf = []

    for ik in range(cats.GetNumMomBins()):

        row = []

        radius = RADIUS_STEP / 2

        while radius < RADIUS_MAX:

            value = cats.EvalWaveFun2(
                ik,
                radius,
                channel
            )

            row.append(value * weight)

            radius += RADIUS_STEP

        wf.append(row)

    return wf


# ------------------------------------------------------------------

def compute_wave_function(
        pdg1=2212,
        pdg2=2212,
        r0=1.25):

    pdg = ROOT.TDatabasePDG.Instance()

    m1 = 1000.0 * pdg.GetParticle(pdg1).Mass()
    m2 = 1000.0 * pdg.GetParticle(pdg2).Mass()

    mu = m1 * m2 / (m1 + m2)

    n_kstar_bins = round(KSTAR_MAX / KSTAR_STEP)

    cats = ROOT.CATS()

    cats.SetMomBins(
        n_kstar_bins,
        0,
        KSTAR_MAX
    )

    cats.SetQ1Q2(1)
    cats.SetQuantumStatistics(True)
    cats.SetRedMass(mu)

    can = ROOT.DLM_CommonAnaFunctions()

    if pdg1 == 2212 and pdg2 == 2212:
        can.SetUpCats_pp(
            cats,
            "AV18",
            "Gauss",
            0,
            0
        )
    else:
        raise RuntimeError(
            "System not implemented"
        )

    cats.SetAnaSource(0, r0)
    cats.KillTheCat()

    gCF = ROOT.TGraph()

    for i in range(n_kstar_bins):

        gCF.SetPoint(
            i,
            cats.GetMomentum(i),
            cats.GetCorrFun(i)
        )

    wf = get_wave_function(cats)

    to_file("wf.dat", wf)

    hWF = table_to_th2d(wf)

    hWF.SetTitle(
        "pp, AV18, |#psi|^{2};r (fm);k* (MeV/c);|#psi|^{2}"
    )

    fout = ROOT.TFile(
        "ppCF_WF_Source.root",
        "RECREATE"
    )

    gCF.Write("gCF")
    hWF.Write("hWF")

    fout.Close()


# ------------------------------------------------------------------

if __name__ == "__main__":
    compute_wave_function()
