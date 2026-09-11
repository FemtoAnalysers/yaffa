import argparse
import os
import re
from dotenv import load_dotenv
from pathlib import Path
import numpy as np

from yaffa import utils

env_path = Path(__file__).resolve().parents[3] / ".env"
print(f"Loading env from {env_path}")
if not load_dotenv(dotenv_path=env_path, verbose=True):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
if not YAFFA_PATH:
    print("\033[33mWARNING: Path to yaffa is empty, something might break!\033[0m")

from ROOT import gInterpreter, TFile, TF1, TCanvas, TLegend, TGraphPainter

gInterpreter.Declare(f'#include "{YAFFA_PATH}/src/cpp/RootFunctions.hxx"')
from ROOT import (
    SourceCountsGauss,
    SourceCountsGaussResonances,
    SourceCountsAAAprrAvg,
    SourceCountsAAApprAvg,
    SourceCountsAAAHypAngle,
    SourceCountsAAApprHypAngle,
    SourceCountsAAAprrHypAngle,
    SourceCountsAAAGaussResonances,
    SourceCountsAAAGaussResonancesHypAngle,
    SourceCountsAAA,
)

utils.style.SetStyle()

# This ensures that TF1s with large Npx are drawn with correct line style (e.g. dashed).
# High granularity TF1s are needed to avoid catching wavefunction fluctuations when using the Koonin-Pratt formula.
TGraphPainter.SetMaxPointsPerLine(1000000)


def main(args):
    oFileBase = Path(args.output).stem

    match = re.search(r"Q3lt(\d+)MeV", Path(args.input).name)
    if not match:
        raise ValueError(f"Could not extract Q3 cut from input file name '{args.input}'")
    Q3cut = int(match.group(1))

    inFile = TFile(args.input)

    if args.fPrim is not None:
        fPrim = args.fPrim
    else:
        nPPP = inFile.Get("triplet/hHypRad_ppp").GetEntries()
        n = inFile.Get("triplet/hHypRad").GetEntries()
        fPrim = (nPPP / n) ** (1. / 3.)

    # Extract rPrim and rSec from the 2B source
    hRStar = inFile.Get("triplet/hRStarVsMt").ProjectionY()

    fSource2B = TF1("fSource2B", SourceCountsGaussResonances, 0, 12, 4)
    fSource2B.SetNpx(100000)
    fSource2B.SetParameter(0, 1.0e5)
    fSource2B.FixParameter(1, fPrim)
    fSource2B.SetParameter(2, 1)
    fSource2B.SetParLimits(2, 0.1, 5)
    fSource2B.SetParameter(3, 1)
    fSource2B.SetParLimits(3, 0, 5)
    hRStar.Fit(fSource2B, "0MRL+")

    rPrim = fSource2B.GetParameter(2)
    rSec = rPrim + fSource2B.GetParameter(3)

    # Draw
    c2B = TCanvas("c2B", "", 600, 600)
    c2B.DrawFrame(0, 0, 20, 1.3 * hRStar.GetMaximum(), ";r* (fm);Counts")
    hRStar.Draw("pe same")

    fSource2B.SetLineColor(2)
    fSource2B.Draw("same")

    fSource2B_pp = TF1("fSource2B_pp", SourceCountsGauss, 0, 20, 2)
    fSource2B_pp.SetNpx(100000)
    fSource2B_pp.SetParameter(0, fPrim * fPrim * fSource2B.GetParameter(0))
    fSource2B_pp.SetParameter(1, rPrim)
    fSource2B_pp.SetLineStyle(7)
    fSource2B_pp.Draw("same")

    fSource2B_ps = TF1("fSource2B_ps", SourceCountsGauss, 0, 20, 2)
    fSource2B_ps.SetNpx(100000)
    fSource2B_ps.SetParameter(0, 2 * fPrim * (1 - fPrim) * fSource2B.GetParameter(0))
    fSource2B_ps.SetParameter(1, np.sqrt((rPrim * rPrim + rSec * rSec) / 2))
    fSource2B_ps.SetLineStyle(8)
    fSource2B_ps.Draw("same")

    fSource2B_ss = TF1("fSource2B_ss", SourceCountsGauss, 0, 20, 2)
    fSource2B_ss.SetNpx(100000)
    fSource2B_ss.SetParameter(0, (1 - fPrim) * (1 - fPrim) * fSource2B.GetParameter(0))
    fSource2B_ss.SetParameter(1, rSec)
    fSource2B_ss.SetLineStyle(9)
    fSource2B_ss.Draw("same")

    leg = TLegend(0.5, 0.5, 0.9, 0.85)
    leg.SetHeader(f"2B, Q_{{3}} < {Q3cut} MeV/c")
    leg.AddEntry(hRStar, f"Total, f_{{prim}} = {fPrim * 100:.2f}%", "pel")
    leg.AddEntry(fSource2B_pp, f"pp: r_{{p}} = {rPrim:.2f} fm", "l")
    leg.AddEntry(
        fSource2B_ps,
        f"ps: r_{{ps}} = {np.sqrt((rPrim * rPrim + rSec * rSec) / 2):.2f} fm",
        "l",
    )
    leg.AddEntry(fSource2B_ss, f"ss: r_{{s}} = {rSec:.2f} fm", "l")

    leg.Draw("same")
    c2B.SaveAs(f"{oFileBase}_2B.pdf")

    # 3B
    hHypRad = inFile.Get("triplet/hHypRad")

    fSource3B = TF1("fSource3B", SourceCountsAAAGaussResonances, 0, 20, 4)
    fSource3B.SetNpx(100000)
    fSource3B.SetParameter(0, hHypRad.GetEntries() / 10)
    fSource3B.FixParameter(1, fPrim)
    fSource3B.SetParameter(2, rPrim)
    fSource3B.SetParameter(3, rSec)

    fSource3B_ppp = TF1("fSource3B_ppp", SourceCountsAAA, 0, 20, 2)
    fSource3B_ppp.SetNpx(100000)
    fSource3B_ppp.SetParameter(0, fPrim**3 * hHypRad.GetEntries() * hHypRad.GetBinWidth(1))
    fSource3B_ppp.SetParameter(1, 2 * rPrim)
    fSource3B_ppp.SetLineStyle(7)

    fSource3B_pps = TF1("fSource3B_pps", SourceCountsAAApprAvg, 0, 20, 3)
    fSource3B_pps.SetNpx(100000)
    fSource3B_pps.SetParameter(0, 3 * fPrim**2 * (1 - fPrim) * hHypRad.GetEntries() * hHypRad.GetBinWidth(1))
    fSource3B_pps.SetParameter(1, rPrim)
    fSource3B_pps.SetParameter(2, rSec)
    fSource3B_pps.SetLineStyle(8)

    fSource3B_pss = TF1("fSource3B_pss", SourceCountsAAAprrAvg, 0, 20, 3)
    fSource3B_pss.SetNpx(100000)
    fSource3B_pss.SetParameter(0, 3 * fPrim * (1 - fPrim) ** 2 * hHypRad.GetEntries() * hHypRad.GetBinWidth(1))
    fSource3B_pss.SetParameter(1, rPrim)
    fSource3B_pss.SetParameter(2, rSec)
    fSource3B_pss.SetLineStyle(9)

    fSource3B_sss = TF1("fSource3B_sss", SourceCountsAAA, 0, 20, 2)
    fSource3B_sss.SetNpx(100000)
    fSource3B_sss.SetParameter(0, (1 - fPrim) ** 3 * hHypRad.GetEntries() * hHypRad.GetBinWidth(1))
    fSource3B_sss.SetParameter(1, 2 * rSec)
    fSource3B_sss.SetLineStyle(10)

    c3B = TCanvas("c3B", "", 600, 600)
    c3B.DrawFrame(0, 0, 20, 1.3 * hHypRad.GetMaximum(), ";#rho (fm);Counts")
    hHypRad.Draw("pe same")

    fSource3B.SetLineColor(2)
    fSource3B.Draw("same")

    fSource3B_ppp.Draw("same")
    fSource3B_pps.Draw("same")
    fSource3B_pss.Draw("same")
    fSource3B_sss.Draw("same")

    leg3B = TLegend(0.55, 0.5, 0.9, 0.9)
    leg3B.SetHeader(f"3B, Q_{{3}} < {Q3cut} MeV/c")
    leg3B.AddEntry(hHypRad, f"Total, f_{{prim}} = {fPrim * 100:.2f}%", "pel")
    leg3B.AddEntry(fSource3B_ppp, f"ppp: #rho_{{ppp}} = {2 * rPrim:.2f} fm", "l")
    leg3B.AddEntry(fSource3B_pps, f"pps", "l")
    leg3B.AddEntry(fSource3B_pss, f"pss", "l")
    leg3B.AddEntry(fSource3B_sss, f"sss: #rho_{{sss}} = {2 * rSec:.2f} fm", "l")

    leg3B.Draw("same")
    c3B.SaveAs(f"{oFileBase}_3B.pdf")

    # 3B hyper-angle: compare with the full ppp/pps/pss/sss mixture, using fPrim, rPrim
    # and rSec from the 2B r* fit above
    hHypAngle = inFile.Get("triplet/hHypAngleVsHypRad").ProjectionY("hHypAngle")
    hHypAngle.Rebin(5)

    fSource3BHypAngle = TF1(
        "fSource3BHypAngle", SourceCountsAAAGaussResonancesHypAngle, 0, np.pi / 2, 4
    )
    fSource3BHypAngle.SetNpx(100000)
    fSource3BHypAngle.SetParameter(0, hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSource3BHypAngle.FixParameter(1, fPrim)
    fSource3BHypAngle.FixParameter(2, rPrim)
    fSource3BHypAngle.FixParameter(3, rSec)

    fSource3BHypAngle_ppp = TF1("fSource3BHypAngle_ppp", SourceCountsAAAHypAngle, 0, np.pi / 2, 1)
    fSource3BHypAngle_ppp.SetNpx(100000)
    fSource3BHypAngle_ppp.SetParameter(0, fPrim**3 * hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSource3BHypAngle_ppp.SetLineStyle(7)

    fSource3BHypAngle_pps = TF1("fSource3BHypAngle_pps", SourceCountsAAApprHypAngle, 0, np.pi / 2, 3)
    fSource3BHypAngle_pps.SetNpx(100000)
    fSource3BHypAngle_pps.SetParameter(0, 3 * fPrim**2 * (1 - fPrim) * hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSource3BHypAngle_pps.SetParameter(1, rPrim)
    fSource3BHypAngle_pps.SetParameter(2, rSec)
    fSource3BHypAngle_pps.SetLineStyle(8)

    fSource3BHypAngle_pss = TF1("fSource3BHypAngle_pss", SourceCountsAAAprrHypAngle, 0, np.pi / 2, 3)
    fSource3BHypAngle_pss.SetNpx(100000)
    fSource3BHypAngle_pss.SetParameter(0, 3 * fPrim * (1 - fPrim) ** 2 * hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSource3BHypAngle_pss.SetParameter(1, rPrim)
    fSource3BHypAngle_pss.SetParameter(2, rSec)
    fSource3BHypAngle_pss.SetLineStyle(9)

    fSource3BHypAngle_sss = TF1("fSource3BHypAngle_sss", SourceCountsAAAHypAngle, 0, np.pi / 2, 1)
    fSource3BHypAngle_sss.SetNpx(100000)
    fSource3BHypAngle_sss.SetParameter(0, (1 - fPrim) ** 3 * hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSource3BHypAngle_sss.SetLineStyle(10)

    cHypAngle = TCanvas("cHypAngle", "", 600, 600)
    cHypAngle.DrawFrame(0, 0, np.pi / 2, 1.6 * hHypAngle.GetMaximum(), ";#varphi (rad);Counts")
    hHypAngle.Draw("pe same")

    fSource3BHypAngle.SetLineColor(2)
    fSource3BHypAngle.Draw("same")

    fSource3BHypAngle_ppp.Draw("same")
    fSource3BHypAngle_pps.Draw("same")
    fSource3BHypAngle_pss.Draw("same")
    fSource3BHypAngle_sss.Draw("same")

    legHypAngle = TLegend(0.22, 0.6, 0.5, 0.9)
    legHypAngle.SetHeader(f"3B, Q_{{3}} < {Q3cut} MeV/c")
    legHypAngle.AddEntry(hHypAngle, f"Total, f_{{prim}} = {fPrim * 100:.2f}%", "pel")
    legHypAngle.AddEntry(fSource3BHypAngle_ppp, "ppp", "l")
    legHypAngle.AddEntry(fSource3BHypAngle_pps, "pps", "l")
    legHypAngle.AddEntry(fSource3BHypAngle_pss, "pss", "l")
    legHypAngle.AddEntry(fSource3BHypAngle_sss, "sss", "l")

    legHypAngle.Draw("same")
    cHypAngle.SaveAs(f"{oFileBase}_HypAngle.pdf")

    # theta12 and theta3_12: polar angles of r12 and r3_12. Since these vectors are
    # isotropically distributed, the expected pdf is 0.5 * sin(theta)
    hTheta12 = inFile.Get("triplet/hTheta12")
    hTheta312 = inFile.Get("triplet/hTheta3_12")

    fTheta12 = TF1("fTheta12", "[0] * 0.5 * sin(x)", 0, np.pi)
    fTheta12.SetNpx(100000)
    fTheta12.SetParameter(0, hTheta12.GetEntries() * hTheta12.GetBinWidth(1))

    fTheta312 = TF1("fTheta312", "[0] * 0.5 * sin(x)", 0, np.pi)
    fTheta312.SetNpx(100000)
    fTheta312.SetParameter(0, hTheta312.GetEntries() * hTheta312.GetBinWidth(1))

    cTheta12 = TCanvas("cTheta12", "", 600, 600)
    cTheta12.DrawFrame(0, 0, np.pi, 1.6 * hTheta12.GetMaximum(), ";#theta_{12} (rad);Counts")
    hTheta12.Draw("pe same")
    fTheta12.SetLineColor(2)
    fTheta12.Draw("same")

    legTheta12 = TLegend(0.55, 0.7, 0.9, 0.9)
    legTheta12.SetHeader(f"3B, Q_{{3}} < {Q3cut} MeV/c")
    legTheta12.AddEntry(hTheta12, f"Total, f_{{prim}} = {fPrim * 100:.2f}%", "pel")
    legTheta12.AddEntry(fTheta12, "0.5 sin(#theta_{12})", "l")
    legTheta12.Draw("same")
    cTheta12.SaveAs(f"{oFileBase}_Theta12.pdf")

    cTheta312 = TCanvas("cTheta312", "", 600, 600)
    cTheta312.DrawFrame(0, 0, np.pi, 1.6 * hTheta312.GetMaximum(), ";#theta_{3,12} (rad);Counts")
    hTheta312.Draw("pe same")
    fTheta312.SetLineColor(2)
    fTheta312.Draw("same")

    legTheta312 = TLegend(0.55, 0.7, 0.9, 0.9)
    legTheta312.SetHeader(f"3B, Q_{{3}} < {Q3cut} MeV/c")
    legTheta312.AddEntry(hTheta312, f"Total, f_{{prim}} = {fPrim * 100:.2f}%", "pel")
    legTheta312.AddEntry(fTheta312, "0.5 sin(#theta_{3,12})", "l")
    legTheta312.Draw("same")
    cTheta312.SaveAs(f"{oFileBase}_Theta312.pdf")

    # Save the fit functions
    oFile = TFile(args.output, "RECREATE")
    fSource2B.Write()
    fSource2B_pp.Write()
    fSource2B_ps.Write()
    fSource2B_ss.Write()
    fSource3B.Write()
    fSource3B_ppp.Write()
    fSource3B_pps.Write()
    fSource3B_pss.Write()
    fSource3B_sss.Write()
    fSource3BHypAngle.Write()
    fSource3BHypAngle_ppp.Write()
    fSource3BHypAngle_pps.Write()
    fSource3BHypAngle_pss.Write()
    fSource3BHypAngle_sss.Write()
    fTheta12.Write()
    fTheta312.Write()
    oFile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to the input root file")
    parser.add_argument("output", help="Path to the output root file")
    parser.add_argument(
        "--fPrim",
        type=float,
        default=None,
        help="Fraction of primary particles. If not given, it is computed as "
        "(N_ppp / N)^(1/3) from triplet/hHypRad_ppp and triplet/hHypRad",
    )
    args = parser.parse_args()

    main(args)
