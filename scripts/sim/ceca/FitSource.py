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

from ROOT import gInterpreter, TFile, TF1, TCanvas, TLegend

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

fPrim = 0.3578


def main(args):
    oFileBase = Path(args.output).stem

    match = re.search(r"Q3lt(\d+)MeV", Path(args.input).name)
    if not match:
        raise ValueError(f"Could not extract Q3 cut from input file name '{args.input}'")
    Q3cut = int(match.group(1))

    inFile = TFile(args.input)

    # Extract rPrim and rSec from the 2B source
    hRStar = inFile.Get("triplet/hRStarVsMt").ProjectionY()

    fSource2B = TF1("fSource2B", SourceCountsGaussResonances, 0, 12, 4)
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
    fSource2B_pp.SetParameter(0, fPrim * fPrim * fSource2B.GetParameter(0))
    fSource2B_pp.SetParameter(1, rPrim)
    fSource2B_pp.SetLineStyle(7)
    fSource2B_pp.Draw("same")

    fSource2B_ps = TF1("fSource2B_ps", SourceCountsGauss, 0, 20, 2)
    fSource2B_ps.SetParameter(0, 2 * fPrim * (1 - fPrim) * fSource2B.GetParameter(0))
    fSource2B_ps.SetParameter(1, np.sqrt((rPrim * rPrim + rSec * rSec) / 2))
    fSource2B_ps.SetLineStyle(8)
    fSource2B_ps.Draw("same")

    fSource2B_ss = TF1("fSource2B_ss", SourceCountsGauss, 0, 20, 2)
    fSource2B_ss.SetParameter(0, (1 - fPrim) * (1 - fPrim) * fSource2B.GetParameter(0))
    fSource2B_ss.SetParameter(1, rSec)
    fSource2B_ss.SetLineStyle(9)
    fSource2B_ss.Draw("same")

    leg = TLegend(0.5, 0.5, 0.9, 0.85)
    leg.SetHeader(f"2B, Q_{{3}} < {Q3cut} MeV/c")
    leg.AddEntry(hRStar, "Total", "pel")
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

    fSourceCountsAAAGaussResonances = TF1("fSourceCountsAAAGaussResonances", SourceCountsAAAGaussResonances, 0, 20, 4)
    fSourceCountsAAAGaussResonances.SetParameter(0, hHypRad.GetEntries() / 10)
    fSourceCountsAAAGaussResonances.FixParameter(1, fPrim)
    fSourceCountsAAAGaussResonances.SetParameter(2, rPrim)
    fSourceCountsAAAGaussResonances.SetParameter(3, rSec)

    fSource3B_ppp = TF1("fSource3B_ppp", SourceCountsAAA, 0, 20, 2)
    fSource3B_ppp.SetParameter(0, fPrim**3 * hHypRad.GetEntries() * hHypRad.GetBinWidth(1))
    fSource3B_ppp.SetParameter(1, 2 * rPrim)
    fSource3B_ppp.SetLineStyle(7)

    fSource3B_pps = TF1("fSource3B_pps", SourceCountsAAApprAvg, 0, 20, 3)
    fSource3B_pps.SetParameter(0, 3 * fPrim**2 * (1 - fPrim) * hHypRad.GetEntries() * hHypRad.GetBinWidth(1))
    fSource3B_pps.SetParameter(1, rPrim)
    fSource3B_pps.SetParameter(2, rSec)
    fSource3B_pps.SetLineStyle(8)

    fSource3B_pss = TF1("fSource3B_pss", SourceCountsAAAprrAvg, 0, 20, 3)
    fSource3B_pss.SetParameter(0, 3 * fPrim * (1 - fPrim) ** 2 * hHypRad.GetEntries() * hHypRad.GetBinWidth(1))
    fSource3B_pss.SetParameter(1, rPrim)
    fSource3B_pss.SetParameter(2, rSec)
    fSource3B_pss.SetLineStyle(9)

    fSource3B_sss = TF1("fSource3B_sss", SourceCountsAAA, 0, 20, 2)
    fSource3B_sss.SetParameter(0, (1 - fPrim) ** 3 * hHypRad.GetEntries() * hHypRad.GetBinWidth(1))
    fSource3B_sss.SetParameter(1, 2 * rSec)
    fSource3B_sss.SetLineStyle(10)

    c3B = TCanvas("c3B", "", 600, 600)
    c3B.DrawFrame(0, 0, 20, 1.3 * hHypRad.GetMaximum(), ";#rho (fm);Counts")
    hHypRad.Draw("pe same")

    fSourceCountsAAAGaussResonances.SetLineColor(2)
    fSourceCountsAAAGaussResonances.Draw("same")

    fSource3B_ppp.Draw("same")
    fSource3B_pps.Draw("same")
    fSource3B_pss.Draw("same")
    fSource3B_sss.Draw("same")

    leg3B = TLegend(0.55, 0.5, 0.9, 0.9)
    leg3B.SetHeader(f"3B, Q_{{3}} < {Q3cut} MeV/c")
    leg3B.AddEntry(hHypRad, "Total", "pel")
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

    fSourceCountsAAAGaussResonancesHypAngle = TF1(
        "fSourceCountsAAAGaussResonancesHypAngle", SourceCountsAAAGaussResonancesHypAngle, 0, np.pi / 2, 4
    )
    fSourceCountsAAAGaussResonancesHypAngle.SetParameter(0, hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSourceCountsAAAGaussResonancesHypAngle.FixParameter(1, fPrim)
    fSourceCountsAAAGaussResonancesHypAngle.FixParameter(2, rPrim)
    fSourceCountsAAAGaussResonancesHypAngle.FixParameter(3, rSec)

    fSourceHypAngle_ppp = TF1("fSourceHypAngle_ppp", SourceCountsAAAHypAngle, 0, np.pi / 2, 1)
    fSourceHypAngle_ppp.SetParameter(0, fPrim**3 * hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSourceHypAngle_ppp.SetLineStyle(7)

    fSourceHypAngle_pps = TF1("fSourceHypAngle_pps", SourceCountsAAApprHypAngle, 0, np.pi / 2, 3)
    fSourceHypAngle_pps.SetParameter(0, 3 * fPrim**2 * (1 - fPrim) * hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSourceHypAngle_pps.SetParameter(1, rPrim)
    fSourceHypAngle_pps.SetParameter(2, rSec)
    fSourceHypAngle_pps.SetLineStyle(8)

    fSourceHypAngle_pss = TF1("fSourceHypAngle_pss", SourceCountsAAAprrHypAngle, 0, np.pi / 2, 3)
    fSourceHypAngle_pss.SetParameter(0, 3 * fPrim * (1 - fPrim) ** 2 * hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSourceHypAngle_pss.SetParameter(1, rPrim)
    fSourceHypAngle_pss.SetParameter(2, rSec)
    fSourceHypAngle_pss.SetLineStyle(9)

    fSourceHypAngle_sss = TF1("fSourceHypAngle_sss", SourceCountsAAAHypAngle, 0, np.pi / 2, 1)
    fSourceHypAngle_sss.SetParameter(0, (1 - fPrim) ** 3 * hHypAngle.GetEntries() * hHypAngle.GetBinWidth(1))
    fSourceHypAngle_sss.SetLineStyle(10)

    cHypAngle = TCanvas("cHypAngle", "", 600, 600)
    cHypAngle.DrawFrame(0, 0, np.pi / 2, 1.6 * hHypAngle.GetMaximum(), ";#varphi (rad);Counts")
    hHypAngle.Draw("pe same")

    fSourceCountsAAAGaussResonancesHypAngle.SetLineColor(2)
    fSourceCountsAAAGaussResonancesHypAngle.Draw("same")

    fSourceHypAngle_ppp.Draw("same")
    fSourceHypAngle_pps.Draw("same")
    fSourceHypAngle_pss.Draw("same")
    fSourceHypAngle_sss.Draw("same")

    legHypAngle = TLegend(0.55, 0.6, 0.9, 0.9)
    legHypAngle.SetHeader(f"3B, Q_{{3}} < {Q3cut} MeV/c")
    legHypAngle.AddEntry(hHypAngle, "Total", "pel")
    legHypAngle.AddEntry(fSourceHypAngle_ppp, "ppp", "l")
    legHypAngle.AddEntry(fSourceHypAngle_pps, "pps", "l")
    legHypAngle.AddEntry(fSourceHypAngle_pss, "pss", "l")
    legHypAngle.AddEntry(fSourceHypAngle_sss, "sss", "l")

    legHypAngle.Draw("same")
    cHypAngle.SaveAs(f"{oFileBase}_HypAngle.pdf")

    # Save the fit functions
    oFile = TFile(args.output, "RECREATE")
    fSource2B.Write()
    fSource2B_pp.Write()
    fSource2B_ps.Write()
    fSource2B_ss.Write()
    fSourceCountsAAAGaussResonances.Write()
    fSource3B_ppp.Write()
    fSource3B_pps.Write()
    fSource3B_pss.Write()
    fSource3B_sss.Write()
    fSourceCountsAAAGaussResonancesHypAngle.Write()
    fSourceHypAngle_ppp.Write()
    fSourceHypAngle_pps.Write()
    fSourceHypAngle_pss.Write()
    fSourceHypAngle_sss.Write()
    oFile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to the input root file")
    parser.add_argument("output", help="Path to the output root file")
    args = parser.parse_args()

    main(args)
