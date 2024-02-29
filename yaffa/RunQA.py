'''
Script to produce the QA plots.
'''

import os
import argparse
import yaml

from ROOT import TFile, TCanvas, EColor, TLegend, gROOT  # pylint: disable=import-error

from yaffa import logger as log
from yaffa.utils.io import Load

parser = argparse.ArgumentParser()
parser.add_argument('cfg', help='Config file')
parser.add_argument('-b', action='store_true', default=False, help='Set batch mode')
args = parser.parse_args()

gROOT.SetBatch(args.b)

# Load yaml file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

suffix = cfg['suffix']
inFile = TFile(cfg['infile'])

if cfg['pair'] == 'dPi':
    part1 = 'Pi'
    part2 = 'd'
    fdnames = {
        0: 'PosPions',
        1: 'NegPions',
        2: 'DeuteronDCA',
        3: 'AntiDeuteronDCA',
    }
    labels = {
        0: '#pi^{+}',
        1: '#pi^{#minus}',
        2: 'd',
        3: '#bar{d}',
    }
else:
    log.error('pair not implemented')


# Event quality
oFileBaseName = 'QA_Evt'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'

oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')
oFile = TFile(oFileName, 'recreate')

cEvt = TCanvas('cEvt', '', 600, 600)
cEvt.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf['))

# Multiplicity
hMult = Load(inFile, f'HMEvtCuts{suffix}/HMEvtCuts{suffix}/after/MultiplicityRef08_after')
hMult.Draw()
cEvt.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# Sphericity
hSphercity = Load(inFile, f'HMEvtCuts{suffix}/HMEvtCuts{suffix}/after/Sphericity_after')
hSphercity.Draw()
cEvt.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

cEvt.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf]'))
# End of event quality plots


# Part1 quality
oFileBaseName = f'QA_{part1}'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'

oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')
oFile = TFile(oFileName, 'recreate')
cPart = TCanvas('cPart1', '', 600, 600)
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf['))

# Pt of Part1
hPt0 = Load(inFile, f'HM{fdnames[0]}{suffix}/HM{fdnames[0]}{suffix}/after/pTDist_after')
hPt0.Draw()
hPt1 = Load(inFile, f'HM{fdnames[1]}{suffix}/HM{fdnames[1]}{suffix}/after/pTDist_after')
hPt1.SetLineColor(EColor.kRed)
hPt1.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hPt0, labels[0])
leg.AddEntry(hPt1, labels[1])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# Eta of Part1
hEta0 = Load(inFile, f'HM{fdnames[0]}{suffix}/HM{fdnames[0]}{suffix}/after/EtaDist_after')
hEta0.Draw()
hEta1 = Load(inFile, f'HM{fdnames[1]}{suffix}/HM{fdnames[1]}{suffix}/after/EtaDist_after')
hEta1.SetLineColor(EColor.kRed)
hEta1.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hEta0, labels[0])
leg.AddEntry(hEta1, labels[1])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# Phi of Part1
hPhi0 = Load(inFile, f'HM{fdnames[0]}{suffix}/HM{fdnames[0]}{suffix}/after/phiDist_after')
hPhi0.Draw()
hPhi1 = Load(inFile, f'HM{fdnames[1]}{suffix}/HM{fdnames[1]}{suffix}/after/phiDist_after')
hPhi1.SetLineColor(EColor.kRed)
hPhi1.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hPhi0, labels[0])
leg.AddEntry(hPhi1, labels[1])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# CrossedRows of Part1
hCrossedRows0 = Load(inFile, f'HM{fdnames[0]}{suffix}/HM{fdnames[0]}{suffix}/after/CrossedRows_after')
hCrossedRows0.Draw('')
hCrossedRows1 = Load(inFile, f'HM{fdnames[1]}{suffix}/HM{fdnames[1]}{suffix}/after/CrossedRows_after')
hCrossedRows1.SetLineColor(EColor.kRed)
hCrossedRows1.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hCrossedRows0, labels[0])
leg.AddEntry(hCrossedRows1, labels[1])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# ClusterOverFindable of Part1
hClusterOverFindable0 = Load(inFile, f'HM{fdnames[0]}{suffix}/HM{fdnames[0]}{suffix}/after/TPCRatio_after')
hClusterOverFindable0.Draw('')
hClusterOverFindable1 = Load(inFile, f'HM{fdnames[1]}{suffix}/HM{fdnames[1]}{suffix}/after/TPCRatio_after')
hClusterOverFindable1.SetLineColor(EColor.kRed)
hClusterOverFindable1.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hClusterOverFindable0, labels[0])
leg.AddEntry(hClusterOverFindable1, labels[1])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TPC dE/dx of Part1+
hTPCdedx = Load(inFile, f'HM{fdnames[0]}{suffix}/HM{fdnames[0]}{suffix}/after/TPCdedx_after')
hTPCdedx.SetTitle(f'TPCdedx_after {labels[0]}')
hTPCdedx.Draw('colz')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TPC dE/dx of Part1-
hTPCdedx = Load(inFile, f'HM{fdnames[1]}{suffix}/HM{fdnames[1]}{suffix}/after/TPCdedx_after')
hTPCdedx.SetTitle(f'TPCdedx_after {labels[1]}')
hTPCdedx.Draw('colz')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TOF beta of Part1+
hTOFbeta = Load(inFile, f'HM{fdnames[0]}{suffix}/HM{fdnames[0]}{suffix}/after/TOFbeta_after')
hTOFbeta.SetTitle(f'TOFbeta_after {labels[0]}')
hTOFbeta.Draw('colz')
cPart.SetLogz()
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TOF beta of Part1-
hTOFbeta = Load(inFile, f'HM{fdnames[1]}{suffix}/HM{fdnames[1]}{suffix}/after/TOFbeta_after')
hTOFbeta.SetTitle(f'TOFbeta_after {labels[1]}')
hTOFbeta.Draw('colz')
cPart.SetLogz()
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# nSigmaTPC of Part1+
hNSigmaTPC = Load(inFile, f'HM{fdnames[0]}{suffix}/HM{fdnames[0]}{suffix}/after/NSigTPC_after')
hNSigmaTPC.SetTitle(f'NSigTPC_after {labels[0]}')
hNSigmaTPC.Draw('colz')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# nSigmaTPC of Part1-
hNSigmaTPC = Load(inFile, f'HM{fdnames[1]}{suffix}/HM{fdnames[1]}{suffix}/after/NSigTPC_after')
hNSigmaTPC.SetTitle(f'NSigTPC_after {labels[1]}')
hNSigmaTPC.Draw('colz')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TOF beta of Part1+
hNSigTOF = Load(inFile, f'HM{fdnames[0]}{suffix}/HM{fdnames[0]}{suffix}/after/NSigTOF_after')
hNSigTOF.SetTitle(f'NSigTOF_after {labels[0]}')
hNSigTOF.Draw('colz')
cPart.SetLogz()
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TOF beta of Part1-
hNSigTOF = Load(inFile, f'HM{fdnames[1]}{suffix}/HM{fdnames[1]}{suffix}/after/NSigTOF_after')
hNSigTOF.SetTitle(f'NSigTOF_after {labels[1]}')
hNSigTOF.Draw('colz')
cPart.SetLogz()
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf]'))
# End of part 1 quality


# Part2 quality
oFileBaseName = f'QA_{part2}'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'

oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')
oFile = TFile(oFileName, 'recreate')
cPart = TCanvas('cPart2', '', 600, 600)
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf['))

# Pt of Part2
hPt2 = Load(inFile, f'HM{fdnames[2]}{suffix}/HM{fdnames[2]}{suffix}/after/pTDist_after')
hPt2.Draw()
hPt3 = Load(inFile, f'HM{fdnames[3]}{suffix}/HM{fdnames[3]}{suffix}/after/pTDist_after')
hPt3.SetLineColor(EColor.kRed)
hPt3.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hPt2, labels[2])
leg.AddEntry(hPt3, labels[3])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# Eta of Part2
hEta2 = Load(inFile, f'HM{fdnames[2]}{suffix}/HM{fdnames[2]}{suffix}/after/EtaDist_after')
hEta2.Draw()
hEta1 = Load(inFile, f'HM{fdnames[3]}{suffix}/HM{fdnames[3]}{suffix}/after/EtaDist_after')
hEta1.SetLineColor(EColor.kRed)
hEta1.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hEta2, labels[2])
leg.AddEntry(hEta1, labels[3])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# Phi of Part2
hPhi2 = Load(inFile, f'HM{fdnames[2]}{suffix}/HM{fdnames[2]}{suffix}/after/phiDist_after')
hPhi2.Draw()
hPhi3 = Load(inFile, f'HM{fdnames[3]}{suffix}/HM{fdnames[3]}{suffix}/after/phiDist_after')
hPhi3.SetLineColor(EColor.kRed)
hPhi3.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hPhi2, labels[2])
leg.AddEntry(hPhi3, labels[3])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# CrossedRows of Part2
hCrossedRows2 = Load(inFile, f'HM{fdnames[2]}{suffix}/HM{fdnames[2]}{suffix}/after/CrossedRows_after')
hCrossedRows2.Draw('')
hCrossedRows3 = Load(inFile, f'HM{fdnames[3]}{suffix}/HM{fdnames[3]}{suffix}/after/CrossedRows_after')
hCrossedRows3.SetLineColor(EColor.kRed)
hCrossedRows3.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hCrossedRows2, labels[2])
leg.AddEntry(hCrossedRows3, labels[3])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# ClusterOverFindable of Part2
hClusterOverFindable2 = Load(inFile, f'HM{fdnames[2]}{suffix}/HM{fdnames[2]}{suffix}/after/TPCRatio_after')
hClusterOverFindable2.Draw('')
hClusterOverFindable3 = Load(inFile, f'HM{fdnames[3]}{suffix}/HM{fdnames[3]}{suffix}/after/TPCRatio_after')
hClusterOverFindable3.SetLineColor(EColor.kRed)
hClusterOverFindable3.Draw('same')
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(hClusterOverFindable2, labels[2])
leg.AddEntry(hClusterOverFindable3, labels[3])
leg.Draw('same')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TPC dE/dx of Part2+
hTPCdedx = Load(inFile, f'HM{fdnames[2]}{suffix}/HM{fdnames[2]}{suffix}/after/TPCdedx_after')
hTPCdedx.SetTitle(f'TPCdedx_after {labels[2]}')
hTPCdedx.Draw('colz')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TPC dE/dx of Part2-
hTPCdedx = Load(inFile, f'HM{fdnames[3]}{suffix}/HM{fdnames[3]}{suffix}/after/TPCdedx_after')
hTPCdedx.SetTitle(f'TPCdedx_after {labels[3]}')
hTPCdedx.Draw('colz')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TOF beta of Part2+
hTOFbeta = Load(inFile, f'HM{fdnames[2]}{suffix}/HM{fdnames[2]}{suffix}/after/TOFbeta_after')
hTOFbeta.SetTitle(f'TOFbeta_after {labels[2]}')
hTOFbeta.Draw('colz')
cPart.SetLogz()
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TOF beta of Part2-
hTOFbeta = Load(inFile, f'HM{fdnames[3]}{suffix}/HM{fdnames[3]}{suffix}/after/TOFbeta_after')
hTOFbeta.SetTitle(f'TOFbeta_after {labels[3]}')
hTOFbeta.Draw('colz')
cPart.SetLogz()
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# nSigmaTPC of Part2+
hNSigmaTPC = Load(inFile, f'HM{fdnames[2]}{suffix}/HM{fdnames[2]}{suffix}/after/NSigTPC_after')
hNSigmaTPC.SetTitle(f'NSigTPC_after {labels[2]}')
hNSigmaTPC.Draw('colz')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# nSigmaTPC of Part2-
hNSigmaTPC = Load(inFile, f'HM{fdnames[3]}{suffix}/HM{fdnames[3]}{suffix}/after/NSigTPC_after')
hNSigmaTPC.SetTitle(f'NSigTPC_after {labels[3]}')
hNSigmaTPC.Draw('colz')
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TOF beta of Part2+
hNSigTOF = Load(inFile, f'HM{fdnames[2]}{suffix}/HM{fdnames[2]}{suffix}/after/NSigTOF_after')
hNSigTOF.SetTitle(f'NSigTOF_after {labels[2]}')
hNSigTOF.Draw('colz')
cPart.SetLogz()
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

# TOF beta of Part2-
hNSigTOF = Load(inFile, f'HM{fdnames[3]}{suffix}/HM{fdnames[3]}{suffix}/after/NSigTOF_after')
hNSigTOF.SetTitle(f'NSigTOF_after {labels[3]}')
hNSigTOF.Draw('colz')
cPart.SetLogz()
cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf'))

cPart.SaveAs(os.path.join(cfg['odir'], oFileBaseName + '.pdf]'))
