'''
Compute the invariant mass of the pairs from the relative momentum k*.
'''

import argparse
import yaml

from ROOT import TFile

from yaffa import utils

parser = argparse.ArgumentParser()
parser.add_argument('cfg')
parser.add_argument('oFileName')
args = parser.parse_args()

with open(args.cfg, 'r') as file:
    cfg = yaml.safe_load(file)

inFileData = TFile(cfg['data']['file'])
hSEData = utils.io.Load(inFileData, cfg['data']['se_path'])
hMEData = utils.io.Load(inFileData, cfg['data']['me_path'])
hSEData = utils.analysis.ChangeUnits(hSEData, 1000)
hMEData = utils.analysis.ChangeUnits(hMEData, 1000)

inFileMC = TFile(cfg['mc']['file'])
hSEMC = utils.io.Load(inFileMC, cfg['mc']['se_path'])
hMEMC = utils.io.Load(inFileMC, cfg['mc']['me_path'])
hSEMC = utils.analysis.ChangeUnits(hSEMC, 1000)
hMEMC = utils.analysis.ChangeUnits(hMEMC, 1000)

firstBin = hSEData.GetXaxis().FindBin(cfg['norm_range'][0] * 1.001)
lastBin = hSEData.GetXaxis().FindBin(cfg['norm_range'][1] * 0.999)
hMEData.Scale(hSEData.Integral(firstBin, lastBin) / hMEData.Integral(firstBin, lastBin))
hMEMC.Scale(hSEMC.Integral(firstBin, lastBin) / hMEMC.Integral(firstBin, lastBin))

oFile = TFile(args.oFileName, 'recreate')

hSEData.Sumw2()
hNum = hSEData - hMEData
hNum.Write()

hDen = hSEMC - hMEMC
hDen.Sumw2()
hDen.Write()

hSEMC.Scale(hSEData.Integral(firstBin, lastBin) / hSEMC.Integral(firstBin, lastBin))
hMEMC.Scale(hSEData.Integral(firstBin, lastBin) / hMEMC.Integral(firstBin, lastBin))
hSEData.Write('hSEData')
hMEData.Write('hMEData')
hSEMC.Write('hSEMC')
hMEMC.Write('hMEMC')

hSEData.Draw()
hMEData.Draw('same')
hSEMC.Draw('same')
hMEMC.Draw('same')

hInvMass = hNum - hDen
hInvMass.Write()

print(f'Output saved in: {args.oFileName}')
