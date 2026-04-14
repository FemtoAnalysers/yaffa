'''
Script to perform the template fits of the DCA distribution.
'''

import os
import argparse
import yaml

from yaffa import logger as log
from yaffa.utils.io import Load
from yaffa.utils import style

from ROOT import gROOT, TFile, TCanvas # pylint: disable=wrong-import-order

parser = argparse.ArgumentParser()
parser.add_argument('cfg', metavar='text', help='yaml configuration file name')
parser.add_argument('-b', action='store_true', default=False, help='Run un batch mode')
args = parser.parse_args()

gROOT.SetBatch(args.b)
style.SetStyle(title=True)

with open(args.cfg, "r", encoding='UTF-8') as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError:
        log.critical('Yaml file not loaded')

# Create output file
oFileName = os.path.join(cfg['odir'], 'TemplateDCA')
if cfg['suffix']:
    oFileName += cfg['suffix']

ptMins = cfg['pt_mins']
ptMaxs = cfg['pt_maxs']

oFile = TFile(oFileName + '.root', 'recreate')

hDcaVsPt = {}
for name, info in cfg['contrib'].items():
    inFile = TFile(info['file'])
    hDcaVsPt[name] = Load(inFile, info['path'])
    hDcaVsPt[name].SetDirectory(0)
    inFile.Close()

cDca = TCanvas('cDca', '', 600, 600)
cDca.SetTopMargin(0.1)
cDca.SaveAs(oFileName + '_distr.pdf[')
for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):

    hDca = {}
    for contrib in cfg['contrib']:
        firstBin = hDcaVsPt[contrib].GetXaxis().FindBin(ptMin * 1.001)
        lastBin = hDcaVsPt[contrib].GetXaxis().FindBin(ptMax * 0.999)
        hDca[contrib] = hDcaVsPt[contrib].ProjectionY(f'hDca_pT{ptMin:1f}_{ptMax:.1f}', firstBin, lastBin)
        hDca[contrib].Rebin(round(cfg['rebin_dca'] / hDca[contrib].GetXaxis().GetBinWidth(0)))
        hDca[contrib].Sumw2()
        hDca[contrib].Scale(1./hDca[contrib].GetEntries())

    # Draw the template distributions
    title = f'{ptMin:.1f} < it{{p}}_{{T}} < {ptMax:.1f} GeV/#it{{c}};DCA'
    cDca.DrawFrame(-1, 0, 1, hDca['data'].GetMaximum() * 1.5, title)
    for contrib in cfg['contrib']:
        hDca[contrib].Draw('same')
    cDca.SaveAs(oFileName + '_distr.pdf')
cDca.SaveAs(oFileName + '_distr.pdf]')
