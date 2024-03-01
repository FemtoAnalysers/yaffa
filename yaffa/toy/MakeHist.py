'''
Make root files with example distributions.
'''

import os

from ROOT import TFile, TH1D, TF1

hGaus = TH1D("hGaus", '', 100, -5, 5)
hGaus.FillRandom('gaus')
hGaus.Scale(1. / hGaus.GetEntries())
hGaus.Draw()

fGaus = TF1('fGaus', 'gausn', -5, 5)
fGaus.SetParameter(0, 0.1)
fGaus.SetParameter(1, 0)
fGaus.SetParameter(2, 1)
fGaus.Draw('same')

oFile = TFile(os.path.join(os.getenv('YAFFA'), 'yaffa/toy/Gaus.root'), 'recreate')
hGaus.Write()
fGaus.Write()
oFile.Close()
