'''
Transformation of random bariables.
'''

import numpy as np

from ROOT import TF1, TH1D, TCanvas, gRandom

hX = TH1D('hX', ';x;Counts', 100, 0, 2*np.pi)
hY = TH1D('hY', ';cos(x);Counts', 100, -1, 1)
hZ = TH1D('hZ', ';cos(x);Counts', 100, 0, 4)
n = 10000
for i in range(n):
    x = gRandom.Uniform(0, 2* np.pi)
    hX.Fill(x)
    hY.Fill(np.cos(x))
    hZ.Fill(x*x/np.pi/np.pi)


fX = TF1('fX', f'{n}*{hX.GetXaxis().GetBinWidth(1)}/2/{np.pi}', 0, 2*np.pi)
fY = TF1('fY', f'{n}*{hY.GetXaxis().GetBinWidth(1)}/{np.pi}/(1 - x*x)^0.5', -1, 1)
fZ = TF1('fZ', f'{n}*{hY.GetXaxis().GetBinWidth(1)}*2/(4*x^0.5)', 0, 4)

c = TCanvas('c', '', 1200, 600)
c.Divide(2, 1)
c.cd(1)
hX.Draw()
fX.Draw('same')
c.cd(2)
hY.Draw()
fY.Draw('same')

c2 = TCanvas('c2', '', 1200, 600)
c2.Divide(2, 1)
c2.cd(1)
hX.Draw()
fX.Draw('same')
c2.cd(2)
hZ.Draw()
fZ.Draw('same')
