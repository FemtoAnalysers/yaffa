'''
Toy MC: sample two points uniformly distributed on the unit sphere and look at
the distribution of the angle they subtend at the origin.
'''

import numpy as np

from ROOT import TH1D, TCanvas, TVector3, gRandom, TF1

def RandomPointOnSphere():
    theta = np.arccos(2 * gRandom.Uniform(0, 1) - 1)
    phi = gRandom.Uniform(0, 2 * np.pi)

    return TVector3(np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta))

n = 100000

hAngle = TH1D('hAngle', ';angle between point 1 and point 2 (rad);Counts', 100, 0, np.pi)
fAngle = TF1('fAngle', f'{n * hAngle.GetBinWidth(1)} * 0.5 * sin(x)', 0, np.pi)

for i in range(n):
    point1 = RandomPointOnSphere()
    point2 = RandomPointOnSphere()
    hAngle.Fill(point1.Angle(point2))

c = TCanvas('c', '', 600, 600)
hAngle.Draw()
fAngle.Draw('same')
