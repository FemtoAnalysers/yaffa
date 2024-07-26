'''
Transformation of random variables.
'''

from ROOT import TF1, TH1D, TCanvas



def kaellen(m1, m2, m3):
    '''Triangular Kaellen function'''
    return m1**2 + m2**2 + m3**2 - 2*m1*m2 - 2*m2*m3 - 2*m1*m3


def BreakUpMomentum(M, m1, m2):
    '''Compute the break-up momentum p of the decay P(M, 0) -> p1(E1, p) + p2(E2, -p)'''
    return kaellen(M**2, m1**2, m2**2)**0.5/2/M



hM = TH1D('hM', ';M;Counts', 1000, 1200, 1500)
hKStar = TH1D('hKStar', ';k*;Counts', 1000, 0, 1600)

n = 10000
mu = 1385
gamma = 16
fM = TF1('fM', f'{n}/((x-{mu})^2 + ({gamma}/2)^2)', 1200,  1600)
for i in range(n):
    x = fM.GetRandom()
    print(x)
    hM.Fill(x)
    hKStar.Fill(BreakUpMomentum(x, 138, 938))

c = TCanvas('c', '', 1200, 600)
c.Divide(2, 1)
c.cd(1)
hM.Draw()
fM.Draw('same')
c.cd(2)
hKStar.Draw()
# fY.Draw('same')

# c2 = TCanvas('c2', '', 1200, 600)
# c2.Divide(2, 1)
# c2.cd(1)
# hKStar.Draw()
# fX.Draw('same')
# c2.cd(2)
# hZ.Draw()
# fZ.Draw('same')
