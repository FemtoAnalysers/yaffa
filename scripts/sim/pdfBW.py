'''
Study the change in shape of a BW when transformed from mass to k*.
'''

from ROOT import TF1, TH1D, TCanvas, gROOT, TLine
gROOT.SetBatch(True)


def kaellen(mass1, mass2, m3):
    '''Triangular Kaellen function'''
    return mass1**2 + mass2**2 + m3**2 - 2*mass1*mass2 - 2*mass2*m3 - 2*mass1*m3


def BreakUpMomentum(M, mass1, mass2):
    '''Compute the break-up momentum p of the decay P(M, 0) -> p1(E1, p) + p2(E2, -p)'''
    k = kaellen(M**2, mass1**2, mass2**2)
    k *= k > 0
    return k**0.5/2/M

def tf1BreakupMomentum(x, pars):
    '''
    Wrapper for ROOT's TF1.
    '''

    M = x[0]
    mass1 = pars[0]
    mass2 = pars[1]
    return BreakUpMomentum(M, mass1, mass2)

hM = TH1D('hM', ';M;Counts', 1000, 1200, 1500)
hKStar = TH1D('hKStar', ';k*;Counts', 1000, 0, 350)

n = 100000
mu = 1385

gamma = 16
m1 = 138
m2 = 1115
fM = TF1('fM', f'{n}/((x-{mu})^2 + ({gamma}/2)^2) * (x > {m1 + m2})', 1000,  2000)
for i in range(n):
    mass = fM.GetRandom()
    hM.Fill(mass)
    hKStar.Fill(BreakUpMomentum(mass, m1, m2))

c = TCanvas('c', '', 1200, 600)
c.Divide(2, 1)
pad1 = c.cd(1)
# pad1.SetLogy()
hM.Draw()
fM.Draw('same')
l2Mu = TLine(mu, 0, mu, n)
l2Mup = TLine(mu+gamma, 0, mu+gamma, n)
l2Mum = TLine(mu-gamma, 0, mu-gamma, n)
l2Mu.Draw('same')
l2Mup.Draw('same')
l2Mum.Draw('same')

pad2 = c.cd(2)
# pad2.SetLogy()
hKStar.Draw()
lMu = TLine(BreakUpMomentum(mu, m1, m2), 0, BreakUpMomentum(mu, m1, m2), n)
lMup = TLine(BreakUpMomentum(mu+gamma, m1, m2), 0, BreakUpMomentum(mu+gamma, m1, m2), n)
lMum = TLine(BreakUpMomentum(mu-gamma, m1, m2), 0, BreakUpMomentum(mu-gamma, m1, m2), n)

l2Mup = TLine(BreakUpMomentum(mu, m1, m2)+gamma, 0, BreakUpMomentum(mu, m1, m2)+gamma, n)
l2Mum = TLine(BreakUpMomentum(mu, m1, m2)-gamma, 0, BreakUpMomentum(mu, m1, m2)-gamma, n)
l2Mup.SetLineColor(2)
l2Mum.SetLineColor(2)

lMu.Draw("same")
lMup.Draw("same")
lMum.Draw("same")
l2Mup.Draw('same')
l2Mum.Draw('same')
c.SaveAs("fig.png")
c.SaveAs("fig.pdf")


###
fKStarVsMass = TF1('fKStarVsMass', tf1BreakupMomentum, m1 + m2 -100, 1500, 2)
fKStarVsMass.SetParameter(0, m1)
fKStarVsMass.SetParameter(1, m2)
fKStarVsMass0 = TF1('fKStarVsMass0', tf1BreakupMomentum, m1 + m2 -100, 1500, 2)
fKStarVsMass0.SetParameter(0, 0)
fKStarVsMass0.SetParameter(1, 0)
cCorr = TCanvas('cCorr', '', 600, 600)
cCorr.DrawFrame(1100, 0, 2000, 1000)

fKStarVsMass.Draw('same')
fKStarVsMass0.Draw('same')
cCorr.SaveAs('cCorr.png')
