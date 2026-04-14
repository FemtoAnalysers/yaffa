'''
Compute the k* distribution of a wide resonance.
'''
import numpy as np

from ROOT import TF1

def BreitWigner(kstar, M, width, m1, m2):
    '''
    Breit wigner formula.
    '''
    gamma = M * (M**2 + width**2) ** 0.5
    norm = 2 * 2**0.5 * M * width * gamma / np.pi / (M**2 + gamma) ** 0.5
    energy2 = m1**2 + m2**2 + 2 * ((m1**2 + kstar**2)**0.5 + (m2**2 + kstar**2)**0.5 + kstar**2)

    return norm/((energy2 - M**2)**2 + M**2 * width**2)

def fBreitWigner(x, pars):
    '''
    Breit wigner formula (ROOT wrapper).
    '''
    kstar = x[0]
    M = pars[0]
    width = pars[1]
    m1 = pars[2]
    m2 = pars[3]
    return BreitWigner(kstar, M, width, m1, m2)

# def kaellen(m1, m2, m3):
#     '''Triangular Kaellen function'''
#     return m1**2 + m2**2 + m3**2 - 2*m1*m2 - 2*m2*m3 - 2*m1*m3


# def BreakUpMomentum(M, m1, m2):
#     '''Compute the break-up momentum p of the decay P(M, 0) -> p1(E1, p) + p2(E2, -p)'''
#     return kaellen(M**2, m1**2, m2**2)**0.5/2/M

# def fBreakUpMomentum(x, pars):
#     kstar = x[0]
#     M = pars[0]
#     width = pars[1]
#     m1 = pars[2]
#     m2 = pars[3]

#     return BreakUpMomentum(M, m1, m2)

fBW = TF1("fBW", fBreitWigner, 0, 1600, 4)
fBW.SetParameter(0, 1385)
fBW.SetParameter(1, 36)
fBW.SetParameter(2, 140)
fBW.SetParameter(3, 1200)
fBW.SetNpx(500)
fBW.Draw()
