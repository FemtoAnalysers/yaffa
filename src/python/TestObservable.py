# pylint: disable-all
 
from ROOT import TFile, gInterpreter
gInterpreter.ProcessLine('#include "Observable.h"')
from ROOT import Observable

from yaffa import utils

inFile = TFile('~/an/LPi/systematics/fit/RawCF_Data_pT017.root')
hObs = utils.io.Load(inFile, 'p02_13/sgn/hCFrew')
hObs.Draw()

oObs = Observable(hObs)
hObs2 = hObs * 1.5
hObs2.Draw('same')
