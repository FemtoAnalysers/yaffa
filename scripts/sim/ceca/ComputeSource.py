import os
import sys
import ctypes
from dotenv import load_dotenv
from pathlib import Path
import numpy as np

env_path = Path(__file__).resolve().parent.parent.parent.parent / ".env"
if not load_dotenv(dotenv_path=env_path):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")

from ROOT import TFile, gInterpreter, TF1, gROOT, TGraphErrors
gInterpreter.Declare(f'#include "{YAFFA_PATH}/src/cpp/RootFunctions.hxx"')
from ROOT import SourceAAA, SourceCountsAAA, SourceCountsGauss

sys.path.append(f'{YAFFA_PATH}/src/python')

from yaffa import utils

# sequence, film/frame, seggiovia, portico
class Chain:
    def __init__(self, chain):
        self._chain = chain
        self.results = {}

    def do(self, func, *args, id=None):
        if id and self.results.get(id):
            raise RuntimeError('results with id="{}" already exists'.format(id))
        
        results = []
        for ring in self._chain:
            results.append(func(ring, *args))

        if id:
            self.results[id] = results

        return self

class FitResult():
    def __init__(self, name, result):
        self.name = name
        self.result = result.Get()
        self.minimizer = None
        self.chi2 = None
        self.ndf = None
        self.pars = None
        self.ncalls = None
        
        self.ready = False
        
    def __compile(self):
        self.minimizer = self.result.MinimizerType()
        if self.result.Status() == 0:
            self.status = '\033[32mOK\033[0m'
        else:
            self.status = '\033[31mFAIL: UNKNOWN REASON\033[0m'

        self.chi2 = self.result.Chi2()
        self.ndf = self.result.Ndf()
        self.ncalls = self.result.NCalls()

        
        if self.chi2 / self.ndf > 5:
            self.status = '\033[31mFAIL: LARGE CHI2/NDF\033[0m'
            
        self.pars = []
        for iPar in range(self.result.NPar()):
            par = self.result.Parameter(iPar)
            unc = self.result.ParError(iPar)
            name = self.result.ParName(iPar)
            
            # Check if parameter is at limit
            lower = ctypes.c_double()
            upper = ctypes.c_double()
            self.result.ParameterBounds(iPar, lower, upper)
            lower = lower.value
            upper = upper.value
            
            if par - lower < (upper - lower) * 1.e-5:
                status = '\033[33m[[AT LOWER LIMIT]]\033[0m'
                self.status = '\033[31mFAIL: PARAMETERS AT LIMIT\033[0m'
            elif upper - par < (upper - lower) * 1.e-5:
                status = '\033[33m[[AT UPPER LIMIT]]\033[0m'
                self.status = '\033[31mFAIL: PARAMETERS AT LIMIT\033[0m'
            else:
                status = ''
            
            self.pars.append([name, par, unc, lower, upper, status])
        self.ready = True

    def __str__(self):
        if not self.ready:
            self.__compile()

        output = f"\nFit to {self.name} with {self.minimizer}: status={self.status} ({self.ncalls} calls) chi2/ndf={self.chi2:.0f}/{self.ndf:.0f}\n"
        for iPar, par in enumerate(self.pars):
            output += f'  {iPar}: {par[0]:<{8}} = {par[1]:>{8}.2e} +/- {par[2]:<{8}.2e}  limits: [{par[3]:.2e}, {par[4]:.2e}]  {par[5]}\n'
        return output


def Fit(obj, func, range, pars, options=''):
    '''
    Create fit function
    '''
    
    fFit = TF1(f'fFit', func, *range, len(pars))
    fFit.SetNpx(10000)
    for iPar, par in enumerate(pars):
        fFit.SetParName(iPar, par[0])
        fFit.SetParameter(iPar, (par[1] + par[2]) / 2)
        fFit.SetParLimits(iPar, par[1], par[2])

    result = FitResult(obj.GetName(), obj.Fit(fFit, 'QLS'))
    if 'Q' not in options:
        print(result)

    return result


def Test2B3BConsistency(hRStarInTriplets, hRho):
    # r* of pairs in triplets
    norm = hRStarInTriplets.GetEntries() * hRStarInTriplets.GetBinWidth(1)
    resultRStar = Fit(hRStarInTriplets, SourceCountsGauss, [0, 8], [['norm', norm / 2, norm * 2] ,['r0', 0.1, 5]])
    avgRStar = hRStarInTriplets.GetMean()
    r0 = resultRStar.pars[1][1]
    expectedAvgRStar = r0 * 4 / np.sqrt(np.pi)
    delta = (expectedAvgRStar - avgRStar) / expectedAvgRStar
    print(f"r0 (fit)   = {r0:.2f} fm")
    print(f"<r*> (exp) = {expectedAvgRStar:.2f} fm")
    print(f"<r*>       = {avgRStar:.2f} fm")
    print(f"delta      = {delta * 100:.2f} %")

    # Hyper-radius of triplets
    norm = hRho.GetEntries() * hRho.GetBinWidth(1)
    resultRho = Fit(hRho, SourceCountsAAA, [0, 8], [['norm', norm / 2, norm * 2] ,['rho0', 0.1, 5]])
    avgRho = hRho.GetMean()
    rho0 = resultRho.pars[1][1]
    expectedAvgRho = rho0 * 15 * np.sqrt(np.pi) / 16
    delta = (expectedAvgRho - avgRho) / expectedAvgRho
    print(f"rho0 (fit)  = {rho0:.2f} fm")
    print(f"<rho> (exp) = {expectedAvgRho:.2f} fm")
    print(f"<rho>       = {avgRho:.2f} fm")
    print(f"delta       = {delta * 100:.2f} %")

    # Compatibility between 2B and 3B
    expectedAvgRho = avgRStar * 15 * np.pi / 32
    deviation = (avgRho - expectedAvgRho) / expectedAvgRho
    print(f"\n<rho> (exp) = {expectedAvgRho:.2f} fm")
    print(f"<rho>       = {avgRho:.2f} fm")
    print(f"delta       = {deviation * 100:.2f} %")

def main(cfg:dict):
    '''
    Compute the source size as a function of mT
    '''
    # Open input file
    inFile = TFile(cfg['infile'])
    oFile = TFile(cfg['ofile'], 'recreate')

    # Load histograms from CECA simulation
    hRho = inFile.Get('hPhiVsRho').ProjectionX()
    hRStarInTriplets = inFile.Get('hRStarInTriplets')
    hFemtoRho = inFile.Get('hFemtoPhiVsRho').ProjectionX()
    hFemtoRStarInTriplets = inFile.Get('hFemtoRStarInTriplets')

    hRho.SetName('hRho')
    hFemtoRho.SetName('hFemtoRho')

    Test2B3BConsistency(hRStarInTriplets, hRho)
    Test2B3BConsistency(hFemtoRStarInTriplets, hFemtoRho)

    # # mT scaling for 3B
    # hRhoVsMt = inFile.Get('hRhoVsMt')
    # chRhos = Chain(utils.analysis.SliceVertically(hRhoVsMt, cfg['mt_bins'], name='hRho_mT'))
     
    # chRhos.do(Fit, SourceCountsAAA, [0, 20], [['norm', 1, 100000], ['rho0', 0, 10]], id='results') \
    #     .do(lambda h : h.Write())
    # [print(r) for r in chRhos.results['results']]

    hRho.Write()
    hRStarInTriplets.Write()
    hFemtoRho.Write()
    hFemtoRStarInTriplets.Write()

    oFile.Close()
    
if __name__ == '__main__':
    import argparse
    import yaml
    
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg')
    args = parser.parse_args()

    with open(args.cfg) as file:
        cfg = yaml.safe_load(file)

    main(cfg)
