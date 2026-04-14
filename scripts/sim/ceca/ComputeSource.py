import sys
import ctypes

from ROOT import TFile, gInterpreter, TF1, gROOT, TGraphErrors
gInterpreter.Declare('#include "../../src/RootFunctions.hxx"')
from ROOT import SourceAAA, SourceCountsAAA

sys.path.append('../../src')

from utils import SliceVertically

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
        self.chi2 = None
        self.ndf = None
        self.pars = None
        self.ncalls = None
        
        self.ready = False
        
    def __compile(self):
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

        output = f"\nFit to {self.name}: status={self.status} ({self.ncalls} calls) chi2/ndf={self.chi2:.0f}/{self.ndf:.0f}\n"
        for iPar, par in enumerate(self.pars):
            output += f'  {iPar}: {par[0]} = {par[1]:>{8}.2f} +/- {par[2]:<{8}.2f}  limits: [{par[3]:.2f}, {par[4]:.2f}]  {par[5]}\n'
        return output


def Fit(obj, func, range, pars):
    '''
    Create fit function
    '''
    
    fFit = TF1(f'fFit', func, *range, len(pars))
    fFit.SetNpx(10000)
    for iPar, par in enumerate(pars):
        fFit.SetParName(iPar, par[0])
        fFit.SetParameter(iPar, (par[1] + par[2]) / 2)
        fFit.SetParLimits(iPar, par[1], par[2])

    result = obj.Fit(fFit, 'QWLLS')
    return FitResult(obj.GetName(), result)


def main(cfg:dict):
    '''
    Compute the source size as a function of mT
    '''
    # Open input file
    inFile = TFile(cfg['infile'])
    oFile = TFile(cfg['ofile'], 'recreate')

    # mT scaling for 3B
    hRhoVsMt = inFile.Get('hRhoVsMt')
    chRhos = Chain(SliceVertically(hRhoVsMt, cfg['mt_bins'], name='hRho_mT'))
     
    chRhos.do(Fit, SourceCountsAAA, [0, 20], [['norm', 1, 100000], ['rho0', 0, 10]], id='results') \
        .do(lambda h : h.Write())
    [print(r) for r in chRhos.results['results']]

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
