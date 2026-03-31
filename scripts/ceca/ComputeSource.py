import sys
import ctypes

from ROOT import TFile, TF1, gSystem, TGraphErrors, RooArgList, RooRealVar, TCanvas, gPad, RooDataHist
from ROOT.RooFit import Range, NormRange, WARNING

# TODO make portable
gSystem.Load('/home/battidan/ph/proj/source3b/ext/yaffa/src/RooFitFunctions_hxx.so')
from ROOT import RooSourceAAA, RooMsgService

msg = RooMsgService.instance()
msg.setGlobalKillBelow(WARNING)


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
            output += f'  {iPar}: {par[0]} = {par[1]:.2f} +/- {par[2]:.2f}  limits: [{par[3]:.2f}, {par[4]:.2f}]  {par[5]}'
        return output

def RooFit(hist, fitrange):
    # Create observables x,y
    radius = RooRealVar("radius", "radius", 0, 15, 'fm')
    rho0 = RooRealVar("rho0", "rho0", 1, 0, 15, 'fm')
    radius.setRange("fitRange", *fitrange)

    pdfSource = RooSourceAAA("aa", "gaussian PDF", radius, rho0)

    hist.Scale(1. / hist.Integral())
    if abs(hist.Integral() - 1) > 1.e-6:
        print("\033{31mERROR: fit function is not normalized!")
        print("  --> Integral: ", hist.Integral())

    dh = RooDataHist("dh", "dh", RooArgList(radius), Import=hist)

    # TODO: clarify which error estimation is better
    result = pdfSource.fitTo(dh, Range("fitRange"), Save=True, PrintLevel=-1, SumW2Error=True) 
    # result.Print()
    # pdfSource.plotOn(frame)
    # dh.plotOn(frame)

    # # Draw plot on canvas
    # c = TCanvas("test", "test", 600, 600)
    # gPad.SetLeftMargin(0.15)
    # frame.GetYaxis().SetTitleOffset(1.6)
    # frame.Draw()
    # c.SaveAs("test.png")
    pdfSource._fit_result = result
    pdfSource._radius = radius
    pdfSource._rho0   = rho0
    pdfSource._data   = dh
    pdfSource._hist   = hist
    pdfSource._fitres = result
    
    return pdfSource


def main(cfg):
    # Open input file
    inFile = TFile(cfg['infile'])
    hRhoVsMt = inFile.Get('hRhoVsMt')
    hRhoVsMt.SetDirectory(0)
    inFile.Close()
    
    oFile = TFile(cfg['ofile'], 'recreate')

    # mT scaling for 3B
    hRho = Chain(SliceVertically(hRhoVsMt, cfg['mt_bins'], name='hRho_mT'))
    hRho.do(RooFit, cfg['fitrange'], id='source').do(lambda h : h.Write())
   
    print(f'Output saved to \'{cfg["ofile"]}\'')

    # def RooDraw(pdf, name):
    #     frame = radius.frame()
        # data.plotOn(frame)
        
    
    results = hRho.results['source']

    radius = RooRealVar("radius", "radius", *cfg['framerange'], 'fm')
    radius.setRange('fitrange', *cfg['fitrange'])

    for iRes, res in enumerate(results):
        print(res, type(res))

        c = TCanvas(f"aaa{iRes}", f"aaa{iRes}", 600, 600)
        gPad.SetLeftMargin(0.15)
        frame = radius.frame()
        res.plotOn(frame, Range('fitrange'), NormRange("fitRange"))
        res._data.plotOn(frame)
        frame.GetYaxis().SetTitleOffset(1.6)
        frame.Draw()
        c.SaveAs(f"aaa{iRes}.png")
        f = res.asTF(RooArgList(res._radius))
        

        npars = res._fitres.floatParsFinal().getSize()
        chi2 = frame.chiSquare(npars)

        print(chi2, npars)

        
        f.Write()
        
        
        # print(res) # Segmentation violation
        # pass
    
    # results.do(lambda pdf : RooDraw(pdf, "aaa"))
        
    oFile.Close()

    # # Draw plot on canvas

def Fit(obj, func, range, pars):
    '''
    Create fit function
    '''
    
    fFit = TF1(f'fFit', func, *range, len(pars))
    fFit.SetNpx(10000)
    for iPar, par in enumerate(pars):
        fFit.SetParName(iPar, par[0])
        fFit.SetParameter(iPar, (par[2] + par[3]) / 2)
        fFit.SetParLimits(iPar, par[2], par[3])

    result = obj.Fit(fFit, 'QWLLS')
    return FitResult(obj.GetName(), result)


def _main(cfg:dict):
    '''
    Compute the source size as a function of mT
    '''
    # Open input file
    inFile = TFile(cfg['infile'])
    hRhoVsMt = inFile.Get('hRhoVsMt')
    hRhoVsMt.SetDirectory(0)
    inFile.Close()
    
    oFile = TFile(cfg['ofile'], 'recreate')

    # mT scaling for 3B
    chRhos = Chain(SliceVertically(hRhoVsMt, cfg['mt_bins'], name='hRho_mT'))
     
    chRhos.do(lambda h : h.Scale(1. / h.Integral() / h.GetBinWidth(1))) \
        .do(Fit, SourceAAA, [0, 20], [['rho0', 1, 0, 10]], id='results') \
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
