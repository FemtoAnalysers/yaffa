import numpy as np
import ctypes

from ROOT import TMath, TF1, gROOT, TLine, TF1
from yaffa import logger as log

gROOT.SetBatch(True)

def Gaus(x, p):
    norm = p[0]
    mean = p[1]
    sigma = p[2]
    return norm / (np.sqrt(2 * np.pi) * sigma) * np.e ** (- 0.5 * (x[0] - mean)**2 / sigma**2)

def Pol0(x, p):
    return p[0]

def Pol1(x, p):
    return Pol0(x, p) + p[1] * x[0]

def Voigt(x, p):
    '''
    Voigtian function

    Parameters
    ----------
    - x: function variable
    - par: function parameters
    '''
    norm = p[0]
    mean = p[1]
    sigma = p[2]
    gamma = p[3]

    return norm * TMath.Voigt(x[0]-mean, sigma, gamma)


# def Pol2(x, p):
#     return Pol1(x, p) + p[2] * pow(x[0], 2)

# def Pol3(x, p):
#     return Pol2(x, p) + p[3] * pow(x[0], 3)

    

class MassFitter:
    def __init__(self, hist, signal, background, range_min, range_max):
        self._hist = hist
        self._signal = signal
        self._background = background
        self._range_min = range_min
        self._range_max = range_max
        self._title = ''
        self._signalWindow = None

        self._hResiduals = None
        self._hSignal = None

        self._fTot = None
        self._fSgn = None
        self._fBkg = None
        self._fPrefit = None
        
        if signal == 'voigt':
            self._sgnFnc = Voigt
            self._nParSgn = 4
        elif signal == 'gaus':
            self._sgnFnc = Gaus
            self._nParSgn = 3
        else:
            log.fatal('Signal function "%s" is not implemented', signal)

        if background == 'pol0':
            self._bkgFnc = Pol0
            self._fPrefit = Pol0
            self._nParBkg = 1
        elif background == 'pol1':
            self._bkgFnc = Pol1
            self._fPrefit = Pol1
            self._nParBkg = 2
        else:
            log.fatal('Background function "%s" is not implemented', background)

    def SetTitle(self, title):
        self._title = title

    def Prefit(self, exclusionWindow):
        
        # Prefit backgorund parameters
        def PrefitFcn(x, p):
            if x[0] > exclusionWindow[0] and x[0] < exclusionWindow[1]:
                TF1.RejectPoint()
            return self._bkgFnc(x, p) # TODO remove hard coded

        self._fPrefit = TF1("fPrefit", PrefitFcn, self._range_min, self._range_max, self._nParBkg)
        self._fPrefit.SetNpx(10000)
        self._hist.Fit(self._fPrefit, 'QLMR0', '')

    def Fit(self):

        # Fit
        def Tot(x, p):
            pBkg = [p[self._nParSgn + i] for i in range(self._nParBkg)]
            return self._sgnFnc(x, p) + self._bkgFnc(x, pBkg)

        self._fTot = TF1("fTot", Tot, self._range_min, self._range_max, self._nParSgn + self._nParBkg)
        self._fTot.SetNpx(10000)

        # TODO remove hard coded parameters
        # Normalization
        self._fTot.SetParameter(0, self._hist.GetMaximum() * 0.004)
        self._fTot.SetParLimits(0, self._hist.GetMaximum() * 0.004 / 2, 0.004 *self._hist.GetMaximum() * 2)

        # Mean
        self._fTot.SetParameter(1, 1.115)
        self._fTot.SetParLimits(1, 1.114, 1.116)

        # Sigma
        self._fTot.SetParameter(2, 0.001)
        self._fTot.SetParLimits(2, 0., 0.003)

        if self._signal == 'voigt':
            self._fTot.SetParameter(3, 0.001)
            self._fTot.SetParLimits(3, 0.0004, 0.002)

        # Background parameters
        for iPar in range(self._nParBkg):
            central = self._fPrefit.GetParameter(iPar)
            lower = self._fPrefit.GetParameter(iPar) / 10
            upper = self._fPrefit.GetParameter(iPar) * 10
            
            self._fTot.SetParameter(self._nParSgn + iPar, central)
            self._fTot.SetParLimits(self._nParSgn + iPar, lower, upper)

        self._hist.Fit(self._fTot, 'QLMR0', '')


        for iPar in range(self._fTot.GetNpar()):
            lowerLim = ctypes.c_double()
            upperLim = ctypes.c_double()
            self._fTot.GetParLimits(iPar, lowerLim, upperLim)
            par = self._fTot.GetParameter(iPar)
            
            lLim = lowerLim.value
            uLim = upperLim.value
            if (par - lLim) / (uLim - lLim) < 1.e-4:
                log.warning('Parameter n. %d hits the lower limit: %.4f', iPar, lLim)
            if (uLim - par) / (uLim - lLim) < 1.e-4:
                log.warning('Parameter n. %d hits the upper limit: %.4f', iPar, lLim)

        # Create Signal function
        self._fSgn = TF1('fSgn', self._sgnFnc, self._range_min, self._range_max, self._nParSgn)
        self._fSgn.SetNpx(10000)
        for iPar in range(self._nParSgn):
            self._fSgn.FixParameter(iPar, self._fTot.GetParameter(iPar))

        # Create Background function
        self._fBkg = TF1('fBkg', self._bkgFnc, self._range_min, self._range_max, self._nParBkg)
        self._fBkg.SetNpx(10000)
        for iPar in range(self._nParBkg):
            self._fBkg.FixParameter(iPar, self._fTot.GetParameter(self._nParSgn + iPar))

    def SetSignalWindow(self, minMass, maxMass):
        self._signalWindow = [minMass, maxMass]


    def Draw(self, pad):
        pad.cd()
        pad.SetLogy()
        frame = pad.DrawFrame(self._range_min, 0.5, self._range_max, self._hist.GetMaximum()*1.3, self._title)
        frame.SetName(f'hframe_{id(self)}')  # make unique

        self._hist.DrawCopy('pesame')
        self._fTot.DrawCopy('same')
        self._fSgn.SetLineColor(3)
        self._fSgn.DrawCopy('same')
        self._fBkg.SetLineColor(4)
        self._fBkg.DrawCopy('same')
        
        self._fPrefit.SetLineColor(6)
        self._fPrefit.DrawCopy('same')
        if self._signalWindow:
            tl = TLine()
            tl.SetLineStyle(7)
            tl.DrawLine(self._signalWindow[0], 0.5, self._signalWindow[0], self._hist.GetMaximum())
            tl.DrawLine(self._signalWindow[1], 0.5, self._signalWindow[1], self._hist.GetMaximum())
                
    def DrawResiduals(self, pad):
        self._fBkg
        pad.cd()

        self._hResiduals = self._hist.Clone('hResiduals')
        self._hResiduals.Reset()

        for iBin in range(self._hist.GetNbinsX()):
            bc = self._hist.GetBinContent(iBin + 1)
            mass = self._hist.GetBinCenter(iBin + 1)
            bkg = self._fTot.Eval(mass)
            self._hResiduals.SetBinContent(iBin + 1, bc - bkg)
        
        pad.DrawFrame(self._range_min, self._hResiduals.GetMinimum() * 1.1, self._range_max, self._hResiduals.GetMaximum()*1.3, self._title)
        self._hResiduals.DrawCopy('same')
        if self._signalWindow:
            tl = TLine()
            tl.SetLineStyle(7)
            tl.DrawLine(self._signalWindow[0], self._hResiduals.GetMinimum() * 1.1, self._signalWindow[0], self._hResiduals.GetMaximum()*1.3)
            tl.DrawLine(self._signalWindow[1], self._hResiduals.GetMinimum() * 1.1, self._signalWindow[1], self._hResiduals.GetMaximum()*1.3)

        
    def DrawSignal(self, pad):
        self._fBkg
        pad.cd()

        self._hSignal = self._hist.Clone('hResiduals')
        self._hSignal.Reset()

        for iBin in range(self._hist.GetNbinsX()):
            bc = self._hist.GetBinContent(iBin + 1)
            mass = self._hist.GetBinCenter(iBin + 1)
            bkg = self._fBkg.Eval(mass)
            self._hSignal.SetBinContent(iBin + 1, bc - bkg)
        
        pad.DrawFrame(self._range_min, self._hSignal.GetMinimum() * 1.1, self._range_max, self._hSignal.GetMaximum()*1.3, self._title)
        self._hSignal.DrawCopy('same')
        if self._signalWindow:
            tl = TLine()
            tl.SetLineStyle(7)
            tl.DrawLine(self._signalWindow[0], self._hSignal.GetMinimum() * 1.1, self._signalWindow[0], self._hSignal.GetMaximum()*1.3)
            tl.DrawLine(self._signalWindow[1], self._hSignal.GetMinimum() * 1.1, self._signalWindow[1], self._hSignal.GetMaximum()*1.3)


    def GetSignal(self):
        firstBin = self._hSignal.FindBin(self._signalWindow[0] * 1.0001)
        lastBin = self._hSignal.FindBin(self._signalWindow[1] * 0.9999)
        integral = self._hSignal.Integral(firstBin, lastBin)
        error = np.sqrt(integral) # TODO think better way
        return integral, error

    def GetBackground(self):
        integral = self._fBkg.Integral(self._signalWindow[0] * 1.0001, self._signalWindow[1] * 0.9999)
        if integral <= 1:
            return 1, 1

        error = np.sqrt(integral) # TODO think better way
        return integral, error

    def GetFitFunc(self):
        return self._fTot
    
    def GetSgnFunc(self):
        if self._fSgn:
            return self._fSgn
    
        log.error('Signal function is None. Did you do the fit?')

    def GetBkgFunc(self):
        if self._fBkg:
            return self._fBkg
    
        log.error('Background function is None. Did you do the fit?')
