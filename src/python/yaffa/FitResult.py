import ctypes

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