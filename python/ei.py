import numpy as np
from pyvar import SimsZhaSVARPrior


class ExternalInstrumentsSVARPrior(SimsZhaSVARPrior):
 
    def __init__(self, *args, **kwargs):
 
        self.gamma = kwargs.pop('gamma', None)
        self.rho = kwargs.pop('rho', None)
        self.signu = kwargs.pop('signu', None)
 
        self.parameterization = 'rho'
        if self.signu is not None:
            self.parameterization = 'signu'
           
        super(ExternalInstrumentsSVARPrior, self).__init__(*args, **kwargs)
 
        self.ny = self.ny+1
       
    def rvs(self, size=1, flatten=True):
        x = super(ExternalInstrumentsSVARPrior, self).rvs(size=size, flatten=flatten)
        gamma = self.gamma.rvs(size=size)
        if self.parameterization == 'rho':
            rho = self.rho.rvs(size=size)
        else:
            rho = self.signu.rvs(size=size)
        x0 = np.c_[x, gamma, rho]
       
        if flatten == True:
            return x0
 
    def para_trans(self, *args, **kwargs):
        dtype=kwargs.pop('dtype', float)
        sqrt = np.sqrt
        if dtype==object:
            import sympy
            sqrt = sympy.sqrt
        if len(args) == 1:
            A0, Aplus = super(ExternalInstrumentsSVARPrior, self).para_trans(args[0][:-2], dtype=dtype)
 
            gam = args[0][-2]
            if self.parameterization == 'signu':
                signu = args[0][-1]
                rho = gam**2 / (gam**2 + signu**2)
            else:
                rho = args[0][-1]
 
            A0t = np.zeros((self.ny, self.ny), dtype=dtype)
            A0t[:-1, :][:, :-1] = A0
            A0t[:-1, -1] = -A0[:, 0] / sqrt( (1-rho)/rho)
            A0t[-1, -1] = 1. / (sqrt( (1-rho)/rho)*gam)
 
            Aplust = np.zeros((self.ny*self.p + self.cons, self.ny), dtype=dtype)
            for i in range(self.p):
                ind0 = i*self.n
                ind0t = i*self.ny
                Aplust[ind0t:ind0t+self.n, :][:, :-1] = Aplus[ind0:ind0+self.n, :]
 
            if self.cons == True:
                Aplust[-1, :-1] = Aplus[-1, :]
 
            Aplust[:, -1] = -Aplust[:, 0] / sqrt( (1-rho)/rho)
            return A0t, Aplust
 
        elif len(args) == 2:
            return args[0], args[1]
   
            
 
    def logpdf(self, x):
        A0, Aplus = super(ExternalInstrumentsSVARPrior, self).para_trans(x[:-2])
        sz = super(ExternalInstrumentsSVARPrior, self).logpdf(A0, Aplus)
 
        if self.parameterization == 'rho':
            pdfx = self.rho.logpdf(x[-1])
        else:
            pdfx = self.signu.logpdf(x[-1])
           
        return sz + self.gamma.logpdf(x[-2]) + pdfx
 
