import numpy as np
import pandas as p 


from pyvar import DummyVarPrior, SimsZhaSVARPrior, BayesianVAR

data_file = '../data/CHdata.txt'
data = p.read_csv(data_file, delim_whitespace=True, index_col='DATES', parse_dates=True)



yy = ['EFFR_LW', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY']

restriction = np.ones((5, 5), dtype=bool)
hyper = np.ones((7,))

svar_model = SimsZhaSVARPrior(data['1990':'1992'][yy], hyper, p=12, 
                              restriction=restriction)

svar_model4 = SimsZhaSVARPrior(data['1990':'1992'][yy[:-1]], hyper, p=12, 
                               restriction=np.ones((4, 4), dtype=bool))

svar_model9 = SimsZhaSVARPrior(np.random.rand(100, 9), hyper, p=12, 
                              restriction=np.ones((9, 9), dtype=bool))

rf_model = BayesianVAR(DummyVarPrior(ny=5, p=12),  data['1993-01':'2007-06'][yy])

models = {'4eq': dict(yy=['EFFR', 'LIPM', 'UNRATE', 'LPPI'], 
                      plot_yy=['EFFR', 'LIPM', 'UNRATE', 'LPPI'], 
                      hyper=[ 0.5, 1, 1, 0.5, 0.5, 1],
                      proxy='MHF', 
                      svar_model=svar_model4),
          '5eq': dict(yy=['EFFR', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'],
                      plot_yy=['EFFR', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'],
                      hyper=[0.5, 3, 1, 0.5, 0.5, 1],
                      proxy='MHF', 
                      svar_model=svar_model),
          '5eq_tight': dict(yy=['EFFR', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'],
                      plot_yy=['EFFR', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'],
                      hyper=[0.5, 3, 1, 0.5, 0.5, 1],
                      proxy='MHF', 
                      svar_model=svar_model),

          '5eq_fin': dict(yy=['EFFR', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'],
                      plot_yy=['EFFR', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'],
                      hyper=[0.5, 3, 1, 0.5, 0.5, 1],
                      proxy='MHF', 
                      svar_model=svar_model),


          '9eq': dict(yy=['EFFR', 'LIPM', 'ADS', 'LPPI', 'BAA10YMOODY', 'UNRATE', 'LPCEND', 'CORELNS', 'LSM'],
                      plot_yy=['EFFR', 'LIPM', 'ADS', 'LPPI', 'BAA10YMOODY', 'UNRATE', 'LPCEND', 'CORELNS', 'LSM'], 
                      hyper=[0.5, 3, 1, 0.5, 0.5, 1],
                      proxy='MHF', 
                      svar_model=svar_model),

          '5eq_cholesky': dict(yy=['LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY', 'EFFR'],
                               plot_yy=['EFFR', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'],
                               hyper=[0.5, 3, 1, 0.5, 0.5, 1],
                               proxy='MHF', 
                               svar_model=svar_model),
          '4eq_cholesky_RR': dict(yy=['LIPM', 'UNRATE', 'LPPI', 'MRR'],
                               plot_yy=['MRR', 'LIPM', 'UNRATE', 'LPPI'],
                               hyper=[0.5, 1, 1, 0.5, 0.5, 1],
                               proxy='MHF', 
                               svar_model=svar_model4), 
          '4eq_cholesky_RRCS': dict(yy=['LIPM', 'UNRATE', 'LPPI', 'MRRCS'],
                               plot_yy=['MRRCS', 'LIPM', 'UNRATE', 'LPPI'],
                               hyper=[0.5, 1, 1, 0.5, 0.5, 1],
                               proxy='MHF', 
                               svar_model=svar_model4), 
          '5eq_cholesky_RR': dict(yy=['LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY', 'MRR'],
                               plot_yy=['MRR', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'], 
                               hyper=[0.5, 3, 1, 0.5, 0.5, 1],
                               proxy='MHF', 
                               svar_model=svar_model), 
          '5eq_cholesky_RRCS': dict(yy=['LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY', 'MRRCS'],
                               plot_yy=['MRRCS', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'], 
                               hyper=[0.5, 3, 1, 0.5, 0.5, 1],
                               proxy='MHF', 
                               svar_model=svar_model), 
}


import scipy.stats as rv

class InvGamma(object):
    
    name = 'inv_gamma'

    def __init__(self, a, b):

        self.a = a
        self.b = b

    def logpdf(self, x):
        a = self.a
        b = self.b
        if x < 0:
            return -1000000000000

        lpdf = (np.log(2) - gammaln(b/2) + b/2*np.log(b*a**2/2)
                -(b+1)/2*np.log(x**2) - b*a**2/(2*x**2))
        return lpdf
        
    
    def rvs(self):
        rn = rv.norm.rvs(size=(int(self.b), 1))
        return np.sqrt(self.b*self.a**2 / np.sum(rn**2, 0))
import numpy as np
from numpy.linalg import svd

def nullspace(A,  atol=1e-13, rtol=0):
    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns

