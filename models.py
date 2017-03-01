import numpy as np
import pandas as p 


from pyvar import DummyVarPrior, SimsZhaSVARPrior, BayesianVAR

data_file = '/if/scratch-m1dxc05/AER_revision/baselineEstimation/vardataJune7.txt'
data = p.read_csv(data_file, delim_whitespace=True, 
                  index_col='DATES', parse_dates=True)

yy = ['FFR_SSR', 'IPM', 'UNRATE', 'PPI_FIN', 'BAA_10YMOODY']
m = ['EGON_KUTTNER_NI']
restriction = np.ones((5, 5), dtype=bool)
hyper = np.ones((7,))
svar_model = SimsZhaSVARPrior(data['1987-8':'1992'][yy], hyper, p=12, 
                              restriction=restriction)

svar_model4 = SimsZhaSVARPrior(data['1987-8':'1992'][yy[:-1]], hyper, p=12, 
                               restriction=np.ones((4, 4), dtype=bool))

rf_model = BayesianVAR(DummyVarPrior(ny=5, p=12),  data['1993-01':'2007-06'][yy])

proxy = data['1994':'2007-06']['EGON_KUTTNER_NI']