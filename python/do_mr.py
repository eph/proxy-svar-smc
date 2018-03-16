import numpy as np 
import pandas as p

from pyvar import DiffusePrior, DummyVarPrior, BayesianVAR, MinnesotaPrior

data_file = '/if/scratch-m1dxc05/AER_revision/baselineEstimation/vardataJune7.txt'
data = p.read_csv(data_file, delim_whitespace=True, 
                  index_col='DATES', parse_dates=True)

yy = ['FFR_SSR', 'IPM', 'UNRATE', 'PPI_FIN', 'BAA_10YMOODY']

presample = data['1987-8':'1993'][yy]

tau, d,  w, lam, mu, root = 0.2, 3.5, 1, 1.0, 1.0, 1
hyper = [tau, d, w, lam, mu, root]
lags = 3
prior = MinnesotaPrior(presample, hyper, p=lags)
prior.ny = len(yy)
bvar = BayesianVAR(prior, data['1993':'2007-06'][yy])

phis, sigmas  = bvar.sample(nsim=5000)

phi = phis.mean(0)
sigma = sigmas.mean(0)

import sys
sys.path.append('/mq/DSGE/research/MPpremia/publicCodes/python/')

from mr import mr_identify 
M = np.atleast_2d(data['1993':'2007-06']['EGON_KUTTNER_NI'].values[lags:]).T
x, y,  z = mr_identify(phi, sigma, M, bvar)

print x

phi, sigma = bvar.mle()
x, y,  z = mr_identify(phi, sigma, M, bvar)

print x