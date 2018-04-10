import numpy as np
import pandas as p

from pyvar import MinnesotaPrior, BayesianVAR


data_file = 'varData.txt'
data = p.read_csv(data_file, delim_whitespace=True, index_col='DATES', parse_dates=True)

yy = ['FFR_SSR', 'IPM', 'UNRATE', 'PCPI', 'BAA_10YMOODY']
presample = data['1990':'1993'][yy]
proxy = data['1994':'2007-06']['EGON_KUTTNER_NI']

signu = 'fixed'
output_dir = '_fortress_tmp'


#------------------------------------------------------------
# Prior as Posterior given VAR data
if 'BAA_10YMOODY' in yy:
    tau, d,  w, lam, mu, root = 0.5, 3, 1, 0.5, 0.5, 1

hyper = [tau, d, w, lam, mu, root]

prior = MinnesotaPrior(presample, hyper, p=2)
prior.ny = len(yy)
bvar = BayesianVAR(prior, data['1993-11':'2007-06'][yy])


yest, xest = bvar.yest, bvar.xest
ydum, xdum = prior.get_pseudo_obs()
yest = np.vstack((ydum, yest))
xest = np.vstack((xdum, xest))

phihatT = np.linalg.solve(xest.T.dot(xest), xest.T.dot(yest))
S = (yest - xest.dot(phihatT)).T.dot(yest - xest.dot(phihatT))
nu = yest.shape[0] - bvar._p * bvar._ny - bvar._cons*1
omega = xest.T.dot(xest)
muphi = phihatT

import pystan

stan_model = pystan.StanModel('proxysvar.stan')

var_data = dict(n=prior.ny,
		p=prior.p,
		cons=int(prior.cons),
		XtXinv=np.linalg.inv(xest.T @ xest),
		phihat=np.ravel(phihatT,order='F'),
		S=S,
		df=nu, 
		y=bvar.yest,
		x=bvar.xest,
		T=bvar.yest.shape[0],
                mm=proxy.values,
                sigmanu=proxy.std())

fit = stan_model.sampling(data=var_data, iter=2000)

print(fit)
