import numpy as np

import pandas as p
from pyvar import MinnesotaPrior, BayesianVAR
from ei import ExternalInstrumentsSVARPrior

from scipy.stats import norm, uniform
from scipy.linalg import qr, cholesky

data_file = '/if/scratch-m1dxc05/AER_revision/baselineEstimation/vardataJune7.txt'
data = p.read_csv(data_file, delim_whitespace=True, 
                  index_col='DATES', parse_dates=True)

yy = ['FFR_SSR', 'IPM', 'UNRATE', 'PPI_FIN', 'BAA_10YMOODY']
#yy = ['FFR_SSR', 'IPM','UNRATE','PPI_FIN','BAA_10YMOODY','YLD1000','LMARKET','PCE']

presample = data['1990':'1993'][yy]

proxy = data['1994':'2007-06']['EGON_KUTTNER_NI']

tau, d,  w, lam, mu, root = 0.5, 3, 1, 0.5, 0.5, 1
#tau, d,  w, lam, mu, root = 0.5, 1, 1, 0.5, 0.5, 1
hyper = [tau, d, w, lam, mu, root]

prior = MinnesotaPrior(presample, hyper, p=12)
prior.ny = len(yy)
bvar = BayesianVAR(prior, data['1993':'2007-06'][yy])

print('Log MDD: {:8.4f}'.format(bvar.logmdd()))


yest, xest = bvar.yest, bvar.xest
ydum, xdum = prior.get_pseudo_obs()
yest = np.vstack((ydum, yest))
xest = np.vstack((xdum, xest))

phihatT = np.linalg.solve(xest.T.dot(xest), xest.T.dot(yest))
S = (yest - xest.dot(phihatT)).T.dot(yest - xest.dot(phihatT))
nu = yest.shape[0] - bvar._p * bvar._ny - bvar._cons*1
omega = xest.T.dot(xest)
muphi = phihatT

np.savetxt('/mq/DSGE/research/MPpremia/missPaper/publicCodesRevision/fortran_estimation/fortran_rf_proxy_fin/phihat.txt', phihatT.flatten(order='F'))
np.savetxt('/mq/DSGE/research/MPpremia/missPaper/publicCodesRevision/fortran_estimation/fortran_rf_proxy_fin/Psi.txt', S)
np.savetxt('/mq/DSGE/research/MPpremia/missPaper/publicCodesRevision/fortran_estimation/fortran_rf_proxy_fin/Omega_inv.txt', np.linalg.inv(omega))
print('nu', nu)


ndraws = 1000
#phis, sigmas = prior.rvs(size=ndraws, flatten_output=False)
phis, sigmas = bvar.sample(nsim=ndraws, flatten_output=False)

npara = phis[0].size + sigmas[0].size + 2
priorsim = np.zeros((ndraws, npara))

from dsge.OtherPriors import InvGamma

from tqdm import tqdm

from scipy.stats import norm
pdf = np.zeros((ndraws))
for i in tqdm(range(ndraws)):

    sigma_tr = cholesky(sigmas[i], lower=True)

    q, r = qr(norm.rvs(size=(bvar._ny, bvar._ny)))

    A0 = np.linalg.inv(sigma_tr).T.dot(q)
    Ap = phis[i].dot(A0)
    beta = 0.1*norm.rvs()
    sigma = InvGamma(0.02, 2.0).rvs()

    priorsim[i] = np.r_[A0.flatten(order='F'), 
                        Ap.flatten(order='F'), 
                        beta, sigma]
    eps = (bvar.yest.dot(A0) - bvar.xest.dot(Ap))[:, 0]
    pdf[i] = norm.logpdf(proxy - beta*eps, scale=sigma).sum()
    

np.savetxt('/mq/DSGE/research/MPpremia/missPaper/publicCodesRevision/fortran_estimation/fortran_rf_proxy_fin/model/priordraws.txt', priorsim)
