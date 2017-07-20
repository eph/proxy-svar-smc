import numpy as np
import pandas as p

from scipy.stats import norm, uniform

from pyvar import MinnesotaPrior, SimsZhaSVARPrior, BayesianVAR
from fortress import make_smc

import sympy 

data_file = 'varData.txt'
data = p.read_csv(data_file, delim_whitespace=True, index_col='DATES', parse_dates=True)

yy = ['FFR_SSR', 'IPM', 'UNRATE', 'PPI_FIN', 'BAA_10YMOODY']
presample = data['1990':'1993'][yy]
proxy = data['1994':'2007-06']['EGON_KUTTNER_NI']

signu = 'fixed'
output_dir = '_fortress_tmp'


#------------------------------------------------------------
# Prior as Posterior given VAR data
if 'BAA_10YMOODY' in yy:
    tau, d,  w, lam, mu, root = 0.5, 3, 1, 0.5, 0.5, 1

hyper = [tau, d, w, lam, mu, root]

prior = MinnesotaPrior(presample, hyper, p=12)
prior.ny = len(yy)
bvar = BayesianVAR(prior, data['1993':'2007-06'][yy])


yest, xest = bvar.yest, bvar.xest
ydum, xdum = prior.get_pseudo_obs()
yest = np.vstack((ydum, yest))
xest = np.vstack((xdum, xest))

phihatT = np.linalg.solve(xest.T.dot(xest), xest.T.dot(yest))
S = (yest - xest.dot(phihatT)).T.dot(yest - xest.dot(phihatT))
nu = yest.shape[0] - bvar._p * bvar._ny - bvar._cons*1
omega = xest.T.dot(xest)
muphi = phihatT

np.savetxt('_fortress_tmp/phistar.txt', phihatT.flatten(order='F'))
np.savetxt('_fortress_tmp/iw_Psi.txt', S)
np.savetxt('_fortress_tmp/Omega_inv.txt', np.linalg.inv(omega))

ndraws = 9600
#phis, sigmas = prior.rvs(size=ndraws, flatten_output=False)
phis, sigmas = bvar.sample(nsim=ndraws, flatten_output=False)

npara = phis[0].size + sigmas[0].size + 1
priorsim = np.zeros((ndraws, npara))

#from dsge.OtherPriors import InvGamma

from tqdm import tqdm
from scipy.stats import norm, uniform
from scipy.linalg import qr, cholesky

phis, sigmas = bvar.sample(nsim=ndraws, flatten_output=False)

for i in tqdm(range(ndraws)):

    sigma_tr = cholesky(sigmas[i], lower=True)
    q, r = qr(norm.rvs(size=(bvar._ny, bvar._ny)))

    A0 = np.linalg.inv(sigma_tr).T.dot(q)
    Ap = phis[i].dot(A0)
    beta = 0.1*norm.rvs()
    #sigma = InvGamma(0.02, 2.0).rvs()

    priorsim[i] = np.r_[A0.flatten(order='F'), 
                        Ap.flatten(order='F'), 
                        beta]


np.savetxt('_fortress_tmp/priordraws.txt', priorsim)

#tau, d,  w, lam, mu, root = 0.5, 1, 1, 0.5, 0.5, 1

#------------------------------------------------------------ 
ny = len(yy)
restriction = np.ones((ny,ny), dtype=bool)
hyper = np.ones((7,))
ei = SimsZhaSVARPrior(data['1990':'1993'][yy], hyper, p=12, restriction=restriction)
ei.nA = restriction.shape[0]**2
ei.nF = ei.ny**2*ei.p + ei.ny
ei.name = 'ei'
ei.T = data['1993':'2007-06'][yy].shape[0]
ei.phistar_file = '_fortress_tmp/phistar.txt'
ei.iw_Psi_file = '_fortress_tmp/iw_Psi.txt'
ei.Omega_inv_file = '_fortress_tmp/Omega_inv.txt'
ei.nu = nu


if signu is 'fixed':
    ei.signu = '0.04_wp'
    ei.nextra_para = 1
    
x = ei.rvs()
x = sympy.symbols(['para({:d})'.format(i+1) for i in range(x.size)], positive=True)
a0, aplus = ei.para_trans(x, dtype=object)
 
mat_str = lambda *x: '    {}({}, {}) = {}'.format(*x)
A0str = [mat_str('A', i+1, j+1, sympy.fcode(value, source_format='free'))
         for (i, j), value in np.ndenumerate(a0) if value > 0]
Apstr = [mat_str('F', i+1, j+1, sympy.fcode(value, source_format='free'))
         for (i, j), value in np.ndenumerate(aplus) if value > 0]
 


varfile = open('svar.f90', 'r').read()
np.savetxt('data.txt',data['1993':'2007-06'][yy])
np.savetxt('proxy.txt',proxy)
varfile = varfile.format(datafile='_fortress_tmp/data.txt', proxyfile='_fortress_tmp/proxy.txt',
                         assign_para = '\n'.join(A0str + Apstr), **vars(ei))

smc = make_smc(varfile, other_files={'data.txt': data['1993':'2007-06'][yy].values,
                                     'proxy.txt': proxy.values})


# Minnesota Prior Hyperparameters





yy = ['FFR_SSR', 'IPM', 'UNRATE', 'PPI_FIN', 'BAA_10YMOODY']
m = ['EGON_KUTTNER_NI']
# ei = ExternalInstrumentsSVARPrior(data['1990':'1993'][yy], hyper, p=12, restriction=restriction,
#                                   gamma=norm, signu=uniform)

#ei.fortran(data['1993':'2007-06'][yy], output_dir=model_dir)
#jfdsklfdsf

#





# import sys

# sys.path.append('/mq/home/m1eph00/projects/dsge-book/code/helper')
# from helper import SMCResults

# smc = SMCResults('svar_proxy-mix', npart=9600, nblocks=25, nphi=2000, lam=2.7)
# parasim = smc.load_draws([0])
# jj = parasim.postsim.argmax()
# A0hat, Aphat = ei.para_trans(parasim[smc['paranames'][:-2]].ix[jj])

# from collections import defaultdict
# irfs = defaultdict(p.DataFrame)
# h = 49
# irfs = p.DataFrame()

# def long_run_coefficients(A0, Ap, equation=0):
#     """
#     Computes the Long Run Coefficients a la Sims-Zha
#     """
#     eqi = equation
#     ny = A0.shape[0]
    
#     alpha = A0[:, eqi]
#     alpha = [np.r_[A0[i, eqi], -Ap[:, eqi][i::ny]] for i in range(ny)]
    
#     delta = -alpha[eqi]
    
#     lr_coeff_diff = np.array([ai.sum()/delta.sum() for ai in alpha])
#     sr_coeff_diff = np.array([ai.sum()/A0[0, 0] for ai in alpha])
#     lr_coeff_lev = np.array([ai.cumsum().sum()/delta.sum() for ai in alpha])

#     return lr_coeff_diff, lr_coeff_lev, sr_coeff_diff


# nsim = 9600

# lrc_diff = np.zeros((nsim, 5))
# lrc_lev = np.zeros((nsim, 5))    
# src_diff = np.zeros((nsim, 5))
# for i in range(9600):
#     A0, Ap = ei.para_trans(parasim[smc['paranames'][:-2]].ix[i])
#     for j in range(5):
#         if np.linalg.inv(A0).dot(A0hat[:, j])[j] > 0:
#             A0[:, j] = -A0[:, j]
#             Ap[:, j] = -Ap[:, j]

#     #irf = ei.structural_irf(A0, Ap, h=h)
#     #irf[0] = irf[0] * np.sign(irf[0][0, 0])

#     #irfs = irfs.append(p.DataFrame(irf[0], columns=yy))

#     lrc_diff[i], lrc_lev[i], src_diff[i] = long_run_coefficients(A0, Ap)


# lrc_diff = p.DataFrame(lrc_diff, columns=yy)
# lrc_lev = p.DataFrame(lrc_lev, columns=yy)
# src_diff = p.DataFrame(src_diff, columns=yy)
# median = irfs.groupby(irfs.index).median()
# q05 = irfs.groupby(irfs.index).quantile(0.05)
# q95 = irfs.groupby(irfs.index).quantile(0.95)


# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(nrows=2, ncols=3)
# axi = ax.reshape(-1)
# for i, name in enumerate(yy):
#     median[name].plot(linewidth=3, ax=axi[i])
#     q05[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')
#     q95[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')

#     axi[i].set_title(name)



    
