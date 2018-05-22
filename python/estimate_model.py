import argparse

import numpy as np
import pandas as p

from scipy.stats import norm, uniform

from pyvar import MinnesotaPrior, SimsZhaSVARPrior, BayesianVAR
from fortress import make_smc

import sympy


from models import models as models, InvGamma, data as data


parser = argparse.ArgumentParser(description='Estimation a model with SMC.')
parser.add_argument('--model', help='the model to estimate', action='store',
                    choices=models.keys(), default='5eq')
parser.add_argument('--nsim', help='# of particles', action='store',
                    type=int, default=9600)
parser.add_argument('--nproc', help='# of processors', action='store',
                    type=int, default=4)


args = parser.parse_args()


model = models[args.model]


presample = data['1990':'1993'][model['yy']]
proxy = data['1994':'2007-06'][model['proxy']]

if 'tight' in args.model:
    signu = 'fixed'
    nextra_para = 2
else:
    signu = 'estimated'
    nextra_para = 2


#------------------------------------------------------------
# Prior as Posterior given VAR data
presample_moments = [presample.mean(), presample.std()]


if 'RR' in args.model:
    data['MRR'] = data['MRR'].cumsum()
    data['MRRCS'] = data['MRRCS'].cumsum()
    presample_moments[0][-1] = 0.0
    presample_moments[1][-1] = 0.25

prior = MinnesotaPrior([], model['hyper'], p=12,
                       presample_moments=presample_moments)


prior.ny = len(model['yy'])
bvar = BayesianVAR(prior, data['1993':'2007-06'][model['yy']].fillna(0.0))


#------------------------------------------------------------
yest, xest = bvar.yest, bvar.xest
ydum, xdum = prior.get_pseudo_obs()
yest = np.vstack((ydum, yest))
xest = np.vstack((xdum, xest))

phihatT = np.linalg.solve(xest.T.dot(xest), xest.T.dot(yest))
S = (yest - xest.dot(phihatT)).T.dot(yest - xest.dot(phihatT))
nu = yest.shape[0] - bvar._p * bvar._ny - bvar._cons * 1
omega = xest.T.dot(xest)
muphi = phihatT

ndraws = args.nsim


phis, sigmas = bvar.sample(nsim=10, flatten_output=False)
npara = phis[0].size + sigmas[0].size + nextra_para
priorsim = np.zeros((ndraws, npara))

from tqdm import tqdm
from scipy.stats import norm, uniform
from scipy.linalg import qr, cholesky

phis, sigmas = bvar.sample(nsim=ndraws, flatten_output=False)


for i in tqdm(range(ndraws)):

    sigma_tr = cholesky(sigmas[i], lower=True)

    if 'cholesky' in args.model:
        q = np.eye(bvar._ny)
    else:
        q, r = qr(norm.rvs(size=(bvar._ny, bvar._ny)))

    A0 = np.linalg.inv(sigma_tr).T.dot(q)
    Ap = phis[i].dot(A0)
    beta = 0.1 * norm.rvs()

    u = bvar.yest @ A0 - bvar.xest @ Ap
    eta = proxy - beta * u[:, 0]
    if signu is 'fixed':
        s = 0.5 * proxy.std()
    else:
        s = InvGamma(0.02, 2).rvs()

    priorsim[i] = np.r_[A0.flatten(order='F'),
                        Ap.flatten(order='F'),
                        beta,
                        s]


if 'cholesky' in args.model:
    results = {'var.%03d' %
               (d + 1): list(priorsim[:, d]) for d in range(priorsim.shape[1])}
    results['weights'] = list(np.ones(priorsim.shape[0]) / priorsim.shape[0])
    import json
    with open('%s.json' % args.model, 'w') as outfile:
        json.dump({'posterior.162': results},  outfile)
else:
    ny = len(model['yy'])

    hyper = np.ones((7,))
    restriction = np.ones((ny, ny), dtype=bool)
    # data['1990':'1993'][model['yy']]
    ei = SimsZhaSVARPrior(np.random.rand(100, ny), hyper, p=12,
                          restriction=restriction)
    ei.nA = restriction.shape[0]**2
    ei.nF = ei.ny**2 * ei.p + ei.ny
    ei.name = 'ei'
    ei.T = data['1993':'2007-06'][model['yy']].shape[0]
    ei.nu = nu

    if signu is 'fixed':
        ei.signu = '%f_wp' % (0.5 * proxy.std())
        print(ei.signu)
        ei.nextra_para = 2
    else:
        ei.signu = 'para(self%nA+self%nF+2)'
        ei.nextra_para = 2
    x = sympy.symbols(['para({:d})'.format(i + 1) for i in range(ei.nA + ei.nF)],
                      positive=True)
    a0, aplus = ei.para_trans(x, dtype=object)

    mat_str = lambda *x: '    {}({}, {}) = {}'.format(*x)
    A0str = [mat_str('A', i + 1, j + 1, sympy.fcode(value, source_format='free'))
             for (i, j), value in np.ndenumerate(a0) if value > 0]
    Apstr = [mat_str('F', i + 1, j + 1, sympy.fcode(value, source_format='free'))
             for (i, j), value in np.ndenumerate(aplus) if value > 0]

    varfile = open('svar.f90', 'r').read()
    varfile = varfile.format(assign_para='\n'.join(A0str + Apstr), **vars(ei))

    other_files = {'data.txt': data['1993':'2007-06'][model['yy']].values,
                   'proxy.txt': proxy.values,
                   'priordraws.txt': priorsim[:args.nsim].T,
                   'phistar.txt': phihatT.flatten(order='F'),
                   'iw_Psi.txt': S,
                   'Omega_inv.txt': np.linalg.inv(omega)}

    output_dir = '_fortress_' + args.model
    smc = make_smc(varfile, other_files=other_files,
                   output_directory=output_dir)

    if 'fin' in args.model:
        print("For financial shock, you need to impose the sign restriction manually in the fortran code!")

    smc.run(npart=args.nsim, nblocks=25, nproc=args.nproc, bend=2.7,
            conditional_covariance=True,
            initial_particles=output_dir + '/priordraws.txt',
            output_file='%s.json' % args.model)
