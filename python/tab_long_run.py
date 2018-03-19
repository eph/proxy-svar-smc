import numpy as np
import pandas as p
import matplotlib.pyplot as plt

from models import svar_model, yy, svar_model4
from simulations import smc_baseline as smc

from tqdm import tqdm

import sys
sys.path.append('/mq/home/m1eph00/projects/dsge-book/code/helper')
from helper import SMCResults

#models = ['svar12_nocs_tight-mix', 'svar12_nocs_medium-mix', 'svar12_nocs_loose-mix', 
#          'svar12_tight-mix', 'svar12_medium-mix', 'svar12_loose-mix']
#models = ['svar12_sign_cpi_loose-mix']

def long_run_coefficients(A0, Ap, equation=0):
    """
    Computes the Long Run Coefficients a la Sims-Zha
    """
    eqi = equation
    ny = A0.shape[0]

    alpha = A0[:, eqi]
    alpha = [np.r_[-A0[i, eqi], Ap[:-1, eqi][i::ny]] for i in range(ny)]
    
    delta = -alpha[eqi]
    
    lr_coeff_diff = np.array([ai.sum()/delta.sum() for ai in alpha])
    sr_coeff_diff = np.array([ai.sum()/A0[0, 0] for ai in alpha])

    sr_coeff_diff[0] = (Ap[:-1, 0][::ny]/A0[0, 0]).sum()
    x = sr_coeff_diff[0].copy()
    sr_coeff_diff = sr_coeff_diff 

    sr_coeff_lev = np.array([ai.cumsum().sum()/A0[0, 0] for ai in alpha])
    sr_coeff_lev = sr_coeff_lev 

    sr_coeff_diff[0] = x
    return sr_coeff_lev, sr_coeff_diff


results = p.DataFrame()

for model in models:
    smc = SMCResults(model, npart=9600, nblocks=25, nphi=500, lam=2.7)

    parasim = smc.load_draws([0])
    nsim = parasim.shape[0]

    coeffs = np.zeros((nsim, 5))

    four_variable = 'nocs' in model
    for i in tqdm(range(nsim)):
        if four_variable:
            A0, Ap = svar_model4.para_trans(parasim[smc['paranames'][:-2]].ix[i])
        else:
            A0, Ap = svar_model.para_trans(parasim[smc['paranames'][:-2]].ix[i])

        src_lev, src_diff = long_run_coefficients(A0, Ap)

        coeffs[i][0] = src_diff[0]
        coeffs[i][1] = src_lev[1]
        coeffs[i][2] = src_diff[2]
        coeffs[i][3] = src_lev[3]
        if four_variable:
            coeffs[i][4] = 0.0
        else:
            coeffs[i][4] = src_diff[4]

    results_model = p.DataFrame(coeffs, columns=yy)
    results_model['model'] = model

    results = results.append(results_model)

tex_labels = ['r', 'y', 'u', '\\pi', 'cs']

by_model = results.groupby(results.model)
nmodels = len(models)

print('&'.join(['                  '] + ['{:^20s}'.format(x[:-4]) for x in models]))
for series, tex in zip(results.columns, tex_labels):
    
    r = [' $\phi_{{ {:5s} }}$ '] + nmodels*['  {:16.2f}  ']
    print('&'.join(r).format(tex, *by_model.mean().ix[models][series]))
    r = ['                  '] + nmodels*[' ({:16.2f}) ']
    print('&'.join(r).format(*by_model.std().ix[models][series]))
