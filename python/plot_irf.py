import argparse

import numpy as np
import pandas as p
import matplotlib.pyplot as plt

import json

from models import svar_model, rf_model, yy,  svar_model4
# from simulations import smc_baseline as smc

h = 49

model = '4eq'
if '4eq' in model: 
    svar_model = svar_model4
    yy = yy[:-1]
if 'cpi' in model:
    yy = ['FFR_SSR', 'IPM', 'UNRATE', 'PCPI', 'BAA_10YMOODY']


results = json.loads(open(model + '.json').read())

parasim = p.DataFrame(results['posterior.162'])

#parasim = smc.load_draws([0])
# ii = parasim.postsim.argmax()
# A0hat, Aphat = svar_model.para_trans(parasim[smc['paranames'][:-2]].ix[ii])

# irf_max = p.DataFrame(-svar_model.structural_irf(A0hat, Aphat, h=h)[0], 
#                       columns=yy)

# fevd_max = p.DataFrame(svar_model.structural_fevd(A0hat, Aphat, h=h)[0], 
#                       columns=yy)

irfs = p.DataFrame()
fevds = p.DataFrame()

from tqdm import tqdm

smc = {}
smc['paranames'] = [name for name in parasim.columns if name.startswith('var')]

nsim = parasim.shape[0]

for i in tqdm(range(nsim)):
    A0, Ap = svar_model.para_trans(parasim[smc['paranames'][:]].iloc[i])
    
    # for j in range(5):
    #     if np.linalg.inv(A0).dot(A0hat[:, j])[j] > 0:
    #         A0[:, j] = -A0[:, j]
    #         Ap[:, j] = -Ap[:, j]

    irf = svar_model.structural_irf(A0, Ap, h=h)
    irf[0] = irf[0] * np.sign(irf[0][0, 0])
    irfs = irfs.append(p.DataFrame(irf[0], columns=yy))
    
    fevd = svar_model.structural_fevd(A0, Ap, h=h)
    fevds = fevds.append(p.DataFrame(fevd[0], columns=yy))

median = irfs.groupby(irfs.index).median()
q05 = irfs.groupby(irfs.index).quantile(0.05)
q95 = irfs.groupby(irfs.index).quantile(0.95)

fig, ax = plt.subplots(nrows=2, ncols=3)
axi = ax.reshape(-1)
for i, name in enumerate(yy):
    #irf_max[name].plot(linewidth=3, ax=axi[i], color='green')
    median[name].plot(linewidth=3, ax=axi[i], color='black')
    q05[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')
    q95[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')

    axi[i].set_title(name)

plt.savefig(model+'_irf.pdf')


fig, ax = plt.subplots(nrows=2, ncols=3)
axi = ax.reshape(-1)
irf0 = irfs[irfs.index==0]
for i, name in enumerate(yy):
    irf0[name].plot(kind='kde', linewidth=3, ax=axi[i])
    axi[i].set_title(name)

plt.tight_layout()
plt.savefig(model+'_irf_impact.pdf')


median = fevds.groupby(irfs.index).median()
q05 = fevds.groupby(irfs.index).quantile(0.05)
q95 = fevds.groupby(irfs.index).quantile(0.95)


fig, ax = plt.subplots(nrows=2, ncols=3)
axi = ax.reshape(-1)
for i, name in enumerate(yy):
    #fevd_max[name].plot(linewidth=3, ax=axi[i], color='green')
    median[name].plot(linewidth=3, ax=axi[i], color='black')
    q05[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')
    q95[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')

    axi[i].set_title(name)
    axi[i].set_ylim(0, 1)

plt.savefig(model+'_fevd.pdf')
