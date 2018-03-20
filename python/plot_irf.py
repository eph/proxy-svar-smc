import argparse

import numpy as np
import pandas as p
import matplotlib.pyplot as plt

from tqdm import tqdm
import json

from collections import defaultdict
from models import yy,  models



parser = argparse.ArgumentParser(description='plot impulse resposes.')
parser.add_argument('--model', help='the model(s) to plot', action='store',
                    choices=models.keys(), default='5eq', nargs='*')
parser.add_argument('--nsim', help='the number of posterior draws to use', 
                    action='store', type=int, default=9600)
parser.add_argument('--horizon', help='horizon of IRF', action='store',
                    type=int, default=49)



args = parser.parse_args()

plot_models = args.model
h = args.horizon
nsim = args.nsim 
nmodels = len(plot_models)

ylims = dict(EFFR=(-0.375, 0.50), 
             LIPM=(-1, 1), 
             UNRATE=(-0.12, 0.1), 
             LPPI=(-1, 0.5), 
             BAA10YMOODY=(-0.1, 0.15))

irfs = defaultdict(p.DataFrame)
fevds = defaultdict(p.DataFrame)


for m in plot_models:
    model = models[m]
    results = json.loads(open(m + '.json').read())

    parasim = p.DataFrame(results['posterior.162'])
    inds = np.random.choice(parasim.shape[0], size=parasim.shape[0], p=parasim.weights)
    parasim = parasim.iloc[inds]

    paranames = [name for name in parasim.columns if name.startswith('var')]

    irf_sel = 0
    if 'cholesky' in m: irf_sel = 4

    for i in tqdm(range(nsim)):
        parai = parasim[paranames[:-2]].iloc[i]
        A0, Ap = model['svar_model'].para_trans(parai)
    
        # for j in range(5):
        #     if np.linalg.inv(A0).dot(A0hat[:, j])[j] > 0:
        #         A0[:, j] = -A0[:, j]
        #         Ap[:, j] = -Ap[:, j]

        irf = model['svar_model'].structural_irf(A0, Ap, h=h)
        irf[irf_sel] = irf[irf_sel] * np.sign(irf[irf_sel][irf_sel, irf_sel])
        irfs[m] = irfs[m].append(p.DataFrame(irf[irf_sel], columns=model['yy']))
         
        fevd = model['svar_model'].structural_fevd(A0, Ap, h=h)
        fevds[m] = fevds[m].append(p.DataFrame(fevd[irf_sel], columns=model['yy']))



fig, ax = plt.subplots(ncols=nmodels, nrows=max([len(models[m]['yy']) for m in plot_models]))
for j, m in enumerate(plot_models):
    median = irfs[m].groupby(irfs[m].index).median()
    q05 = irfs[m].groupby(irfs[m].index).quantile(0.05)
    q95 = irfs[m].groupby(irfs[m].index).quantile(0.95)


    if nmodels > 1:
        axi = ax[:, j].reshape(-1)
    else:
        axi = ax.reshape(-1)


    for i, name in enumerate(models[m]['yy']):
        median[name].plot(linewidth=3, ax=axi[i], color='black')
        q05[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')
        q95[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')

        #axi[i].set_ylim(ylims[name]);
        axi[i].set_title(name)

fig.set_size_inches(nmodels*6, 15)
plt.tight_layout()
plt.savefig('_'.join(plot_models)+'__irf.pdf')


# fig, ax = plt.subplots(nrows=nmodels, ncols=max([len(m.yy) for m in plot_models]))
# for j in range(nmodels):
#     axi = ax[:, j].reshape(-1)
#     irf0 = irfs[irfs.index==0]
#     for i, name in enumerate(yy):
#         irf0[name].plot(kind='kde', linewidth=3, ax=axi[i])
#         axi[i].set_title(name)

# fig.set_size_inches(nmodels*3, 15)
# plt.tight_layout()
# plt.savefig(model+'_irf_impact.pdf')


# median = fevds.groupby(irfs.index).median()
# q05 = fevds.groupby(irfs.index).quantile(0.05)
# q95 = fevds.groupby(irfs.index).quantile(0.95)


# fig, ax = plt.subplots(nrows=2, ncols=3)
# axi = ax.reshape(-1)
# for i, name in enumerate(yy):
#     median[name].plot(linewidth=3, ax=axi[i], color='black')
#     q05[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')
#     q95[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')

#     axi[i].set_title(name)
#     axi[i].set_ylim(0, 1)

# plt.savefig(model+'_fevd.pdf')
