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
parser.add_argument('--sim-dir', help='directory of simulations', action='store',
                    type=str, default='./')
parser.add_argument('--output-file', help='file name', action='store', 
                    type=str, default=None)



args = parser.parse_args()

plot_models = args.model
h = args.horizon
nsim = args.nsim 
nmodels = len(plot_models)


fevds = defaultdict(p.DataFrame)


for m in plot_models:
    model = models[m]
    results = json.loads(open(args.sim_dir + '/' + m + '.json').read())

    parasim = p.DataFrame(results['posterior.162'])
    inds = np.random.choice(parasim.shape[0], size=parasim.shape[0], p=parasim.weights)
    parasim = parasim.iloc[inds]

    paranames = [name for name in parasim.columns if name.startswith('var')]

    irf_sel = 0
    if 'cholesky' in m: irf_sel = 4

    for i in tqdm(range(nsim)):
        parai = parasim[paranames[:-2]].iloc[i]
        A0, Ap = model['svar_model'].para_trans(parai)
         
        fevd = model['svar_model'].structural_fevd(A0, Ap, h=h)
        fevds[m] = fevds[m].append(p.DataFrame(fevd[irf_sel], columns=model['yy']))



fig, ax = plt.subplots(ncols=nmodels, nrows=max([len(models[m]['yy']) for m in plot_models]))
for j, m in enumerate(plot_models):
    median = fevds[m].groupby(fevds[m].index).median()
    q05 = fevds[m].groupby(fevds[m].index).quantile(0.05)
    q95 = fevds[m].groupby(fevds[m].index).quantile(0.95)


    if nmodels > 1:
        axi = ax[:, j].reshape(-1)
    else:
        axi = ax.reshape(-1)


    for i, name in enumerate(models[m]['plot_yy']):
        median[name].plot(linewidth=3, ax=axi[i], color='black')
        axi[i].fill_between(q05[name].index, q05[name], q95[name], color='orange', alpha=0.5)
        # q05[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')
        # q95[name].plot(linewidth=3, ax=axi[i], color='red', linestyle='dashed')
        axi[i].set_ylim(0, 1)
        axi[i].set_title(name)

if args.output_file is None:
    output_file = args.sim_dir + '_'.join(plot_models)+'__irf.pdf'
else:
    output_file = args.sim_dir + '/' + args.output_file

fig.set_size_inches(nmodels*6, 15)
plt.tight_layout()
plt.savefig(output_file)



