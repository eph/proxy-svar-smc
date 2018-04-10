import argparse

import numpy as np
import pandas as p
import matplotlib.pyplot as plt

from tqdm import tqdm
import json

from collections import defaultdict
from models import yy,  models

import os

parser = argparse.ArgumentParser(description='plot impulse resposes.')
parser.add_argument('--model', help='the model(s) to plot', action='store',
                    choices=models.keys(), default='5eq', nargs='*')
parser.add_argument('--nsim', help='the number of posterior draws to use',
                    action='store', type=int, default=9600)
parser.add_argument('--horizon', help='horizon of IRF', action='store',
                    type=int, default=49)
parser.add_argument('--sim-dir', help='directory of simulations', action='store',
                    type=str, default='.')
parser.add_argument('--output-file', help='file name', action='store',
                    type=str, default=None)


args = parser.parse_args()


nsim = args.nsim


def long_run_coefficients(A0, Ap, equation=0):
    """
    Computes the Long Run Coefficients a la Sims-Zha
    """
    eqi = equation
    ny = A0.shape[0]

    alpha = A0[:, eqi]
    alpha = [np.r_[-A0[i, eqi], Ap[:-1, eqi][i::ny]] for i in range(ny)]

    delta = -alpha[eqi]

    lr_coeff_diff = np.array([ai.sum() / delta.sum() for ai in alpha])
    sr_coeff_diff = np.array([ai.sum() / A0[0, 0] for ai in alpha])

    sr_coeff_diff[0] = (Ap[:-1, 0][::ny] / A0[0, 0]).sum()
    x = sr_coeff_diff[0].copy()
    sr_coeff_diff = sr_coeff_diff

    sr_coeff_lev = np.array([ai.cumsum().sum() / A0[0, 0] for ai in alpha])
    sr_coeff_lev = sr_coeff_lev

    sr_coeff_diff[0] = x
    return sr_coeff_lev, sr_coeff_diff


cum_elasticities = p.DataFrame()
con_elasticities = p.DataFrame()
for m in args.model:
    model = models[m]
    results = json.loads(open(os.path.join(args.sim_dir, m + '.json')).read())

    parasim = p.DataFrame(results['posterior.162'])
    inds = np.random.choice(
        parasim.shape[0], size=parasim.shape[0], p=parasim.weights)
    parasim = parasim.iloc[inds]

    paranames = [name for name in parasim.columns if name.startswith('var')]

    cumulative_elasticities = np.zeros((nsim, len(model['yy'])))
    contemp_elasticities = np.zeros((nsim, len(model['yy']) - 1))

    four_variable = '4' in m

    for i in tqdm(range(nsim)):
        parai = parasim[paranames[:-2]].iloc[i]
        A0, Ap = model['svar_model'].para_trans(parai)

        src_lev, src_diff = long_run_coefficients(A0, Ap)

        cumulative_elasticities[i][0] = src_diff[0]
        cumulative_elasticities[i][1] = src_lev[1]
        cumulative_elasticities[i][2] = src_diff[2]
        cumulative_elasticities[i][3] = src_lev[3]
        if not four_variable:
            cumulative_elasticities[i][4] = src_diff[4]

        contemp_elasticities[i] = -A0[1:, 0] / A0[0, 0]

    elasticities_m = p.DataFrame(cumulative_elasticities, columns=model['yy'])
    elasticities_m['model'] = m
    cum_elasticities = cum_elasticities.append(elasticities_m)

    elasticities_m = p.DataFrame(contemp_elasticities, columns=model['yy'][1:])
    elasticities_m['model'] = m
    con_elasticities = con_elasticities.append(elasticities_m)

if args.output_file is None:
    output_file = args.sim_dir + '_'.join(args.model) + '__elasticities.tex'
else:
    output_file = args.sim_dir + '/' + args.output_file


import sys
sys.stdout = open(output_file,  'w')


def q05(x): return np.percentile(x, 5)


def q95(x): return np.percentile(x, 95)


results = con_elasticities.groupby('model').aggregate([np.median, q05, q95])

rows = ['BAA10YMOODY', 'LPPI', 'LIPM', 'UNRATE']

nmodels = len(args.model)
for r in rows:
    row_template = '{:20s} & ' + ' & '.join(nmodels * ['    {: 5.2f}     '])
    print(row_template.format(r, *results[r]['median']))
    row_template = '{:20s} & ' + ' & '.join(nmodels * ['[{: 5.2f}, {: 5.2f}]'])
    print(row_template.format('', *[results[r].ix[m][q]
                                    for m in args.model for q in ['q05', 'q95']]))

print('')
results = cum_elasticities.groupby('model').aggregate([np.median, q05, q95])

rows = ['BAA10YMOODY', 'LPPI', 'LIPM', 'UNRATE', 'EFFR']

nmodels = len(args.model)
for r in rows:
    row_template = '{:20s} & ' + ' & '.join(nmodels * ['    {: 5.2f}     '])
    print(row_template.format(r, *results[r]['median']))
    row_template = '{:20s} & ' + ' & '.join(nmodels * ['[{: 5.2f}, {: 5.2f}]'])
    print(row_template.format('', *[results[r].ix[m][q]
                                    for m in args.model for q in ['q05', 'q95']]))
