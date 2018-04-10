import argparse

import pandas as p
import statsmodels.formula.api as smf


from models import data as var_data

parser = argparse.ArgumentParser(description='table of local projections')
parser.add_argument('--output-file', help='file name', action='store',
                    type=str, default=None)
args = parser.parse_args()


key = dict(FFR='EFFR_LW',
           IP='LIPM',
           U='UNRATE',
           P='LPPI',
           CS='BAA10YMOODY',
           DIP='DIP',
           DP='DP')



var_data['DP'] = var_data[key['P']].diff()
var_data['DIP'] = var_data[key['IP']].diff()

max_lag = 48
to_lag = ['CS', 'DIP', 'DP', 'U', 'FFR']
for pr in to_lag:
    for i in range(1, max_lag + 1):
        var_data[key[pr] + '_L%d' % i] = var_data[key[pr]].shift(i)

c = 12
for h in range(49):
    # growth rates
    var_data['DIP_H' + str(h)] = c * (var_data[key['IP']] -
                                      var_data[key['IP']].shift(h + 1)).shift(-(h + 1)).shift(1) / (h + 1)
    var_data['DP_H' + str(h)] = c * (var_data[key['P']] -
                                     var_data[key['P']].shift(h + 1)).shift(-(h + 1)).shift(1) / (h + 1)

    # levels
    var_data['U_H' + str(h)] = var_data[key['U']].shift(-h)
    var_data['CS_H' + str(h)] = var_data[key['CS']].shift(-h)
    var_data['FFR_H' + str(h)] = var_data[key['FFR']].shift(-h)


def uniform_weights(h): return np.ones(h + 1)


def egon_predictability_regression(h=0, p=12, d0='1994-01', d1='2007-06',
                                   shock='MRR', controls=['U', 'DIP', 'FFR', 'DP'], y='DLOGIP_H'):
    y = y + str(h)

    controls_vec = []

    for c in controls:
        controls_vec += [key[c] + '_L%d' % i for i in range(1, p + 1)]

    lags = '+'.join(controls_vec)

    formula = y + '~' + lags + "+" + shock

    model = smf.ols(formula=formula, data=var_data[d0:d1])
    res = model.fit().get_robustcov_results(cov_type='HAC', maxlags=h - 1)

    return res

# if you want to mess with model selection


def max_aic_reg(pmax=12, *args, **kwargs):
    aic0 = 10e10
    preferred_model = None
    for pp in range(1, pmax + 1):
        res = egon_predictability_regression(p=pp, *args, **kwargs)
        if res.aic < aic0:
            aic0 = res.aic
            preferred_model = res
            preferred_model.p = pp

    return preferred_model


import numpy as np
hmax = 19
rescs = {}
res = {}
rescs_cs = {}
res_cs = {}

to_run = egon_predictability_regression
for y in ['FFR_H', 'DIP_H', 'U_H', 'DP_H', 'CS_H']:

    rrcs_pred = np.zeros((hmax, 2))
    rr_pred = np.zeros((hmax, 2))

    rrcs_pred_cs = np.zeros((hmax, 2))
    rr_pred_cs = np.zeros((hmax, 2))

    for h in range(1, hmax + 1):

        controls = ['U', 'DIP', 'FFR', 'DP']

        d = to_run(p=12, controls=controls, h=h, shock='MRR', y=y)
        rr_pred[h - 1] = [d.params[-1], d.bse[-1]]

        d = to_run(p=12, controls=controls, h=h, shock='MRRCS', y=y)
        rrcs_pred[h - 1] = [d.params[-1], d.bse[-1]]

        controls.append('CS')

        d = to_run(p=12, controls=controls, h=h, shock='MRR', y=y)
        rr_pred_cs[h - 1] = [d.params[-1], d.bse[-1]]
        
        d = to_run(p=12, controls=controls, h=h, shock='MRRCS', y=y)
        rrcs_pred_cs[h - 1] = [d.params[-1], d.bse[-1]]



    rescs[y] = p.DataFrame(rrcs_pred, columns=['mn', 'se'])
    res[y] = p.DataFrame(rr_pred, columns=['mn', 'se'])

    rescs_cs[y] = p.DataFrame(rrcs_pred_cs, columns=['mn', 'se'])
    res_cs[y] = p.DataFrame(rr_pred_cs, columns=['mn', 'se'])


if args.output_file is None:
    pass
else:
    import sys
    sys.stdout = open(args.output_file,  'w')

h_to_print = [0, 18]
line = '{:25s} & {: 5.2f} % {: 5.2f} \\\\'.format
for h in h_to_print:
    print('          h = %d' % h)
    for y in ['FFR_H', 'DIP_H', 'U_H', 'DP_H', 'CS_H']:
        print(line(y, res[y].iloc[h].mn, rescs_cs[y].iloc[h].mn))
        print(line('', res[y].iloc[h].se, rescs_cs[y].iloc[h].se))

    print()


import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=5, ncols=2)

names = ['FFR_H', 'DIP_H', 'U_H', 'DP_H', 'CS_H']
for i, n in enumerate(names):
    for j, v in enumerate([res, rescs_cs]):
        v[n].mn.plot(ax=ax[i, j])
        ax[i, j].fill_between(v[n].mn.index,
                              v[n].mn - 1.96 * v[n].se,
                              v[n].mn + 1.96 * v[n].se,
                              alpha=0.3)
fig.tight_layout()
plt.savefig('local_projections.pdf')
