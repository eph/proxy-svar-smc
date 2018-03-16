import pandas as p
import statsmodels.formula.api as smf


data_file = '/home/eherbst/Dropbox/Bayesian-Proxy-SVAR/code/data/vardataJune7.txt'
var_data = p.read_csv(data_file, delim_whitespace=True)
var_data.index = p.period_range(freq='M', start='1973-01', periods=var_data.shape[0])

import statsmodels.api as sm


use_real_IP = False


var_data['DPPCE'] = var_data.PPCE.diff()
var_data['DLMRET'] = var_data.LMRET.diff()
var_data['DUNRATE'] = var_data['UNRATE'].diff()
var_data['DPPI_FIN'] = var_data['PPI_FIN'].diff()

if use_real_IP:
    var_data['DLOGIP'] = var_data['USIP'].diff()
else:
    var_data['DLOGIP'] = var_data['IPM'].diff()

pro = ['BAA_10YMOODY','DLOGIP','DPPI_FIN', 'UNRATE', 'FFR_SSR','IPM','PPI_FIN']
lag = 48
for pr in pro:
    for i in range(1, lag+1):
        var_data[pr+'LAG{:02d}'.format(i)] = var_data[pr].shift(i)

c = 12

for h in range(49):
    if use_real_IP:
        var_data['DLOGIP_H'+str(h)] = c*(var_data.USIP - var_data.USIP.shift(h+1)).shift(-(h+1)).shift(1)  / (h+1)
    else:
        var_data['DLOGIP_H'+str(h)] = c*(var_data.IPM - var_data.IPM.shift(h+1)).shift(-(h+1)).shift(1)  / (h+1)
    var_data['DPPI_FIN_H'+str(h)] = c*(var_data.PPI_FIN - var_data.PPI_FIN.shift(h+1)).shift(-(h+1)).shift(1)  / (h+1)
    var_data['UNRATE_H'+str(h)] = var_data.UNRATE.shift(-h)
    var_data['BAA_10YMOODY_H'+str(h)] = var_data.BAA_10YMOODY.shift(-h)
    var_data['FFR_SSR_H'+str(h)] = var_data.FFR_SSR.shift(-h)
    var_data['IPM_H'+str(h)] = var_data.IPM.shift(-h)
    var_data['PPI_FIN_H'+str(h)] = var_data.PPI_FIN.shift(-h)

def egon_predictability_regression(h=0, p=12, d0='1994-01', d1='2007-06', shock='EBP_OA', use_EBP=False,
controls=['DLOGIP','FFR_SSR','DPPCE'],y='DLOGIP_H'):
    y = y + str(h) 

    controls_vec = []
    #controls = ['DLOGIP_', 'FFR_EFF', 'DPPCE']
    for c in controls:
        controls_vec += [c + 'LAG{:02d}'.format(i) for i in range(1, p+1)]

    lags = '+'.join(controls_vec) #['DLOGIP_LAG{:02d}'.format(i) for i in range(1, p+1)])
    
    if use_EBP:
        formula = y + '~' + lags + '+ RRATE + EBP_OALAG01 + EBP_OALAG02 + EBP_OALAG03 +' + shock 
    else:
        formula = y + '~' + lags + "+" + shock

    #res = smf.ols(formula=formula, data=var_data[d0:d1]).fit().get_robustcov_results(cov_type='HAC', maxlags=h-1, kernel='uniform')
    res = smf.ols(formula=formula, data=var_data[d0:d1]).fit().get_robustcov_results()
    return res

def max_aic_reg(pmax=12, *args, **kwargs):
    aic0 = 10e10
    preferred_model = None
    for pp in range(1, pmax+1):
        res = egon_predictability_regression(p=pp, *args, **kwargs)
        if res.aic < aic0:
            aic0 = res.aic
            preferred_model = res
            preferred_model.p = pp
            
    return preferred_model



import numpy as np
hmax = 48
rescs = {}
res = {}
rescs_cs = {}
res_cs = {}

to_run = egon_predictability_regression 
for y in ['FFR_SSR_H','DLOGIP_H','UNRATE_H', 'DPPI_FIN_H', 'BAA_10YMOODY_H']:

    rrcs_pred = np.zeros((hmax,3))
    rr_pred = np.zeros((hmax,3))

    rrcs_pred_cs = np.zeros((hmax,3))
    rr_pred_cs = np.zeros((hmax,3))

    for h in range(1,hmax+1):

        controls=['UNRATE','DLOGIP','FFR_SSR','DPPI_FIN']

        d = to_run(p=12,controls=controls,h=h,shock='RR_SHOCK_SPREAD',y=y)
        rrcs_pred[h-1] = [d.params[-1], d.params[-1] - 1.96*d.bse[-1], d.params[-1] + 1.96*d.bse[-1]]

        d = to_run(p=12,controls=controls,h=h,shock='RR_SHOCK_NOSPREAD',y=y)
        rr_pred[h-1] = [d.params[-1], d.params[-1] - 1.96*d.bse[-1], d.params[-1] + 1.96*d.bse[-1]]

        controls=['UNRATE','DLOGIP','FFR_SSR','DPPI_FIN','BAA_10YMOODY']

        d = to_run(p=12,controls=controls,h=h,shock='RR_SHOCK_SPREAD',y=y)
        rrcs_pred_cs[h-1] = [d.params[-1], d.params[-1] - 1.96*d.bse[-1], d.params[-1] + 1.96*d.bse[-1]]

        d = to_run(p=12,controls=controls,h=h,shock='RR_SHOCK_NOSPREAD',y=y)
        rr_pred_cs[h-1] = [d.params[-1], d.params[-1] - 1.96*d.bse[-1], d.params[-1] + 1.96*d.bse[-1]]


    rescs[y] = p.DataFrame(rrcs_pred, columns=['mn','lb','ub'])
    res[y] = p.DataFrame(rr_pred, columns=['mn','lb','ub'])

    
    rescs_cs[y] = p.DataFrame(rrcs_pred_cs, columns=['mn','lb','ub'])
    res_cs[y] = p.DataFrame(rr_pred_cs, columns=['mn','lb','ub'])
