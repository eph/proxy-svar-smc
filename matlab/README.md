#Replication code for "Monetary Policy, Credit Spreads, and Business Cycle Fluctuations"
### by Dario Caldara and Ed Herbst

This MATLAB code runs an MCMC algorithm to sample the posteriors of the
Proxy-SVARs.

The file `main_MCMC_mixture_proposal_priors.m` simulates from the posterior and
creates figures.

## Figures 1 + 2
+ run main_MCMC_mixture_proposal_priors with the following settings: 
  model_vec = {'m2levNoSp', 'm2lev'}
  instrList = {'EGON_KUTTNER_NI'}

  
