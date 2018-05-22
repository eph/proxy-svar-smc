# Replication code for "Monetary Policy, Credit Spreads, and Business Cycle Fluctuations"
### by Dario Caldara and Ed Herbst

You can get the paper [here](https://www.aeaweb.org/articles?id=10.1257/mac.20170294).

Email ed.herbst@gmail.com with comments/questions.


### Data

Data used in the paper:

* `data/CHdata.txt` contains the macroeconomic series used to estimate the
  proxy SVARs and the local projections.
  * The federal funds rate is the average effective rate over the last week of the month
  * The unemployment rate, and the producer
    price index for final goods are taken from Coibion (2012) "Are the effects
    of monetary policy shocks big or small?"
  * We use manufacturing industrial production index (NAICS);
  * The BAA spread is Moody's Seasoned Baa Corporate Bond Yield Relative to
    Yield on 10-Year Treasury Constant Maturity, downloaded from the St. Louis
    Fed FRED database.
  * The Arouba, Diebold, Scotti business condition index is downloaded from the
    Federal Reserve Bank of Philadelphia
    [website](https://www.philadelphiafed.org/research-and-data/real-time-center/business-conditions-index). 
  * We use real personal consumption expenditures in nondurable goods, deflated
    using its own price deflator (downloaded from FRED);
  * We use the value-weighted total stock market index; Core loans are the sum
	of loans to households and businesses. 
  * Business loans include commercial and industrial (C&I) loans and business
	loans secured by commercial real estate; household loans include residential
	mortgages, credit card loans, and other consumer loans. All series were
	obtained from the Federal Reserve's H.8 Statistical Release.
  * mhf is the series of monetary policy surprises constructed using high frequency data
  * mrr and mrrcs are the series of Romer and Romer shocks constructed following Equation (16);
  * mcgcs is the series of monetary policy shocks constructed following Equation (18)
		  
* `data/RR-CG-data.dta` contains the data used to estimate the Romer and Romer
   (2004) and the Coibion Gorodnichenko regressions (equations (16) and (18) in
   the paper). The naming of the variables is consistent with these two papers.



### Estimating Models

For most models, you can either use Matlab (Gibbs Sampler) or Python/Fortran
(Sequential Monte Carlo) to estimate the models.  See the anaconda environment
file (proxy-svar.yaml) for the python packages necessary.

| Model                       | Matlab                       | Python                                             |
|-----------------------------|------------------------------|----------------------------------------------------|
| 4 Equation BP-SVAR          | matlab main_BPSVAR.m         | python estimate_model.py --model 4eq               |
| 5 Equation BP-SVAR          | matlab main_BPSVAR.m         | python estimate_model.py --model 5eq               |
| 5 Equation Cholesky         | matlab main_cholesky.m       | python estimate_model.py --model 5eq_cholesky      |
| 5 Equation BP-SVAR, w/fin   | n/a                          | python estimate_model.py --model 5eq_fin           |
| 4 Equation Hybrid VAR, RRCS | matlab main_cholesky.m       | python estimate_model.py --model 4eq_cholesky_RRCS |
| 4 Equation Hybrid VAR, RR   | matlab main_cholesky.m       | python estimate_model.py --model 4eq_cholesky_RR   |
| 5 Equation Hybrid VAR, RRCS | matlab main_cholesky.m       | python estimate_model.py --model 5eq_cholesky_RRCS |
| 5 Equation Hybrid VAR, RR   | matlab main_cholesky.m       | python estimate_model.py --model 5eq_cholesky_RR   |
| 9 Equation BP-SVAR          | n/a                          | python estimate_model.py --model 9eq               |
| Frequentist Estimation      | matlab main_proxy_bootstra.m | n/a                                                |

### Generating Tables and Figures


1. Figure 1: Impulse Reponse to a Monteray Policy Shock

   ```sh 
   cd python
   python compare_irfs.py --model 4eq 5eq 
   ```

2. Figure 2: Contribution of Monetary Policy Shocks to FEVD
   
   ```sh 
   cd python
   python compare_fevds.py --model 4eq 5eq
   ```

3. Table 1: Coefficients in the Monetary Policy Equation

   ```sh 
   cd python
   python compare_elasticities.py --model 4eq 5eq
   ```
	
4. Figure 3: Impulse Reponse to a Monetary Policy Shock

	Python code:
	```sh
	cd python
    python compare_irfs.py --model 5eq_cholesky --overlay 5eq
	```
	
5. Figure 4: Macroeconomic Implications of Financial Shocks

   Python code:
	```sh
	cd python
    python compare_irfs.py --model 5eq_cholesky --overlay 5eq_fin
	```

6. Table 2: Coefficients in the Monetary Policy Equation

	```sh 
	cd python
    python compare_elasticities.py --model 5eq 5eq_tight
	```

7. Figure 5: Impulse Responses to a Monetary Policy Shock

	```sh
	cd python
    python compare_irfs.py --model 5eq_tight
    ```

8. Table 3: Determinants of the Change in the Intended Fed Funds Rate

   ```sh
   cd stata
   stata -b do RR-CG-reg.do
   ```

9. Figure 6: Impulse Response to a Monetary Policy Shock 

   ```sh
   python compare_irfs.py --model 4eq_cholesky_RRCS 5eq_cholesky_RRCS --overlay 4eq_cholesky_RR 5eq_cholesky_RR 
   ```

10. Table 4: Local Projections

    ```sh
	cd python
    python local_projections.py
    ```

11. Figure 7: Impulse Responses to a Monetary Policy Shock

    The figure is created automatically after the matlab estimation.
    
12. Table 5: Determinants of the Federal Funds Rate

    ```sh
    cd stata
    stata -b do RR-CG-reg.do
    ```
