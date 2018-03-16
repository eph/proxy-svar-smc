# Replication code for "Monetary Policy, Credit Spreads, and Business Cycle Fluctuations"
## by Dario Caldara and Ed Herbst

You can get the paper here: [paper](https://www.aeaweb.org/articles?id=10.1257/mac.20170294)

Email ed.herbst@gmail.com with comments/questions.


1. Figure 1: Impulse Response to a Monetary Policy Shock

	Python code:
	```sh 
	cd python
	python estimate_model.py --model 4eq --method smc
	python estimate_model.py --model 5eq --method smc
	python plot_figure_1.py
	```

	Matlab code:
	```sh
	cd matlab 
	matlab figure_1.m
	```


2. Figure 2: Contribution of Monetary Policy Shocks to FEVD
   
   Python code:
	```sh 
	cd python
	python estimate_model.py --model 4eq --method smc
	python estimate_model.py --model 5eq --method smc
	python plot_figure_2.pyy
	```

3. Table 1: Coefficients in the Monetary Policy Equation

   Python code:
	```sh 
	cd python
	python estimate_model.py --model 4eq --method smc
	python estimate_model.py --model 5eq --method smc
	python plot_figure_2.pyy
	```
	
4. Figure 3: Impulse Reponse to a Monetary Policy Shock

	Python code:

	```sh
	cd python
	python estimate_model.py --model 4eq --method smc
	python estimate_model.py --model 5eq --method smc
	python plot_figure_2.pyy
	```
	
	
5. Figure 4: Macroeconomic Implications of Financial Shocks

6. Table 2: Coefficients in the Monetary Policy Equation

7. Figure 5: Impulse Responses to a Monetary Policy Shock

8. Table 3: Determinants of the Change in the Intended Fed Funds Rate

9. Figure 6: Impulse Response to a Monetary Policy Shock 

10. Table 4: Local Projections

11. Figure 7: Impulse Responses to a Monetary Policy Shock

12. Figure 8: Impulse Responses to a Monetary Policy Shock

13. Table 5: Determinants of the Federal Funds Rate
