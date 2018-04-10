OD=/msu/scratch-m1eph00/proxy-svar/final-results/new-observables

# figures
python compare_irfs.py --model 4eq 5eq --output-file figure_1.pdf --sim-dir $OD 
python compare_fevds.py --model 4eq 5eq --output-file figure_2.pdf --sim-dir $OD
python compare_irfs.py --model 5eq_cholesky --overlay 5eq --output-file figure_3.pdf --sim-dir $OD 
python compare_irfs.py --model 5eq_cholesky --overlay 5eq_fin --output-file figure_4.pdf --sim-dir $OD
python compare_irfs.py --model 5eq_tight --output-file figure_5.pdf --sim-dir $OD 
python compare_irfs.py --model 4eq_cholesky_RRCS 5eq_cholesky_RRCS --overlay 4eq_cholesky_RR 5eq_cholesky_RR --output-file figure_6.pdf --sim-dir $OD 

# tables
python compare_elasticities.py --model 4eq 5eq  --output-file table_1.tex --sim-dir $OD
python compare_elasticities.py --model 5eq 5eq_tight --output-file table_2.tex --sim-dir $OD
python local_projections.py --output-file $OD/table_4.tex
