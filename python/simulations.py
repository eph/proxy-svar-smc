import sys
sys.path.append('/msu/home/m1eph00/projects/dsge-book/code/helper')
from helper import SMCResults

smc_baseline = SMCResults('svar12_loose-mix', npart=9600, nblocks=50, nphi=500, lam=2.7)

smc_tight = SMCResults('svar12_tight-mix', npart=9600, nblocks=25, nphi=2000, lam=2.7)

