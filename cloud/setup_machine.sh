# TODO: miniconda + yes installers
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
cd ../python
conda env create -f proxy-svar.yaml
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
python estimate_model.py --npart 9600 --nproc 28 --model 4eq
python 
