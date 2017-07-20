sudo apt-get update
sudo apt-get install gcc
sudo apt-get install mpich
sudo apt-get install make
sudo apt-get install awscli
# TODO: miniconda + yes installers
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install python=3.5
conda install fortress -c eherbst
conda install pyvar -c eherbst
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
python create_models.py
