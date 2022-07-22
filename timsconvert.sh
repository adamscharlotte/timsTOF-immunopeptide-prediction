# Install TIMSCONVERT following instructions by https://github.com/gtluu/timsconvert

wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
bash Anaconda3-2022.05-Linux-x86_64.sh -u
export PATH=$PATH:anaconda3/bin

conda create -n timsconvert python=3.7
source activate timsconvert
conda install -c bioconda nextflow

git clone https://www.github.com/gtluu/timsconvert
pip install -r timsconvert/requirements.txt

pip install git+https://github.com/gtluu/pyimzML

# Test the workflow
cd timsconvert/test
make download_test
make run_test
make run_nextflow_test

# Run TIMSCONVERT
ssh cadams@10.152.135.57

source activate timsconvert
nextflow run timsconvert/nextflow.nf --input /media/c2m/5_user_files/cadams/raw/E1_50fmol_bruker_S1-E2_1_6339.d

python timsconvert/bin/run.py --input /media/c2m/5_user_files/cadams/raw/E1_50fmol_bruker_S1-E2_1_6339.d

timsconvert/test/massive.ucsd.edu/MSV000088438/updates/2022-05-05_gbass_c9522c16/raw/timsconvert_raw_data_4/mini_tdf_1_1_1_1565.d


