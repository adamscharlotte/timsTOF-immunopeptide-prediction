# conda create -n proslsit-annotate python=3.8
conda activate prosit-annotate

cd /Users/adams/Projects/300K
cd /Users/adams/Projects/300K/prosit_io
pip install -e .
cd /Users/adams/Projects/300K/fundamentals
pip install -e .
pip install git+https://github.com/wilhelm-lab/PROSPECT
# pip install -e git+https://gitlab.lrz.de/proteomics/prosit_tools/prosit_grpc.git@v1.2.0#egg=prosit_grpc

# Copy files 
# scp -P 4024 full-truncated-qc.hdf5 root@ac922a.ucc.in.tum.de:/scratch/data/
# scp -P 4024 /Users/adams/Projects/300K/2022-library-run/Annotation/full-truncated-qc/annotated/full-truncated-qc-fixed_*.hdf5 root@ac922a.ucc.in.tum.de:/scratch/data/
scp -P 4024 /Users/adams/Projects/300K/2022-library-run/Annotation/full-truncated-qc/annotated/full-truncated-qc-un-callibrated_*.hdf5 root@ac922a.ucc.in.tum.de:/scratch/data/
scp -P 4024 /Users/adams/Projects/300K/prosit_non_tryptic_fragmentation/hla_hcd/*.hdf5 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof/training
scp -P 4024 /Users/adams/Projects/300K/2022-library-run/Annotation/precursor-consensus/annotated-hdf5/full-truncated-qc-un-callibrated-precursor-consensus_*.hdf5 root@ac922a.ucc.in.tum.de:/scratch/data/
scp -P 4024 /Users/adams/Projects/300K/2022-library-run/Annotation/precursor-consensus/calibrated-hdf5/precursor-20ppm-calibrated_*.hdf5 root@ac922a.ucc.in.tum.de:/scratch/data/

t3_W1xi_MhG3maK9_a$

ssh -p 4024 root@ac922a.ucc.in.tum.de
t3_W1xi_MhG3maK9_a$

# change config file
cd /scratch/vitor/training/tims_tof

# Open screen
screen -S tims_tof_train
screen -r tims_tof_train
screen -r -d tims_train

# Remove log
cd /scratch/vitor/training/tims_tof/training

# Train model
cd /scratch/vitor/prosit/prosit
# python training.py -w "/scratch/vitor/training/tims_tof/training/weight_163_0.11385.hdf5" /scratch/vitor/training/tims_tof

# With weights (refinement)
# python training.py -r /scratch/vitor/training/tims_tof
python training.py -w /scratch/vitor/training/tims_tof/training/weights_163_0.11385.hdf5 /scratch/vitor/training/tims_tof

# From scratch
python training.py /scratch/vitor/training/tims_tof

cd /scratch/vitor/training/tims_tof

# -----------------------------------------------------------------------------

cd /Users/adams/Projects/300K/2022-library-run/Annotation/log/
scp -P 4024 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof/training-1/training/training_log.csv .
scp -P 4024 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof/training/training_log.csv .


cd /Users/adams/Code/timsTOF-immunopeptide-prediction
# -----------------------------------------------------------------------------

ssh -p 4024 root@ac922a.ucc.in.tum.de
t3_W1xi_MhG3maK9_a$

git clone https://gitlab.lrz.de/compmass/prosit/tools/fundamentals.git
charlotte.adams
proteomics

git clone https://gitlab.lrz.de/compmass/prosit/tools/prosit_io.git
git clone https://gitlab.lrz.de/proteomics/prosit_tools/prosit_grpc.git

conda create -n prosit-annotate python=3.8
conda activate prosit-annotate
cd fundamentals
pip install .
cd prosit_io
pip install .

cd prosit_grpc
pip install Cython
git checkout oktoberfest_revamp
# go to pyproject.toml and remove line containing fundamentals.
pip install .





# --------
cd /scratch/vitor/training/tims_tof