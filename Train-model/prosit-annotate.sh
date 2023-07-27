# --------------------------------- INSTALL ANNOTATION CODE ---------------------------------
# conda create -n proslsit-annotate python=3.8
# conda activate prosit-annotate

# cd /Users/adams/Projects/300K
# cd /Users/adams/Projects/300K/prosit_io
# pip install -e .
# cd /Users/adams/Projects/300K/fundamentals
# pip install -e .
# pip install git+https://github.com/wilhelm-lab/PROSPECT
# pip install -e git+https://gitlab.lrz.de/proteomics/prosit_tools/prosit_grpc.git@v1.2.0#egg=prosit_grpc

# --------------------------------------- TRAIN MODEL ---------------------------------------

# Copy files
# scp -P 4024 /Users/adams/Projects/300K/2022-library-run/Annotation/precursor-consensus/calibrated-hdf5/precursor-20-ppm-calibrated_*.hdf5 root@ac922a.ucc.in.tum.de:/scratch/data/
# scp -P 4024 /Users/adams/Projects/300K/2022-library-run/Annotation/scan-consensus/calibrated-hdf5/scan-40-ppm-calibrated-mapped_*.hdf5 root@ac922a.ucc.in.tum.de:/scratch/data/
scp -P 4024 /Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/hdf5/scan-40-ppm-calibrated_*.hdf5 root@ac922a.ucc.in.tum.de:/scratch/data/

# Change config file
cd /scratch/vitor/training/tims_tof
cd /scratch/vitor/training/tims_tof_frozen_ce

# Open screen
# screen -S tims_tof
screen -r tims_tof_train
screen -r tims_train

# Remove log
cd /scratch/vitor/training/tims_tof/training
cd /scratch/vitor/training/tims_tof_frozen_ce/training

# Train model
cd /scratch/vitor/prosit/prosit

# With weights (refinement)
# python training.py -r /scratch/vitor/training/tims_tof
python training.py -w /scratch/vitor/training/tims_tof/training/weights_163_0.11385.hdf5 /scratch/vitor/training/tims_tof

python training_wandb_tims.py -w /scratch/vitor/training/tims_tof/training/weights_163_0.11385.hdf5 /scratch/vitor/training/tims_tof
# 
# python training.py -w /scratch/vitor/training/tims_tof_frozen_ce/training/weights_163_0.11385.hdf5 /scratch/vitor/training/tims_tof_frozen_ce
python training_wandb_tims.py -w /scratch/vitor/training/tims_tof_frozen_ce/training/weights_163_0.11385.hdf5 /scratch/vitor/training/tims_tof_frozen_ce

# From scratch
# python training.py /scratch/vitor/training/tims_tof
python training_wandb_tims.py /scratch/vitor/training/tims_tof

# ------------------------------------ PREDICT TEST SET -------------------------------------
cd /scratch/vitor/prosit/prosit

python prediction.py /scratch/vitor/training/tims_tof

#  --weights_file /scratch/vitor/training/tims_tof/training/weights_28_0.143392.hdf5
#  -d /scratch/data/scan-40-ppm-calibrated_test.hdf5 -n buker_test

python prediction.py --help

cd /scratch/vitor/training/tims_tof

scp -r -P 4024 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof/prediction .
scp -P 4024 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof/prediction/* .

# -------------------------------------- DOWNLOAD LOG ---------------------------------------

cd /Users/adams/Projects/300K/2022-library-run/Annotation/log/
scp -P 4024 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof/training/training_log.csv .
scp -P 4024 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof_frozen_ce/training/training_log.csv .


scp -P 4024 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof/config.yml .


# ------------------------------------ CHANGE THE MODEL -------------------------------------

scp -P 4024 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof_frozen_ce/model.yml .

scp -P 4024 root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof_frozen_ce/model.yml .


scp -P 4024 /Users/adams/Downloads/model.yml root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof_frozen_ce/
scp -P 4024 /Users/adams/Downloads/Frozen-me/model.yml root@ac922a.ucc.in.tum.de:/scratch/vitor/training/tims_tof_frozen_ce/


# ------------------------------------------ CODE? ------------------------------------------

# git clone https://gitlab.lrz.de/compmass/prosit/tools/fundamentals.git
# charlotte.adams
# proteomics

# git clone https://gitlab.lrz.de/compmass/prosit/tools/prosit_io.git
# git clone https://gitlab.lrz.de/proteomics/prosit_tools/prosit_grpc.git

# conda create -n prosit-annotate python=3.8
# conda activate prosit-annotate
# cd fundamentals
# pip install .
# cd prosit_io
# pip install .

# cd prosit_grpc
# pip install Cython
# git checkout oktoberfest_revamp
# # go to pyproject.toml and remove line containing fundamentals.
# pip install .
