# /Users/adams/opt/miniconda3/envs/ann_solo/bin/python3
# ssh cadams@10.152.135.57
# /home/cadams/anaconda3/envs/prosit-annotate/bin/python3

# from pickle import FRAME
from re import I
import pandas as pd
import numpy as np
import os
from pyteomics import mgf

import prosit_grpc
from prosit_grpc.predictPROSIT import PROSITpredictor

# import argparse, pathlib
# parser = argparse.ArgumentParser()
# parser.add_argument('csv_path', type=str)					# Filename
# parser.add_argument('mgf_path', type=str)					# Filename
# parser.add_argument('charge_name', type=str)					# Filename
# args = parser.parse_args()

# csv_path = args.csv_path
# mgf_path = args.mgf_path
# charge_name = args.charge_name

mgf_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/mgf/HLA_133_tims_pred.mgf"
# model = "Prosit_2020_intensity_hcd"
# model = "Prosit_2020_intensity_cid"
# collision_energy = 30
model = "Prosit_2023_intensity_tims"
collision_energy = 31

path = "/home/cadams/anaconda3/envs/oktoberfest-0_1/lib/python3.8/site-packages/oktoberfest/certificates"
predictor = PROSITpredictor(server="10.152.135.57:8500",
            path_to_ca_certificate=os.path.join(path, "Proteomicsdb-Prosit-v2.crt"),
            path_to_certificate=os.path.join(path, "oktoberfest-production.crt"),
            path_to_key_certificate=os.path.join(path, "oktoberfest-production.key"),)

data = {
    "SCAN_NUMBER": "14565",
    "MASS": 1194.5385,
    "RTs": 1275,
    "MODIFIED_SEQUENCE": ["GVDAANSAAQQY"],
    "PRECURSOR_CHARGE": [1],
    "ORIG_COLLISION_ENERGY": [collision_energy]
}

df_133 = pd.DataFrame(data, columns=["SCAN_NUMBER", "MASS", "RTs", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE", "ORIG_COLLISION_ENERGY"])

predictions = predictor.predict(sequences=df_133['MODIFIED_SEQUENCE'].values.tolist(),
                                charges=df_133["PRECURSOR_CHARGE"].values.tolist(),
                                collision_energies=df_133["ORIG_COLLISION_ENERGY"].values/100.0,
                                                    models=[model],
                                                    disable_progress_bar=True)

for key, value in predictions[model].items():
     print(key)

df_133['PREDICTED_INTENSITY'] = predictions[model]['intensity'].tolist()
df_133['MZ'] = predictions[model]['fragmentmz'].tolist()
df_133.columns

spectra_list = []
for i in range (len(df_133)):
    mz_list = df_133.iloc[i]["MZ"]
    mz_array = np.array(mz_list)
    intensity_list = df_133.iloc[i]["PREDICTED_INTENSITY"]
    intensity_array = np.array(intensity_list)
    scan = df_133.iloc[i]["SCAN_NUMBER"]
    sequence = df_133.iloc[i]["MODIFIED_SEQUENCE"]
    title = scan
    mass = df_133.iloc[i]["MASS"]
    retentiontime = df_133.iloc[i]["RTs"]
    charge = df_133.iloc[i]["PRECURSOR_CHARGE"]
    mgf_params = dict({'title' : title, 'pepmass' : mass, 'rtinseconds' : retentiontime, 'charge' : charge})
    spectra_dict = dict({'m/z array' : mz_array, 'intensity array' : intensity_array, 'params' : mgf_params})
    spectra_list.append(spectra_dict)

mgf.write(spectra = spectra_list, output = mgf_path)

# scp -r cadams@10.152.135.57:/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/mgf/ .