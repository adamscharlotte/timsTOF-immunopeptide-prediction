# /Users/adams/opt/miniconda3/envs/ann_solo/bin/python3

# from pickle import FRAME
from re import I
import pandas as pd
import numpy as np
import os
from pyteomics import mgf

# import argparse, pathlib
# parser = argparse.ArgumentParser()
# parser.add_argument('pool', type=str)					# Filename
# args = parser.parse_args()

# pool = args.pool
pool = "TUM_HLA2_88"

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
# csv_path = base_path + "full-truncated-qc/un-annotated/" + pool + ".csv"
# mgf_path = base_path + "full-truncated-qc/un-annotated-mgf/" + pool + ".mgf"
csv_path = base_path + "full-truncated-qc/annotated/" + pool + ".csv"
mgf_path = base_path + "full-truncated-qc/annotated-mgf/" + pool + ".mgf"

csv_df = pd.read_csv(csv_path)

csv_df.INTENSITIES = csv_df.INTENSITIES.str.split(';').apply(lambda s: [float(x) for x in s])
csv_df.MZ = csv_df.MZ.str.split(';').apply(lambda s: [float(x) for x in s])
csv_df.columns

spectra_list = []
for i in range (len(csv_df)):
    mz_list = csv_df.iloc[i]["MZ"]
    mz_array = np.array(mz_list)
    intensity_list = csv_df.iloc[i]["INTENSITIES"]
    intensity_array = np.array(intensity_list).astype(int)
    frame = csv_df.iloc[i]["FRAME"].astype(str)
    precursor = csv_df.iloc[i]["PRECURSOR"].astype(str)
    title = frame + "_" + precursor
    mass = csv_df.iloc[i]["MASS"]
    retentiontime = csv_df.iloc[i]["RETENTION_TIME"]
    # charge = csv_df.iloc[i]["CHARGE"]
    charge = csv_df.iloc[i]["PRECURSOR_CHARGE"]
    mgf_params = dict({'title' : title, 'pepmass' : mass, 'rtinseconds' : retentiontime, 'charge' : charge})
    spectra_dict = dict({'m/z array' : mz_array, 'intensity array' : intensity_array, 'params' : mgf_params})
    spectra_list.append(spectra_dict)

mgf.write(spectra = spectra_list, output = mgf_path)

# ----------------------------------- PRECURSOR -----------------------------------
# csv_path = base_path + "precursor-consensus/un-annotated/" + pool + ".csv"
# mgf_path = base_path + "precursor-consensus/un-annotated-mgf/" + pool + ".mgf"
csv_path = base_path + "precursor-consensus/annotated/" + pool + ".csv"
mgf_path = base_path + "precursor-consensus/annotated-mgf/" + pool + ".mgf"

csv_df = pd.read_csv(csv_path)

# csv_df["INTENSITIES"] = csv_df.combined_INTENSITIES.str.split(';').apply(lambda s: [int(x) for x in s])
# csv_df["MZ"] = csv_df.combined_MZ.str.split(';').apply(lambda s: [float(x) for x in s])

csv_df.INTENSITIES = csv_df.INTENSITIES.str.split(';').apply(lambda s: [float(x) for x in s])
csv_df.MZ = csv_df.MZ.str.split(';').apply(lambda s: [float(x) for x in s])

spectra_list = []
for i in range (len(csv_df)):
    mz_list = csv_df.iloc[i]["MZ"]
    mz_array = np.array(mz_list)
    intensity_list = csv_df.iloc[i]["INTENSITIES"]
    intensity_array = np.array(intensity_list).astype(int)
    precursor = csv_df.iloc[i]["PRECURSOR"].astype(str)
    title = precursor
    mass = csv_df.iloc[i]["MASS"]
    retentiontime = csv_df.iloc[i]["RETENTION_TIME"]
    # charge = csv_df.iloc[i]["CHARGE"]
    charge = csv_df.iloc[i]["PRECURSOR_CHARGE"]
    mgf_params = dict({'title' : title, 'pepmass' : mass, 'rtinseconds' : retentiontime, 'charge' : charge})
    spectra_dict = dict({'m/z array' : mz_array, 'intensity array' : intensity_array, 'params' : mgf_params})
    spectra_list.append(spectra_dict)

mgf.write(spectra = spectra_list, output = mgf_path)