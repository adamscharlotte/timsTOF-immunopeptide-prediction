# ssh cadams@10.152.135.57
# /home/cadams/anaconda3/envs/prosit-annotate/bin/python3

from pickle import TRUE
import pandas as pd
import numpy as np
import os

from sqlalchemy import true

from fundamentals import constants
from fundamentals.fragments import initialize_peaks
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import h5py
import argparse, pathlib
from fundamentals.mod_string import parse_modstrings, maxquant_to_internal
from fundamentals.constants import ALPHABET
import prosit_grpc
from prosit_grpc.predictPROSIT import PROSITpredictor

def int_to_onehot(charge):
    precursor_charge = np.full((6), False)
    precursor_charge[charge-1] = True
    return precursor_charge

def get_spectral_angle(intensities):
    pred= np.array(intensities[0])
    true = np.array(intensities[1])
    epsilon = 1e-7
    list_1 = np.argwhere(true>0)
    list_2 = np.argwhere(pred>0)
    indices = np.union1d(list_1,list_2)
    pred_masked = pred[indices]
    true_masked = true[indices]
    true_masked += epsilon
    pred_masked += epsilon
    true_norm = true_masked*(1/np.sqrt(np.sum(np.square(true_masked), axis=0)))
    pred_norm = pred_masked*(1/np.sqrt(np.sum(np.square(pred_masked), axis=0)))
    product = np.sum(true_norm*pred_norm, axis=0)
    arccos = np.arccos(product)
    return 1-2*arccos/np.pi

parser = argparse.ArgumentParser()
parser.add_argument('pool', type=str)					# Filename
args = parser.parse_args()

# pool = "TUM_HLA2_7"
pool = args.pool
base_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/" # nolint
annot_path = base_path + "total-scan-consensus/annotated-40-ppm/" + pool + ".csv"
ce_sa_path = base_path + "total-scan-consensus/spectral-angle-comparison/" + pool

full_df = pd.read_csv(annot_path)
full_df.INTENSITIES = full_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
full_df.MZ = full_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
full_df.SEQUENCE_INT = full_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

full_df.rename(columns = {"median_CE": "ORIG_COLLISION_ENERGY"}, inplace=True)
full_df.columns

# Filter based on score
col_filter = ['PRECURSOR_CHARGE']
full_df[col_filter] = full_df[full_df[col_filter] <= 3][col_filter]
filtered_annot_df = full_df.dropna()
print("How many PSMs with charge 4?", len(full_df.index) - len(filtered_annot_df.index))

# -----------------------------------------------------------------------------

charges = filtered_annot_df['PRECURSOR_CHARGE'].unique()
number_of_charges = charges.__len__()

filtered_annot_df['PRECURSOR_CHARGE'].value_counts()
filtered_annot_df['SCORE'].unique()
filtered_annot_df['ORIG_COLLISION_ENERGY'].mean()

# predictor = PROSITpredictor(server="131.159.152.7:8500")
predictor = PROSITpredictor(server="10.152.135.57:8500")

appended_data = []

predictions = predictor.predict(sequences=filtered_annot_df['MODIFIED_SEQUENCE'].values.tolist(),
                                charges=filtered_annot_df["PRECURSOR_CHARGE"].values.tolist(),
                                collision_energies=filtered_annot_df["ORIG_COLLISION_ENERGY"].values/100.0,
                                                    models=['Prosit_2020_intensity_hcd'],
                                                    disable_progress_bar=True)
filtered_annot_df['PREDICTED_INTENSITY'] = predictions['Prosit_2020_intensity_hcd']['intensity'].tolist()
filtered_annot_df["SPECTRAL_ANGLE"] = filtered_annot_df[['INTENSITIES','PREDICTED_INTENSITY']].apply(lambda x : get_spectral_angle(x), axis=1)
filtered_annot_df["SPECTRAL_ANGLE"].fillna(0, inplace=True)
filtered_annot_df.to_csv(ce_sa_path + '.csv')

filtered_annot_df['SPECTRAL_ANGLE'].mean()

