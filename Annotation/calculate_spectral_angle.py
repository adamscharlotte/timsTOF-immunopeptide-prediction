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
    indices = np.union1d(list_1,list_1)
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

pool = args.pool
# pool = "TUM_HLA2_7"

base_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/" # nolint

calibrated_unsummed_path = base_path + "full-truncated-qc/calibrated-20ppm/" + pool + ".csv"
calibrated_sum_path = base_path + "precursor-consensus/calibrated-20ppm/" + pool + ".csv"
calibrated_scan_path = base_path + "precursor-consensus/calibrated-20ppm/" + pool + ".csv"
sa_path = base_path + "spectral-angle/"

unsummed_df = pd.read_csv(calibrated_unsummed_path)
sum_df = pd.read_csv(calibrated_sum_path)
scan_df = pd.read_csv(calibrated_scan_path)

unsummed_df.INTENSITIES = unsummed_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
unsummed_df.MZ = unsummed_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
unsummed_df.SEQUENCE_INT = unsummed_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

sum_df.INTENSITIES = sum_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
sum_df.MZ = sum_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
sum_df.SEQUENCE_INT = sum_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

scan_df.INTENSITIES = scan_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
scan_df.MZ = scan_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
scan_df.SEQUENCE_INT = scan_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

# ------------------------------------------- FRAME ---------------------------------------------

col_filter = ['SCORE']
unsummed_df[col_filter] = unsummed_df[unsummed_df[col_filter] >= 70][col_filter]
filtered_unsummed_df = unsummed_df.dropna()

filtered_unsummed_df.columns

filtered_unsummed_df['PRECURSOR_CHARGE'].value_counts()
filtered_unsummed_df['SCORE'].unique()
filtered_unsummed_df['ORIG_COLLISION_ENERGY'].min()
filtered_unsummed_df['aligned_collision_energy'].min()

grouped_charge_df = filtered_unsummed_df.groupby('PRECURSOR_CHARGE')

predictor = PROSITpredictor(server="10.152.135.57:8500")

for charge, df_charge in grouped_charge_df:
    random_10_precursors = df_charge['PRECURSOR'].drop_duplicates().sample(n = 10, random_state = 43)
    df_random_10 = df_charge[df_charge['PRECURSOR'].isin(random_10_precursors)]
    df_random_10['PRECURSOR'].drop_duplicates()
    predictions = predictor.predict(sequences=df_random_10['MODIFIED_SEQUENCE'].values.tolist(),
                                    charges=df_random_10["PRECURSOR_CHARGE"].values.tolist(),
                                    collision_energies=df_random_10["ORIG_COLLISION_ENERGY"].values/100.0,
                                                        models=['Prosit_2020_intensity_hcd'],
                                                        disable_progress_bar=True)
    df_random_10['PREDICTED_INTENSITY'] = predictions['Prosit_2020_intensity_hcd']['intensity'].tolist()
    df_random_10["SPECTRAL_ANGLE"] = df_random_10[['INTENSITIES','PREDICTED_INTENSITY']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    predictions_cal = predictor.predict(sequences=df_random_10['MODIFIED_SEQUENCE'].values.tolist(),
                                    charges=df_random_10["PRECURSOR_CHARGE"].values.tolist(),
                                    collision_energies=df_random_10["aligned_collision_energy"].values/100.0,
                                                        models=['Prosit_2020_intensity_hcd'],
                                                        disable_progress_bar=True)
    df_random_10['PREDICTED_INTENSITY_CAL'] = predictions_cal['Prosit_2020_intensity_hcd']['intensity'].tolist()
    df_random_10["SPECTRAL_ANGLE_CAL"] = df_random_10[['INTENSITIES','PREDICTED_INTENSITY_CAL']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10["SPECTRAL_ANGLE_CAL"].fillna(0, inplace=True)
    df_random_10.to_csv(sa_path + "frames/20-ppm/" + pool + '_' + str(charge) + '.csv')

# ------------------------------------------ PRECURSOR ------------------------------------------

col_filter = ['SCORE']
sum_df[col_filter] = sum_df[sum_df[col_filter] >= 70][col_filter]
filtered_sum_df = sum_df.dropna()

filtered_sum_df.columns

filtered_sum_df['PRECURSOR_CHARGE'].value_counts()
filtered_sum_df['SCORE'].unique()
filtered_sum_df['ORIG_COLLISION_ENERGY'].min()
filtered_sum_df['aligned_collision_energy'].min()

grouped_charge_df = filtered_sum_df.groupby('PRECURSOR_CHARGE')

predictor = PROSITpredictor(server="10.152.135.57:8500")

for charge, df_charge in grouped_charge_df:
    random_10_precursors = df_charge['SEQUENCE'].drop_duplicates().sample(n = 10, random_state = 43)
    df_random_10 = df_charge[df_charge['SEQUENCE'].isin(random_10_precursors)]
    df_random_10['SEQUENCE'].drop_duplicates()
    predictions = predictor.predict(sequences=df_random_10['MODIFIED_SEQUENCE'].values.tolist(),
                                    charges=df_random_10["PRECURSOR_CHARGE"].values.tolist(),
                                    collision_energies=df_random_10["ORIG_COLLISION_ENERGY"].values/100.0,
                                                        models=['Prosit_2020_intensity_hcd'],
                                                        disable_progress_bar=True)
    df_random_10['PREDICTED_INTENSITY'] = predictions['Prosit_2020_intensity_hcd']['intensity'].tolist()
    df_random_10["SPECTRAL_ANGLE"] = df_random_10[['INTENSITIES','PREDICTED_INTENSITY']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    predictions_cal = predictor.predict(sequences=df_random_10['MODIFIED_SEQUENCE'].values.tolist(),
                                    charges=df_random_10["PRECURSOR_CHARGE"].values.tolist(),
                                    collision_energies=df_random_10["aligned_collision_energy"].values/100.0,
                                                        models=['Prosit_2020_intensity_hcd'],
                                                        disable_progress_bar=True)
    df_random_10['PREDICTED_INTENSITY_CAL'] = predictions_cal['Prosit_2020_intensity_hcd']['intensity'].tolist()
    df_random_10["SPECTRAL_ANGLE_CAL"] = df_random_10[['INTENSITIES','PREDICTED_INTENSITY_CAL']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10["SPECTRAL_ANGLE_CAL"].fillna(0, inplace=True)
    df_random_10.to_csv(sa_path + "precursors/20-ppm/" + pool + '_' + str(charge) + '.csv')

# -------------------------------------------- SCAN ---------------------------------------------

col_filter = ['SCORE']
scan_df[col_filter] = scan_df[scan_df[col_filter] >= 70][col_filter]
filtered_scan_df = scan_df.dropna()

filtered_scan_df.columns

filtered_scan_df['PRECURSOR_CHARGE'].value_counts()
filtered_scan_df['SCORE'].unique()
filtered_scan_df['ORIG_COLLISION_ENERGY'].min()
filtered_scan_df['aligned_collision_energy'].min()

grouped_charge_df = filtered_scan_df.groupby('PRECURSOR_CHARGE')

predictor = PROSITpredictor(server="10.152.135.57:8500")

for charge, df_charge in grouped_charge_df:
    random_10_precursors = df_charge['SEQUENCE'].drop_duplicates().sample(n = 10, random_state = 43)
    df_random_10 = df_charge[df_charge['SEQUENCE'].isin(random_10_precursors)]
    df_random_10['SEQUENCE'].drop_duplicates()
    predictions = predictor.predict(sequences=df_random_10['MODIFIED_SEQUENCE'].values.tolist(),
                                    charges=df_random_10["PRECURSOR_CHARGE"].values.tolist(),
                                    collision_energies=df_random_10["ORIG_COLLISION_ENERGY"].values/100.0,
                                                        models=['Prosit_2020_intensity_hcd'],
                                                        disable_progress_bar=True)
    df_random_10['PREDICTED_INTENSITY'] = predictions['Prosit_2020_intensity_hcd']['intensity'].tolist()
    df_random_10["SPECTRAL_ANGLE"] = df_random_10[['INTENSITIES','PREDICTED_INTENSITY']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    predictions_cal = predictor.predict(sequences=df_random_10['MODIFIED_SEQUENCE'].values.tolist(),
                                    charges=df_random_10["PRECURSOR_CHARGE"].values.tolist(),
                                    collision_energies=df_random_10["aligned_collision_energy"].values/100.0,
                                                        models=['Prosit_2020_intensity_hcd'],
                                                        disable_progress_bar=True)
    df_random_10['PREDICTED_INTENSITY_CAL'] = predictions_cal['Prosit_2020_intensity_hcd']['intensity'].tolist()
    df_random_10["SPECTRAL_ANGLE_CAL"] = df_random_10[['INTENSITIES','PREDICTED_INTENSITY_CAL']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10["SPECTRAL_ANGLE_CAL"].fillna(0, inplace=True)
    df_random_10.to_csv(sa_path + "scans/20-ppm/" + pool + '_' + str(charge) + '.csv')

