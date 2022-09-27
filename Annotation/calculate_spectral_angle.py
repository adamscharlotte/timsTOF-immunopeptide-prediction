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

# parser = argparse.ArgumentParser()
# parser.add_argument('pool', type=str)					# Filename
# args = parser.parse_args()

pool = "TUM_HLA2_7"

base_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/" # nolint
# annot_path = base_path + "precursor-consensus/annotated/" + args.pool + ".csv"
# sum_path = base_path + "precursor-consensus/summed/" + args.pool + ".csv"
sum_path = base_path + "precursor-consensus/annotated/" + pool + ".csv"
# sum_path = base_path + "precursor-consensus/summed/" + pool + ".csv"
unsummed_path = base_path + "full-truncated-qc/annotated/" + pool + ".csv"
calibrated_path = base_path + "full-truncated-qc/calibrated/" + pool + ".csv"
calibrated_sum_path = base_path + "precursor-consensus/calibrated/" + pool + ".csv"
sa_path = base_path + "spectral-angle/"

full_df = pd.read_csv(annot_path)
sum_df = pd.read_csv(sum_path)
unsummed_df = pd.read_csv(unsummed_path)
cal_sum_df = pd.read_csv(calibrated_sum_path)

full_df.INTENSITIES = full_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
full_df.MZ = full_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
full_df.SEQUENCE_INT = full_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

# un_annot_df_combined.INTENSITIES = un_annot_df_combined.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
# un_annot_df_combined.MZ = un_annot_df_combined.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
# un_annot_df_combined.SEQUENCE_INT = un_annot_df_combined.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

sum_df.INTENSITIES = sum_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
sum_df.MZ = sum_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
sum_df.SEQUENCE_INT = sum_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

cal_sum_df.INTENSITIES = cal_sum_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
cal_sum_df.MZ = cal_sum_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
cal_sum_df.SEQUENCE_INT = cal_sum_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

# ------------------------------------------- SUMMED --------------------------------------------

col_filter = ['SCORE']
cal_sum_df[col_filter] = cal_sum_df[cal_sum_df[col_filter] >= 70][col_filter]
filtered_cal_df = cal_sum_df.dropna()

filtered_cal_df.columns

# cal_df = filtered_cal_df[['RAW_FILE', 'SCAN_NUMBER', 'MODIFIED_SEQUENCE', 'PRECURSOR_CHARGE', 'PRECURSOR', 'aligned_collision_energy']].drop_duplicates()
# un_cal_df = filtered_un_cal_df[['RAW_FILE', 'SCAN_NUMBER', 'MODIFIED_SEQUENCE', 'PRECURSOR_CHARGE', 'PRECURSOR', 'COLLISION_ENERGY']].drop_duplicates()
# result = pd.concat([cal_df, un_cal_df], axis=1)
# result_na = result[result.isna().any(axis=1)]
# filtered_annot_df['PRECURSOR_CHARGE'].value_counts()
# filtered_annot_df['SCORE'].unique()
# filtered_annot_df['COLLISION_ENERGY'].min()

grouped_charge_df = filtered_cal_df.groupby('PRECURSOR_CHARGE')

predictor = PROSITpredictor(server="10.152.135.57:8500")

for charge, df_charge in grouped_charge_df:
    random_10_precursors = df_charge['SEQUENCE'].drop_duplicates().sample(n = 10, random_state = 43)
    df_random_10 = df_charge[df_charge['SEQUENCE'].isin(random_10_precursors)]
    df_random_10['SEQUENCE'].drop_duplicates()
    predictions = predictor.predict(sequences=df_random_10['MODIFIED_SEQUENCE'].values.tolist(),
                                    charges=df_random_10["PRECURSOR_CHARGE"].values.tolist(),
                                    collision_energies=df_random_10["COLLISION_ENERGY"].values/100.0,
                                                        models=['Prosit_2020_intensity_hcd'],
                                                        disable_progress_bar=True)
    df_random_10['PREDICTED_INTENSITY'] = predictions['Prosit_2020_intensity_hcd']['intensity'].tolist()
    df_random_10["SPECTRAL_ANGLE"] = df_random_10[['INTENSITIES','PREDICTED_INTENSITY']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    df_random_10.to_csv(sa_path + "precursors/un-calibrated/" + pool + '_' + str(charge) + '.csv')

# ----------------------------------- SUMMED CALIBRATED -----------------------------------

col_filter = ['SCORE']
cal_sum_df[col_filter] = cal_sum_df[cal_sum_df[col_filter] >= 70][col_filter]
filtered_annot_df = cal_sum_df.dropna()

# filtered_annot_df.columns
# filtered_annot_df['PRECURSOR_CHARGE'].value_counts()
# filtered_annot_df['SCORE'].unique()
# # filtered_annot_df['COLLISION_ENERGY'].min()
# filtered_annot_df['aligned_collision_energy'].min()

# col_filter = ['PRECURSOR_CHARGE']
# filtered_annot_df[col_filter] = filtered_annot_df[filtered_annot_df[col_filter] <= 3][col_filter]
# filfiltered_annot_df = filtered_annot_df.dropna()
# grouped_charge_df = filfiltered_annot_df.groupby('PRECURSOR_CHARGE')
grouped_charge_df = filtered_annot_df.groupby('PRECURSOR_CHARGE')

predictor = PROSITpredictor(server="10.152.135.57:8500")

for charge, df_charge in grouped_charge_df:
    random_10_precursors = df_charge['SEQUENCE'].drop_duplicates().sample(n = 10, random_state = 43)
    df_random_10 = df_charge[df_charge['SEQUENCE'].isin(random_10_precursors)]
    df_random_10["SEQUENCE"].drop_duplicates()
    predictions = predictor.predict(sequences=df_random_10['MODIFIED_SEQUENCE'].values.tolist(),
                                    charges=df_random_10["PRECURSOR_CHARGE"].values.tolist(),
                                    collision_energies=df_random_10["aligned_collision_energy"].values/100.0,
                                                        models=['Prosit_2020_intensity_hcd'],
                                                        disable_progress_bar=True)
    df_random_10['PREDICTED_INTENSITY'] = predictions['Prosit_2020_intensity_hcd']['intensity'].tolist()
    df_random_10["SPECTRAL_ANGLE"] = df_random_10[['INTENSITIES','PREDICTED_INTENSITY']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    df_random_10.to_csv(sa_path + "precursors/calibrated/" + pool + '_' + str(charge) + '.csv')
