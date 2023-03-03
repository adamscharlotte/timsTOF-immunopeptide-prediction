# ssh cadams@10.152.135.57
# /home/cadams/anaconda3/envs/prosit-annotate/bin/python3

from pickle import TRUE
import pandas as pd
import numpy as np
import os

from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import HuberRegressor
import matplotlib.pyplot as plt
from sqlalchemy import true
import prosit_grpc
from prosit_grpc.predictPROSIT import PROSITpredictor
from plotnine import *

from fundamentals import constants
from fundamentals.fragments import initialize_peaks
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import argparse, pathlib
from fundamentals.mod_string import parse_modstrings, maxquant_to_internal
from fundamentals.constants import ALPHABET

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

base_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/" # nolint
train_path = base_path + "total-scan-consensus/split/non-tryptic/annotated-40-ppm-train.csv"
ce_sa_path = base_path + "total-scan-consensus/ce_calibration/non-tryptic-train"
cali_path = base_path + "total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-train-non-tryptic.csv"

annot_df = pd.read_csv(train_path)
annot_df.INTENSITIES = annot_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
annot_df.MZ = annot_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
annot_df.SEQUENCE_INT = annot_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

annot_df.rename(columns = {"median_CE": "ORIG_COLLISION_ENERGY"}, inplace=True)
# full_df.columns

# -----------------------------------------------------------------------------

charges = annot_df['PRECURSOR_CHARGE'].unique()
number_of_charges = charges.__len__()

annot_df['PRECURSOR_CHARGE'].value_counts()
annot_df['SCORE'].unique()
annot_df['ORIG_COLLISION_ENERGY'].min()

grouped_charge_df = annot_df.groupby('PRECURSOR_CHARGE')

# predictor = PROSITpredictor(server="131.159.152.7:8500")
predictor = PROSITpredictor(server="10.152.135.57:8500")

CE_RANGE = range(5, 45)
appended_data = []

for charge, df_charge in grouped_charge_df:
    val_10p = round(len(df_charge)/10)
    top_10p_df = df_charge.sort_values(['SCORE'], ascending=False).head(val_10p)
    nrow = len(top_10p_df)
    top_10p_df = pd.concat([top_10p_df for _ in CE_RANGE], axis=0)
    top_10p_df["COLLISION_ENERGY"] = np.repeat(CE_RANGE, nrow)
    top_10p_df.reset_index(inplace=True)
    predictions = predictor.predict(sequences=top_10p_df['MODIFIED_SEQUENCE'].values.tolist(),
                                charges=top_10p_df["PRECURSOR_CHARGE"].values.tolist(),
                                collision_energies=top_10p_df["COLLISION_ENERGY"].values/100.0,
                                                    models=['Prosit_2020_intensity_hcd'],
                                                    disable_progress_bar=True)
    top_10p_df['PREDICTED_INTENSITY'] = predictions['Prosit_2020_intensity_hcd']['intensity'].tolist()
    top_10p_df["SPECTRAL_ANGLE"] = top_10p_df[['INTENSITIES','PREDICTED_INTENSITY']].apply(lambda x : get_spectral_angle(x), axis=1)
    top_10p_df["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    top_10p_df.to_csv(ce_sa_path + '_' + str(charge) + '.csv')
    groups = top_10p_df.groupby(by=['ORIG_COLLISION_ENERGY', "COLLISION_ENERGY", "MASS"])["SPECTRAL_ANGLE"].mean()
    groups_2 = groups.reset_index()
    ids = groups_2.groupby(['ORIG_COLLISION_ENERGY'])['SPECTRAL_ANGLE'].transform(max) == groups_2['SPECTRAL_ANGLE']
    calib_group = groups_2[ids]
    calib_group['delta_collision_energy'] = calib_group['COLLISION_ENERGY'] - calib_group['ORIG_COLLISION_ENERGY']
    # Fit a linear model
    X = calib_group[['MASS']] # input feature
    y = calib_group['delta_collision_energy'] # target variable
    model = HuberRegressor().fit(X, y)
    # Use the model to predict delta_CE for a range of masses
    min_mass = calib_group['MASS'].min()
    max_mass = calib_group['MASS'].max()
    mass_range = np.linspace(min_mass, max_mass, 100).reshape(-1, 1)
    predicted_delta_CE = model.predict(mass_range)
    r_squared = model.score(X, y)
    # Plot the model and data points
    plt.scatter(calib_group['MASS'], calib_group['delta_collision_energy'], label='Data')
    plt.plot(mass_range, predicted_delta_CE, color='red', label='Linear Model')
    plt.xlabel('Mass')
    plt.ylabel('Delta_CE')
    plt.legend()
    plt.suptitle('CE calibration for non-tryptic with charge ' + str(charge))
    plt.title(" HuberRegressor with R2:" + str(r_squared))
    plt.savefig('/home/cadams/Figures/model_non-tryptic_train' + str(charge) + '.png')
    plt.cla()
    # df_charge['delta_ce'] = model.predict(df_charge[['MASS']])
    # df_charge['aligned_collision_energy'] = df_charge['ORIG_COLLISION_ENERGY'] + df_charge['delta_ce']
    # appended_data.append(df_charge)


p = (ggplot(calib_group, aes('MASS', 'delta_collision_energy')) # , color='SPECTRAL_ANGLE'
 + geom_point()
 + geom_smooth(method='lm')
 + labs(x='MASS', y='delta_collision_energy', title='Scatter Plot with Linear Regression Line')
 )

calibrated_annot_df = pd.concat(appended_data)

calibrated_annot_df['aligned_collision_energy'] = calibrated_annot_df['aligned_collision_energy'].apply(lambda x : float(round(x)))

calibrated_annot_df_save = calibrated_annot_df
calibrated_annot_df_save["MZ"] = [';'.join(map(str, l)) for l in calibrated_annot_df_save['MZ']]
calibrated_annot_df_save["INTENSITIES"] = [';'.join(map(str, l)) for l in calibrated_annot_df_save['INTENSITIES']]
calibrated_annot_df_save.to_csv(cali_path)
