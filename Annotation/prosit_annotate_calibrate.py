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

base_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/" # nolint
un_annot_path = base_path + "full-truncated-qc/un-annotated-ce/" + args.pool + ".csv"
# un_annot_path = base_path + "full-truncated-qc/un-annotated-ce/" + pool + ".csv"
ce_sa_path = base_path + "full-truncated-qc/spectral-angle/" + args.pool 
# ce_sa_path = base_path + "full-truncated-qc/spectral-angle/" + pool

un_annot_df = pd.read_csv(un_annot_path)
file_path = base_path + "full-truncated-qc/annotated/full-truncated-qc-calibrated.hdf5"

un_annot_df.INTENSITIES = un_annot_df.INTENSITIES.str.split(';').apply(lambda s: [float(x) for x in s])
un_annot_df.MZ = un_annot_df.MZ.str.split(';').apply(lambda s: [float(x) for x in s])

un_annot_df.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)
un_annot_df.rename(columns = {"COLLISION_ENERGY": "ORIG_COLLISION_ENERGY"}, inplace=True)
un_annot_df["REVERSE"].fillna(False, inplace=True)
un_annot_df["REVERSE"].replace("+", True, inplace=True)
un_annot_df["MODIFIED_SEQUENCE"] = maxquant_to_internal(un_annot_df["MODIFIED_SEQUENCE"].to_numpy())
un_annot_df["SEQUENCE"] = internal_without_mods(un_annot_df["MODIFIED_SEQUENCE"])
un_annot_df['PEPTIDE_LENGTH'] = un_annot_df["SEQUENCE"].apply(lambda x: len(x))
un_annot_df["PRECURSOR_CHARGE_ONEHOT"] = un_annot_df.apply(lambda row: int_to_onehot(row.PRECURSOR_CHARGE), axis = 1)

generator_sequence_numeric = parse_modstrings(list(un_annot_df['MODIFIED_SEQUENCE'].values), ALPHABET, translate=True)
enum_gen_seq_num = enumerate(generator_sequence_numeric)
array = np.zeros((len(list(un_annot_df['MODIFIED_SEQUENCE'].values)),30), dtype=np.uint8)
for i, sequence_numeric in enum_gen_seq_num:
    if len(sequence_numeric) > 30:
        if filter:
            pass # don't overwrite 0 in the array that is how we can differentiate
        else:
            raise Exception(f"The Sequence {sequence_numeric}, has {len(sequence_numeric)} Amino Acids."
                        f"The maximum number of amino acids allowed is {C.SEQ_LEN}")
    else:
        array[i, 0:len(sequence_numeric)] = sequence_numeric

un_annot_df['SEQUENCE_INT'] = array.tolist()

annot_df = annotate_spectra(un_annot_df)

full_df = pd.concat([un_annot_df.drop(columns = ["INTENSITIES", "MZ"]), annot_df], axis=1)

# Filter based on score
col_filter = ['SCORE']
full_df[col_filter] = full_df[un_annot_df[col_filter] >= 70][col_filter]
filtered_annot_df = full_df.dropna()

# -----------------------------------------------------------------------------

charges = filtered_annot_df['PRECURSOR_CHARGE'].unique()
number_of_charges = charges.__len__()

filtered_annot_df['PRECURSOR_CHARGE'].value_counts()
filtered_annot_df['SCORE'].unique()
filtered_annot_df['ORIG_COLLISION_ENERGY'].min()

grouped_charge_df = filtered_annot_df.groupby('PRECURSOR_CHARGE')

predictor = PROSITpredictor(server="10.152.135.57:8500")
CE_RANGE = range(5, 45)
appended_data = []

for charge, df_charge in grouped_charge_df:
    top_100_df = df_charge.sort_values(['SCORE'], ascending=False).head(100)
    nrow = len(top_100_df)
    top_100_df = pd.concat([top_100_df for _ in CE_RANGE], axis=0)
    top_100_df["COLLISION_ENERGY"] = np.repeat(CE_RANGE, nrow)
    top_100_df.reset_index(inplace=True)
    
    predictions = predictor.predict(sequences=top_100_df['MODIFIED_SEQUENCE'].values.tolist(),
                                    charges=top_100_df["PRECURSOR_CHARGE"].values.tolist(),
                                    collision_energies=top_100_df["COLLISION_ENERGY"].values/100.0,
                                                        models=['Prosit_2020_intensity_hcd'],
                                                        disable_progress_bar=True)
    
    top_100_df['PREDICTED_INTENSITY'] = predictions['Prosit_2020_intensity_hcd']['intensity'].tolist()
    top_100_df["SPECTRAL_ANGLE"] = top_100_df[['INTENSITIES','PREDICTED_INTENSITY']].apply(lambda x : get_spectral_angle(x), axis=1)
    top_100_df["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    top_100_df.to_csv(ce_sa_path + '_' + str(charge) + '.csv')
    groups = top_100_df.groupby(by=['ORIG_COLLISION_ENERGY', "COLLISION_ENERGY"])["SPECTRAL_ANGLE"].mean()
    groups_2 = groups.reset_index()
    ids = groups_2.groupby(['ORIG_COLLISION_ENERGY'])['SPECTRAL_ANGLE'].transform(max) == groups_2['SPECTRAL_ANGLE']
    calib_group = groups_2[ids]
    calib_group['delta_collision_energy'] = calib_group['COLLISION_ENERGY'] - calib_group['ORIG_COLLISION_ENERGY']
    ce_alignment = calib_group['delta_collision_energy']
    best_ce = ce_alignment.median()
    best_ce
    df_charge['aligned_collision_energy'] = df_charge['ORIG_COLLISION_ENERGY'] + best_ce
    appended_data.append(df_charge)

calibrated_annot_df = pd.concat(appended_data)

# -----------------------------------------------------------------------------

calibrated_annot_df['aligned_collision_energy'] = calibrated_annot_df['aligned_collision_energy'].apply(lambda x : float(round(x)))

collision_energies = np.asarray(calibrated_annot_df.ORIG_COLLISION_ENERGY.values)
adjusted_collision_energies = np.asarray(calibrated_annot_df.aligned_collision_energy.values/100.0)
collision_energies_normalized = np.asarray(calibrated_annot_df.ORIG_COLLISION_ENERGY.values/100.0)
intensities = np.asarray(calibrated_annot_df.INTENSITIES.values.tolist())
masses = np.asarray(calibrated_annot_df.MZ.values.tolist())
methods = np.asarray(calibrated_annot_df.FRAGMENTATION.values,dtype ='S3')
precursor_charges = np.asarray(calibrated_annot_df.PRECURSOR_CHARGE_ONEHOT.values.tolist())
raw_files = np.asarray(calibrated_annot_df.RAW_FILE.values,dtype='S120')
scan_numbers = np.asarray(calibrated_annot_df.SCAN_NUMBER.values)
scores = np.asarray(calibrated_annot_df.SCORE.values)
sequence_integers = np.asarray(calibrated_annot_df.SEQUENCE_INT.values.tolist())
retention_time = np.asarray(calibrated_annot_df.RETENTION_TIME.values)

if os.path.exists(file_path):
    with h5py.File(file_path, 'a') as hf:
        hf["collision_energy"].resize((hf["collision_energy"].shape[0] + collision_energies.shape[0]), axis = 0)
        hf["collision_energy"][-collision_energies.shape[0]:] = collision_energies
        
        hf["collision_energy_aligned_normed"].resize((hf["collision_energy_aligned_normed"].shape[0] + adjusted_collision_energies.shape[0]), axis = 0)
        hf["collision_energy_aligned_normed"][-adjusted_collision_energies.shape[0]:] = adjusted_collision_energies
        
        hf["collision_energy_normed"].resize((hf["collision_energy_normed"].shape[0] + collision_energies_normalized.shape[0]), axis = 0)
        hf["collision_energy_normed"][-collision_energies_normalized.shape[0]:] = collision_energies_normalized
        
        hf["intensities_raw"].resize((hf["intensities_raw"].shape[0] + intensities.shape[0]), axis = 0)
        hf["intensities_raw"][-intensities.shape[0]:] = intensities
        
        hf["masses_raw"].resize((hf["masses_raw"].shape[0] + masses.shape[0]), axis = 0)
        hf["masses_raw"][-masses.shape[0]:] = masses
        
        hf["method"].resize((hf["method"].shape[0] + methods.shape[0]), axis = 0)
        hf["method"][-methods.shape[0]:] = methods
        
        hf["precursor_charge_onehot"].resize((hf["precursor_charge_onehot"].shape[0] + precursor_charges.shape[0]), axis = 0)
        hf["precursor_charge_onehot"][-precursor_charges.shape[0]:] = precursor_charges
        
        hf["rawfile"].resize((hf["rawfile"].shape[0] + raw_files.shape[0]), axis = 0)
        hf["rawfile"][-raw_files.shape[0]:] = raw_files
        
        hf["scan_number"].resize((hf["scan_number"].shape[0] + scan_numbers.shape[0]), axis = 0)
        hf["scan_number"][-scan_numbers.shape[0]:] = scan_numbers
        
        hf["score"].resize((hf["score"].shape[0] + scores.shape[0]), axis = 0)
        hf["score"][-scores.shape[0]:] = scores
        
        hf["sequence_integer"].resize((hf["sequence_integer"].shape[0] + sequence_integers.shape[0]), axis = 0)
        hf["sequence_integer"][-sequence_integers.shape[0]:] = sequence_integers
        
        hf["retention_time"].resize((hf["retention_time"].shape[0] + retention_time.shape[0]), axis = 0)
        hf["retention_time"][-retention_time.shape[0]:] = retention_time
else:
    h5f = h5py.File(file_path,'w')
    h5f.create_dataset('collision_energy', data=collision_energies, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset('collision_energy_aligned_normed', data=adjusted_collision_energies, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset('collision_energy_normed', data=collision_energies_normalized, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset('intensities_raw', data=intensities, compression="gzip", chunks=True, maxshape=(None,174)) 
    h5f.create_dataset('masses_raw', data=masses, compression="gzip", chunks=True, maxshape=(None,174)) 
    h5f.create_dataset('method', data=methods, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset('precursor_charge_onehot', data=precursor_charges, compression="gzip", chunks=True, maxshape=(None,6)) 
    h5f.create_dataset('rawfile', data=raw_files, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset('scan_number', data=scan_numbers, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset('score', data=scores, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset('sequence_integer', data=sequence_integers, compression="gzip", chunks=True, maxshape=(None,32))
    h5f.create_dataset('retention_time', data=retention_time, compression="gzip", chunks=True, maxshape=(None,))
    h5f.close()
