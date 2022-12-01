# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

import pandas as pd
import numpy as np
import os

from fundamentals import constants
from fundamentals.fragments import initialize_peaks
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import h5py
import argparse, pathlib
from fundamentals.mod_string import parse_modstrings, maxquant_to_internal
from fundamentals.constants import ALPHABET
from mgf_filter.util import timeStamped
from mgf_filter.masterSpectrum import MasterSpectrum

def int_to_onehot(charge):
    precursor_charge = np.full((6), False)
    precursor_charge[charge-1] = True
    return precursor_charge

parser = argparse.ArgumentParser()
parser.add_argument("pool", type=str)					# Filename
args = parser.parse_args()
pool = args.pool
# pool = "TUM_lysn_33"

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
calibrated_path = base_path + "scan-consensus/calibrated-40-ppm/" + pool + ".csv"
file_path = base_path + "scan-consensus/calibrated-hdf5/scan-40-ppm-calibrated-mapped.hdf5"
map_path = base_path + "full-length-map/" + pool + ".csv"

map_df = pd.read_csv(map_path)

calibrated_annot_df = pd.read_csv(calibrated_path)
calibrated_annot_df = calibrated_annot_df.rename({"SEQUENCE": "OBS_SEQUENCE", "MODIFIED_SEQUENCE": "OBS_MODIFIED_SEQUENCE", "SEQUENCE_INT": "OBS_SEQUENCE_INT"}, axis='columns')
calibrated_annot_df_merged = pd.merge(left=calibrated_annot_df, right=map_df, how='left', left_on='OBS_SEQUENCE', right_on='OBS_SEQUENCE')
calibrated_annot_df_merged.INTENSITIES = calibrated_annot_df_merged.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
calibrated_annot_df_merged.MZ = calibrated_annot_df_merged.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
calibrated_annot_df_merged.OBS_SEQUENCE_INT = calibrated_annot_df_merged.OBS_SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

calibrated_annot_df_merged = calibrated_annot_df_merged.dropna()
calibrated_annot_df_merged["MODIFIED_SEQUENCE"] = '_' + calibrated_annot_df_merged['SEQUENCE'].astype(str) + '_'
# calibrated_annot_df_merged["MODIFIED_SEQUENCE"] = calibrated_annot_df_merged['SEQUENCE'].astype(str)
calibrated_annot_df_merged["MODIFIED_SEQUENCE"] = maxquant_to_internal(calibrated_annot_df_merged["MODIFIED_SEQUENCE"].to_numpy())

generator_sequence_numeric = parse_modstrings(list(calibrated_annot_df_merged["MODIFIED_SEQUENCE"].values), ALPHABET, translate=True)
enum_gen_seq_num = enumerate(generator_sequence_numeric)
array = np.zeros((len(list(calibrated_annot_df_merged["MODIFIED_SEQUENCE"].values)),30), dtype=np.uint8)
for i, sequence_numeric in enum_gen_seq_num:
    if len(sequence_numeric) > 30:
        if filter:
            pass # don"t overwrite 0 in the array that is how we can differentiate
        else:
            raise Exception(f"The Sequence {sequence_numeric}, has {len(sequence_numeric)} Amino Acids."
                        f"The maximum number of amino acids allowed is {C.SEQ_LEN}")
    else:
        array[i, 0:len(sequence_numeric)] = sequence_numeric

calibrated_annot_df_merged["SEQUENCE_INT"] = array.tolist()

calibrated_annot_df_merged["PRECURSOR_CHARGE_ONEHOT"] = calibrated_annot_df_merged.apply(lambda row: int_to_onehot(row.PRECURSOR_CHARGE), axis = 1)

collision_energies = np.asarray(calibrated_annot_df_merged.ORIG_COLLISION_ENERGY.values)
adjusted_collision_energies = np.asarray(calibrated_annot_df_merged.aligned_collision_energy.values/100.0)
collision_energies_normalized = np.asarray(calibrated_annot_df_merged.ORIG_COLLISION_ENERGY.values/100.0)
intensities = np.asarray(calibrated_annot_df_merged.INTENSITIES.values.tolist())
masses = np.asarray(calibrated_annot_df_merged.MZ.values.tolist())
methods = np.asarray(calibrated_annot_df_merged.FRAGMENTATION.values,dtype ='S3')
precursor_charges = np.asarray(calibrated_annot_df_merged.PRECURSOR_CHARGE_ONEHOT.values.tolist())
raw_files = np.asarray(calibrated_annot_df_merged.RAW_FILE.values,dtype='S120')
scan_numbers = np.asarray(calibrated_annot_df_merged.SCAN_NUMBER.values)
scores = np.asarray(calibrated_annot_df_merged.SCORE.values)
sequence_integers = np.asarray(calibrated_annot_df_merged.SEQUENCE_INT.values.tolist())
obs_sequence_integers = np.asarray(calibrated_annot_df_merged.OBS_SEQUENCE_INT.values.tolist())
retention_time = np.asarray(calibrated_annot_df_merged.RETENTION_TIME.values)

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

        hf["obs_sequence_integer"].resize((hf["obs_sequence_integer"].shape[0] + obs_sequence_integers.shape[0]), axis = 0)
        hf["obs_sequence_integer"][-sequence_integers.shape[0]:] = obs_sequence_integers

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
    h5f.create_dataset('obs_sequence_integer', data=obs_sequence_integers, compression="gzip", chunks=True, maxshape=(None,32))
    h5f.create_dataset('retention_time', data=retention_time, compression="gzip", chunks=True, maxshape=(None,))
    h5f.close()
