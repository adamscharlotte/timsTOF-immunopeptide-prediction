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

parser = argparse.ArgumentParser()
parser.add_argument("pool", type=str)					# Filename
args = parser.parse_args()

# pool = "TUM_HLA_16"

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
file_path = base_path + "precursor-consensus/annotated-hdf5/full-truncated-qc-un-callibrated-precursor-consensus.hdf5"
sum_path = base_path + "precursor-consensus/summed/" + args.pool + ".csv"
annot_path = base_path + "precursor-consensus/annotated/" + args.pool + ".csv"

# sum_path = base_path + "precursor-consensus/summed/" + pool + ".csv"
# annot_path = base_path + "precursor-consensus/annotated/" + pool + ".csv"

full_df = pd.read_csv(annot_path)
un_annot_df_combined = pd.read_csv(sum_path)

full_df.INTENSITIES = full_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
full_df.MZ = full_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
full_df.SEQUENCE_INT = full_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

un_annot_df_combined.INTENSITIES = un_annot_df_combined.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
un_annot_df_combined.MZ = un_annot_df_combined.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
un_annot_df_combined.SEQUENCE_INT = un_annot_df_combined.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

# ------------------------ Filter based on score ------------------------

col_filter = ["SCORE"]
full_df[col_filter] = full_df[un_annot_df_combined[col_filter] >= 70][col_filter]
filtered_annot_df = full_df.dropna()

# ---------------------------- one-hot charge ---------------------------

def int_to_onehot(charge):
    precursor_charge = np.full((6), False)
    precursor_charge[charge-1] = True
    return precursor_charge

filtered_annot_df["PRECURSOR_CHARGE_ONEHOT"] = filtered_annot_df.apply(lambda row: int_to_onehot(row.PRECURSOR_CHARGE), axis = 1)

collision_energies = np.asarray(filtered_annot_df.COLLISION_ENERGY.values)
collision_energies_normalized = np.asarray(filtered_annot_df.COLLISION_ENERGY.values/100.0)
intensities = np.asarray(filtered_annot_df.INTENSITIES.values.tolist())
masses = np.asarray(filtered_annot_df.MZ.values.tolist())
methods = np.asarray(filtered_annot_df.FRAGMENTATION.values,dtype ="S3")
precursor_charges = np.asarray(filtered_annot_df.PRECURSOR_CHARGE_ONEHOT.values.tolist())
raw_files = np.asarray(filtered_annot_df.RAW_FILE.values,dtype="S120")
scan_numbers = np.asarray(filtered_annot_df.SCAN_NUMBER.values)
scores = np.asarray(filtered_annot_df.SCORE.values)
sequence_integers = np.asarray(filtered_annot_df.SEQUENCE_INT.values.tolist())
retention_time = np.asarray(filtered_annot_df.RETENTION_TIME.values)

if os.path.exists(file_path):
    with h5py.File(file_path, "a") as hf:
        hf["collision_energy"].resize((hf["collision_energy"].shape[0] + collision_energies.shape[0]), axis = 0)
        hf["collision_energy"][-collision_energies.shape[0]:] = collision_energies
        
        # hf["collision_energy_aligned_normed"].resize((hf["collision_energy_aligned_normed"].shape[0] + adjusted_collision_energies.shape[0]), axis = 0)
        # hf["collision_energy_aligned_normed"][-adjusted_collision_energies.shape[0]:] = adjusted_collision_energies
        
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
    h5f = h5py.File(file_path,"w")
    h5f.create_dataset("collision_energy", data=collision_energies, compression="gzip", chunks=True, maxshape=(None,)) 
    # h5f.create_dataset("collision_energy_aligned_normed", data=adjusted_collision_energies, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset("collision_energy_normed", data=collision_energies_normalized, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset("intensities_raw", data=intensities, compression="gzip", chunks=True, maxshape=(None,174)) 
    h5f.create_dataset("masses_raw", data=masses, compression="gzip", chunks=True, maxshape=(None,174)) 
    h5f.create_dataset("method", data=methods, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset("precursor_charge_onehot", data=precursor_charges, compression="gzip", chunks=True, maxshape=(None,6)) 
    h5f.create_dataset("rawfile", data=raw_files, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset("scan_number", data=scan_numbers, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset("score", data=scores, compression="gzip", chunks=True, maxshape=(None,)) 
    h5f.create_dataset("sequence_integer", data=sequence_integers, compression="gzip", chunks=True, maxshape=(None,32))
    h5f.create_dataset("retention_time", data=retention_time, compression="gzip", chunks=True, maxshape=(None,))
    h5f.close()
