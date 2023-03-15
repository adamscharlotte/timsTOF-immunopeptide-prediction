# ssh cadams@10.152.135.57
# /home/cadams/anaconda3/envs/prosit-annotate/bin/python3

import pandas as pd
import numpy as np
import os
import h5py

def int_to_onehot(charge):
    precursor_charge = np.full((6), False)
    precursor_charge[charge-1] = True
    return precursor_charge

base_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/" # nolint
tryp_val_path = base_path + "total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-validation-tryptic.csv"
non_tryp_val_path = base_path + "total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-validation-non-tryptic.csv"
tryp_train_path = base_path + "total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-train-tryptic.csv"
non_tryp_train_path = base_path + "total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-train-non-tryptic.csv"
tryp_test_path = base_path + "total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-test-tryptic.csv"
non_tryp_test_path = base_path + "total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-test-non-tryptic.csv"

hdf_val_path = base_path + "total-scan-consensus/hdf5/scan-40-ppm-calibrated_validation.hdf5"
hdf_train_path = base_path + "total-scan-consensus/hdf5/scan-40-ppm-calibrated_train.hdf5"
hdf_test_path = base_path + "total-scan-consensus/hdf5/scan-40-ppm-calibrated_test.hdf5"

calibrated_tryp_val_df = pd.read_csv(tryp_val_path)
calibrated_non_val_df = pd.read_csv(non_tryp_val_path)
calibrated_tryp_train_df = pd.read_csv(tryp_train_path)
calibrated_non_train_df = pd.read_csv(non_tryp_train_path)
calibrated_tryp_test_df = pd.read_csv(tryp_test_path)
calibrated_non_test_df = pd.read_csv(non_tryp_test_path)

validation_df = pd.concat([calibrated_tryp_val_df, calibrated_non_val_df])
train_df = pd.concat([calibrated_tryp_train_df, calibrated_non_train_df])
test_df = pd.concat([calibrated_tryp_test_df, calibrated_non_test_df])

def generate_hdf5(set_df, hdf5_path):
    set_df.INTENSITIES = set_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
    set_df.MZ = set_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
    set_df.OBS_SEQUENCE_INT = set_df.OBS_SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])
    set_df.SEQUENCE_INT = set_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])
    set_df['PRECURSOR_CHARGE'] = set_df['PRECURSOR_CHARGE'].astype(int)
    
    set_df["PRECURSOR_CHARGE_ONEHOT"] = set_df.apply(lambda row: int_to_onehot(row.PRECURSOR_CHARGE), axis = 1)

    collision_energies = np.asarray(set_df.ORIG_COLLISION_ENERGY.values)
    adjusted_collision_energies = np.asarray(set_df.aligned_collision_energy.values/100.0)
    collision_energies_normalized = np.asarray(set_df.ORIG_COLLISION_ENERGY.values/100.0)
    intensities = np.asarray(set_df.INTENSITIES.values.tolist())
    masses = np.asarray(set_df.MZ.values.tolist())
    methods = np.asarray(set_df.FRAGMENTATION.values,dtype ='S3')
    precursor_charges = np.asarray(set_df.PRECURSOR_CHARGE_ONEHOT.values.tolist())
    raw_files = np.asarray(set_df.RAW_FILE.values,dtype='S120')
    scan_numbers = np.asarray(set_df.SCAN_NUMBER.values)
    scores = np.asarray(set_df.SCORE.values)
    sequence_integers = np.asarray(set_df.SEQUENCE_INT.values.tolist())
    obs_sequence_integers = np.asarray(set_df.OBS_SEQUENCE_INT.values.tolist())
    retention_time = np.asarray(set_df.median_RETENTION_TIME.values)
    if os.path.exists(hdf5_path):
        with h5py.File(hdf5_path, 'a') as hf:
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
        h5f = h5py.File(hdf5_path,'w')
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

generate_hdf5(validation_df, hdf_val_path)
generate_hdf5(train_df, hdf_train_path)
generate_hdf5(test_df, hdf_test_path)
