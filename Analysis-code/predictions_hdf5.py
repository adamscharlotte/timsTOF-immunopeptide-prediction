# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction
20230315104417_prediction.hdf5  20230315111004_prediction.hdf5

import h5py
import pandas as pd
import numpy as np

file_name_1 = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/20230315104417_prediction.hdf5"
file_name_2 = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/20230315111004_prediction.hdf5"
file_name_3 = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/20230316105826_prediction.hdf5"
file_name_4 = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/20230316105504_prediction.hdf5"

def read_hdf5_to_dataframe(file_path, keys=None):
    """
    Reads data from an HDF5 file and returns a Pandas DataFrame.
    Parameters:
    file_path (str): The path to the HDF5 file.
    keys (list or None): A list of keys to include as columns. If None, include all keys.
    Returns:
    A Pandas DataFrame containing the data.
    """
    with h5py.File(file_path, 'r') as f:
        # Create an empty dictionary to hold the data
        data_dict = {}
        # Get the length of the first dataset
        length = f[list(f.keys())[0]].shape[0]
        # Loop through the keys in the HDF5 file and add them as columns to the dictionary
        for key in f.keys():
            # Check if the dataset has the same length as the first dataset
            if f[key].shape[0] == length:
                # Check if the key is in the list of keys to include
                if keys is None or key in keys:
                    # Add the dataset to the dictionary
                    data_dict[key] = f[key][:]
    # Convert the dictionary to a Pandas DataFrame
    df = pd.DataFrame(data_dict)
    return df

df_1 = read_hdf5_to_dataframe(file_name_1, ['spectral_angle', 'scan_number', 'collision_energy_aligned_normed', 'score',
    'rawfile', 'retention_time', 'collision_energy_normed', 'collision_energy'])

df_2 = read_hdf5_to_dataframe(file_name_2, ['spectral_angle', 'scan_number', 'collision_energy_aligned_normed', 'score',
    'rawfile', 'retention_time', 'collision_energy_normed', 'collision_energy'])

df_3 = read_hdf5_to_dataframe(file_name_3, ['spectral_angle', 'scan_number', 'collision_energy_aligned_normed', 'score',
    'rawfile', 'retention_time', 'collision_energy_normed', 'collision_energy'])

df_4 = read_hdf5_to_dataframe(file_name_4, ['spectral_angle', 'scan_number', 'collision_energy_aligned_normed', 'score',
    'rawfile', 'retention_time', 'collision_energy_normed', 'collision_energy'])

file_out_1 = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/20230315104417_prediction.csv"
file_out_2 = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/20230315111004_prediction.csv"
file_out_3 = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/20230316105826_prediction.csv"
file_out_4 = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/20230316105504_prediction.csv"

df_1.to_csv(file_out_1, index=False)
df_2.to_csv(file_out_2, index=False)
df_3.to_csv(file_out_3, index=False)
df_4.to_csv(file_out_4, index=False)