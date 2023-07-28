# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import pandas as pd
import seaborn as sns
import os
import h5py

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
    # df = pd.DataFrame(data_dict)
    return data_dict

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

# ---------------------- Import data ----------------------
base_path = "/Users/adams/Projects/300K/2022-library-run/training-sets/"  # nolint
# hcd_files = [base_path+"prediction_hcd_train.hdf5", base_path+"prediction_hcd_val.hdf5", base_path+"prediction_hcd_ho.hdf5"]
hcd_files = [base_path+"prediction_hcd_ho.hdf5"]
tof_files = [base_path+"prediction_tof_train.hdf5", base_path+"prediction_tof_validation.hdf5", base_path+"prediction_tof_test.hdf5"]
# tof_files = [base_path+"prediction_tof_test.hdf5"]
f = h5py.File(base_path+"prediction_tof_train.hdf5", 'r')
[key for key in f.keys()]

tofid_dict_list = []
for file in tof_files:
    temp_df = read_hdf5_to_dataframe(file, ['precursor_charge_onehot', 'sequence_integer'])
    tofid_dict_list.append(temp_df)

tofid_dict = {k:[] for k in tofid_dict_list[0].keys()}

for d in tofid_dict_list:
    for k, v in d.items():
        tofid_dict[k].append(v)

for k, v in tofid_dict.items():
    tofid_dict[k] = np.concatenate(v, axis=0)

filter_array = np.unique(tofid_dict["sequence_integer"], axis = 0)
len(filter_array)

hcd_dict_list = []
for file in hcd_files:
    temp_df = read_hdf5_to_dataframe(file, ['precursor_charge_onehot', 'sequence_integer',
        # 'collision_energy_aligned_normed', 'score', 'rawfile', 'collision_energy_normed',
        'collision_energy', 'intensities_raw', 'masses_raw'])
    hcd_dict_list.append(temp_df)

filtered_hcd_dict = {k:[] for k in hcd_dict_list[0].keys()}
len(hcd_dict_list[0]['sequence_integer'])

for d in hcd_dict_list:
    for x in zip(*d.values()):
        if x[-1] in filter_array:
            for k, xi in zip(hcd_dict_list[0].keys(), x):
                filtered_hcd_dict[k].append(xi)

filter_array = np.unique(filtered_hcd_dict["sequence_integer"], axis = 0)
len(filter_array)

tof_dict_list = []
for file in tof_files:
    temp_df = read_hdf5_to_dataframe(file, ['precursor_charge_onehot', 'sequence_integer',
        # 'collision_energy_aligned_normed', 'score', 'rawfile', 'collision_energy_normed',
        'collision_energy', 'intensities_raw', 'masses_raw'])
    tof_dict_list.append(temp_df)

filtered_tof_dict = {k:[] for k in tof_dict_list[0].keys()}

for d in tof_dict_list:
    for x in zip(*d.values()):
        if x[-1] in filter_array:
            for k, xi in zip(tof_dict_list[0].keys(), x):
                filtered_tof_dict[k].append(xi)

filtered_tof_dict["precursor_charge_onehot"] = filtered_tof_dict["precursor_charge_onehot"].tolist()
filtered_tof_dict["sequence_integer"] = filtered_tof_dict["sequence_integer"].tolist()
filtered_tof_dict["intensities_raw"] = filtered_tof_dict["intensities_raw"].tolist()
filtered_tof_dict["masses_raw"] = filtered_tof_dict["masses_raw"].tolist()
tof_df = pd.DataFrame(filtered_tof_dict)

filtered_hcd_dict["precursor_charge_onehot"] = filtered_hcd_dict["precursor_charge_onehot"].tolist()
filtered_hcd_dict["sequence_integer"] = filtered_hcd_dict["sequence_integer"].tolist()
filtered_hcd_dict["intensities_raw"] = filtered_hcd_dict["intensities_raw"].tolist()
filtered_hcd_dict["masses_raw"] = filtered_hcd_dict["masses_raw"].tolist()
hcd_df = pd.DataFrame(filtered_hcd_dict)

hcd_df.rename(columns = {'masses_raw':'MZ_HCD', 'intensities_raw': 'INTENSITY_HCD', 'collision_energy':'CE_HCD'}, inplace = True)
tof_df.rename(columns = {'masses_raw':'MZ_TOF', 'intensities_raw': 'INTENSITY_TOF', 'collision_energy':'CE_TOF'}, inplace = True)

hcd_df["precursor_charge_onehot"] = hcd_df["precursor_charge_onehot"].apply(str)
hcd_df["sequence_integer"] = hcd_df["sequence_integer"].apply(str)
hcd_df["CE_HCD"] = hcd_df["CE_HCD"].apply(int)
tof_df["precursor_charge_onehot"] = tof_df["precursor_charge_onehot"].apply(str)
tof_df["sequence_integer"] = tof_df["sequence_integer"].apply(str)
tof_df["CE_TOF"] = tof_df["CE_TOF"].apply(int)

hcd_tof_df = pd.merge(hcd_df, tof_df, how="inner", on=["precursor_charge_onehot", "sequence_integer"])
hcd_tof_df["SPECTRAL_ANGLE"] = hcd_tof_df[['INTENSITY_TOF','INTENSITY_HCD']].apply(lambda x : get_spectral_angle(x), axis=1)

mean_df = hcd_tof_df[["CE_HCD", "CE_TOF", "SPECTRAL_ANGLE"]].groupby(["CE_HCD", "CE_TOF"]).median("SPECTRAL_ANGLE").reset_index()
glue = mean_df.pivot("CE_HCD", "CE_TOF", "SPECTRAL_ANGLE")

cm = 1/2.54  # centimeters in inches
width = 18*cm
height = 7*cm
fig, ax = plt.subplots(figsize=(width, height))

sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_context("paper")

cmap0 = mpl.colors.LinearSegmentedColormap.from_list(
        'green2red', ['#0E1C36' ,'#7a8db3', '#C4E6B3'])

# plotting the heatmap
sns.heatmap(data=glue, cmap=cmap0, annot=False,
            vmin=0, vmax=1,
            cbar_kws={'label': 'Median spectral angle'})

ax.set_xlabel("Colission energy - timsTOF")
ax.set_ylabel("Colission energy - Orbitrap")

plot_path = "/Users/adams/Projects/300K/Results/Figures/paper-fig-1b.png"
plt.savefig(plot_path, dpi=600, bbox_inches='tight')