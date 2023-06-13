# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import matplotlib.pyplot as plt
import statistics
import pandas as pd
import numpy as np
import seaborn as sns
import glob
import os
from prosit_grpc.predictPROSIT import PROSITpredictor
import h5py
from fundamentals.annotation.annotation import annotate_spectra

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

def load_files_to_dfs(root_folder, file_pattern):
    # Initialize an empty list to store the file paths
    file_paths = []
    # Traverse through each sample folder and find the result file
    for sample_folder in os.listdir(root_folder):
        sample_folder_path = os.path.join(root_folder, sample_folder)
        result_folder_path = os.path.join(sample_folder_path, "results/percolator")
        # Check if the result folder exists
        if os.path.isdir(result_folder_path):
            # Find the files matching the file pattern within the result folder
            for file in os.listdir(result_folder_path):
                if file.endswith(file_pattern):
                    file_path = os.path.join(result_folder_path, file)
                    # Add the file path to the list
                    file_paths.append(file_path)
    # Initialize an empty list to store individual DataFrames
    dfs = []
    # Read each file as a DataFrame and append it to the list
    for file_path in file_paths:
        df = pd.read_csv(file_path, sep="\t")  # Modify the read function based on the file type
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in combined_df["PSMId"]]
    combined_df["SCAN_NUMBER"] = [item.split("-")[2] for item in combined_df["PSMId"]]
    combined_df["GRPC_SEQUENCE"] = [item.split("-")[3] for item in combined_df["PSMId"]]
    combined_df["SCAN_NUMBER"] = combined_df["SCAN_NUMBER"].astype(int)
    return combined_df

def load_hdf5_to_dfs(root_folder, file_pattern):
    # Initialize an empty list to store the file paths
    file_paths = []
    # Traverse through each sample folder and find the result file
    for sample_folder in os.listdir(root_folder):
        sample_folder_path = os.path.join(root_folder, sample_folder)
        result_folder_path = os.path.join(sample_folder_path, "mzML")
        # Check if the result folder exists
        if os.path.isdir(result_folder_path):
            # Find the files matching the file pattern within the result folder
            for file in os.listdir(result_folder_path):
                if file.endswith(file_pattern):
                    file_path = os.path.join(result_folder_path, file)
                    # Add the file path to the list
                    file_paths.append(file_path)
    # Initialize an empty list to store individual DataFrames
    dfs = []
    # Read each file as a DataFrame and append it to the list
    for file_path in file_paths:
        df = pd.read_hdf(file_path, sep="\t")  # Modify the read function based on the file type
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def load_pkl_to_dfs(root_folder, file_pattern):
    # Initialize an empty list to store the file paths
    file_paths = []
    # Traverse through each sample folder and find the result file
    for sample_folder in os.listdir(root_folder):
        sample_folder_path = os.path.join(root_folder, sample_folder)
        # Check if the result folder exists
        if os.path.isdir(sample_folder_path):
            # Find the files matching the file pattern within the result folder
            for file in os.listdir(sample_folder_path):
                if file.endswith(file_pattern):
                    file_path = os.path.join(sample_folder_path, file)
                    # Add the file path to the list
                    file_paths.append(file_path)
    # Initialize an empty list to store individual DataFrames
    dfs = []
    # Read each file as a DataFrame and append it to the list
    for file_path in file_paths:
        df = pd.read_pickle(file_path)  # Modify the read function based on the file type
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def calculate_pred_sa(df, predictor, model):
    predictions = predictor.predict(sequences=df["MODIFIED_SEQUENCE"].values.tolist(),
                                charges=df["PRECURSOR_CHARGE"].values.tolist(),
                                collision_energies=df["COLLISION_ENERGY"].values/100.0,
                                                    models=[model],
                                                    disable_progress_bar=True)
    df["PREDICTED_INTENSITY"] = predictions[model]["intensity"].tolist()
    df["SPECTRAL_ANGLE"] = df[["INTENSITIES","PREDICTED_INTENSITY"]].apply(lambda x : get_spectral_angle(x), axis=1)
    df["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    return df


base_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/"  # nolint
tims_path = base_path + "reresults/d-tims"
cid_path = base_path + "reresults/d-cid"
hcd_path = base_path + "reresults/d-hcd"

psm_tims = load_files_to_dfs(tims_path, "rescore_target.psms")
psm_cid = load_files_to_dfs(cid_path, "rescore_target.psms")
psm_hcd = load_files_to_dfs(hcd_path, "rescore_target.psms")

hdf_tims = load_hdf5_to_dfs(tims_path, ".mzML_pred.hdf5")
hdf_cid = load_hdf5_to_dfs(cid_path, ".mzML_pred.hdf5")
hdf_hcd = load_hdf5_to_dfs(hcd_path, ".mzML_pred.hdf5")

pkl_df = load_pkl_to_dfs(tims_path, ".pkl")

psm_tims_1 = psm_tims[psm_tims["q-value"] < 0.01]
psm_cid_1 = psm_cid[psm_cid["q-value"] < 0.01]
psm_hcd_1 = psm_hcd[psm_hcd["q-value"] < 0.01]

# psm_tims_1["SCAN_NUMBER"] = psm_tims_1["SCAN_NUMBER"].astype(int)
# psm_cid_1["SCAN_NUMBER"] = psm_cid_1["SCAN_NUMBER"].astype(int)
# psm_hcd_1["SCAN_NUMBER"] = psm_hcd_1["SCAN_NUMBER"].astype(int)

tims_ce_df = pd.merge(psm_tims_1, hdf_tims, how="inner", on=["RAW_FILE","SCAN_NUMBER", "GRPC_SEQUENCE"])
cid_ce_df = pd.merge(psm_cid_1, hdf_cid, how="inner", on=["RAW_FILE","SCAN_NUMBER", "GRPC_SEQUENCE"])
hcd_ce_df = pd.merge(psm_hcd_1, hdf_hcd, how="inner", on=["RAW_FILE","SCAN_NUMBER", "GRPC_SEQUENCE"])

pkl_df = pkl_df.drop(columns=["COLLISION_ENERGY", "MZ_RANGE", "MASS_ANALYZER", "FRAGMENTATION"])
tims_df = pd.merge(tims_ce_df, pkl_df, how="inner", on=["RAW_FILE","SCAN_NUMBER"])
cid_df = pd.merge(cid_ce_df, pkl_df, how="inner", on=["RAW_FILE","SCAN_NUMBER"])
hcd_df = pd.merge(hcd_ce_df, pkl_df, how="inner", on=["RAW_FILE","SCAN_NUMBER"])

annot_tims_df = annotate_spectra(tims_df)
annot_cid_df = annotate_spectra(cid_df)
annot_hcd_df = annotate_spectra(hcd_df)

tims_df = tims_df.drop(columns=["CALCULATED_MASS"])
tims_full_df = pd.concat([tims_df.drop(columns = ["INTENSITIES", "MZ"]), annot_tims_df], axis=1)
cid_df = cid_df.drop(columns=["CALCULATED_MASS"])
cid_full_df = pd.concat([cid_df.drop(columns = ["INTENSITIES", "MZ"]), annot_cid_df], axis=1)
hcd_df = hcd_df.drop(columns=["CALCULATED_MASS"])
hcd_full_df = pd.concat([hcd_df.drop(columns = ["INTENSITIES", "MZ"]), annot_hcd_df], axis=1)

tims_full_df["label"] = "HCD TIMS Prosit 2023"
cid_full_df["label"] = "CID Prosit 2020"
hcd_full_df["label"] = "HCD Prosit 2020"

path = "/home/cadams/anaconda3/envs/oktoberfest-0_1/lib/python3.8/site-packages/oktoberfest/certificates"
predictor = PROSITpredictor(server="10.152.135.57:8500",
            path_to_ca_certificate=os.path.join(path, "Proteomicsdb-Prosit-v2.crt"),
            path_to_certificate=os.path.join(path, "oktoberfest-production.crt"),
            path_to_key_certificate=os.path.join(path, "oktoberfest-production.key"),)

tims_sa_df = calculate_pred_sa(tims_full_df, predictor, "Prosit_2023_intensity_tims")
cid_sa_df = calculate_pred_sa(cid_full_df, predictor, "Prosit_2020_intensity_cid")
hcd_sa_df = calculate_pred_sa(hcd_full_df, predictor, "Prosit_2020_intensity_hcd")

tims_vs_cid = tims_sa_df.append(cid_sa_df)
tims_vs_cid["type"] = "test"

tims_vs_hcd = tims_sa_df.append(hcd_sa_df)
tims_vs_hcd["type"] = "test"

plt.rcParams["figure.figsize"] = [8, 5]

fig = plt.figure()
ax = fig.gca()
ax = sns.violinplot(x="type", y="SPECTRAL_ANGLE", hue="label", inner=None,
                    edgecolor=None, linewidth=0, palette=["#0E1C36", "#CDEAC0"],
                    hue_order=["HCD TIMS Prosit 2023", "HCD Prosit 2020"],
                    data=tims_vs_hcd, split=True)
ax.legend(loc="lower left", bbox_to_anchor=(1.0, 0.1), frameon=False)
#sns.move_legend(ax, "lower right")
ax.set_ylabel("Spectral angle", color="black")
ax.set(xlabel=None)
ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")
ax.spines["left"].set_color("none")
ax.spines["bottom"].set_color("none")
sns.set_style("whitegrid")
ax.yaxis.grid(True, linewidth=0.5, which="major", color="lightgrey", alpha=0.5)
plt.yticks([i/10 for i in range(0, 11)])
ax.tick_params(axis="x", colors="#464646")
ax.tick_params(axis="y", colors="#464646")
ax_x = plt.gca().xaxis
ax_x.set_tick_params(pad=-5)
plot_path = base_path + "/Figures/Violin/tims-vs-hcd.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
