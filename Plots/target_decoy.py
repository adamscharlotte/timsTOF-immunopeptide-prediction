import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

def plot_target_decoy(target: pd.DataFrame, decoy: pd.DataFrame, type: str, search_type: str, directory: str):
    """Generate target-decoy distribution of the score."""
    plt.figure(figsize=(6, 4))
    bins = np.linspace(-3, 2, 45)
    plt.hist(target.score, bins, label="Targets", rwidth=0.8, color="#7a8db3")
    plt.hist(decoy.score, bins, label="Decoys", rwidth=0.8, color="#0e1c36")
    plt.axvline(color='#f33b16')
    plt.xlabel("Score", size=14)
    plt.ylim(0,140000)
    # plt.title(f"{search_type} Target vs Decoys {type}", size=14)
    plt.title(f"", size=14)
    plt.legend(loc="upper right", frameon=False)
    plt.savefig(directory + f"/{search_type}_Target_vs_Decoys_{type}_bins.png", dpi=300)

def load_psms_to_df(root_folder, file_pattern):
    file_paths = []
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
        df = pd.read_csv(file_path, sep="\t")  # Modify the read function based on the file Sample
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df["RAW_FILE"] = ["-".join(item.split("-")[0:1]) for item in combined_df["PSMId"]]
    return combined_df

def load_psm_to_df(root_folder, file_pattern):
    file_paths = []
    for file in os.listdir(root_folder):
        if file.endswith(file_pattern):
            file_path = os.path.join(root_folder, file)
            # Add the file path to the list
            file_paths.append(file_path)
    # Initialize an empty list to store individual DataFrames
    dfs = []
    # Read each file as a DataFrame and append it to the list
    for file_path in file_paths:
        df = pd.read_csv(file_path, sep="\t")  # Modify the read function based on the file Sample
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df["RAW_FILE"] = ["-".join(item.split("-")[0:1]) for item in combined_df["PSMId"]]
    return combined_df

prosit_psms_target = load_psms_to_df(root_folder, "rescore_target.psms")
prosit_psms_decoy = load_psms_to_df(root_folder, "rescore_decoy.psms")

andromeda_psms_target = load_psms_to_df(root_folder, "original_target.psms")
andromeda_psms_decoy = load_psms_to_df(root_folder, "original_decoy.psms")

# prosit_psm_target = load_psm_to_df(root_folder, "rescore_target.psms")
# prosit_psm_decoy = load_psm_to_df(root_folder, "rescore_decoy.psms")

# andromeda_psm_target = load_psm_to_df(root_folder, "original_target.psms")
# andromeda_psm_decoy = load_psm_to_df(root_folder, "original_decoy.psms")

# plot_target_decoy(prosit_pep_target, prosit_pep_decoy, "Peptides", "Rescore", percolator_path)
plot_target_decoy(prosit_psms_target, prosit_psms_decoy, "PSMs", "After Rescoring", figure_path)
# plot_target_decoy(andromeda_pep_target, andromeda_pep_decoy, "Peptides", "Original", percolator_path)
plot_target_decoy(andromeda_psms_target, andromeda_psms_decoy, "PSMs", "Before Rescoring", figure_path)
