import pandas as pd
import numpy as np
import logging
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Set up input and output file names

def load_scp_to_dfs(root_folder, file_pattern, dir_list):
    # Initialize an empty list to store the file paths
    file_paths = []
    # Traverse through each sample folder and find the result file
    for sample_folder in dir_list:
        sample_folder_path = os.path.join(root_folder, sample_folder)
        result_folder_path = os.path.join(sample_folder_path, f"results/percolator")
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
        df = pd.read_csv(file_path, sep="\t")  # Modify the read function based on the file Sample
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df["RAW_FILE"] = ["-".join(item.split("-")[0:1]) for item in combined_df["PSMId"]]
    return combined_df


root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/A375/reresults/A375-low-input-HLAI"
# dir_1e6 = ["E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep1_Slot1-06_1_3498",
#             "E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep3_Slot2-1_1_3508",
#             "E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep4_Slot2-1_1_3509"]
# dir_5e6 = ["E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep1_Slot1-37_1_3502",
#             "E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep2_Slot2-2_1_3510",
#             "E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep3_Slot2-2_1_3511"]
# dir_1e7 = ["E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep1_Slot2-3_1_3512",
#            "E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep2_Slot2-3_1_3513",
#            "E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep3_Slot2-3_1_3514"]
# dir_2e7 = ["E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep1_Slot2-4_1_3515",
#            "E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep3_Slot2-4_1_3518",
#            "E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep4_Slot2-4_1_3519"]
# dir_4e7 = ["E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496",
#            "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep2_Slot2-5_1_3520",
#            "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep3_Slot2-5_1_3521",
#            "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep4_Slot2-5_1_3522"]

dir_list = ["E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep1_Slot1-06_1_3498",
            "E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep3_Slot2-1_1_3508",
            "E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep4_Slot2-1_1_3509",
            "E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep1_Slot1-37_1_3502",
            "E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep2_Slot2-2_1_3510",
            "E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep3_Slot2-2_1_3511",
            "E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep1_Slot2-3_1_3512",
            "E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep2_Slot2-3_1_3513",
            "E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep3_Slot2-3_1_3514",
            "E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep1_Slot2-4_1_3515",
            "E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep3_Slot2-4_1_3518",
            "E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep4_Slot2-4_1_3519",
        #    "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496",
             "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep2_Slot2-5_1_3520",
           "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep3_Slot2-5_1_3521",
           "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep4_Slot2-5_1_3522"]

output_dir = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/A375/Figures"

scp_prosit_pep = load_scp_to_dfs(root_folder, "rescore_target.peptides", dir_list)
scp_adromeda_pep = load_scp_to_dfs(root_folder, "original_target.peptides", dir_list)

scp_prosit_psms = load_scp_to_dfs(root_folder, "rescore_target.psms", dir_list)
scp_adromeda_psms = load_scp_to_dfs(root_folder, "original_target.psms", dir_list)

join_col = ["Sample","peptide"]
join_col = ["PSMId"]

plot_name = "scp-psm-hlai"

def plot_gain_loss(prosit_target: pd.DataFrame, andromeda_target: pd.DataFrame, join_col, directory: str, plot_name: str, type: str):
    prosit_target["Sample"] = prosit_target["RAW_FILE"].str[28:-30]
    andromeda_target["Sample"] = andromeda_target["RAW_FILE"].str[28:-30]
    prosit_df = prosit_target[prosit_target["q-value"] < 0.01]
    adromeda_df = andromeda_target[andromeda_target["q-value"] < 0.01]
    merged_df = prosit_df.merge(adromeda_df, how="inner", on=join_col)
    adromeda_df = adromeda_df.drop_duplicates(subset=join_col)
    prosit_df = prosit_df.drop_duplicates(subset=join_col)
    merged_df = merged_df.drop_duplicates(subset=join_col)
    len(prosit_df)/len(adromeda_df) # 50% increase
    merged_df["Value"] = merged_df.groupby('Sample')[type].transform("count")
    prosit_df["Value"] = prosit_df.groupby('Sample')[type].transform("count")
    adromeda_df["Value"] = adromeda_df.groupby('Sample')[type].transform("count")
    merged_df["Label"] = "Shared"
    adromeda_df["Label"] = "Lost"
    prosit_df["Label"] = "Gained"
    merged = merged_df[["Value", "Sample", "Label"]].drop_duplicates().reset_index()[["Value", "Sample", "Label"]].sort_values("Sample")
    andromeda = adromeda_df[["Value", "Sample", "Label"]].drop_duplicates().reset_index()[["Value", "Sample", "Label"]].sort_values("Sample")
    prosit = prosit_df[["Value", "Sample", "Label"]].drop_duplicates().reset_index()[["Value", "Sample", "Label"]].sort_values("Sample")
    andromeda["Value"] = andromeda["Value"] - merged["Value"]
    prosit = prosit.sort_values("Sample")
    merged = merged.sort_values("Sample")
    andromeda = andromeda.sort_values("Sample")
    full_df = prosit.append(merged).append(andromeda)
    full_df['Value'] = np.where(full_df['Label'] == 'Lost', -full_df['Value'], full_df['Value'])
    plt.figure()
    sns.set_context("talk")
    # width = 4.5
    # height = 6.0
    width = 4.5
    height = 6
    fig, ax = plt.subplots(figsize=(width, height))
    ax = sns.barplot(data=full_df, x='Sample', y='Value', hue='Label', dodge=False, palette=["#cdeac0", "#7a8db3", "#f33b16"], order = ["1e6", "5e6", "1e7", "2e7", "4e7"])
    ax.axhline(0, lw=2, color='black')
    ax.legend(loc="center right", bbox_to_anchor=(1.8, 0.5))#,fancybox=True, ncol=1)
    # plt.ylim(0, max_value)
    # plt.xticks(rotation = 45) # Rotates X-Axis Ticks by 45-degrees
    ylabels = ["{:,.0f}".format(x) + "K" for x in ax.get_yticks()/1000]
    ax.set_yticklabels(ylabels)
    plt.ylabel('unique HLA-I ligands')
    plt.xlabel('Number of cells')
    sns.despine()
    plt.savefig(directory + f"/{plot_name}_1%_FDR.png", dpi=300, bbox_inches="tight")

plot_gain_loss(scp_prosit_psms, scp_adromeda_psms, ["Sample", "PSMId"], output_dir, "scp-psm-hlai", "PSMId")
plot_gain_loss(scp_prosit_pep, scp_adromeda_pep, ["Sample", "peptide"], output_dir, "scp-pep-hlai", "peptide")

def plot_mean(prosit_target: pd.DataFrame, andromeda_target: pd.DataFrame, join_col, directory: str, plot_name: str, type: str, y_axis_name):
    prosit_df = prosit_target[prosit_target["q-value"] < 0.01]
    adromeda_df = andromeda_target[andromeda_target["q-value"] < 0.01]
    adromeda_df = adromeda_df.drop_duplicates(subset=join_col)
    prosit_df = prosit_df.drop_duplicates(subset=join_col)
    prosit_df["Value"] = prosit_df.groupby('RAW_FILE')[type].transform("count")
    adromeda_df["Value"] = adromeda_df.groupby('RAW_FILE')[type].transform("count")
    adromeda_df["Label"] = "MaxQuant"
    prosit_df["Label"] = "Rescored"
    andromeda = adromeda_df[["Value", "RAW_FILE", "Label"]].drop_duplicates().reset_index()[["Value", "RAW_FILE", "Label"]].sort_values("RAW_FILE")
    prosit = prosit_df[["Value", "RAW_FILE", "Label"]].drop_duplicates().reset_index()[["Value", "RAW_FILE", "Label"]].sort_values("RAW_FILE")
    prosit = prosit.sort_values("RAW_FILE")
    andromeda = andromeda.sort_values("RAW_FILE")
    full_df = andromeda.append(prosit)
    full_df["Sample"] = full_df["RAW_FILE"].str[28:-30]
    plt.figure()
    sns.set_context("talk")
    # width = 4.5
    # height = 6.0
    width = 4.5
    height = 5
    order_x = ["1e6", "5e6", "1e7", "2e7", "4e7"]
    fig, ax = plt.subplots(figsize=(width, height))
    # ax = sns.barplot(data=full_df, x='Sample', y='Value', hue='Label', errorbar="sd",capsize=.4,
    #     errwidth=1, palette=["#0e1c36", "#cdeac0"], order=order_x)
    ax = sns.barplot(data=full_df, x='Sample', y='Value', hue='Label', errorbar=None,capsize=.4,
        palette=["#0e1c36", "#cdeac0"], order=order_x)
    ax.axhline(0, lw=2, color='black')
    ax.legend(loc="center right", bbox_to_anchor=(2, 0.5))#,fancybox=True, ncol=1)
    # plt.ylim(0, max_value)
    # plt.xticks(rotation = 45) # Rotates X-Axis Ticks by 45-degrees
    ylabels = ["{:,.0f}".format(x) + "K" for x in ax.get_yticks()/1000]
    ax.set_yticklabels(ylabels)
    plt.ylabel(y_axis_name)
    plt.xlabel('Cell input')
    sns.despine()
    plt.savefig(directory + f"/{plot_name}_1%_FDR.png", dpi=300, bbox_inches="tight")

# plot_mean(scp_prosit_psms, scp_adromeda_psms, ["RAW_FILE", "PSMId"], output_dir, "mean-psm-err", "PSMId", "PSMs")
# plot_mean(scp_prosit_pep, scp_adromeda_pep, ["RAW_FILE", "peptide"], output_dir, "mean-pep-err", "peptide", "unique HLA-I ligands")

plot_mean(scp_prosit_psms, scp_adromeda_psms, ["RAW_FILE", "PSMId"], output_dir, "mean-psm", "PSMId", "PSMs")
plot_mean(scp_prosit_pep, scp_adromeda_pep, ["RAW_FILE", "peptide"], output_dir, "mean-pep", "peptide", "unique HLA-I ligands")

prosit_df = scp_prosit_pep[scp_prosit_pep["q-value"] < 0.01]
adromeda_df = scp_adromeda_pep[scp_adromeda_pep["q-value"] < 0.01]
adromeda_df = adromeda_df.drop_duplicates(subset=["RAW_FILE", "peptide"])
prosit_df = prosit_df.drop_duplicates(subset=["RAW_FILE", "peptide"])
prosit_df["Value"] = prosit_df.groupby('RAW_FILE')["peptide"].transform("count")
adromeda_df["Value"] = adromeda_df.groupby('RAW_FILE')["peptide"].transform("count")
adromeda_df["Label"] = "MaxQuant"
prosit_df["Label"] = "MaxQuant-Rescored"
andromeda = adromeda_df[["Value", "RAW_FILE", "Label"]].drop_duplicates().reset_index()[["Value", "RAW_FILE", "Label"]].sort_values("RAW_FILE")
prosit = prosit_df[["Value", "RAW_FILE", "Label"]].drop_duplicates().reset_index()[["Value", "RAW_FILE", "Label"]].sort_values("RAW_FILE")
prosit = prosit.sort_values("RAW_FILE")
andromeda = andromeda.sort_values("RAW_FILE")
full_df = andromeda.append(prosit)
full_df["Sample"] = full_df["RAW_FILE"].str[28:-30]
