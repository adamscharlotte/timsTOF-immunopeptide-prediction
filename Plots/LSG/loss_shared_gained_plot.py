import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import re
from tqdm import tqdm
import matplotlib.ticker as mtick
import logging

sns.set_context("paper")

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
    return combined_df

def load_jurkat_to_dfs(root_folder, file_pattern, model, dir_list):
    # Initialize an empty list to store the file paths
    file_paths = []
    # Traverse through each sample folder and find the result file
    for sample_folder in dir_list:
        sample_folder_path = os.path.join(root_folder, sample_folder)
        result_folder_path = os.path.join(sample_folder_path, f"results-{model}/percolator")
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
    combined_df["RAW_FILE"] = combined_df["PSMId"].str.extract(r"(.+)-\d+-\D+-\d+--\d+")[0]
    return combined_df

def plot_gain_loss(prosit_target: pd.DataFrame, andromeda_target: pd.DataFrame, type: str, directory: str, plot_name: str):
    """Generate gain-loss plot (peptides/PSMs 1% FDR)."""
    if type == "Peptides":
        join_col = ["RAW_FILE","peptide"]
        andromeda_target = andromeda_target[andromeda_target["q-value"] < 0.01]
        prosit_target = prosit_target[prosit_target["q-value"] < 0.01]
        merged_df = prosit_target.merge(andromeda_target, how="inner", on=join_col)
        andromeda_target = andromeda_target.drop_duplicates(subset=["RAW_FILE", "peptide"])
        prosit_target = prosit_target.drop_duplicates(subset=["RAW_FILE", "peptide"])
        merged_df = merged_df.drop_duplicates(subset=["RAW_FILE", "peptide"])
    else:
        join_col = "PSMId"
        andromeda_target = andromeda_target[andromeda_target["q-value"] < 0.01]
        prosit_target = prosit_target[prosit_target["q-value"] < 0.01]
        merged_df = prosit_target.merge(andromeda_target, how="inner", on=join_col)
    shared = len(merged_df.index)
    gained = len(prosit_target.index) - shared
    lost = len(andromeda_target.index) - shared
    fig, ax = plt.subplots(1, figsize=(1.2, 9))
    labels = [""]
    ax1 = ax.bar(labels, shared, width=0.5, color="#7a8db3")
    ax2 = ax.bar(labels, gained, width=0.5, bottom=shared, color="#cdeac0")
    ax3 = ax.bar(labels, -lost, color="#f33b16", width=0.5)
    for r1, r2, r3, v1, v2, v3 in zip(ax1, ax2, ax3, [shared], [gained], [lost]):
        h1 = r1.get_height()
        h2 = r2.get_height()
        plt.text(
            r1.get_x() + r1.get_width() / 2.0,
            h1 / 2.0,
            "%d" % v1,
            ha="center",
            va="bottom",
            color="black",
            # fontsize=14,
        )
        plt.text(
            r2.get_x() + r2.get_width() / 2.0,
            h1 + h2 / 2,
            "%d" % v2,
            ha="center",
            va="bottom",
            color="black",
            # fontsize=14,
        )
        plt.text(
            r3.get_x() + r3.get_width() / 2.0,
            -0.025 * gained,
            "%d" % -v3,
            ha="center",
            va="bottom",
            color="black",
            # fontsize=14,
        )
    plt.ylim(-lost - 100, h1 + h2 + 30)
    # remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    # grid
    # ax.set_ylabel("Percentage")#, fontsize=14)
    ax.set_axisbelow(True)
    # ax.yaxis.grid(color="gray")
    # ax.tick_params(axis="y", which="major", labelsize=13)
    ax.set_xticks([])
    # for minor ticks
    # ax.set_yticks(np.arange(0, max(prosit_psms_table)+1, 100000))
    # ylabels = ["{:,.0f}".format(x) + "K" for x in ax.get_yticks()/1000]
    # ax.set_yticklabels(ylabels)
    ax.set_xticks([], minor=True)
    legend_label = ["Common", "Gained", "Lost"]
    sns.set_context("talk")
    plt.legend(legend_label, ncol=1, bbox_to_anchor=([1.2, 0.5, 0, 0]), frameon=False)
    plt.title(f"{type} 1% FDR\n")#, fontsize=14)
    plt.savefig(directory + f"/{plot_name}_{type}_1%_FDR.png", dpi=300, bbox_inches="tight")

def prep_multiple_gain_loss(root_folder, reorderlist):
    prosit_pep_target = load_files_to_dfs(root_folder, "rescore_target.peptides")
    andromeda_pep_target = load_files_to_dfs(root_folder, "original_target.peptides")
    prosit_pep_target["Sample"] = prosit_pep_target["RAW_FILE"].str.replace(r'(-\d+$)', "")
    andromeda_pep_target["Sample"] = andromeda_pep_target["RAW_FILE"].str.replace(r'(-\d+$)', "")
    join_col = ["Sample","peptide"]
    andromeda_target = andromeda_pep_target[andromeda_pep_target["q-value"] < 0.01]
    prosit_target = prosit_pep_target[prosit_pep_target["q-value"] < 0.01]
    merged_df = prosit_target.merge(andromeda_target, how="inner", on=join_col)
    len_prosit = len(prosit_target)
    len_adromeda = len(andromeda_target)
    gain_fold = len(prosit_target)/len(andromeda_target)
    logging.info(f"There are {len_adromeda} PSMs identified with MaxQuant.")
    logging.info(f"There are {len_prosit} PSMs identified after rescoring.")
    logging.info(f"There is a gain of {gain_fold} fold.")
    print(f"There are {len_adromeda} PSMs identified with MaxQuant.")
    print(f"There are {len_prosit} PSMs identified after rescoring.")
    print(f"There is a gain of {gain_fold} fold.")
    andromeda_target = andromeda_target.drop_duplicates(subset=["Sample", "peptide"])
    prosit_target = prosit_target.drop_duplicates(subset=["Sample", "peptide"])
    merged_df = merged_df.drop_duplicates(subset=["Sample", "peptide"])
    merged_df["Value"] = merged_df.groupby("Sample")['peptide'].transform("count")
    andromeda_target["Value"] = andromeda_target.groupby("Sample")['peptide'].transform("count")
    prosit_target["Value"] = prosit_target.groupby("Sample")['peptide'].transform("count")
    merged_df["Label"] = "Shared"
    andromeda_target["Label"] = "Lost"
    prosit_target["Label"] = "Gained"
    merged = merged_df[["Value", "Sample", "Label"]].drop_duplicates().reset_index()[["Value", "Sample", "Label"]].sort_values("Sample")
    andromeda = andromeda_target[["Value", "Sample", "Label"]].drop_duplicates().reset_index()[["Value", "Sample", "Label"]].sort_values("Sample")
    prosit = prosit_target[["Value", "Sample", "Label"]].drop_duplicates().reset_index()[["Value", "Sample", "Label"]].sort_values("Sample")
    andromeda["Value"] = andromeda["Value"] - merged["Value"]
    prosit = prosit.sort_values("Sample")
    merged = merged.sort_values("Sample")
    andromeda = andromeda.sort_values("Sample")
    full_df = prosit.append(merged).append(andromeda)
    full_df["Sample"] = full_df["Sample"].str.replace(r'(_Tue39L243_\d+%_orbitrap_DDA_Rep\d)', "")
    full_df["Sample"] = full_df["Sample"].str.replace(r'(_Tue39L243_\d+%_DDA_Rep\d)', "")
    full_df["Sample"] = full_df["Sample"].str.replace(r'(_Tue39L243_\d+%_Rep\d)', "")
    full_df["Sample"] = full_df["Sample"].str.replace(r'(^\d+_[A-Z]+_[a-z]+_)', "")
    full_df["Sample"] = full_df["Sample"].str.replace(r'(^UDN\d+_)', "")
    # pd.set_option('display.max_colwidth', None)
    # full_df.drop_duplicates("Sample")
    full_df["Sample"] = full_df["Sample"].str.replace(r'(_W6-\d+_\d+%_[a-z]+_DDA_Rep\d)', "")
    full_df["Sample"] = full_df["Sample"].str.replace(r'(_W6-\d+_\d+%_[A-Z]+_Rep\d)', "")
    full_df["Sample"] = full_df["Sample"].str.replace(r'(_W6-32_17%_Rep\d)', "")
    full_df['Value'] = np.where(full_df['Label'] == 'Lost', -full_df['Value'], full_df['Value'])
    return full_df

root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/d-tims-ii"
# root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/d-tims"
output_dir = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/Figures/LSG"
prosit_pep_target = load_files_to_dfs(root_folder, "rescore_target.peptides")
andromeda_pep_target = load_files_to_dfs(root_folder, "original_target.peptides")
prosit_psms_target = load_files_to_dfs(root_folder, "rescore_target.psms")
andromeda_psms_target = load_files_to_dfs(root_folder, "original_target.psms")

plot_name = "HLA-II-tims"
plot_gain_loss(prosit_pep_target, andromeda_pep_target, "Peptides", output_dir, plot_name)
plot_gain_loss(prosit_psms_target, andromeda_psms_target, "PSMs", output_dir, plot_name)

len(prosit_psms_target.RAW_FILE.drop_duplicates())
len(prosit_pep_target.RAW_FILE.drop_duplicates())

directory = output_dir
join_col = "PSMId"
andromeda_target = andromeda_psms_target[andromeda_psms_target["q-value"] < 0.01]
prosit_target = prosit_psms_target[prosit_psms_target["q-value"] < 0.01]
merged_df = prosit_target.merge(andromeda_target, how="inner", on=join_col)

shared = len(merged_df.index)
gained = len(prosit_target.index) - shared
lost = len(andromeda_target.index) - shared

# --------------------- Combined LSG plots ---------------------

# root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/Jurkat-A549/reresults"
# dir_list = ["IPX67_HLAI_A_S1-A2_1_10830"]#, "IPX67_HLAI_B_S1-A8_1_10836"]
# output_dir = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/Jurkat-A549/Figures/LSG"

# A549_prosit_pep_tims = load_jurkat_to_dfs(root_folder, "rescore_target.peptides", "tims", dir_list)
# A549_adromeda_pep_tims = load_jurkat_to_dfs(root_folder, "original_target.peptides", "tims", dir_list)

# A549_prosit_psms_tims = load_jurkat_to_dfs(root_folder, "rescore_target.psms", "tims", dir_list)
# A549_adromeda_psms_tims = load_jurkat_to_dfs(root_folder, "original_target.psms", "tims", dir_list)

# plot_name = "Jurkat-A-tims"

# plot_gain_loss(A549_prosit_pep_tims, A549_adromeda_pep_tims, "Peptides", output_dir, plot_name)
# plot_gain_loss(A549_prosit_psms_tims, A549_adromeda_psms_tims, "PSMs", output_dir, plot_name)

# --------------------- Multiple LSG plots ---------------------

# tims root folder
tims_folder_1 = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/d-tims"
orbi_folder_1 = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/raw-cid"
tims_folder_2 = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/d-tims-ii"
orbi_folder_2 = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/raw"

reorderlist = ["Brain", "Spleen", "Lung", "Thymus", "Liver", "RCC", "HNSCC", "PBMC", "CLL_01", "CLL_02"]

orbi_df_1 = prep_multiple_gain_loss(orbi_folder_1, reorderlist)
tims_df_1 = prep_multiple_gain_loss(tims_folder_1, reorderlist)
orbi_df_2 = prep_multiple_gain_loss(orbi_folder_2, reorderlist)
tims_df_2 = prep_multiple_gain_loss(tims_folder_2, reorderlist)

# Plot figure
plt.figure()
sns.set_context("paper")
cm = 1/2.54  # centimeters in inches
width = 18*cm
height = 14*cm

fig, axes = plt.subplots(2, 2, figsize=(width, height))
axes = np.ravel(axes)

sns.barplot(data=orbi_df_1, x="Sample", y='Value', hue='Label', dodge=False, palette=["#cdeac0", "#7a8db3", "#f33b16"], order = reorderlist, errorbar=None, ax = axes[0])
sns.barplot(data=tims_df_1, x="Sample", y='Value', hue='Label', dodge=False, palette=["#cdeac0", "#7a8db3", "#f33b16"], order = reorderlist, errorbar=None, ax = axes[1])
sns.barplot(data=orbi_df_2, x="Sample", y='Value', hue='Label', dodge=False, palette=["#cdeac0", "#7a8db3", "#f33b16"], order = reorderlist, errorbar=None, ax = axes[2])
sns.barplot(data=tims_df_2, x="Sample", y='Value', hue='Label', dodge=False, palette=["#cdeac0", "#7a8db3", "#f33b16"], order = reorderlist, errorbar=None, ax = axes[3])

axes[0].set_ylim(-250, 12000)
axes[1].set_ylim(-250, 12000)
axes[2].set_ylim(-250, 12000)
axes[3].set_ylim(-250, 12000)

axes[0].set_ylabel('unique HLA-I ligands')
axes[2].set_ylabel('unique HLA-II ligands')
axes[1].set_ylabel('')
axes[3].set_ylabel('')

axes[0].set_xlabel('')
axes[1].set_xlabel('')

for ax in axes:
    ax.set_xticklabels(reorderlist, rotation = 45)
    # ax.set_xlabel("Sample")
    ylabels = ["{:,.0f}".format(x) + "K" for x in ax.get_yticks()/1000]
    ax.set_yticklabels(ylabels)

axes[1].get_legend().remove()
axes[2].get_legend().remove()
axes[3].get_legend().remove()

# axes[1].axes.get_yaxis().set_visible(False)
axes[0].legend(loc='upper left')#, bbox_to_anchor=(1, 1.15)

for i, (ax, c) in enumerate(zip(axes, 'abcd')):
    ax.annotate(c, xy=(-0.18, 1.05), xycoords='axes fraction',
                fontsize='xx-large', weight='bold')

sns.despine()
fig.tight_layout()

plot_name = "paper-orbi-vs-tims-samples"
plt.savefig(directory + f"/{plot_name}_1%_FDR.png", dpi=300, bbox_inches="tight")

# --------------------- Multiple orbi vs tims plots ---------------------
directory = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/Figures/LSG"

root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/raw-cid"
orbi_pep_target_rescore = load_files_to_dfs(root_folder, "rescore_target.peptides")
orbi_pep_target = load_files_to_dfs(root_folder, "original_target.peptides")
root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/d-tims"
tims_pep_target_rescore = load_files_to_dfs(root_folder, "rescore_target.peptides")
tims_pep_target = load_files_to_dfs(root_folder, "original_target.peptides")

orbi_pep_target_rescore["Sample"] = orbi_pep_target_rescore["RAW_FILE"].str[:-4]
tims_pep_target_rescore["Sample"] = tims_pep_target_rescore["RAW_FILE"].str[:-4]

tims_pep_target_rescore["type"] = tims_pep_target_rescore["Sample"].str.replace(r'(^\d+_[A-Z]+_[a-z]+_)', "")
tims_pep_target_rescore["type"] = tims_pep_target_rescore["type"].str.replace(r'(_W6-\d+_\d+%_[A-Z]+_$)', "")
tims_pep_target_rescore["type"] = tims_pep_target_rescore["type"].str.replace(r'(^UDN\d+_)', "")
tims_pep_target_rescore["type"] = tims_pep_target_rescore["type"].str.replace(r'(_W6-32_17%_$)', "")

orbi_pep_target_rescore["type"] = orbi_pep_target_rescore["Sample"].str.replace(r'(^\d+_[A-Z]+_[a-z]+_)', "")
orbi_pep_target_rescore["type"] = orbi_pep_target_rescore["type"].str.replace(r'(_W6-\d+_\d+%_[a-z]+_DDA_$)', "")
orbi_pep_target_rescore["type"] = orbi_pep_target_rescore["type"].str.replace(r'(^UDN\d+_)', "")

join_col = ["type","peptide"]
orbi_target_rescore = orbi_pep_target_rescore[orbi_pep_target_rescore["q-value"] < 0.01]
tims_target_rescore = tims_pep_target_rescore[tims_pep_target_rescore["q-value"] < 0.01]
merged_df_rescore = tims_target_rescore.merge(orbi_target_rescore, how="inner", on=join_col)

orbi_target_rescore = orbi_target_rescore.drop_duplicates(subset=["type", "peptide"])
tims_target_rescore = tims_target_rescore.drop_duplicates(subset=["type", "peptide"])
merged_df_rescore = merged_df_rescore.drop_duplicates(subset=["type", "peptide"])

len(orbi_target_rescore)
len(tims_target_rescore)
len(tims_target_rescore)/len(orbi_target_rescore)

merged_df_rescore["Value"] = merged_df_rescore.groupby('type')['peptide'].transform("count")
orbi_target_rescore["Value"] = orbi_target_rescore.groupby('type')['peptide'].transform("count")
tims_target_rescore["Value"] = tims_target_rescore.groupby('type')['peptide'].transform("count")

merged_df_rescore["Label"] = "Shared"
orbi_target_rescore["Label"] = "Orbitrap"
tims_target_rescore["Label"] = "timsTOF"

merged_rescore = merged_df_rescore[["Value", "type", "Label"]].drop_duplicates().reset_index()[["Value", "type", "Label"]].sort_values("type")
orbi_rescore = orbi_target_rescore[["Value", "type", "Label"]].drop_duplicates().reset_index()[["Value", "type", "Label"]].sort_values("type")
tims_rescore = tims_target_rescore[["Value", "type", "Label"]].drop_duplicates().reset_index()[["Value", "type", "Label"]].sort_values("type")

tims_rescore["Value"] = tims_rescore["Value"] + orbi_rescore["Value"]
merged_rescore["Value"] = merged_rescore["Value"] + orbi_rescore["Value"]
# prosit["Value"] = prosit["Value"]

tims_rescore = tims_rescore.sort_values("type")
merged_rescore = merged_rescore.sort_values("type")
orbi_rescore = orbi_rescore.sort_values("type")

max_value = tims_rescore["Value"].max()

reorderlist = ["Brain", "Spleen", "Lung", "Thymus", "Liver", "RCC", "HNSCC", "PBMC", "CLL_01", "CLL_02"]
full_df_rescore = tims_rescore.append(merged_rescore).append(orbi_rescore)


tims_pep_target["Sample"] = tims_pep_target["RAW_FILE"].str[:-4]
orbi_pep_target["Sample"] = orbi_pep_target["RAW_FILE"].str[:-4]

tims_pep_target["type"] = tims_pep_target["Sample"].str.replace(r'(^\d+_[A-Z]+_[a-z]+_)', "")
tims_pep_target["type"] = tims_pep_target["type"].str.replace(r'(_W6-\d+_\d+%_[A-Z]+_$)', "")
tims_pep_target["type"] = tims_pep_target["type"].str.replace(r'(^UDN\d+_)', "")
tims_pep_target["type"] = tims_pep_target["type"].str.replace(r'(_W6-32_17%_$)', "")

orbi_pep_target["type"] = orbi_pep_target["Sample"].str.replace(r'(^\d+_[A-Z]+_[a-z]+_)', "")
orbi_pep_target["type"] = orbi_pep_target["type"].str.replace(r'(_W6-\d+_\d+%_[a-z]+_DDA_$)', "")
orbi_pep_target["type"] = orbi_pep_target["type"].str.replace(r'(^UDN\d+_)', "")

join_col = ["type","peptide"]

orbi_target = orbi_pep_target[orbi_pep_target["q-value"] < 0.01]
tims_target = tims_pep_target[tims_pep_target["q-value"] < 0.01]
merged_df = tims_target.merge(orbi_target, how="inner", on=join_col)

orbi_target = orbi_target.drop_duplicates(subset=["type", "peptide"])
tims_target = tims_target.drop_duplicates(subset=["type", "peptide"])
merged_df = merged_df.drop_duplicates(subset=["type", "peptide"])

len(orbi_target)
len(tims_target)
len(tims_target)/len(orbi_target)

merged_df["Value"] = merged_df.groupby('type')['peptide'].transform("count")
orbi_target["Value"] = orbi_target.groupby('type')['peptide'].transform("count")
tims_target["Value"] = tims_target.groupby('type')['peptide'].transform("count")

merged_df["Label"] = "Shared"
orbi_target["Label"] = "Orbitrap"
tims_target["Label"] = "timsTOF"

merged = merged_df[["Value", "type", "Label"]].drop_duplicates().reset_index()[["Value", "type", "Label"]].sort_values("type")
orbi = orbi_target[["Value", "type", "Label"]].drop_duplicates().reset_index()[["Value", "type", "Label"]].sort_values("type")
tims = tims_target[["Value", "type", "Label"]].drop_duplicates().reset_index()[["Value", "type", "Label"]].sort_values("type")

tims["Value"] = tims["Value"] + orbi["Value"]
merged["Value"] = merged["Value"] + orbi["Value"]
# prosit["Value"] = prosit["Value"]

tims = tims.sort_values("type")
merged = merged.sort_values("type")
orbi = orbi.sort_values("type")

reorderlist = ["Brain", "Spleen", "Lung", "Thymus", "Liver", "RCC", "HNSCC", "PBMC", "CLL_01", "CLL_02"]
full_df = tims.append(merged).append(orbi)
# full_df = full_df.sort_values(by=["Sample", "Label"], key=[sorter, sorter_label])
# full_df['Value'] = np.where(full_df['Label'] == 'Lost', -full_df['Value'], full_df['Value'])

plt.figure()
sns.set_context("talk")
# width = 4.5
# height = 6.0
width = 4.5
height = 4.8

fig, ax = plt.subplots(figsize=(width, height))
# ax = sns.barplot(data=full_df, x='type', y='Value', hue='Label', dodge=False, palette=["#cdeac0", "#7a8db3", "#0e1c36"], order = reorderlist)
ax = sns.barplot(data=full_df_rescore, x='type', y='Value', hue='Label', dodge=False, palette=["#cdeac0", "#7a8db3", "#0e1c36"], order = reorderlist)
ax.axhline(0, lw=2, color='black')
ax.legend(loc="center right", bbox_to_anchor=(1.8, 0.5))#,fancybox=True, ncol=1)
plt.ylim(0, max_value)
plt.xticks(rotation = 45) # Rotates X-Axis Ticks by 45-degrees
ylabels = ["{:,.0f}".format(x) + "K" for x in ax.get_yticks()/1000]
ax.set_yticklabels(ylabels)
plt.ylabel('unique HLA-I ligands')
plt.xlabel("Sample")

sns.despine()
# plot_name = "Samples-MaxQuant"
plot_name = "Samples-Rescore"
plt.savefig(directory + f"/{plot_name}_1%_FDR.png", dpi=300, bbox_inches="tight")

# --------------------- Q-value ---------------------
directory = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/Figures/q-value"

root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/raw-cid"
orbi_psm_target_rescore = load_files_to_dfs(root_folder, "rescore_target.psms")
orbi_psm_target = load_files_to_dfs(root_folder, "original_target.psms")
plot_name = "orbitrap"

len(orbi_psm_target_rescore[orbi_psm_target_rescore["q-value"] < 0.01])
len(orbi_psm_target[orbi_psm_target["q-value"] < 0.01])

x = []
prosit_psms_table = []
andromeda_psms_table = []

fdr = 0.01
for i in range(0,51):
    fdr = round((i/1000), 3)
    x.append(fdr)
    prosit_psms_table.append(len(orbi_psm_target_rescore[orbi_psm_target_rescore["q-value"]<=fdr]))
    andromeda_psms_table.append(len(orbi_psm_target[orbi_psm_target["q-value"]<=fdr]))

def plot_q_value(prosit_psms_table: pd.DataFrame, andromeda_psms_table: pd.DataFrame, directory: str, plot_name: str):
    plt.figure()
    plt.style.use(["seaborn-white", "seaborn-paper"])
    sns.set_context("talk")
    width = 8.2
    height = 5.5
    # height = width / 1.618
    fig, ax = plt.subplots(figsize=(width, height))
    ax.plot(x, prosit_psms_table,"-", label = "MaxQuant-Rescored", color="#cdeac0")
    ax.plot(x, andromeda_psms_table,"-", label = "MaxQuant", color="#0e1c36")
    ax.set_xlabel("q-value")
    ax.set_ylabel("Number of PSMs")
    # ax.set_ylim([0, max(prosit_psms_table)])
    plt.xlim([0, max(x)+0.002])
    plt.ylim([0, 300000])
    ax.legend(title_fontsize=10,loc="center right", bbox_to_anchor=(1, 0.2),fancybox=True, ncol=1)
    # ax.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol="%", is_latex=False))
    ax.set_yticks(np.arange(0, max(prosit_psms_table)+1, 100000))
    ylabels = ["{:,.0f}".format(x) + "K" for x in ax.get_yticks()/1000]
    ax.set_yticklabels(ylabels)
    sns.despine()
    plt.savefig(directory + f"/{plot_name}_q-value.png", dpi=300, bbox_inches="tight")
    # plt.savefig(directory + f"/{plot_name}_q-value-MaxQuant.png", dpi=300, bbox_inches="tight")
    plt.close()

plot_q_value(prosit_psms_table, andromeda_psms_table, directory, plot_name)

# --------------------- Q-value ---------------------
directory = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/Figures/q-value"
x = []
prosit_psms_table = []
andromeda_psms_table = []

fdr = 0.01
for i in range(0,51):
    fdr = round((i/1000), 3)
    x.append(fdr)
    prosit_psms_table.append(len(prosit_psms_target[prosit_psms_target["q-value"]<=fdr]))
    andromeda_psms_table.append(len(andromeda_psms_target[andromeda_psms_target["q-value"]<=fdr]))

plot_q_value(prosit_psms_table, andromeda_psms_table, directory, plot_name)

# --------------------- FDR Plot ---------------------
output_dir = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/Figures/q-value"

root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/d-tims"
tims_psms_target = load_files_to_dfs(root_folder, "rescore_target.psms")
andromeda_psms_target = load_files_to_dfs(root_folder, "original_target.psms")

root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/d-cid"
cid_psms_target = load_files_to_dfs(root_folder, "rescore_target.psms")

root_folder = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/d-hcd"
hcd_psms_target = load_files_to_dfs(root_folder, "rescore_target.psms")

tims_target = tims_psms_target[tims_psms_target["q-value"] < 0.10]
cid_target = cid_psms_target[cid_psms_target["q-value"] < 0.10]
hcd_target = hcd_psms_target[hcd_psms_target["q-value"] < 0.10]

plot_name = "compare_models_FDR"

directory = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/Figures/q-value"
x = []
tims_psms_table = []
cid_psms_table = []
hcd_psms_table = []

fdr = 0.00001
for i in range(0, 100001):
    fdr = round((i/1000000), 6)
    x.append(fdr)
    tims_psms_table.append(len(tims_target[tims_target["q-value"]<=fdr]))
    cid_psms_table.append(len(cid_target[cid_target["q-value"]<=fdr]))
    hcd_psms_table.append(len(hcd_target[hcd_target["q-value"]<=fdr]))

tims_psms_table = tims_psms_table[::-1]
cid_psms_table = cid_psms_table[::-1]
hcd_psms_table = hcd_psms_table[::-1]

x[::-1]
plt.figure()

plt.style.use(["seaborn-white", "seaborn-paper"])
sns.set_context("talk")

width = 9
height = width / 1.618
fig, ax = plt.subplots(figsize=(width, height))

ax.plot(x, tims_psms_table,"-", label = "TIMS", color="#0e1c36", linewidth=1)
ax.plot(x, cid_psms_table,"-", label = "CID", color="#7a8db3", linewidth=1)
ax.plot(x, hcd_psms_table,"--", label = "HCD", color="#7a8db3", linewidth=1)
ax.set_xlabel("q-value")
ax.set_ylabel("Number of PSMs")
ax.set_ylim([0, max(prosit_psms_table)])
ax.legend(title_fontsize=10,loc="center right", bbox_to_anchor=(0.4, 0.2),fancybox=True, ncol=1)
# ax.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol="%", is_latex=False))
ax.set_yticks(np.arange(0, max(tims_psms_table)+10000, 100000))
ylabels = ["{:,.0f}".format(x) + "K" for x in ax.get_yticks()/1000]
ax.set_yticklabels(ylabels)
sns.despine()

plt.savefig(directory + f"/{plot_name}_q-value.png", dpi=300, bbox_inches="tight")
# plt.savefig(directory + f"/{plot_name}_q-value-MaxQuant.png", dpi=300, bbox_inches="tight")
plt.close()

