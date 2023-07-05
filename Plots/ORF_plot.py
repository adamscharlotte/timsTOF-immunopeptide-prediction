# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
import seaborn as sns
import os

from Bio import SeqIO
import itertools
from pandas.api.types import CategoricalDtype

def read_file_list(file_list, merge_col):
    gained_shared_lost = []
    maxquant = []
    rescored = []
    for file in file_list:
        psm_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/rescore-tims/" + file + ".psms"
        max_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/andromeda/" + file + ".psms"
        msms_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/msms/" + file + ".txt"
        df = pd.read_csv(psm_file, sep="\t")
        df_max = pd.read_csv(max_file, sep="\t")
        df["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in df["PSMId"]]
        df_max["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in df_max["PSMId"]]
        df["unmod_peptides"] = df["proteinIds"].str.replace(r"\[UNIMOD:\d+\]", "", regex=True)
        df_max["unmod_peptides"] = df_max["proteinIds"].str.replace(r"\[UNIMOD:\d+\]", "", regex=True)
        df["Scan number"] = df["PSMId"].str.extract(r"-(\d+)-")
        df_max["Scan number"] = df_max["PSMId"].str.extract(r"-(\d+)-")
        psm_df = df[df["q-value"] < 0.01]
        max_df = df_max[df_max["q-value"] < 0.01]
        # psm_df = df[df["q-value"] < 0.01].drop_duplicates(subset=["RAW_FILE", "unmod_peptides"])
        # max_df = df_max[df_max["q-value"] < 0.01].drop_duplicates(subset=["RAW_FILE", "unmod_peptides"])
        msms_df = pd.read_csv(msms_file, sep="\t")
        msms_df["Scan number"] = msms_df["Scan number"].astype(str)
        psm_msms = pd.merge(psm_df, msms_df, on="Scan number", how="left").drop_duplicates(subset=merge_col)
        max_msms = pd.merge(max_df, msms_df, on="Scan number", how="left").drop_duplicates(subset=merge_col)
        merged_df = pd.merge(psm_msms, max_msms, on=merge_col, how="inner")
        indicator_df = pd.merge(psm_msms, max_msms, on=merge_col, how='outer', indicator=True)
        gained_df = indicator_df.loc[indicator_df._merge == 'left_only', merge_col]
        lost_df = indicator_df.loc[indicator_df._merge == 'right_only', merge_col]
        merged_df["Label"] = "Shared"
        gained_df["Label"] = "Gained"
        lost_df["Label"] = "Lost"
        gained_shared_lost.append(merged_df)
        gained_shared_lost.append(gained_df)
        gained_shared_lost.append(lost_df)
        maxquant.append(max_msms)
        rescored.append(psm_msms)
    lsg_df = pd.concat(gained_shared_lost, ignore_index=True)
    maxquant_df = pd.concat(maxquant, ignore_index=True)
    rescored_df = pd.concat(rescored, ignore_index=True)
    return lsg_df, maxquant_df, rescored_df

def map_file_list(file_list):
    orf = []
    for file in file_list:
        psm_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/rescore-tims/" + file + ".psms"
        msms_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/msms/" + file + ".txt"
        df = pd.read_csv(psm_file, sep="\t")
        psm_df = df[df["q-value"] < 0.01]
        psm_df["unmod_peptides"] = psm_df["proteinIds"].str.replace(r"\[UNIMOD:\d+\]", "", regex=True)
        psm_df["Scan number"] = psm_df["PSMId"].str.extract(r"-(\d+)-")
        msms_df = pd.read_csv(msms_file, sep="\t")
        msms_df["Scan number"] = msms_df["Scan number"].astype(str)
        merged_df = pd.merge(psm_df, msms_df, on="Scan number", how="left")
        filtered_df = merged_df[~merged_df['Proteins'].str.contains(';', na=True)]
        # Map full fasta headers to the protein names
        filtered_df["full_prot"] = filtered_df.Proteins.map(header_dict)
        orf_df = filtered_df[filtered_df["full_prot"].str.contains("nuORFdb_human", na=False)]
        # orf_count = len(orf_df)
        # logging.info(f"There are {orf_count} nuORFs identified.")
        # print(f"There are {orf_count} nuORFs identified.")
        # pd.options.display.max_colwidth = 100
        # print(orf_df[["q-value", "full_prot"]])
        orf.append(orf_df)
    full_orf_df = pd.concat(orf, ignore_index=True)
    return full_orf_df

# ---------------------- Import data ----------------------
# Data needed for the ORF Plot

file_list = ["E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep1_Slot1-06_1_3498",
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
             "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496",
             "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep2_Slot2-5_1_3520",
             "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep3_Slot2-5_1_3521",
             "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep4_Slot2-5_1_3522"]

headers = []
with open("/Users/adams/Projects/300K/MSV000091456-SCP/fasta/Gencode_v34_3nr.602contams.2043smorfs.nuORFv1.fasta", "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        headers.append(record.description)

headers_split_0 = [v.split(" ", 1)[0] for v in headers]
headers_split_1 = [v.split(" ", 1)[1:] for v in headers]
headers_split_1 = list(itertools.chain(*headers_split_1))
header_dict = dict(zip(headers_split_0, headers_split_1))

full_orf_df = map_file_list(file_list)
full_orf_df["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in full_orf_df["PSMId"]]

full_orf_df.value_counts("full_prot")

lsg_orf_df, maxquant_df, rescored_df = read_file_list(file_list, ["RAW_FILE", "Proteins"])
lsg_orf_df[["RAW_FILE", "Proteins", "Label"]].drop_duplicates()

maxquant_df.value_counts("Label")
lsg_orf_df.value_counts("Label")

def merge_orf(full_orf_df, lsg_orf_df):
    orf_merged = pd.merge(full_orf_df, lsg_orf_df, on = ["RAW_FILE", "Proteins"])
    orf_merged["Cell_input"] = orf_merged["RAW_FILE"].str.replace(r'_directIP_titration_rep\d_Slot\d+-\d+_\d+_\d+$','')
    orf_merged["Cell_input"] = orf_merged["Cell_input"].str.replace(r'E_\d+_NO\d+_\d+nL_HLAc1_','')
    orf_merged["full_prot"] = orf_merged["full_prot"].str.replace('nuORFdb_human','')
    orf_merged["full_prot"] = orf_merged["full_prot"].str.replace('|','')
    orf_merged_small = orf_merged[~orf_merged['full_prot'].isin(["Sense Overlapping","TUCP","TEC","Other"])]
    return orf_merged_small

orf_merged = merge_orf(full_orf_df, lsg_orf_df)
orf_max = merge_orf(full_orf_df, maxquant_df)
orf_rescored = merge_orf(full_orf_df, rescored_df)

orf_rescored.value_counts("full_prot")
orf_max.value_counts("full_prot")

custom_order = ['1e6', '5e6', '1e7', '2e7', '4e7']
orf_rescored['Cell_input'] = orf_rescored['Cell_input'].astype(CategoricalDtype(categories=custom_order, ordered=True))
orf_max['Cell_input'] = orf_max['Cell_input'].astype(CategoricalDtype(categories=custom_order, ordered=True))
sorted_full_prot = ["3' dORF", "3' Overlap dORF","5' Overlap uORF", "5' uORF", 'Antisense', 'lincRNA', 'ncRNA Processed Transcript', 'ncRNA Retained Intron', 'Out-of-Frame', 'Pseudogene']
sorted_full_prot.reverse()
orf_rescored['full_prot'] = orf_rescored['full_prot'].astype(CategoricalDtype(categories=sorted_full_prot, ordered=True))
orf_max['full_prot'] = orf_max['full_prot'].astype(CategoricalDtype(categories=sorted_full_prot, ordered=True))

stacked_rescored = orf_rescored.groupby(['Cell_input', 'full_prot']).size().unstack()
stacked_max = orf_max.groupby(['Cell_input', 'full_prot']).size().unstack()

stacked_rescored.iloc[0:4] = stacked_rescored.iloc[0:4].div(3)
stacked_rescored.iloc[4] = stacked_rescored.iloc[4].div(4)
stacked_max.iloc[0:4] = stacked_max.iloc[0:4].div(3)
stacked_max.iloc[4] = stacked_max.iloc[4].div(4)

# ---------- Plot ----------

sns.set_context("paper")
sns.set_style("whitegrid", {'axes.grid' : False})

cm = 1/2.54  # centimeters in inches
width = 18*cm
height = 18*cm
fig = plt.figure(figsize=(width, height))
axes = fig.subplot_mosaic([
     ['a', 'b']])

# ORF Plot
color_list = ['#0e1c36', '#DCD9F4', '#6DED80', '#2F2963', '#CDEAC0', '#746EA6', '#4B7F52', '#FFA69E', '#eee3ca', '#f33b16']
color_list.reverse()
stacked_max.plot.bar(stacked=True, color=color_list, ax = axes["a"])
axes["a"].set_ylabel('Unique nuORF source proteins')
axes['a'].set_xlabel('Cell input')
axes["a"].set_xticklabels(['1e6', '5e6', '1e7', '2e7', '4e7'], rotation=0)

# reordering the labels
handles, labels = axes["a"].get_legend_handles_labels()
order = [9,8,7,6,5,4,3,2,1,0]
axes["a"].legend([handles[i] for i in order], [labels[i] for i in order],
                 loc='upper left', frameon=False, ncol = 1,
                 bbox_to_anchor=(0, 1.1), fontsize="7",
                 title="nuORF category")

stacked_rescored.plot.bar(stacked=True, color=color_list, ax = axes["b"])
axes["b"].set_ylabel('Unique nuORF source proteins')
axes['b'].set_xlabel('Cell input')
axes["b"].set_xticklabels(['1e6', '5e6', '1e7', '2e7', '4e7'], rotation=0)

# reordering the labels
handles, labels = axes["b"].get_legend_handles_labels()
order = [9,8,7,6,5,4,3,2,1,0]
axes["b"].legend([handles[i] for i in order], [labels[i] for i in order],
                 loc='upper left', frameon=False, ncol = 1,
                 bbox_to_anchor=(0, 1.1), fontsize="7",
                 title="nuORF category")

sns.despine()
# fig.tight_layout(pad=0.4, w_pad=-10, h_pad=0.2)
# fig.tight_layout(pad=0.5, w_pad=1, h_pad=0.2)
fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

for label, ax in axes.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    ax.text(0.0, 1, label, transform=ax.transAxes + trans,
            fontsize='x-large', va='bottom', weight='bold')

plot_path = "/Users/adams/Projects/300K/Results/Figures/paper-fig-orf.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
