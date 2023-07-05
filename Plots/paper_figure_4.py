# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
import seaborn as sns
import os
import logomaker
from logomaker.src.matrix import transform_matrix
from Bio import SeqIO
import itertools
from pandas.api.types import CategoricalDtype

def load_scp_to_dfs(root_folder, file_pattern):
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

def plot_mean(prosit_target: pd.DataFrame, andromeda_target: pd.DataFrame, join_col, plot_name: str, type: str, ax_1):
    prosit_df = scp_prosit_pep[scp_prosit_pep["q-value"] < 0.01]
    adromeda_df = scp_adromeda_pep[scp_adromeda_pep["q-value"] < 0.01]
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
    order_x = ["1e6", "5e6", "1e7", "2e7", "4e7"]
    sns.barplot(data=full_df, x='Sample', y='Value', hue='Label', capsize=.4, errorbar=('ci', 68), errwidth=1,
        palette=["#0e1c36", "#cdeac0"], order=order_x, ax = ax_1)
    sns.despine()

def read_file_list(file_list, merge_col):
    gained_shared_lost = []
    # maxquant = []
    # rescored = []
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
        # maxquant.append(max_msms)
        # rescored.append(psm_msms)
    lsg_df = pd.concat(gained_shared_lost, ignore_index=True)
    # maxquant_df = pd.concat(maxquant, ignore_index=True)
    # rescored_df = pd.concat(rescored, ignore_index=True)
    return lsg_df#, maxquant_df, rescored_df

def read_hla_list(hla_list):
    hla_dfs = pd.DataFrame()
    for hla in hla_list:
        hla_file = "/Users/adams/Code/timsTOF-immunopeptide-prediction/NetMHCpan/" + hla
        hla_df = pd.read_csv(hla_file, sep=",")
        hla_df = hla_df.drop(columns=['Core', 'Of', 'Gp', 'Gl', 'Ip', 'Il', 'Icore',
        'Identity', 'Score_EL'])
        hla_dfs = pd.concat([hla_dfs.reset_index(drop=True), hla_df], axis=1)
    hla_dfs = hla_dfs.T.drop_duplicates().T
    return hla_dfs

def merge_hla(hla_dfs, lsg_df):
    hla_merged = pd.merge(hla_dfs, lsg_df, left_on="peptide", right_on="unmod_peptides")
    hla_merged["min_binding_score"] = hla_merged[['HLA-C-16_02', 'HLA-A-01_01', 'HLA-A-02_02', 'HLA-B-44_03',
           'HLA-B-57_01', 'HLA-C-06_02']].min(axis=1)
    hla_merged["log_min_binding_score"] = np.log10(hla_merged[["min_binding_score"]])
    return(hla_merged)

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

# Data needed for the mean SD plot
root_folder = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/peptides/"
scp_prosit_pep = load_scp_to_dfs(root_folder, "rescore_target.peptides")
scp_adromeda_pep = load_scp_to_dfs(root_folder, "original_target.peptides")

prosit_df = scp_prosit_pep[scp_prosit_pep["q-value"] < 0.01]
adromeda_df = scp_adromeda_pep[scp_adromeda_pep["q-value"] < 0.01]
adromeda_df = adromeda_df.drop_duplicates(subset=["RAW_FILE", "peptide"])
prosit_df = prosit_df.drop_duplicates(subset=["RAW_FILE", "peptide"])
prosit_df["Value"] = prosit_df.groupby('RAW_FILE')["peptide"].transform("count")
adromeda_df["Value"] = adromeda_df.groupby('RAW_FILE')["peptide"].transform("count")
adromeda_df["Label"] = "MaxQuant"
prosit_df["Label"] = "Rescored"
andromeda = adromeda_df[["Value", "RAW_FILE", "Label"]].drop_duplicates().reset_index()[["Value", "RAW_FILE", "Label"]].sort_values("RAW_FILE")
prosit = prosit_df[["Value", "RAW_FILE", "Label"]].drop_duplicates().reset_index()[["Value", "RAW_FILE", "Label"]].sort_values("RAW_FILE")
prosit = prosit.sort_values("RAW_FILE")
andromeda = andromeda.sort_values("RAW_FILE")
full_df = andromeda.append(prosit)
full_df["Sample"] = full_df["RAW_FILE"].str[28:-30]
prosit["Sample"] = prosit["RAW_FILE"].str[28:-30]
andromeda["Sample"] = andromeda["RAW_FILE"].str[28:-30]
prosit.groupby(['Sample'])['Value'].mean()/andromeda.groupby(['Sample'])['Value'].mean()

# Lists needed to generate the logos

# file_list = ["E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep1_Slot1-06_1_3498",
#              "E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep3_Slot2-1_1_3508",
#              "E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep4_Slot2-1_1_3509",
#              "E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep1_Slot1-37_1_3502",
#              "E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep2_Slot2-2_1_3510",
#              "E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep3_Slot2-2_1_3511",
#              "E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep1_Slot2-3_1_3512",
#              "E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep2_Slot2-3_1_3513",
#              "E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep3_Slot2-3_1_3514",
#              "E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep1_Slot2-4_1_3515",
#              "E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep3_Slot2-4_1_3518",
#              "E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep4_Slot2-4_1_3519",
#              "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496",
#              "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep2_Slot2-5_1_3520",
#              "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep3_Slot2-5_1_3521",
#              "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep4_Slot2-5_1_3522"]

# lsg_df = read_file_list(file_list, ["RAW_FILE", "unmod_peptides"])[["unmod_peptides", "Label"]].drop_duplicates()
# lsg_df[lsg_df["Label"] == "Shared"].to_csv(sep=',', columns=["unmod_peptides"],header = False, index = False, path_or_buf="/Users/adams/Projects/300K/MSV000091456-SCP/shared_peptides.txt")
# lsg_df[lsg_df["Label"] == "Gained"].to_csv(sep=',', columns=["unmod_peptides"],header = False, index = False, path_or_buf="/Users/adams/Projects/300K/MSV000091456-SCP/gained_peptides.txt")
# lsg_df[lsg_df["Label"] == "Lost"].to_csv(sep=',', columns=["unmod_peptides"],header = False, index = False, path_or_buf="/Users/adams/Projects/300K/MSV000091456-SCP/lost_peptides.txt")

# Data needed for the logo plots

gained_psfm = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_gained.txt", sep=" ")
shared_psfm = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_shared.txt", sep=" ")
lost_psfm = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_lost.txt", sep=" ")
psfm_0101 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_0101.txt", sep=" ")

inf_lost_matrix = transform_matrix(lost_psfm,
                 from_type = "probability",
                 to_type = "information")
inf_gained_matrix = transform_matrix(gained_psfm,
                 from_type = "probability",
                 to_type = "information")
inf_shared_matrix = transform_matrix(shared_psfm,
                 from_type = "probability",
                 to_type = "information")
inf_0101 = transform_matrix(psfm_0101,
                 from_type = "probability",
                 to_type = "information")

# Data needed for the CDF plot

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

hla_list = ["all-HLA-C16_02.txt",
            "all-HLA-A01_01.txt",
            "all-HLA-A02_02.txt",
            "all-HLA-B44_03.txt",
            "all-HLA-B57_01.txt",
            "all-HLA-C06_02.txt"]

lsg_df = read_file_list(file_list, ["RAW_FILE", "unmod_peptides"])
hla_dfs = read_hla_list(hla_list)
hla_merged = merge_hla(hla_dfs, lsg_df)
counts = hla_merged.value_counts("Label")
counts_wb = hla_merged[hla_merged["min_binding_score"] < 2].value_counts("Label")
hla_merged[hla_merged["log_min_binding_score"] < np.log10(2)].value_counts("Label")
counts_sb = hla_merged[hla_merged["min_binding_score"] < 0.5].value_counts("Label")

counts_wb/counts
counts_sb/counts
orf_merged.unmod_peptides
hla_merged.columns
hla_orf = pd.merge(orf_merged, hla_merged, on = ["RAW_FILE", "unmod_peptides"])
counts = hla_orf.value_counts("Label_y")
orf_wb = hla_orf[hla_orf["min_binding_score"] <= 2].value_counts("Label_y")
orf_sb = hla_orf[hla_orf["min_binding_score"] <= 0.5].value_counts("Label_y")
orf_wb/counts
orf_sb/counts
len(hla_orf[hla_orf["min_binding_score"] < 2])/len(hla_orf)
len(hla_orf[hla_orf["min_binding_score"] < 0.5])/len(hla_orf)

# Data needed for the ORF Plot
headers = []
with open("/Users/adams/Projects/300K/MSV000091456-SCP/fasta/Gencode_v34_3nr.602contams.2043smorfs.nuORFv1.fasta", "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        headers.append(record.description)

headers_split_0 = [v.split(" ", 1)[0] for v in headers]
headers_split_1 = [v.split(" ", 1)[1:] for v in headers]
headers_split_1 = list(itertools.chain(*headers_split_1))
header_dict = dict(zip(headers_split_0, headers_split_1))

full_orf_df = map_file_list(file_list)
full_orf_df.value_counts("full_prot")
full_orf_df["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in full_orf_df["PSMId"]]

lsg_orf_df = read_file_list(file_list, ["RAW_FILE", "Proteins"])[["RAW_FILE", "Proteins", "Label"]].drop_duplicates()
lsg_orf_df.value_counts("Label")

orf_merged = pd.merge(full_orf_df, lsg_orf_df, on = ["RAW_FILE", "Proteins"])
orf_merged.value_counts("Label")
orf_merged.value_counts("full_prot")
orf_merged.value_counts("Cell_input")

orf_merged["Cell_input"] = orf_merged["Raw file"].str.replace(r'_directIP_titration_rep\d_Slot\d+-\d+_\d+_\d+$','')
orf_merged["Cell_input"] = orf_merged["Cell_input"].str.replace(r'E_\d+_NO\d+_\d+nL_HLAc1_','')
orf_merged["full_prot"] = orf_merged["full_prot"].str.replace('nuORFdb_human','')
orf_merged["full_prot"] = orf_merged["full_prot"].str.replace('|','')
orf_merged_small = orf_merged[~orf_merged['full_prot'].isin(["Sense Overlapping","TUCP","TEC","Other"])]

custom_order = ['1e6', '5e6', '1e7', '2e7', '4e7']
orf_merged_small['Cell_input'] = orf_merged_small['Cell_input'].astype(CategoricalDtype(categories=custom_order, ordered=True))
sorted_full_prot = ["3' dORF", "3' Overlap dORF","5' Overlap uORF", "5' uORF", 'Antisense', 'lincRNA', 'ncRNA Processed Transcript', 'ncRNA Retained Intron', 'Out-of-Frame', 'Pseudogene']
sorted_full_prot.reverse()
orf_merged_small['full_prot'] = orf_merged_small['full_prot'].astype(CategoricalDtype(categories=sorted_full_prot, ordered=True))

stacked_df = orf_merged_small.groupby(['Cell_input', 'full_prot']).size().unstack()
div_stacked_df = stacked_df
div_stacked_df.iloc[0:4] = div_stacked_df.iloc[0:4].div(3)
div_stacked_df.iloc[4] = div_stacked_df.iloc[4].div(4)

# ---------------------- Plot ----------------------
sns.set_context("paper")
sns.set_style("whitegrid", {'axes.grid' : False})

cm = 1/2.54  # centimeters in inches
width = 18*cm
height = 24*cm
# constrained_layout=True, 
fig = plt.figure(figsize=(width, height))
# axes = fig.subplot_mosaic([
#      ['a', 'b'],
#      ['a', ''],
#      ['c', ' '],
#      ['c', '  '],
#      ['d', 'BLANK'],
#      ['d', 'BLANK']], empty_sentinel="BLANK")
axes = fig.subplot_mosaic([
     ['a', 'b'],
     ['a', 'b'],
     ['a', 'b'],
     ['a', ''],
     ['c', ''],
     ['c', ''],
     ['c', ' '],
     ['c', ' '],
     ['d', ' '],
     ['d', '  '],
     ['d', '  '],
     ['d', '  ']])

# Mean SD Plot
plot_mean(scp_prosit_pep, scp_adromeda_pep, ["RAW_FILE", "peptide"],
          "mean-pep", "peptide", axes['a'])

axes['a'].legend(frameon=False).set_title('')
axes['a'].axhline(0, lw=2, color='black')
ylabels = ["{:,.0f}".format(x) + "K" for x in axes['a'].get_yticks()/1000]
axes['a'].set_yticklabels(ylabels)
axes['a'].set_ylabel("Unique HLA-I ligands")
axes['a'].set_xlabel('Cell input')

# Logo Plots
color_dict =  {
        'A': '#0e1c36',
        'C': '#cdeac0',
        'D': '#f33b16',
        'E': '#f33b16',
        'F': '#0e1c36',
        'G': '#cdeac0',
        'H': '#7a8db3',
        'I': '#0e1c36',
        'K': '#7a8db3',
        'L': '#0e1c36',
        'M': '#0e1c36',
        'N': '#b8b3e9',
        'P': '#0e1c36',
        'Q': '#b8b3e9',
        'R': '#7a8db3',
        'S': '#cdeac0',
        'T': '#cdeac0',
        'V': '#0e1c36',
        'W': '#0e1c36',
        'Y': '#cdeac0'
    }

logomaker.Logo(inf_0101,
                color_scheme=color_dict,
                center_values=True,
                vpad=.1,
                width=.8,
                ax = axes['b'])

logomaker.Logo(inf_gained_matrix,
                color_scheme=color_dict,
                center_values=True,
                vpad=.1,
                width=.8,
                ax = axes[''])
logomaker.Logo(inf_shared_matrix,
                color_scheme=color_dict,
                center_values=True,
                vpad=.1,
                width=.8,
                ax = axes[' '])
logomaker.Logo(inf_lost_matrix,
                color_scheme=color_dict,
                center_values=True,
                vpad=.1,
                width=.8,
                ax = axes['  '])

motif_axes = [axes['b'], axes[''], axes[' '], axes['  ']]
for ax in motif_axes:
    ax.set(xticks=np.arange(1, 10))
    ax.set_ylabel("Bits", labelpad=-0.5)

axes['b'].set_title("HLA-A*0101", fontweight = "bold", y=0.9, pad=-1)
axes[''].set_title("Gained", fontweight = "bold", y=0.9, pad=-1)
axes[' '].set_title("Shared", fontweight = "bold", y=0.9, pad=-1)
axes['  '].set_title("Lost", fontweight = "bold", y=0.9, pad=-1)
axes['  '].set_ylabel("Bits", labelpad=-2)
axes['  '].set_xlabel("Position")

list_color = ["#f33b16", "#7a8db3", "#0e1c36", "#b8b3e9", "#cdeac0"]
list_lab = ['Acidic','Basic','Hydrophobic', 'Neutral', 'Polar']

handles, labels = axes["b"].get_legend_handles_labels()
for col, lab in zip(list_color, list_lab):
    patch = mpatches.Patch(color=col, label=lab)
    handles.append(patch) 

# plot the legend
axes['b'].legend(handles=handles, loc='upper left', ncol=5,
                handlelength=1, handleheight=1, handletextpad=0.4,
                bbox_to_anchor=(0.05, 1.18),
                fontsize="7", frameon=False, columnspacing = 0.7)

# CDF affinity Plot
sns.ecdfplot(data=hla_merged, x="log_min_binding_score", hue="Label",
             palette=["#7A8DB3", "#CDEAC0", "#0e1c36"], ax=axes["c"], legend=True)
sns.move_legend(axes["c"], loc="upper left", frameon=False, title = None)
axes["c"].set_xlabel('Log affinity score')
axes["c"].grid(False, which="both")
axes["c"].axvline(np.log10(2), color="#f33b16") #WB
axes["c"].axvline(np.log10(0.5), color="#f33b16") #SB
axes["c"].text(np.log10(2)-0.2, 1.03, "WB", color="#f33b16")
axes["c"].text(np.log10(0.5)-0.2,1.03, "SB", color="#f33b16")

# ORF Plot
color_list = ['#0e1c36', '#DCD9F4', '#6DED80', '#2F2963', '#CDEAC0', '#746EA6', '#4B7F52', '#FFA69E', '#eee3ca', '#f33b16']
color_list.reverse()
div_stacked_df.plot.bar(stacked=True, color=color_list, ax = axes["d"])
axes["d"].set_ylabel('Unique nuORF source proteins')
axes['d'].set_xlabel('Cell input')
axes["d"].set_xticklabels(['1e6', '5e6', '1e7', '2e7', '4e7'], rotation=0)

# reordering the labels
handles, labels = axes["d"].get_legend_handles_labels()
order = [9,8,7,6,5,4,3,2,1,0]
leg = axes["d"].legend([handles[i] for i in order], [labels[i] for i in order],
                 loc='upper left', frameon=False, ncol = 1,
                 bbox_to_anchor=(0, 1), fontsize="7",
                 title="nuORF category")#, alignment = "left")

leg._legend_box.align = "left"
sns.despine()
# fig.tight_layout(pad=0.4, w_pad=-10, h_pad=0.2)
# fig.tight_layout(pad=0.5, w_pad=1, h_pad=0.2)
fig.tight_layout(pad=0.5)

for label, ax in axes.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    ax.text(0.0, 1, label, transform=ax.transAxes + trans,
            fontsize='x-large', va='bottom', weight='bold')

plot_path = "/Users/adams/Projects/300K/Results/Figures/paper-fig-4.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")




# lost_matrix = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/lost_matrix.txt", sep=" ")
# gained_matrix = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/gained_matrix.txt", sep=" ")
# shared_matrix = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/shared_matrix.txt", sep=" ")

# lsg_df, maxquant_df, rescored_df = read_file_list(file_list, ["RAW_FILE", "Proteins"])
# orf_max = pd.merge(full_orf_df, maxquant_df, on = ["RAW_FILE", "Proteins"])
# orf_rescored = pd.merge(full_orf_df, rescored_df, on = ["RAW_FILE", "Proteins"])

# orf_max.value_counts("full_prot")
# orf_rescored.value_counts("full_prot")


# def kl(p, q):
#     """Kullback-Leibler divergence D(P || Q) for discrete distributions
#     Parameters
#     ----------
#     p, q : array-like, dtype=float, shape=n
#     Discrete probability distributions.
#     """
#     p = np.asarray(p, dtype=np.float)
#     q = np.asarray(q, dtype=np.float)
#     return np.sum(np.where(p != 0, p * np.log(p / q), 0))

# kl(inf_gained_matrix, inf_0101)
# kl(inf_0101, inf_shared_matrix)
# kl(inf_shared_matrix, inf_0101)
# kl(inf_lost_matrix, inf_0101)
# kl(inf_0101, inf_lost_matrix)
