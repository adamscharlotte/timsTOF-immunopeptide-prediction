# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import pandas as pd
import logging
import re
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker

from Bio import SeqIO
import itertools
from pandas.api.types import CategoricalDtype


# Set up input and output file names
file = "E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep4_Slot2-1_1_3509"
psm_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/rescore-tims/" + file + ".psms"
msms_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/msms/" + file + ".txt"
max_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/andromeda/" + file + ".psms"

# Read fasta headers
headers = []
with open("/Users/adams/Projects/300K/MSV000091456-SCP/fasta/Gencode_v34_3nr.602contams.2043smorfs.nuORFv1.fasta", "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        headers.append(record.description)

headers_split_0 = [v.split(" ", 1)[0] for v in headers]
headers_split_1 = [v.split(" ", 1)[1:] for v in headers]
headers_split_1 = list(itertools.chain(*headers_split_1))
header_dict = dict(zip(headers_split_0, headers_split_1))

def map_to_msms(psm_file, msms_file):
    df = pd.read_csv(psm_file, sep="\t")
    psm_df = df[df["q-value"] < 0.01]
    psm_df["unmod_peptides"] = psm_df["proteinIds"].str.replace(r"\[UNIMOD:\d+\]", "", regex=True)
    psm_df["Scan number"] = psm_df["PSMId"].str.extract(r"-(\d+)-")
    msms_df = pd.read_csv(msms_file, sep="\t")
    msms_df["Scan number"] = msms_df["Scan number"].astype(str)
    merged_df = pd.merge(psm_df, msms_df, on="Scan number", how="left")
    filtered_df = merged_df.dropna(subset=['Proteins'])[merged_df["Proteins"].str.startswith('SNP')]
    snp_df = filtered_df[~filtered_df['Proteins'].str.contains(';')]
    print(filtered_df[["Proteins", "Sequence"]])
    snp_count = len(snp_df)
    logging.info(f"There are {snp_count} missense SNPs identified.")
    print(f"There are {snp_count} missense SNPs identified.")
    pd.options.display.max_colwidth = 100
    print(snp_df[["q-value", "score", "PSMId", "Proteins"]])

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

def read_file_list(file_list):
    gained_shared_lost = []
    for file in file_list:
        # psm_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/small-fasta/rescore-tims/" + file + ".psms"
        # max_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/small-fasta/andromeda/" + file + ".psms"
        psm_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/rescore-tims/" + file + ".psms"
        max_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/andromeda/" + file + ".psms"
        df = pd.read_csv(psm_file, sep="\t")
        df_max = pd.read_csv(max_file, sep="\t")
        df["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in df["PSMId"]]
        df_max["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in df_max["PSMId"]]
        df["unmod_peptides"] = df["proteinIds"].str.replace(r"\[UNIMOD:\d+\]", "", regex=True)
        df_max["unmod_peptides"] = df_max["proteinIds"].str.replace(r"\[UNIMOD:\d+\]", "", regex=True)
        psm_df = df[df["q-value"] < 0.01].drop_duplicates(subset=["RAW_FILE", "unmod_peptides"])
        max_df = df_max[df_max["q-value"] < 0.01].drop_duplicates(subset=["RAW_FILE", "unmod_peptides"])
        merged_df = pd.merge(psm_df, max_df, on=["RAW_FILE", "unmod_peptides"], how="inner")
        indicator_df = pd.merge(psm_df, max_df, on=["RAW_FILE", "unmod_peptides"], how='outer', indicator=True)
        gained_df = indicator_df.loc[indicator_df._merge == 'left_only', ["RAW_FILE", "unmod_peptides"]]
        lost_df = indicator_df.loc[indicator_df._merge == 'right_only', ["RAW_FILE", "unmod_peptides"]]
        merged_df["Label"] = "Shared"
        gained_df["Label"] = "Gained"
        lost_df["Label"] = "Lost"
        gained_shared_lost.append(merged_df)
        gained_shared_lost.append(gained_df)
        gained_shared_lost.append(lost_df)
    lsg_df = pd.concat(gained_shared_lost, ignore_index=True)
    return lsg_df

def read_max_rescore(file_list):
    max_rescore = []
    for file in file_list:
        psm_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/rescore-tims/" + file + ".psms"
        max_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/andromeda/" + file + ".psms"
        df = pd.read_csv(psm_file, sep="\t")
        df_max = pd.read_csv(max_file, sep="\t")
        df["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in df["PSMId"]]
        df_max["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in df_max["PSMId"]]
        df["unmod_peptides"] = df["proteinIds"].str.replace(r"\[UNIMOD:\d+\]", "", regex=True)
        df_max["unmod_peptides"] = df_max["proteinIds"].str.replace(r"\[UNIMOD:\d+\]", "", regex=True)
        psm_df = df[df["q-value"] < 0.01].drop_duplicates(subset=["RAW_FILE", "unmod_peptides"])
        max_df = df_max[df_max["q-value"] < 0.01].drop_duplicates(subset=["RAW_FILE", "unmod_peptides"])
        psm_df["Label"] = "Rescored"
        max_df["Label"] = "MaxQuant"
        max_rescore.append(psm_df)
        max_rescore.append(max_df)
        max_rescore_df = pd.concat(max_rescore, ignore_index=True)
    return max_rescore_df

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

# file_list = ["E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496",
#              "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep2_Slot2-5_1_3520",
#              "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep3_Slot2-5_1_3521",
#              "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep4_Slot2-5_1_3522"]

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

full_orf_df = map_file_list(file_list)
full_orf_df.value_counts("full_prot")

lsg_df = read_file_list(file_list)[["RAW_FILE", "unmod_peptides", "Label"]].drop_duplicates()
lsg_df.value_counts("RAW_FILE", sort=True, ascending=False, dropna=False)
lsg_df.value_counts("Label")

full_orf_df["RAW_FILE"] = ["-".join(item.split("-")[0:2]) for item in full_orf_df["PSMId"]]
lgs_orf = pd.merge(full_orf_df, lsg_df, on=["RAW_FILE", "unmod_peptides"])
unique_lgs_orf = lgs_orf[["RAW_FILE", "Label", "Proteins", "unmod_peptides"]].drop_duplicates()
unique_lgs = lsg_df[["RAW_FILE", "Label", "unmod_peptides"]].drop_duplicates()

# 2% of immunopeptides could be mapped to nuORFs
2251/102539

unique_lgs_orf = lgs_orf[["RAW_FILE", "Label", "Proteins"]].drop_duplicates()
unique_lgs_orf.groupby(["Label"])["Label"].count()
unique_lgs_orf.value_counts("Proteins")
unique_lgs_orf[unique_lgs_orf["Proteins"] == "ENST00000617998.4_2_17:7388053-7417318:+|GN=POLR2A"]

full_orf_df["Cell_input"] = full_orf_df["Raw file"].str.replace(r'_directIP_titration_rep\d_Slot\d+-\d+_\d+_\d+$','')
full_orf_df["Cell_input"] = full_orf_df["Cell_input"].str.replace(r'E_\d+_NO\d+_\d+nL_HLAc1_','')
full_orf_df["full_prot"] = full_orf_df["full_prot"].str.replace('nuORFdb_human','')
full_orf_df["full_prot"] = full_orf_df["full_prot"].str.replace('|','')

prot_orf_df = full_orf_df[["Cell_input", "full_prot", "Proteins"]].drop_duplicates()
prot_orf_df.value_counts("Cell_input")
prot_orf_df.value_counts("full_prot")
custom_order = ['1e6', '5e6', '1e7', '2e7', '4e7']
prot_orf_df['Cell_input'] = prot_orf_df['Cell_input'].astype(CategoricalDtype(categories=custom_order, ordered=True))
# unique_full_prot = prot_orf_df['full_prot'].unique()
# sorted_full_prot = sorted(unique_full_prot, key=lambda x: x[0].lower())
sorted_full_prot = ["3' dORF", "3' Overlap dORF", "5' Overlap uORF", "5' uORF", 'Antisense', 'lincRNA', 'ncRNA Processed Transcript', 'ncRNA Retained Intron', 'Out-of-Frame', 'Pseudogene', 'Sense Overlapping', 'TUCP', 'TEC', 'Other']
prot_orf_df['full_prot'] = prot_orf_df['full_prot'].astype(CategoricalDtype(categories=sorted_full_prot, ordered=True))

stacked_df = prot_orf_df.groupby(['Cell_input', 'full_prot']).size().unstack()
#37123C
#2E294E
#5E2BFF
#0df040


width = 7
height = width / 1.618
fig, axes = plt.subplots(1, figsize=(width * 1, height * 1))
axes = np.ravel(axes)
#F97862
axes[0] = stacked_df.plot.bar(stacked=True, color=['#0e1c36', '#454372', '#b8b3e9', '#4B7F52', '#cdeac0', '#6DED80', '#F7C967', '#eee3ca', '#F4AE96','#f33b16', '#FFAE03', '#9DB18E', '#8CB589'])
fig.tight_layout()

plt.savefig('/Users/adams/Projects/300K/MSV000091456-SCP/Figures/4e7-max-vs-rescore.pdf', dpi=300, bbox_inches='tight')
# plt.show()
plt.close()








lsg_df = read_file_list(file_list)
lsg_df.value_counts("Label")

lsg_df[lsg_df["Label"] == "Shared"].to_csv(sep=',', columns=["unmod_peptides"],header = False, index = False, path_or_buf="/Users/adams/Projects/300K/MSV000091456-SCP/shared_peptides.txt")
lsg_df[lsg_df["Label"] == "Gained"].to_csv(sep=',', columns=["unmod_peptides"],header = False, index = False, path_or_buf="/Users/adams/Projects/300K/MSV000091456-SCP/gained_peptides.txt")
lsg_df[lsg_df["Label"] == "Lost"].to_csv(sep=',', columns=["unmod_peptides"],header = False, index = False, path_or_buf="/Users/adams/Projects/300K/MSV000091456-SCP/lost_peptides.txt")

max_rescore_df = read_max_rescore(file_list)
max_rescore_df.value_counts("Label")
max_rescore_df.value_counts("RAW_FILE")

# -------------------- Generate a list of peptides to predict with NetMHCpan --------------------

hla_list = ["all-HLA-C16_02.txt",
            "all-HLA-A01_01.txt",
            "all-HLA-A02_02.txt",
            "all-HLA-B44_03.txt",
            "all-HLA-B57_01.txt",
            "all-HLA-C06_02.txt"]

hla_dfs = read_hla_list(hla_list)
hla_merged = merge_hla(hla_dfs, lsg_df)
hla_max = merge_hla(hla_dfs, max_rescore_df)

lsg_df.value_counts("Label")
hla_merged.value_counts("Label")

# peptide_list = list(lsg_df.unmod_peptides.drop_duplicates())
# peptide_list_done = list(hla_dfs.peptide.drop_duplicates())
 
# todo_pep_list = [i for i in peptide_list if i not in peptide_list_done]
# len(todo_pep_list)
# with open('/Users/adams/Projects/300K/MSV000091456-SCP/reresults/long_fasta_todo.txt', mode='wt', encoding='utf-8') as myfile:
#     myfile.write('\n'.join(todo_pep_list))

axes[1] = sns.ecdfplot(data=hla_max, x="log_min_binding_score", hue="Label", palette=["#CDEAC0", "#0e1c36"])
axes[2] = sns.ecdfplot(data=hla_merged, x="log_min_binding_score", hue="Label", palette=["#7A8DB3", "#CDEAC0", "#0e1c36"])

# axes[0].set_xscale('log')
# axes[1].set_xscale('log')

for i, (ax, c) in enumerate(zip(axes, 'abcd')):
    ax.annotate(c, xy=(-0.17, 1.1), xycoords='axes fraction',
                fontsize='xx-large', weight='bold')

for ax in axes:
    ax.set_xlabel('Affinity score')
    ax.grid(False, which="both")
    ax.legend(loc='lower right', frameon=False)
    sns.despine(ax=ax)

fig.tight_layout()

plt.savefig('/Users/adams/Projects/300K/MSV000091456-SCP/Figures/4e7-max-vs-rescore.pdf', dpi=300, bbox_inches='tight')
# plt.show()
plt.close()






width = 6.5
height = 6.0
# width = 4
# height = 5
# penguins.species.drop_duplicates()
fig, ax = plt.subplots(figsize=(width, height))

# ax = sns.ecdfplot(data=hla_merged, x="log_min_binding_score", hue="Label", palette=["#7A8DB3", "#CDEAC0", "#0e1c36"])
ax = sns.ecdfplot(data=hla_max, x="log_min_binding_score", hue="Label", palette=["#CDEAC0", "#0e1c36"])
plt.legend(loc="center right", bbox_to_anchor=(1.8, 0.5))#,fancybox=True, ncol=1)
sns.despine()
ax.grid(False, which="both")
plt.xlabel('Affinity score')
plt.axvline(np.log10(2), color="#f33b16") #WB
plt.axvline(np.log10(0.5), color="#f33b16") #SB
# plt.axvline(np.log10(0.036), color="#f33b16") #neo-epitope
# plt.axvline(np.log10(0.007), color="#f33b16") #neo-epitope

plt.savefig("/Users/adams/Projects/300K/MSV000091456-SCP/Figures/4e7-max-vs-rescore.png", dpi=300, bbox_inches="tight")

hla_merged.value_counts("Label")
30058 + 25149
hla_merged[hla_merged["log_min_binding_score"]<np.log10(2)].value_counts("Label")
26001 + 21696
47697/55207
# psm_df["Scan number"] = psm_df["PSMId"].str.extract(r"-(\d+)-")
# msms_df = pd.read_csv(msms_file, sep="\t")
# msms_df["Scan number"] = msms_df["Scan number"].astype(str)
# merged_df = pd.merge(psm_df, msms_df, on="Scan number", how="left")
# filtered_df = merged_df.dropna(subset=['Proteins'])[merged_df["Proteins"].str.startswith('SNP')]
# snp_df = filtered_df[~filtered_df['Proteins'].str.contains(';')]
# gained_shared_lost.append()