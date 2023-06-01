import pandas as pd
import logging
import re


# Set up input and output file names
file = "E_20221201_NO30_400nL_HLAc1_2e7_directIP_titration_rep3_Slot2-4_1_3518"
psm_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/rescore-tims/" + file + ".psms"
msms_file = "/Users/adams/Projects/300K/MSV000091456-SCP/reresults/msms/" + file + ".txt"

# Load the psm file
df = pd.read_csv(psm_file, sep="\t")
psm_df = df[df["q-value"] < 0.01]
psm_df["unmod_peptides"] = psm_df["proteinIds"].str.replace(r"\[UNIMOD:\d+\]", "", regex=True)
psm_df["Scan number"] = psm_df["PSMId"].str.extract(r"-(\d+)-")

psm_file
psm_df.columns
msms_df = pd.read_csv(msms_file, sep="\t")
msms_df["Scan number"] = msms_df["Scan number"].astype(str)
merged_df = pd.merge(psm_df, msms_df, on="Scan number", how="left")
filtered_df = merged_df.dropna(subset=['Proteins'])[merged_df["Proteins"].str.startswith('SNP')]
snp_df = filtered_df[~filtered_df['Proteins'].str.contains(';')]
filtered_df[["Proteins", "Sequence"]]

snp_count = len(snp_df)
logging.info(f"There are {snp_count} missense SNPs identified.")
print(f"There are {snp_count} missense SNPs identified.")

pd.options.display.max_colwidth = 100
snp_df[["q-value", "score", "PSMId", "Proteins"]]


# # df[df["PSMId"] == "IPX72_HLAI_ElA_S3-H2_1_11987-14452-RSLTRHLKY-3--1"]

# psm_max_file = "/Users/adams/Projects/300K/Jurkat-A549/reresults/rescore-tims/IPX67_HLAI_A_S1-A2_1_10830-andromeda.psms"
# max_df = pd.read_csv(psm_max_file, sep="\t")
# max_df.columns

# df_12168 = max_df[max_df["PSMId"].str.contains("12168")]

# df_12168["q-value"]