# ssh cadams@10.152.135.57
# /home/cadams/anaconda3/envs/prosit-annotate/bin/python3

import pandas as pd
import numpy as np
import os

from fundamentals import constants
from fundamentals.fragments import initialize_peaks
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import argparse, pathlib
from fundamentals.mod_string import parse_modstrings, maxquant_to_internal
from fundamentals.constants import ALPHABET
from mgf_filter.util import timeStamped
from mgf_filter.masterSpectrum import MasterSpectrum
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("csv_path", type=str)
args = parser.parse_args()
csv_path = args.csv_path

# csv_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/test-d/211113_SS_malignant_HNSCC_Tue39L243_20%_DDA_Rep1.csv" # nolint
base_path = os.path.dirname(csv_path)
csv_df = pd.read_csv(csv_path)

csv_df.combined_INTENSITIES = csv_df.combined_INTENSITIES.str.split(";")
csv_df.combined_INTENSITIES = csv_df.combined_INTENSITIES.apply(lambda s: [float(x) for x in s])
csv_df.combined_MZ = csv_df.combined_MZ.str.split(";")
csv_df.combined_MZ = csv_df.combined_MZ.apply(lambda s: [float(x) for x in s])

def binning(inp, ignoreCharges):
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignoreCharges)
    ms.export_to_csv("/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/test-d/temp.txt")
    comb_ms = pd.read_csv("/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/test-d/temp.txt")
    scan = inp["SCAN_NUMBER"]
    comb_ms["SCAN_NUMBER"] = scan
    comb_ms = comb_ms.drop(columns=["counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"])
    return comb_ms

bin_result_df = pd.DataFrame()
for index, line in csv_df.iterrows():
    bin_result = binning(line, True)
    bin_result_df = bin_result_df.append(bin_result)

bin_result_df_collapsed = bin_result_df.groupby("SCAN_NUMBER").agg(list)
csv_df_combined = pd.merge(csv_df, bin_result_df_collapsed, on="SCAN_NUMBER")

csv_df_comb = csv_df_combined.drop(columns = ["combined_INTENSITIES", "combined_MZ"]).rename(columns={"intensity":"INTENSITIES", "mz":"MZ"})
type(csv_df_comb["MZ"].iloc[0])
csv_df_comb['INTENSITIES'] = csv_df_comb['INTENSITIES'].apply(lambda x: np.array(x))
csv_df_comb['MZ'] = csv_df_comb['MZ'].apply(lambda x: np.array(x))

for i in range (len(csv_df_comb)):
    zipped_list = zip(csv_df_comb.iloc[i]["MZ"], csv_df_comb.iloc[i]["INTENSITIES"])
    sorted_pair = sorted(zipped_list)
    tuples = zip(*sorted_pair)
    list_1, list_2 = [np.array(tuple) for tuple in tuples]
    csv_df_comb.at[i, "MZ"] = list_1
    csv_df_comb.at[i, "INTENSITIES"] = list_2

csv_df_comb["MZ_RANGE"] = "0"

file_name = csv_df["RAW_FILE"][0]

csv_df_comb.to_pickle(base_path + "/" + file_name + ".pkl")