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

def binning(inp, ignoreCharges):
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignoreCharges)
    ms.export_to_csv(base_path + "/temp_3.txt")
    comb_ms = pd.read_csv(base_path + "/temp_3.txt")
    scan = inp["SCAN_NUMBER"]
    comb_ms["SCAN_NUMBER"] = scan
    comb_ms = comb_ms.drop(columns=["counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"])
    return comb_ms


parser = argparse.ArgumentParser()
parser.add_argument("csv_path", type=str)
args = parser.parse_args()
csv_path = args.csv_path

# csv_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/Jurkat-A549/reresults/IPX_A549_A_S1-A2_1_11144/IPX_A549_A_S1-A2_1_11144.csv" # nolint
# csv_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/d/220404_NHG_benign_UDN10_Liver_Tue39L243_17%_DDA_Rep2/220404_NHG_benign_UDN10_Liver_Tue39L243_17%_DDA_Rep2.csv" # nolint
# csv_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/A375/reresults/A375-low-input-HLAI/E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496/E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496.csv" # nolint

base_path = os.path.dirname(csv_path)

total_list = []
chunksize = 1000
with pd.read_csv(csv_path, chunksize=chunksize) as reader:
    for i, chunk in enumerate(reader):
        # if i == 2:
        #     break
        print("Starting chunk" + str(i))
        chunk.combined_INTENSITIES = chunk.combined_INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
        chunk.combined_MZ = chunk.combined_MZ.str.split(";").apply(lambda s: [float(x) for x in s])
        bin_result_list = []
        for index, line in chunk.iterrows():
            bin_result = binning(line, True)
            bin_result_list.append(bin_result)
        bin_result_df = pd.concat(bin_result_list)
        bin_result_df_collapsed = bin_result_df.groupby("SCAN_NUMBER").agg(list)
        chunk_combined = pd.merge(chunk, bin_result_df_collapsed, on="SCAN_NUMBER")
        chunk_comb = chunk_combined.drop(columns = ["combined_INTENSITIES", "combined_MZ"]).rename(columns={"intensity":"INTENSITIES", "mz":"MZ", "median_CE":"COLLISION_ENERGY"})
        chunk_comb['INTENSITIES'] = chunk_comb['INTENSITIES'].apply(lambda x: np.array(x))
        chunk_comb['MZ'] = chunk_comb['MZ'].apply(lambda x: np.array(x))
        for i in range (len(chunk_comb)):
            zipped_list = zip(chunk_comb.iloc[i]["MZ"], chunk_comb.iloc[i]["INTENSITIES"])
            sorted_pair = sorted(zipped_list)
            tuples = zip(*sorted_pair)
            list_1, list_2 = [np.array(tuple) for tuple in tuples]
            chunk_comb.at[i, "MZ"] = list_1
            chunk_comb.at[i, "INTENSITIES"] = list_2
        total_list.append(chunk_comb)

total_df = pd.concat(total_list, ignore_index=True)
total_df["MZ_RANGE"] = "0"
file_name = total_df["RAW_FILE"][0]

total_df.to_pickle(base_path + "/" + file_name + ".pkl")

# file_name = "E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep3_Slot2-5_1_3521"
# total_df = pd.read_pickle(base_path + "/" + file_name + ".pkl")

# total_df
# df_raw_drop = total_df[["RAW_FILE", "SCAN_NUMBER", "INTENSITIES", "MZ", "MZ_RANGE","COLLISION_ENERGY"]]

total_df.columns