# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

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

parser = argparse.ArgumentParser()
parser.add_argument("pool", type=str)					# Filename
args = parser.parse_args()

# pool = "TUM_first_pool_5"
pool = args.pool

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
un_annot_path = base_path + "total-scan-consensus/filtered/" + pool + ".csv"
sum_path = base_path + "total-scan-consensus/summed-40-ppm/" + pool + ".csv"

un_annot_df = pd.read_csv(un_annot_path)
un_annot_df["combined_INTENSITIES"]

un_annot_df.combined_INTENSITIES = un_annot_df.combined_INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
un_annot_df.combined_MZ = un_annot_df.combined_MZ.str.split(";").apply(lambda s: [float(x) for x in s])

def binning(inp, ignoreCharges):
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignoreCharges)
    ms.export_to_csv(base_path + "scan-consensus/tmp-MasterSpectrum/24102022.csv")
    comb_ms = pd.read_csv(base_path + "scan-consensus/tmp-MasterSpectrum/24102022.csv")
    scan = inp["SCAN_NUMBER"]
    comb_ms["SCAN_NUMBER"] = scan
    comb_ms = comb_ms.drop(columns=["counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"])
    return comb_ms

bin_result_df = pd.DataFrame()
for index, line in un_annot_df.iterrows():
    bin_result = binning(line, True)
    bin_result_df = bin_result_df.append(bin_result)

bin_result_df_collapsed = bin_result_df.groupby("SCAN_NUMBER").agg(list)
un_annot_df_combined = pd.merge(un_annot_df, bin_result_df_collapsed, on="SCAN_NUMBER")

un_annot_df_combined.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)
un_annot_df_combined["REVERSE"].fillna(False, inplace=True)
un_annot_df_combined["REVERSE"].replace("+", True, inplace=True)
un_annot_df_combined["MODIFIED_SEQUENCE"] = maxquant_to_internal(un_annot_df_combined["MODIFIED_SEQUENCE"].to_numpy())
un_annot_df_combined["SEQUENCE"] = internal_without_mods(un_annot_df_combined["MODIFIED_SEQUENCE"])
un_annot_df_combined["PEPTIDE_LENGTH"] = un_annot_df_combined["SEQUENCE"].apply(lambda x: len(x))
un_annot_df_combined = un_annot_df_combined.drop(columns=["combined_INTENSITIES", "combined_MZ"])
un_annot_df_combined.rename(columns = {"mz": "MZ"}, inplace=True)
un_annot_df_combined.rename(columns = {"intensity": "INTENSITIES"}, inplace=True)

generator_sequence_numeric = parse_modstrings(list(un_annot_df_combined["MODIFIED_SEQUENCE"].values), ALPHABET, translate=True)
enum_gen_seq_num = enumerate(generator_sequence_numeric)
array = np.zeros((len(list(un_annot_df_combined["MODIFIED_SEQUENCE"].values)),30), dtype=np.uint8)
for i, sequence_numeric in enum_gen_seq_num:
    if len(sequence_numeric) > 30:
        if filter:
            pass # don"t overwrite 0 in the array that is how we can differentiate
        else:
            raise Exception(f"The Sequence {sequence_numeric}, has {len(sequence_numeric)} Amino Acids."
                        f"The maximum number of amino acids allowed is {C.SEQ_LEN}")
    else:
        array[i, 0:len(sequence_numeric)] = sequence_numeric

un_annot_df_combined["SEQUENCE_INT"] = array.tolist()

# annot_df = annotate_spectra(un_annot_df_combined)
# full_df = pd.concat([un_annot_df_combined.drop(columns = ["INTENSITIES", "MZ"]), annot_df], axis=1)

un_annot_df_combined["MZ"] = [';'.join(map(str, l)) for l in un_annot_df_combined['MZ']]
un_annot_df_combined["INTENSITIES"] = [';'.join(map(str, l)) for l in un_annot_df_combined['INTENSITIES']]

un_annot_df_combined.to_csv(sum_path)

print("Are spectra lost?", len(un_annot_df.index) > len(un_annot_df_combined.index))

# full_df["MZ"] = [';'.join(map(str, l)) for l in full_df['MZ']]
# full_df["INTENSITIES"] = [';'.join(map(str, l)) for l in full_df['INTENSITIES']]
# full_df.to_csv(annot_path)