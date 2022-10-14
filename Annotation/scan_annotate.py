# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

from operator import concat
import pandas as pd
import numpy as np
import os

from fundamentals import constants
from fundamentals.fragments import initialize_peaks
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import h5py
import argparse, pathlib
from fundamentals.mod_string import parse_modstrings, maxquant_to_internal
from fundamentals.constants import ALPHABET
from mgf_filter.util import timeStamped
from mgf_filter.masterSpectrum import MasterSpectrum

parser = argparse.ArgumentParser()
parser.add_argument("pool", type=str)					# Filename
args = parser.parse_args()

# pool = "TUM_HLA2_88"
pool = args.pool

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
sum_path = base_path + "scan-consensus/summed-20-ppm/" + pool + ".csv"
annot_path = base_path + "scan-consensus/annotated-20-ppm/" + pool + ".csv"

un_annot_df_combined = pd.read_csv(sum_path)

un_annot_df_combined.MZ = un_annot_df_combined.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
un_annot_df_combined.INTENSITIES = un_annot_df_combined.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])

for i in range (len(un_annot_df_combined)):
    zipped_list = zip(un_annot_df_combined.iloc[i]["MZ"], un_annot_df_combined.iloc[i]["INTENSITIES"])
    sorted_pair = sorted(zipped_list)
    tuples = zip(*sorted_pair)
    list_1, list_2 = [ list(tuple) for tuple in  tuples]
    un_annot_df_combined.at[i, "MZ"] = list_1
    un_annot_df_combined.at[i, "INTENSITIES"] = list_2

annot_df = annotate_spectra(un_annot_df_combined)
full_df = pd.concat([un_annot_df_combined.drop(columns = ["INTENSITIES", "MZ"]), annot_df], axis=1)

full_df["MZ"] = [';'.join(map(str, l)) for l in full_df['MZ']]
full_df["INTENSITIES"] = [';'.join(map(str, l)) for l in full_df['INTENSITIES']]
full_df.to_csv(annot_path)

# full_df[full_df["PRECURSOR"] == 8200][["INTENSITIES", "MZ"]].values
