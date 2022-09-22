# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

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

# pool = "TUM_HLA_16"

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
sum_path = base_path + "precursor-consensus/summed/" + args.pool + ".csv"
annot_path = base_path + "precursor-consensus/annotated/" + args.pool + ".csv"

# sum_path = base_path + "precursor-consensus/summed/" + pool + ".csv"
# annot_path = base_path + "precursor-consensus/annotated/" + pool + ".csv"

un_annot_df_combined = pd.read_csv(sum_path)

un_annot_df_combined.MZ = un_annot_df_combined.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
un_annot_df_combined.INTENSITIES = un_annot_df_combined.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])

annot_df = annotate_spectra(un_annot_df_combined)
full_df = pd.concat([un_annot_df_combined.drop(columns = ["INTENSITIES", "MZ"]), annot_df], axis=1)

full_df["MZ"] = [';'.join(map(str, l)) for l in full_df['MZ']]
full_df["INTENSITIES"] = [';'.join(map(str, l)) for l in full_df['INTENSITIES']]
full_df.to_csv(annot_path)