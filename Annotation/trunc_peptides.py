# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

import pandas as pd
import numpy as np
import os

# from fundamentals import constants
# from fundamentals.fragments import initialize_peaks
# from fundamentals.annotation.annotation import annotate_spectra
# from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
# import argparse, pathlib
# from fundamentals.mod_string import parse_modstrings, maxquant_to_internal
# from fundamentals.constants import ALPHABET

parser = argparse.ArgumentParser()
parser.add_argument("pool", type=str)					# Filename
args = parser.parse_args()
pool = args.pool
# pool = "TUM_lysn_33"

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
calibrated_path = base_path + "scan-consensus/calibrated-40-ppm/" + pool + ".csv"
file_path = base_path + "scan-consensus/calibrated-mapped/" + pool + ".csv"
map_path = base_path + "full-length-map/" + pool + ".csv"

map_df = pd.read_csv(map_path)

calibrated_annot_df = pd.read_csv(calibrated_path)
calibrated_annot_df = calibrated_annot_df.rename({"SEQUENCE": "OBS_SEQUENCE", "MODIFIED_SEQUENCE": "OBS_MODIFIED_SEQUENCE", "SEQUENCE_INT": "OBS_SEQUENCE_INT"}, axis='columns')
calibrated_annot_df_merged = pd.merge(left=calibrated_annot_df, right=map_df, how='left', left_on='SCAN_NUMBER', right_on='SCAN_NUMBER')

calibrated_annot_df_merged.to_csv(file_path)