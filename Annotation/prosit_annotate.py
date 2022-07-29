# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

# import logging
# from tkinter.font import names
# from matplotlib.colors import cnames

import pandas as pd
import numpy as np
# from psims import INTENSITY_ARRAY

from fundamentals import constants
from fundamentals.fragments import initialize_peaks
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import h5py
# import fundamentals.constants as C

un_annot_df = pd.read_csv("/Users/adams/Projects/300K/TUM_aspn_10.csv")

un_annot_df.INTENSITIES = un_annot_df.INTENSITIES.str.split(';').apply(lambda s: [float(x) for x in s])
un_annot_df.MZ = un_annot_df.MZ.str.split(';').apply(lambda s: [float(x) for x in s])

list(un_annot_df.columns)

un_annot_df.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)
un_annot_df["REVERSE"].fillna(False, inplace=True)
un_annot_df["REVERSE"].replace("+", True, inplace=True)
un_annot_df["MODIFIED_SEQUENCE"] = maxquant_to_internal(un_annot_df["MODIFIED_SEQUENCE"].to_numpy())
un_annot_df["SEQUENCE"] = internal_without_mods(un_annot_df["MODIFIED_SEQUENCE"])
un_annot_df['PEPTIDE_LENGTH'] = un_annot_df["SEQUENCE"].apply(lambda x: len(x))

annot_df = annotate_spectra(un_annot_df)

full_df = pd.concat([un_annot_df, annot_df], axis=1)