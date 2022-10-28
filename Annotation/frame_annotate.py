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

parser = argparse.ArgumentParser()
parser.add_argument('pool', type=str)					# Filename
args = parser.parse_args()

# pool = "TUM_HLA2_7"
pool = args.pool

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
un_annot_path = base_path + "full-truncated-qc/un-annotated/" + pool + ".csv"
un_annot_df = pd.read_csv(un_annot_path)
annot_path = base_path + "full-truncated-qc/annotated-40-ppm/" + pool + ".csv"

un_annot_df.INTENSITIES = un_annot_df.INTENSITIES.str.split(';').apply(lambda s: [float(x) for x in s])
un_annot_df.MZ = un_annot_df.MZ.str.split(';').apply(lambda s: [float(x) for x in s])

un_annot_df.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)
un_annot_df["REVERSE"].fillna(False, inplace=True)
un_annot_df["REVERSE"].replace("+", True, inplace=True)
un_annot_df["MODIFIED_SEQUENCE"] = maxquant_to_internal(un_annot_df["MODIFIED_SEQUENCE"].to_numpy())
un_annot_df["SEQUENCE"] = internal_without_mods(un_annot_df["MODIFIED_SEQUENCE"])
un_annot_df['PEPTIDE_LENGTH'] = un_annot_df["SEQUENCE"].apply(lambda x: len(x))

generator_sequence_numeric = parse_modstrings(list(un_annot_df['MODIFIED_SEQUENCE'].values), ALPHABET, translate=True)
enum_gen_seq_num = enumerate(generator_sequence_numeric)
array = np.zeros((len(list(un_annot_df['MODIFIED_SEQUENCE'].values)),30), dtype=np.uint8)
for i, sequence_numeric in enum_gen_seq_num:
    if len(sequence_numeric) > 30:
        if filter:
            pass # don't overwrite 0 in the array that is how we can differentiate
        else:
            raise Exception(f"The Sequence {sequence_numeric}, has {len(sequence_numeric)} Amino Acids."
                        f"The maximum number of amino acids allowed is {C.SEQ_LEN}")
    else:
        array[i, 0:len(sequence_numeric)] = sequence_numeric

un_annot_df['SEQUENCE_INT'] = array.tolist()

annot_df = annotate_spectra(un_annot_df)

full_df = pd.concat([un_annot_df.drop(columns = ["INTENSITIES", "MZ"]), annot_df], axis=1)
full_df["MZ"] = [';'.join(map(str, l)) for l in full_df['MZ']]
full_df["INTENSITIES"] = [';'.join(map(str, l)) for l in full_df['INTENSITIES']]
full_df.to_csv(annot_path)
