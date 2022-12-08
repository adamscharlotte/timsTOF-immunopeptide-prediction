# ssh cadams@10.152.135.57
/home/cadams/anaconda3/envs/prosit-tims/bin/python3
# /Users/adams/opt/miniconda3/envs/prosit-tims/bin/python3
# screen -r tims_tof
# conda activate prosit-tims
# pip install python-dev-tools
# pip install opentims_bruker_bridge
# pip install opentimspy
# pip uninstall numpy
# pip install numpy==1.19.3
# pip install git+https://github.com/michalsta/opentims

import pandas as pd
import numpy as np
import os

from timspy.df import TimsPyDF
from mgf_filter.masterSpectrum import MasterSpectrum
from fundamentals import constants
from fundamentals.fragments import initialize_peaks
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import argparse, pathlib
from fundamentals.mod_string import parse_modstrings, maxquant_to_internal
from fundamentals.constants import ALPHABET

try:
    import opentims_bruker_bridge
    all_columns = ('frame','scan','tof','intensity','mz','inv_ion_mobility','retention_time')
except ModuleNotFoundError:
    print("Without Bruker proprietary code we cannot yet perform tof-mz and scan-dt transformations.")
    print("Download 'opentims_bruker_bridge' if you are on Linux or Windows.")
    print("Otherwise, you will be able to use only these columns:")
    all_columns = ('frame','scan','tof','intensity','retention_time')

# --------------------------------------- IMPORT DATA ---------------------------------------

file_name = "HLAI_1_96_p2-A3_S2-A3_1_6928"
# path to MaxQuant txt output folder
maxquant_folder = "/Users/adams/Projects/300K/2022-library-run/MaxQuant/TUM_HLA_3-txt/"
SCAN_NUMBERS = [727, 876, 3810]
# base_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/"
base_path = "/Users/adams/Projects/300K/2022-library-run/"
# path to the bruker .d folder
d_folder_path = base_path + "raw-folders/" + file_name + ".d"
# output path
csv_path = base_path + "Annotation/annotated-tims/" + file_name + ".csv"
temp_path = base_path + "Annotation/temp-tims.csv"

# Get the wanted precursors
path_precursors = maxquant_folder + "accumulatedMsmsScans.txt"
df_precursors = pd.read_csv(path_precursors, sep = "\t")
df_precursors.columns = df_precursors.columns.str.replace(" ", "")

df_precursors_filtered = df_precursors[df_precursors.Scannumber.isin(SCAN_NUMBERS)]
df_precursors_filtered.PASEFprecursorIDs = df_precursors_filtered.PASEFprecursorIDs.str.split(";").apply(lambda s: [int(x) for x in s])
df_precursors_filtered = df_precursors_filtered.explode("PASEFprecursorIDs")
df_precursors_filtered.columns
df_precursors_filtered[["Scannumber", "PASEFprecursorIDs"]]
# Get the frames
path_pasef = maxquant_folder + "pasefMsmsScans.txt"
df_pasef = pd.read_csv(path_pasef, sep = "\t")
df_pasef.columns = df_pasef.columns.str.replace(" ", "")

df_pasef_filtered = df_pasef[df_pasef.Precursor.isin(df_precursors_filtered["PASEFprecursorIDs"])]
frames = df_pasef_filtered.Frame.to_numpy()
df_pasef_filtered.columns
df_pasef_filtered[["Precursor", "Frame", "ScanNumBegin", "ScanNumEnd"]]

# Load the timsTOF folder
D = TimsPyDF(d_folder_path)
print(D)
print(len(D))
print(all_columns)

df_tims = D.query(frames=frames, columns=all_columns)
df_raw_pasef_all_scans = df_tims.merge(df_pasef_filtered, left_on='frame', right_on='Frame')
df_raw_pasef = df_raw_pasef_all_scans[df_raw_pasef_all_scans.ScanNumBegin <= df_raw_pasef_all_scans.scan]
df_raw_pasef = df_raw_pasef[df_raw_pasef.scan <= df_raw_pasef.ScanNumEnd]
df_raw_pasef.columns
df_raw_pasef[["Precursor", "Frame"]].drop_duplicates()
df_raw_pasef_mapped = df_raw_pasef.drop(columns=["scan"]).merge(df_precursors_filtered, left_on='Precursor', right_on='PASEFprecursorIDs')
df_raw_pasef_mapped.columns
df_intensity = df_raw_pasef_mapped.groupby(['Scannumber'])['intensity'].apply(list).reset_index(name='combined_INTENSITIES')
df_tof = df_raw_pasef_mapped.groupby(['Scannumber'])['tof'].apply(list).reset_index(name='combined_MZ')
df_raw_pasef_mapped.rename(columns = {"Charge": "CHARGE"}, inplace=True)

df_ms_values = df_raw_pasef_mapped[['Scannumber', 'Sequence', 'Modifiedsequence', 'CHARGE']].drop_duplicates().merge(df_tof, left_on=['Scannumber'], right_on=['Scannumber']).merge(df_intensity, left_on=['Scannumber'], right_on=['Scannumber'])
df_ms_values.columns
# .drop_duplicates()

# -------------------------------------- COMBINE SCANS --------------------------------------

def binning(inp, ignoreCharges):
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignoreCharges)
    ms.export_to_csv(temp_path)
    comb_ms = pd.read_csv(temp_path)
    scan = inp["Scannumber"]
    comb_ms["Scannumber"] = scan
    comb_ms = comb_ms.drop(columns=["counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"])
    return comb_ms

bin_result_df = pd.DataFrame()
for index, line in df_ms_values.iterrows():
    bin_result = binning(line, True)
    bin_result_df = bin_result_df.append(bin_result)

bin_result_df_collapsed = bin_result_df.groupby("Scannumber").agg(list)
un_annot_df_combined = pd.merge(df_ms_values, bin_result_df_collapsed, on="Scannumber")

un_annot_df_combined["Modifiedsequence"] = maxquant_to_internal(un_annot_df_combined["Modifiedsequence"].to_numpy())
generator_sequence_numeric = parse_modstrings(list(un_annot_df_combined["Modifiedsequence"].values), ALPHABET, translate=True)
enum_gen_seq_num = enumerate(generator_sequence_numeric)
array = np.zeros((len(list(un_annot_df_combined["Modifiedsequence"].values)),30), dtype=np.uint8)
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


def int_to_onehot(charge):
    precursor_charge = np.full((6), False)
    precursor_charge[charge-1] = True
    return precursor_charge

un_annot_df_combined["PRECURSOR_CHARGE_ONEHOT"] = un_annot_df_combined.apply(lambda row: int_to_onehot(row.CHARGE), axis = 1)

for i in range (len(un_annot_df_combined)):
    zipped_list = zip(un_annot_df_combined.iloc[i]["combined_MZ"], un_annot_df_combined.iloc[i]["combined_INTENSITIES"])
    sorted_pair = sorted(zipped_list)
    tuples = zip(*sorted_pair)
    list_1, list_2 = [ list(tuple) for tuple in  tuples]
    un_annot_df_combined.at[i, "combined_MZ"] = list_1
    un_annot_df_combined.at[i, "combined_INTENSITIES"] = list_2

un_annot_df_combined
un_annot_df_combined.rename(columns = {"Modifiedsequence": "MODIFIED_SEQUENCE"}, inplace=True)
un_annot_df_combined.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)
un_annot_df_combined.rename(columns = {"combined_MZ": "MZ"}, inplace=True)
un_annot_df_combined.rename(columns = {"combined_INTENSITIES": "INTENSITIES"}, inplace=True)
un_annot_df_combined["MASS_ANALYZER"] = "TOF"
annot_df = annotate_spectra(un_annot_df_combined)
