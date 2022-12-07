# /Users/adams/opt/miniconda3/envs/prosit-tims/bin/python3

import sqlite3
import pandas as pd
import numpy as np
from mgf_filter.masterSpectrum import MasterSpectrum
from pyteomics import mgf

# import argparse, pathlib
# parser = argparse.ArgumentParser()
# parser.add_argument("folder_path", type=str)
# parser.add_argument("mgf_path", type=str)
# args = parser.parse_args()
# folder_path = args.folder_path
# mgf_path = args.mgf_path
folder_path = "/Users/adams/Projects/300K/2022-library-run/raw-folders/HLAI_1_96_p2-A3_S2-A3_1_6928.d"
csv_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/precursor-mapped/TUM_HLA_3.csv"
mgf_path = "/Users/adams/Projects/300K/2022-library-run/precursor-mgf/test/TUM_HLA_3.mgf"
temp_path = "/Users/adams/Projects/300K/2022-library-run/temp-tims.csv"

def binning(inp, ignoreCharges):
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignoreCharges)
    ms.export_to_csv(temp_path)
    comb_ms = pd.read_csv(temp_path)
    scan = inp["PRECURSOR"]
    comb_ms["PRECURSOR"] = scan
    comb_ms = comb_ms.drop(columns=["counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"])
    return comb_ms

# Load the pasef
tdf_path = folder_path + "/analysis.tdf"
con = sqlite3.connect(tdf_path)
cur = con.cursor()
pasef_df = pd.read_sql_query("SELECT * from PasefFrameMsMsInfo", con)
con.close()

# Load the csv
df_csv = pd.read_csv(csv_path)
columns = df_csv.columns.str.replace(' ', '_')
df_csv.columns = [x.upper() for x in columns]
df_csv.INTENSITIES = df_csv.INTENSITIES.str.split(';').apply(lambda s: [float(x) for x in s])
df_csv.MZ = df_csv.MZ.str.split(';').apply(lambda s: [float(x) for x in s])

# Merge all scans for the same precursor
df_intensity = df_csv[['PRECURSOR', 'INTENSITIES']].explode('INTENSITIES').groupby(['PRECURSOR'])['INTENSITIES'].apply(list).reset_index(name='combined_INTENSITIES')
df_mz = df_csv[['PRECURSOR', 'MZ']].explode('MZ').groupby(['PRECURSOR'])['MZ'].apply(list).reset_index(name='combined_MZ')

df_ms_values = df_csv[['PRECURSOR']].drop_duplicates().merge(df_mz, left_on=['PRECURSOR'], right_on=['PRECURSOR']).merge(df_intensity, left_on=['PRECURSOR'], right_on=['PRECURSOR'])
df_ms_values["CHARGE"] = np.nan

df_bin_result = pd.DataFrame()
for index, line in df_ms_values.iterrows():
    bin_result = binning(line, True)
    df_bin_result = df_bin_result.append(bin_result)

df_bin_collapsed = df_bin_result.groupby("PRECURSOR").agg(list)
df_precursor_scans = pd.merge(df_bin_collapsed, pasef_df[['Precursor', 'IsolationMz', 'CollisionEnergy']].drop_duplicates(), left_on="PRECURSOR", right_on="Precursor")
df_rt = df_csv[['PRECURSOR', 'RETENTION_TIME']].groupby("PRECURSOR").mean()
df_precursor_scans_rt = pd.merge(df_rt, df_precursor_scans, left_on="PRECURSOR", right_on="Precursor")

# Sort the m/Z
for i in range (len(df_precursor_scans_rt)):
    zipped_list = zip(df_precursor_scans_rt.iloc[i]["mz"], df_precursor_scans_rt.iloc[i]["intensity"])
    sorted_pair = sorted(zipped_list)
    tuples = zip(*sorted_pair)
    list_1, list_2 = [ list(tuple) for tuple in  tuples]
    df_precursor_scans_rt.at[i, "mz"] = list_1
    df_precursor_scans_rt.at[i, "intensity"] = list_2

# Write to an mgf
df_precursor_scans_rt["CHARGE"] = 0
spectra_list = []
for i in range (len(df_precursor_scans_rt)):
    mz_list = df_precursor_scans_rt.iloc[i]["mz"]
    mz_array = np.array(mz_list)
    intensity_list = df_precursor_scans_rt.iloc[i]["intensity"]
    intensity_array = np.array(intensity_list).astype(float)
    # frame = df_precursor_scans.iloc[i]["FRAME"].astype(str)
    precursor = df_precursor_scans_rt.iloc[i]["Precursor"].astype(str)
    # title = frame + "_" + precursor
    title = "test." + precursor + "." + precursor + ".0"
    isolationmz = df_precursor_scans_rt.iloc[i]["IsolationMz"]
    collisionenergy = df_precursor_scans_rt.iloc[i]["CollisionEnergy"]
    # mass = df_precursor_scans.iloc[i]["MASS"]
    retentiontime = df_precursor_scans_rt.iloc[i]["RETENTION_TIME"]
    charge = df_precursor_scans_rt.iloc[i]["CHARGE"]
    # mgf_params = dict({'title' : title, 'pepmass' : mass, 'rtinseconds' : retentiontime, 'charge' : charge})
    mgf_params = dict({'title' : title, 'pepmass' : isolationmz, 'rtinseconds' : retentiontime, 'charge' : charge})
    spectra_dict = dict({'m/z array' : mz_array, 'intensity array' : intensity_array, 'params' : mgf_params})
    spectra_list.append(spectra_dict)

mgf.write(spectra = spectra_list, output = mgf_path)

pasef_df['Precursor'].nunique()
df_precursor_scans_rt['Precursor'].nunique()
