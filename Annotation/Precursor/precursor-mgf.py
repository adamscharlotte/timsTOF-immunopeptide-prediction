# ssh vsc20709@login-leibniz.hpc.uantwerpen.be
ssh cadams@10.152.135.57
cd /home/cadams/anaconda3/envs/
conda activate prosit-tims
conda install -c anaconda python=3.8
conda install -c conda-forge python-devtools 
pip install opentims_bruker_bridge
pip install opentimspy

import sqlite3
import pandas as pd
from timspy.df import TimsPyDF
from mgf_filter.masterSpectrum import MasterSpectrum
import numpy as np
from pyteomics import mgf

# import argparse, pathlib
# parser = argparse.ArgumentParser()
# parser.add_argument("folder_path", type=str)
# parser.add_argument("mgf_path", type=str)
# args = parser.parse_args()
# folder_path = args.folder_path
# mgf_path = args.mgf_path
folder_path = "/Users/adams/Projects/300K/2022-library-run/raw-folders/HLAI_1_96_p2-A3_S2-A3_1_6928.d"
mgf_path = "/Users/adams/Projects/300K/2022-library-run/precursor-mgf/HLAI_1_96_p2-A3_S2-A3_1_6928.mgf"
temp_path = "/Users/adams/Projects/300K/2022-library-run/temp-tims.csv"

try:
    import opentims_bruker_bridge
    all_columns = ('frame','scan','tof','intensity','mz','inv_ion_mobility','retention_time')
except ModuleNotFoundError:
    print("Without Bruker proprietary code we cannot yet perform tof-mz and scan-dt transformations.")
    print("Download 'opentims_bruker_bridge' if you are on Linux or Windows.")
    print("Otherwise, you will be able to use only these columns:")
    all_columns = ('frame','scan','tof','intensity','retention_time')

def binning(inp, ignoreCharges):
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignoreCharges)
    ms.export_to_csv(temp_path)
    comb_ms = pd.read_csv(temp_path)
    scan = inp["Precursor"]
    comb_ms["Precursor"] = scan
    comb_ms = comb_ms.drop(columns=["counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"])
    return comb_ms

# Get the pasef information
tdf_path = folder_path + "/analysis.tdf"
con = sqlite3.connect(tdf_path)
cur = con.cursor()
pasef_df = pd.read_sql_query("SELECT * from PasefFrameMsMsInfo", con)
con.close()

# Get all frames
frames = pasef_df.Frame.to_numpy()

# Get the raw data
D = TimsPyDF(folder_path)
df_tims = D.query(frames=frames, columns=all_columns)

# Add precursor ID
df_tims_pasef = df_tims.merge(pasef_df, left_on='frame', right_on='Frame')

# Get the scans per percursor frame
df_frame_pasef = df_tims_pasef[df_tims_pasef.ScanNumBegin <= df_tims_pasef.scan]
df_frame_pasef = df_frame_pasef[df_frame_pasef.scan <= df_frame_pasef.ScanNumEnd]
# df_frame_pasef[["Precursor", "Frame"]].drop_duplicates().len = df_frame_pasef.len

# Write intensities and m/Z in lists
df_intensity = df_frame_pasef.groupby(['Precursor'])['intensity'].apply(list).reset_index(name='combined_INTENSITIES')
df_tof = df_frame_pasef.groupby(['Precursor'])['mz'].apply(list).reset_index(name='combined_MZ')
df_precursor_values = df_frame_pasef[['Precursor', 'CollisionEnergy']].drop_duplicates().merge(df_tof, left_on=['Precursor'], right_on=['Precursor']).merge(df_intensity, left_on=['Precursor'], right_on=['Precursor'])

# Merge all scans for the same precursor
df_bin_result = pd.DataFrame()
for index, line in df_precursor_values.iterrows():
    bin_result = binning(line, True)
    df_bin_result = df_bin_result.append(bin_result)

df_bin_collapsed = df_bin_result.groupby("Precursor").agg(list)
df_precursor_scans = pd.merge(df_precursor_values, df_bin_collapsed, on="Precursor")

# Write to an mgf
spectra_list = []
for i in range (len(df_precursor_scans)):
    mz_list = df_precursor_scans.iloc[i]["MZ"]
    mz_array = np.array(mz_list)
    intensity_list = df_precursor_scans.iloc[i]["INTENSITIES"]
    intensity_array = np.array(intensity_list).astype(int)
    frame = df_precursor_scans.iloc[i]["FRAME"].astype(str)
    precursor = df_precursor_scans.iloc[i]["PRECURSOR"].astype(str)
    title = frame + "_" + precursor
    mass = df_precursor_scans.iloc[i]["MASS"]
    retentiontime = df_precursor_scans.iloc[i]["RETENTION_TIME"]
    # charge = df_precursor_scans.iloc[i][charge_name]
    mgf_params = dict({'title' : title, 'pepmass' : mass, 'rtinseconds' : retentiontime, 'charge' : charge})
    spectra_dict = dict({'m/z array' : mz_array, 'intensity array' : intensity_array, 'params' : mgf_params})
    spectra_list.append(spectra_dict)

mgf.write(spectra = spectra_list, output = mgf_path)


mgf_path