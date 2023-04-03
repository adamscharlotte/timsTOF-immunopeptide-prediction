# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

import os
import pandas as pd
import numpy as np
from pyteomics import mgf
from mgf_filter.masterSpectrum import MasterSpectrum
import time

def binning(inp, ignoreCharges):
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignoreCharges)
    ms.export_to_csv(base_path + "total-scan-consensus/tmp-MasterSpectrum/" + time.strftime("%Y%m%d") + ".csv")
    comb_ms = pd.read_csv(base_path + "total-scan-consensus/tmp-MasterSpectrum/" + time.strftime("%Y%m%d") + ".csv")
    scan = inp["Scannumber"]
    comb_ms["Scannumber"] = scan
    comb_ms = comb_ms.drop(columns=["counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"])
    return comb_ms

# Import extracted_psm
base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
csv_path = base_path + "precursor-mapped/TUM_HLA_133.csv"
mgf_path = base_path + "mapped-summed-mgf/TUM_HLA_133-2.mgf"

csv_df = pd.read_csv(csv_path)
csv_df.columns = csv_df.columns.str.replace(" ", "")
csv_df['combined_MZ'] = csv_df.groupby('Scannumber')['mz'].transform(lambda x: ';'.join(map(str, x)))
csv_df['combined_INTENSITIES'] = csv_df.groupby('Scannumber')['intensities'].transform(lambda x: ';'.join(map(str, x)))
csv_df['rt_combined'] = csv_df.groupby('Scannumber')['retention_time'].transform(lambda x: ';'.join(map(str, x)))

un_annot_df = csv_df[['Precursor', 'Scannumber', 'Rawfile', 'Sequence', 'Score',
    'combined_MZ', 'rt_combined', 'combined_INTENSITIES', 'Charge', 'Mass']]
un_annot_df.combined_INTENSITIES = un_annot_df.combined_INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
un_annot_df.combined_MZ = un_annot_df.combined_MZ.str.split(";").apply(lambda s: [float(x) for x in s])
un_annot_df.rename(columns = {"Charge": "CHARGE"}, inplace=True)

# Sum scans
bin_result_df = pd.DataFrame()
for index, line in un_annot_df.iterrows():
    bin_result = binning(line, True)
    bin_result_df = bin_result_df.append(bin_result)

bin_result_df_collapsed = bin_result_df.groupby("Scannumber").agg(list)
un_annot_df_combined = pd.merge(un_annot_df, bin_result_df_collapsed, on="Scannumber")
un_annot_df_combined.columns
un_annot_df_combined['rt_combined'] = un_annot_df_combined['rt_combined'].str.split(';')
un_annot_df_combined['rt_median'] = un_annot_df_combined['rt_combined'].apply(lambda x: pd.Series(x).median())

# Write mgf
spectra_list = []
for i in range (len(un_annot_df_combined)):
    mz_list = un_annot_df_combined.iloc[i]["mz"]
    mz_array = np.array(mz_list)
    intensity_list = un_annot_df_combined.iloc[i]["intensity"]
    intensity_array = np.array(intensity_list)
    scan = un_annot_df_combined.iloc[i]["Scannumber"].astype(str)
    sequence = un_annot_df_combined.iloc[i]["Sequence"]
    title = scan + "_" + sequence
    mass = un_annot_df_combined.iloc[i]["Mass"]
    retentiontime = un_annot_df_combined.iloc[i]["rt_median"]
    charge = un_annot_df_combined.iloc[i]["CHARGE"]
    mgf_params = dict({'title' : title, 'pepmass' : mass, 'rtinseconds' : retentiontime, 'charge' : charge})
    spectra_dict = dict({'m/z array' : mz_array, 'intensity array' : intensity_array, 'params' : mgf_params})
    spectra_list.append(spectra_dict)

mgf.write(spectra = spectra_list, output = mgf_path)
