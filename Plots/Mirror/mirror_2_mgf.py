# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import os
import logging

import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from spectrum_utils.spectrum import MsmsSpectrum
from typing import Iterator
from typing import List
from pyteomics import mgf
import pandas as pd
import numpy as np
import itertools
import seaborn as sns

from fundamentals import constants
from fundamentals.fragments import initialize_peaks
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.mod_string import maxquant_to_internal

def get_spectral_angle(intensities):
    pred= np.array(intensities[0])
    true = np.array(intensities[1])
    epsilon = 1e-7
    list_1 = np.argwhere(true>0)
    list_2 = np.argwhere(pred>0)
    indices = np.union1d(list_1,list_2)
    pred_masked = pred[indices]
    true_masked = true[indices]
    true_masked += epsilon
    pred_masked += epsilon
    true_norm = true_masked*(1/np.sqrt(np.sum(np.square(true_masked), axis=0)))
    pred_norm = pred_masked*(1/np.sqrt(np.sum(np.square(pred_masked), axis=0)))
    product = np.sum(true_norm*pred_norm, axis=0)
    arccos = np.arccos(product)
    return 1-2*arccos/np.pi

def read_mgf_df(filename: str):
    identifier = []
    precursor_mz = []
    retention_time = []
    charge = []
    mz = []
    intensity = []
    for i, mgf_spectrum in enumerate(mgf.read(filename)):
        identifier.append(mgf_spectrum['params']['title'])
        precursor_mz.append(float(mgf_spectrum['params']['pepmass'][0]))
        retention_time.append(float(mgf_spectrum['params']['rtinseconds']))
        charge.append(int(mgf_spectrum['params']['charge'][0]))
        mz.append(mgf_spectrum['m/z array'])
        intensity.append(mgf_spectrum['intensity array'])
    spectrum_df = pd.DataFrame({"identifier":identifier, "precursor_mz":precursor_mz,
        "retention_time":retention_time, "charge":charge, "mz":mz, "intensity":intensity})
    return spectrum_df

def verify_extension(supported_extensions: List[str], filename: str) -> None:
    """
    Check that the given file name has a supported extension.

    Parameters
    ----------
    supported_extensions : List[str]
        A list of supported file extensions.
    filename : str
        The file name to be checked.

    Raises
    ------
    FileNotFoundError
        If the file name does not have one of the supported extensions.
    """
    _, ext = os.path.splitext(os.path.basename(filename))
    if ext.lower() not in supported_extensions:
        logging.error('Unrecognized file format: %s', filename)
        raise FileNotFoundError(f'Unrecognized file format (supported file '
                                f'formats: {", ".join(supported_extensions)})')
    elif not os.path.isfile(filename):
        logging.error('File not found: %s', filename)
        raise FileNotFoundError(f'File {filename} does not exist')

def read_mgf(filename: str) -> Iterator[MsmsSpectrum]:
    """
    Read all spectra from the given mgf file.
    Parameters
    ----------
    filename: str
        The mgf file name from which to read the spectra.
    Returns
    -------
    Iterator[Spectrum]
        An iterator of spectra in the given mgf file.
    """
    # Test if the given file is an mgf file.
    verify_extension(['.mgf'], filename)
    # Get all query spectra.
    for i, mgf_spectrum in enumerate(mgf.read(filename)):
        # Create spectrum.
        identifier = mgf_spectrum['params']['title']
        precursor_mz = float(mgf_spectrum['params']['pepmass'][0])
        retention_time = float(mgf_spectrum['params']['rtinseconds'])
        if 'charge' in mgf_spectrum['params']:
            precursor_charge = int(mgf_spectrum['params']['charge'][0])
        else:
            precursor_charge = None
        spectrum = MsmsSpectrum(identifier, precursor_mz, precursor_charge,
                                mgf_spectrum['m/z array'],
                                mgf_spectrum['intensity array'],
                                retention_time=retention_time)
        spectrum.index = i
        spectrum.is_processed = False
        yield spectrum

def order_peaks(df):
    for i in range (len(df)):
        zipped_list = zip(df.iloc[i]["mz"], df.iloc[i]["intensity"])
        sorted_pair = sorted(zipped_list)
        tuples = zip(*sorted_pair)
        list_1, list_2 = [ list(tuple) for tuple in  tuples]
        df.at[i, "mz"] = list_1
        df.at[i, "intensity"] = list_2

def merge_df(df_mgf, msms):
    merged_df = pd.merge(df_mgf.drop('charge', axis=1), msms.drop('Intensities', axis=1), left_on='id', right_on='Scan number', how='inner')
    new_columns = {col: col.replace(' ', '_').upper() for col in merged_df.columns}
    merged_df.rename(columns=new_columns, inplace=True)
    merged_df.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)
    merged_df.rename(columns = {"INTENSITY": "INTENSITIES"}, inplace=True)
    merged_df['PRECURSOR_CHARGE'] = merged_df['PRECURSOR_CHARGE'].astype(int)
    merged_df["MODIFIED_SEQUENCE"] = maxquant_to_internal(merged_df["MODIFIED_SEQUENCE"].to_numpy())
    return merged_df

def concat_and_filter(merged_df, annot_df):
    full_df = pd.concat([merged_df.drop(columns = ["INTENSITIES", "MZ"]), annot_df], axis=1)
    filtered_df = full_df[full_df["SEQUENCE"] == peptide][full_df["PRECURSOR_CHARGE"] == 1]
    return filtered_df

# -----------------------------------------------------------------------------
charge = "1"
peptide = "GVDAANSAAQQY"

# mgf1_filename="/Users/adams/Downloads/02445d_BB3-TUM_HLA_15_01_01-3xHCD-1h-R4.mgf"
mgf1_filename="/Users/adams/Downloads/02446d_GD1-TUM_HLA_133_01_01-3xHCD-1h-R4.mgf"
mgf2_filename="/Users/adams/Projects/300K/2022-library-run/Annotation/mapped-summed-mgf/TUM_HLA_133-3.mgf"
# mgf2_filename="/Users/adams/Projects/300K/PXD038782-comparison/mgf/HLA_133_orbi_pred.mgf"
# msms1_path="/Users/adams/Downloads/TUM_HLA_15_01_01_3xHCD-1h-R4-unspecific/msms.txt"
msms1_path="/Users/adams/Downloads/TUM_HLA_133_01_01_3xHCD-1h-R4-unspecific/msms.txt"
msms2_path="/Users/adams/Projects/300K/2022-library-run/msms-txt/TUM_HLA_133.txt"
# mgf2_id = "3035"
# mgf2_id = "14565"
mgf2_id = "3378"

# -----------------------------------------------------------------------------
df_mgf1 = read_mgf_df(mgf1_filename)
df_mgf2 = read_mgf_df(mgf2_filename)

order_peaks(df_mgf1)
order_peaks(df_mgf2)

df_mgf1['id'] = df_mgf1['identifier'].str[41:].astype(int)
df_mgf2['id'] = df_mgf2['identifier'].astype(int)

msms1 = pd.read_csv(msms1_path, sep='\t')
msms2 = pd.read_csv(msms2_path, sep='\t')

merged_df1 = merge_df(df_mgf1, msms1)
merged_df2 = merge_df(df_mgf2, msms2)
# merged_df2 = merge_df(df_mgf2, msms1)

annot_df1 = annotate_spectra(merged_df1)
annot_df2 = annotate_spectra(merged_df2)

filtered_df1 = concat_and_filter(merged_df1, annot_df1)
filtered_df2 = concat_and_filter(merged_df2, annot_df2)

filtered_df1.rename(columns = {"INTENSITIES": "INTENSITIES_1", "MZ":"MZ_1"}, inplace=True)
filtered_df2.rename(columns = {"INTENSITIES": "INTENSITIES_2", "MZ":"MZ_2"}, inplace=True)

result = pd.merge(filtered_df1.assign(dummy=1), filtered_df2.assign(dummy=1), on='dummy', how='outer').drop('dummy', axis=1)
result["SPECTRAL_ANGLE"] = result[['INTENSITIES_1','INTENSITIES_2']].apply(lambda x : get_spectral_angle(x), axis=1)
result["SPECTRAL_ANGLE"].fillna(0, inplace=True)

# -----------------------------------------------------------------------------
filtered_df = filtered_df1

# Generate all possible combinations of elements
combinations = list(itertools.product(filtered_df['IDENTIFIER'], filtered_df['IDENTIFIER']))

# Create a new dataframe that contains these combinations and their corresponding values
combination_values = [(combination[0], combination[1], filtered_df.loc[filtered_df['IDENTIFIER']==combination[0], 'INTENSITIES'].iloc[0], filtered_df.loc[filtered_df['IDENTIFIER']==combination[1], 'INTENSITIES'].iloc[0]) for combination in combinations]
df_combinations = pd.DataFrame(combination_values, columns=['A', 'B', 'valueA', 'valueB'])
df_combinations_filter = df_combinations[df_combinations["A"]!=df_combinations["B"]]

df_combinations_filter["SPECTRAL_ANGLE"] = df_combinations_filter[['valueA','valueB']].apply(lambda x : get_spectral_angle(x), axis=1)
df_combinations_filter["SPECTRAL_ANGLE"].fillna(0, inplace=True)

# Compute the average similarity score for each element in column 'A'
avg_similarity = df_combinations_filter.groupby('A')['SPECTRAL_ANGLE'].mean().reset_index()

# Select the element with the lowest average similarity score
most_similar = avg_similarity.loc[avg_similarity['SPECTRAL_ANGLE'].idxmax(), 'A']

# -----------------------------------------------------------------------------
result[result["IDENTIFIER_x"] == most_similar]["SPECTRAL_ANGLE"]
# result[result["IDENTIFIER_x"] == "controllerType=0 controllerNumber=1 scan=14565"]["SPECTRAL_ANGLE"]

# -----------------------------------------------------------------------------
mgf1_spectrum = None
for spec in read_mgf(mgf1_filename):
    if spec.identifier == "controllerType=0 controllerNumber=1 scan=14565":
    # if spec.identifier == most_similar:
        mgf1_spectrum = spec
        break

if mgf1_spectrum is None:
    raise ValueError('Could not find the specified query spectrum in mgf 1')

mgf2_spectrum = None
for spec in read_mgf(mgf2_filename):
    if spec.identifier == mgf2_id:
        mgf2_spectrum = spec
        mgf2_spectrum.precursor_charge = mgf1_spectrum.precursor_charge
        break

if mgf2_spectrum is None:
    raise ValueError('Could not find the specified query spectrum in mgf 2')

mgf1_spectrum.annotate_proforma(peptide, 10, "ppm", ion_types="by")
mgf2_spectrum.annotate_proforma(peptide, 10, "ppm", ion_types="by")

# fig, ax = plt.subplots(figsize=(10, 5))
cm = 1/2.54  # centimeters in inches

width = 18*cm
height = 7*cm
# height = 11*cm
# width = 14
# height = 6
fig, ax = plt.subplots(figsize=(width, height))

plt.style.use(["seaborn-white", "seaborn-paper"])
sns.set_context("paper")

sup.mirror(mgf1_spectrum, mgf2_spectrum)
sns.despine()
ax.grid(False, which="both")
# ax.right_ax.grid(False)
name_plot = "paper-orbitrap-vs-orbitrap-HLA_133-" + peptide

plot_path = "/Users/adams/Projects/300K/Results/Figures/Mirror/" + name_plot + ".png"
plt.savefig(plot_path, dpi=600, bbox_inches='tight')
plt.close()

