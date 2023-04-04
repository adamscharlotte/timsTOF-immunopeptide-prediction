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

# def annotate_ion_type(annotation, ion_types="aby"):
#     if annotation.ion_type[0] in ion_types:
#         if abs(annotation.isotope) == 1:
#             iso = "+i" if annotation.isotope > 0 else "-i"
#         elif annotation.isotope != 0:
#             iso = f"{annotation.isotope:+}i"
#         else:
#             iso = ""
#         nl = {"-NH3": "*", "-H2O": "o"}.get(annotation.neutral_loss, "")
#         return f"{annotation.ion_type}{iso}{'+' * annotation.charge}{nl}"
#     else:
#         return ""

# def annotate_ion_type(annotation, ion_types="by"):
#     if annotation.ion_type[0] in ion_types:
#         if abs(annotation.isotope) == 1:
#             iso = "+i" if annotation.isotope > 0 else "-i"
#         elif annotation.isotope != 0:
#             iso = f"{annotation.isotope:+}i"
#         else:
#             iso = ""
#         return f"{annotation.ion_type}{iso}{'+' * annotation.charge}{nl}"
#     else:
#         return ""

# -----------------------------------------------------------------------------
charge = "1"
peptide = "GVDAANSAAQQY"
name_plot = "tims-vs-orbitrap"
spectral_angle = "xx"

mgf1_filename="/Users/adams/Downloads/02446d_GD1-TUM_HLA_133_01_01-3xHCD-1h-R4.mgf"
mgf2_filename="/Users/adams/Projects/300K/2022-library-run/Annotation/mapped-summed-mgf/TUM_HLA_133-3.mgf"
# msms1_path="/Users/adams/Downloads/TUM_HLA_15_01_01_3xHCD-1h-R4-unspecific/msms.txt"
msms1_path="/Users/adams/Downloads/TUM_HLA_133_01_01_3xHCD-1h-R4-unspecific/msms.txt"
msms2_path="/Users/adams/Projects/300K/2022-library-run/msms-txt/TUM_HLA_133.txt"

mgf1_id = "controllerType=0 controllerNumber=1 scan=" + "14634"
mgf2_id = "3378"
# -----------------------------------------------------------------------------
df_mgf1 = read_mgf_df(mgf1_filename)
df_mgf2 = read_mgf_df(mgf2_filename)

def order_peaks(df):
    for i in range (len(df)):
        zipped_list = zip(df.iloc[i]["mz"], df.iloc[i]["intensity"])
        sorted_pair = sorted(zipped_list)
        tuples = zip(*sorted_pair)
        list_1, list_2 = [ list(tuple) for tuple in  tuples]
        df.at[i, "mz"] = list_1
        df.at[i, "intensity"] = list_2

order_peaks(df_mgf1)
order_peaks(df_mgf2)

df_mgf1['id'] = df_mgf1['identifier'].str[41:].astype(int)
df_mgf2['id'] = df_mgf2['identifier'].astype(int)

df_mgf2.columns

msms1 = pd.read_csv(msms1_path, sep='\t')
msms2 = pd.read_csv(msms2_path, sep='\t')

def merge_df(df_mgf, msms):
    merged_df = pd.merge(df_mgf.drop('charge', axis=1), msms.drop('Intensities', axis=1), left_on='id', right_on='Scan number', how='inner')
    new_columns = {col: col.replace(' ', '_').upper() for col in merged_df.columns}
    merged_df.rename(columns=new_columns, inplace=True)
    merged_df.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)
    merged_df.rename(columns = {"INTENSITY": "INTENSITIES"}, inplace=True)
    merged_df['PRECURSOR_CHARGE'] = merged_df['PRECURSOR_CHARGE'].astype(int)
    merged_df["MODIFIED_SEQUENCE"] = maxquant_to_internal(merged_df["MODIFIED_SEQUENCE"].to_numpy())
    return merged_df

merged_df1 = merge_df(df_mgf1, msms1)
merged_df2 = merge_df(df_mgf2, msms2)

annot_df1 = annotate_spectra(merged_df1)
annot_df2 = annotate_spectra(merged_df2)

full_df1 = pd.concat([merged_df1.drop(columns = ["INTENSITIES", "MZ"]), annot_df1], axis=1)
full_df2 = pd.concat([merged_df2.drop(columns = ["INTENSITIES", "MZ"]), annot_df2], axis=1)

filtered_df1 = full_df1[full_df1["SEQUENCE"] == peptide][full_df1["PRECURSOR_CHARGE"] == 1]
filtered_df1.rename(columns = {"INTENSITIES": "INTENSITIES_1", "MZ":"MZ_1"}, inplace=True)
filtered_df2 = full_df2[full_df2["SEQUENCE"] == peptide]
filtered_df2.rename(columns = {"INTENSITIES": "INTENSITIES_2", "MZ":"MZ_2"}, inplace=True)

result = pd.merge(filtered_df1.assign(dummy=1), filtered_df2.assign(dummy=1), on='dummy', how='outer').drop('dummy', axis=1)
result["SPECTRAL_ANGLE"] = result[['INTENSITIES_1','INTENSITIES_2']].apply(lambda x : get_spectral_angle(x), axis=1)
result["SPECTRAL_ANGLE"].fillna(0, inplace=True)

highest_row = result.nlargest(1, 'SPECTRAL_ANGLE')
mgf1_id = highest_row["IDENTIFIER_x"].item()

# -----------------------------------------------------------------------------

mgf1_spectrum = None
for spec in read_mgf(mgf1_filename):
    if spec.identifier == mgf1_id:
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
mgf2_spectrum.annotate_proforma(peptide, 20, "ppm", ion_types="by")

fig, ax = plt.subplots(figsize=(12, 6))
sup.spectrum(mgf1_spectrum, grid=False, ax=ax)
ax.set_title(peptide, fontdict={"fontsize": "xx-large"})
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.savefig("proforma_ex1.png", bbox_inches="tight", dpi=300, transparent=True)
plt.close()

fig, ax = plt.subplots(figsize=(12, 6))
sup.spectrum(mgf2_spectrum, grid=False, ax=ax)
ax.set_title(peptide, fontdict={"fontsize": "xx-large"})
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.savefig("proforma_ex2.png", bbox_inches="tight", dpi=300, transparent=True)
plt.close()

fig, ax = plt.subplots(figsize=(12, 6))
sup.mirror(mgf1_spectrum, mgf2_spectrum)
plt.savefig("mirror.png", dpi=300, bbox_inches="tight", transparent=True)
plt.close()


# --------
mgf1_spectrum.annotate_proforma(
    peptide,
    fragment_tol_mass=10,
    fragment_tol_mode="ppm",
    ion_types="by",
    max_ion_charge=3
)

fig, ax = plt.subplots(figsize=(12, 6))
sup.spectrum(mgf1_spectrum, annot_fmt=annotate_ion_type, grid=False, ax=ax)
ax.set_title(peptide, fontdict={"fontsize": "xx-large"})
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.savefig("annot_fmt.png", dpi=300, bbox_inches="tight", transparent=True)
plt.close()

if mgf2_spectrum is None:
    raise ValueError('Could not find the specified query spectrum in mgf 2')

set_matching_peaks(mgf1_spectrum, mgf2_spectrum)

peak_matches = spectrum_match.get_best_match(
    mgf1_spectrum, [mgf2_spectrum],
    config.fragment_mz_tolerance, config.allow_peak_shifts)[2]

plot.colors[None] = '#757575'
fig, ax = plt.subplots(figsize=(15, 7))

plot.mirror(mgf1_spectrum, mgf2_spectrum,
            {'color_ions': True, 'annotate_ions': False}, ax)


ax.set_ylim(-1.1, 1.05)
ax.text(0.5, 1.06, f'{mgf1_usi}',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='x-large', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 1.02, f'Precursor (top): {mgf1} {mgf1_id}, '
                f'Precursor (bottom): {mgf2} {mgf2_id}, '
                f'Charge: {charge}, '
                f'SA: {spectral_angle}',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='large', transform=ax.transAxes)

plot_path = "/Users/adams/Projects/300K/Results/Figures/Mirror/" + name_plot + ".png"
plt.savefig(plot_path, dpi=300, bbox_inches='tight')
plt.close()
