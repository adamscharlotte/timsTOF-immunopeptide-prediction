# /Users/adams/anaconda3/envs/ann_solo/bin/python3
# /Users/adams/opt/miniconda3/envs/ann_solo/bin/python3

from typing import Sequence
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
# import spectrum_utils
from spectrum_utils import plot
from spectrum_utils.spectrum import PeptideFragmentAnnotation

from ann_solo import reader
from ann_solo import spectrum_match
from ann_solo.config import config
from ann_solo.spectrum import process_spectrum

import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument('Filename', type=str)					# Filename
parser.add_argument('Scan', type=str)
parser.add_argument('Charge', type=str)
parser.add_argument('Sequence', type=str)
args = parser.parse_args()

# pool = args.pool
# scan = args.Scan
# charge = args.Charge
# sequence = args.Sequence

pool = "TUM_HLA2_12"
precursor = "13045"
frame_1 = "4264"
frame_2 = "4269"
charge = "1"
sequence = "PLAGQEAVVDL"
# query_filename = "/Users/adams/Projects/SARS-CoV-2/Workspace/mgf/qx017092.mgf"
# library_filename = "/Users/adams/Projects/SARS-CoV-2/Workspace/tmp/massivekb_sars_cov_2_targetdecoy.splib"
base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
query_filename = base_path + "full-truncated-qc/un-annotated-mgf/" + pool + ".mgf"

def set_matching_peaks(annot_spectrum, top_spectrum):
    peak_matches = spectrum_match.get_best_match(
        top_spectrum, [annot_spectrum],
        config.fragment_mz_tolerance, config.allow_peak_shifts)[2]
    top_spectrum.annotation = np.full_like(top_spectrum.mz, None, object)
    for peak_match in peak_matches:
        library_annotation = annot_spectrum.annotation[peak_match[1]]
        if library_annotation is not None and library_annotation.ion_type in 'by':
            top_spectrum.annotation[peak_match[0]] = library_annotation
        # else:
        #     fragment_annotation = PeptideFragmentAnnotation(1, 1, 'z', 0)
        #     fragment_annotation.ion_type = 'unknown'
        #     top_spectrum.annotation[peak_match[0]] =\
        #         annot_spectrum.annotation[peak_match[1]] =\
        #         fragment_annotation


# # Read the mzTab file.
# metadata = {}
# mztabFile = "/Users/adams/Projects/SARS-CoV-2/Workspace/ANN-SoLo/mztab/" + filename + ".mztab"
# mztabFile = "/Users/adams/Projects/SARS-CoV-2/Workspace/ANN-SoLo/mztab_selection/" + filename + ".mztab"
# with open(mztabFile) as f_mztab:
#     for line in f_mztab:
#         line_split = line.strip().split('\t')
#         if line_split[0] == 'MTD':
#             metadata[line_split[1]] = line_split[2]
#         else:
#             break   # Metadata lines should be on top.

# ssms = reader.read_mztab_ssms(mztabFile)
# # make sure the SSM ids are strings.
# ssms.index = ssms.index.map(str)
# ssms.columns

# # Recreate the search configuration.
# settings = []
# # Search settings.
# for key in metadata:
#     if 'software[1]-setting' in key:
#         param = metadata[key][: metadata[key].find(' ')]
#         value = metadata[key][metadata[key].rfind(' ') + 1:]
#         if value != 'None':
#             if value != 'False':
#                 settings.append('--{}'.format(param))
#             if value not in ('False', 'True'):
#                 settings.append(value)

# # File names.
# settings.append('dummy_spectral_library_filename')
# settings.append('dummy_query_filename')
# settings.append('dummy_output_filename')
# config.parse(' '.join(settings))

# Retrieve information on the requested query.
top_id = frame_1 + "_" + precursor
annot_id = frame_2 + "_" + precursor
query_usi = pool + ":" + precursor + ":" + frame_1 + ":" + frame_2 + ":" + charge + ":" + sequence #[UNIMOD:121]    [UNIMOD:21]    [UNIMOD:21#g1]	[-28.031300]
txt = 'm/z'
top_spectrum_number = top_id
annot_spectrum_number = annot_id
# ssm = ssms.loc[top_id]
# library_id = ssm['accession']
# score = ssm['search_engine_score[1]']


# Read library and query spectrum.
# with reader.SpectralLibraryReader(query_filename) as lib_reader:
#     annot_spectrum = lib_reader.get_spectrum(bottom_id, False)

annot_spectrum = None
for spec in reader.read_mgf(query_filename):
    if spec.identifier == annot_spectrum_number:
        annot_spectrum = spec
        break

top_spectrum = None
for spec in reader.read_mgf(query_filename):
    if spec.identifier == top_spectrum_number:
        top_spectrum = spec
        top_spectrum.precursor_charge = annot_spectrum.precursor_charge
        break

# verify that the query spectrum was found
if top_spectrum is None:
    raise ValueError('Could not find the specified query spectrum')

# Set the matching peaks in the query spectrum to correctly color them.
# set_matching_peaks(annot_spectrum, top_spectrum)
# Modify the colors to differentiate non-matching peaks.
plot.colors[None] = '#757575'

# Plot the match.
fig, ax = plt.subplots(figsize=(15, 7))
# Plot with annotations.
# plot.mirror(top_spectrum, annot_spectrum,
# 			{'color_ions': True, 'annotate_ions': True}, ax)

plot.mirror(top_spectrum, annot_spectrum,
            {'color_ions': True, 'annotate_ions': False}, ax)
# Add annotations to the library spectrum.
max_intensity = annot_spectrum.intensity.max()
# for i, annotation in enumerate(annot_spectrum.annotation):
#     if annotation is not None and annotation.ion_type != 'unknown':
#         x = annot_spectrum.mz[i]
#         y = - annot_spectrum.intensity[i] / max_intensity
#         ax.text(x, y, str(annotation),
#                 color=plot.colors[annotation.ion_type], zorder=5,
#                 horizontalalignment='right', verticalalignment='center',
#                 rotation=90, rotation_mode='anchor')

# Add annotations to the query spectrum.
max_intensity = top_spectrum.intensity.max()
# for i, annotation in enumerate(top_spectrum.annotation):
#     if annotation is not None and annotation.ion_type != 'unknown':
#         x = top_spectrum.mz[i]
#         y = top_spectrum.intensity[i] / max_intensity + 0.14
#         ax.text(x, y, str(annotation),
#                 color=plot.colors[annotation.ion_type], zorder=5,
#                 horizontalalignment='right', verticalalignment='center',
#                 rotation=90, rotation_mode='anchor')

ax.set_ylim(-1.1, 1.05)
ax.text(0.5, 1.06, f'{query_usi}',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='x-large', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 1.02,  f'Precursor ${txt}$ (top): {top_spectrum.precursor_mz:.4f}, '
                    f'Library ${txt}$ (bottom): {annot_spectrum.precursor_mz:.4f}, '
                    f'Charge: {top_spectrum.precursor_charge}',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='large', transform=ax.transAxes)

plot_path = "/Users/adams/Projects/300K/Results/Figures/Mirror/" + query_usi + ".png"
plt.savefig(plot_path, dpi=300, bbox_inches='tight')
plt.close()
