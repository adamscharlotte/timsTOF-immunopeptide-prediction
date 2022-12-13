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

pool = "TUM_HLA2_88"
precursor_1 = "9947"
precursor_2 = "9998"
charge = "1"
sequence = "LDGFSVK"
spectral_angle = "0.6741211572229355"
# query_filename = "/Users/adams/Projects/SARS-CoV-2/Workspace/mgf/qx017092.mgf"
# library_filename = "/Users/adams/Projects/SARS-CoV-2/Workspace/tmp/massivekb_sars_cov_2_targetdecoy.splib"
base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
query_filename = base_path + "precursor-consensus/annotated-20ppm-mgf/" + pool + ".mgf"
# query_filename = base_path + "precursor-consensus/un-annotated-mgf/" + pool + ".mgf"

# Retrieve information on the requested query.
top_id = precursor_1
annot_id = precursor_2
query_usi = pool + "_" + charge + "_" + sequence
name_plot = pool + "_" + precursor_1 + "_" + precursor_2 + "_" + sequence #[UNIMOD:121]    [UNIMOD:21]    [UNIMOD:21#g1]	[-28.031300]
txt = 'm/z'
top_spectrum_number = top_id
annot_spectrum_number = annot_id
# ssm = ssms.loc[top_id]
# library_id = ssm['accession']
# score = ssm['search_engine_score[1]']

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

# Add annotations to the query spectrum.
max_intensity = top_spectrum.intensity.max()

ax.set_ylim(-1.1, 1.05)
ax.text(0.5, 1.06, f'{query_usi}',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='x-large', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 1.02, f'Precursor (top): {top_id}, '
                f'Precursor (bottom): {annot_id}, '
                f'Charge: {charge}, '
                f'SA: {spectral_angle}',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='large', transform=ax.transAxes)

# plot_path = "/Users/adams/Projects/300K/Results/Figures/Mirror/precursor/un-annotated/" + name_plot + ".png"
plot_path = "/Users/adams/Projects/300K/Results/Figures/Mirror/precursor/annotated/" + name_plot + ".png"
plt.savefig(plot_path, dpi=300, bbox_inches='tight')
plt.close()
