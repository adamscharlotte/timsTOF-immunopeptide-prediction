# /Users/adams/opt/miniconda3/envs/ann_solo/bin/python3
import os
import numpy as np
import matplotlib.pyplot as plt
from spectrum_utils import plot
from spectrum_utils.spectrum import PeptideFragmentAnnotation
from ann_solo import reader
from ann_solo import spectrum_match
from ann_solo.config import config

# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('Filename', type=str)					# Filename
# parser.add_argument('Scan', type=str)
# parser.add_argument('Charge', type=str)
# parser.add_argument('Sequence', type=str)
# args = parser.parse_args()

def set_matching_peaks(mgf_1, mgf_2):
    peak_matches = spectrum_match.get_best_match(
        mgf_1, [mgf_2],
        config.fragment_mz_tolerance, config.allow_peak_shifts)[2]
    mgf_1.annotation = np.full_like(mgf_1.mz, None, object)
    for peak_match in peak_matches:
        library_annotation = mgf_2.annotation[peak_match[1]]
        if library_annotation is not None and library_annotation.ion_type in 'by':
            mgf_1.annotation[peak_match[0]] = library_annotation



# pool = "TUM_HLA2_88"
# precursor_1 = "9947"
# precursor_2 = "9998"
# charge = "1"
# sequence = "LDGFSVK"
# spectral_angle = "0.6741211572229355"
mgf1 = "TUM_HLA2_12"
mgf2 = "TUM_HLA2_88"
charge = "1"
sequence = "GVDAANSAAQQY"
name_plot = "Test 1"
spectral_angle = "xx"

# base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
# mgf1_filename = base_path + "precursor-consensus/un-annotated-mgf/TUM_HLA2_12.mgf"
# mgf2_filename = base_path + "precursor-consensus/un-annotated-mgf/TUM_HLA2_88.mgf"

mgf1_filename="/Users/adams/Downloads/02446d_GD1-TUM_HLA_133_01_01-3xHCD-1h-R4.mgf"
mgf2_filename="/Users/adams/Downloads/02446d_GD1-TUM_HLA_133_01_01-3xHCD-1h-R4.mgf"

mgf1_id = "controllerType=0 controllerNumber=1 scan=" + "14634"
mgf2_id = "controllerType=0 controllerNumber=1 scan=" + "12615"
mgf1_usi = mgf1 + "_" + charge + "_" + sequence
mgf2_usi = mgf2 + "_" + charge + "_" + sequence

txt = 'm/z'

mgf1_spectrum = None
for spec in reader.read_mgf(mgf1_filename):
    if spec.identifier == mgf1_id:
        mgf1_spectrum = spec
        break

mgf2_spectrum = None
for spec in reader.read_mgf(mgf2_filename):
    if spec.identifier == mgf2_id:
        mgf2_spectrum = spec
        mgf2_spectrum.precursor_charge = mgf1_spectrum.precursor_charge
        break

if mgf1_spectrum is None:
    raise ValueError('Could not find the specified query spectrum in mgf 1')

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

# Add annotations to the mgf 2 spectrum.
max_intensity = mgf2_spectrum.intensity.max()
for i, annotation in enumerate(mgf2_spectrum.annotation):
    if annotation is not None and annotation.ion_type != 'unknown':
        x = mgf2_spectrum.mz[i]
        y = - mgf2_spectrum.intensity[i] / max_intensity
        ax.text(x, y, str(annotation),
                color=plot.colors[annotation.ion_type], zorder=5,
                horizontalalignment='right', verticalalignment='center',
                rotation=90, rotation_mode='anchor')

# Add annotations to the mgf 1 spectrum.
max_intensity = mgf1_spectrum.intensity.max()
for i, annotation in enumerate(mgf1_spectrum.annotation):
	if annotation is not None and annotation.ion_type != 'unknown':
		x = mgf1_spectrum.mz[i]
		y = mgf1_spectrum.intensity[i] / max_intensity + 0.14
		ax.text(x, y, str(annotation),
				color=plot.colors[annotation.ion_type], zorder=5,
				horizontalalignment='right', verticalalignment='center',
				rotation=90, rotation_mode='anchor')

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
