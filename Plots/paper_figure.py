# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
import matplotlib.transforms as mtransforms
import os
from fundamentals.annotation.annotation import annotate_spectra
import spectrum_utils.plot as sup
from typing import Iterator
from typing import List
from spectrum_utils.spectrum import MsmsSpectrum
from pyteomics import mgf

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


# ---------------------- Import data ----------------------
# Data for bar plots

data_tryptic = {'Category': ['Tryptic', 'Tryptic'],
        'Charge': [2, 3],
        'Training': [125075, 28734],
        'Validation': [13316, 3167],
        'Test': [10794, 3468]}

data_nontryptic = {'Category': ['Non-Tryptic', 'Non-Tryptic', 'Non-Tryptic'],
        'Charge': [1, 2, 3],
        'Training': [20652, 45263, 11662],
        'Validation': [1745, 4420, 1613],
        'Test': [1943, 4623, 1306]}

df_tryptic = pd.DataFrame(data_tryptic)
df_nontryptic = pd.DataFrame(data_nontryptic)
tryptic_melt = pd.melt(df_tryptic, id_vars=['Category', 'Charge'], var_name='Condition', value_name='Peptide')
nontryptic_melt = pd.melt(df_nontryptic, id_vars=['Category', 'Charge'], var_name='Condition', value_name='Peptide')

# Data for Violin Plot

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/"  # nolint
HCD_path = base_path + "20230315104417_prediction.csv"
TIMS_path = base_path + "20230315111004_prediction.csv"

base = "/Users/adams/Projects/300K/2022-library-run/"
test_non_path = base + "Annotation/total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-test-non-tryptic.csv"
train_non_path = base + "Annotation/total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-train-non-tryptic.csv"
test_tryp_path = base + "Annotation/total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-test-tryptic.csv"
train_tryp_path = base + "Annotation/total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-train-tryptic.csv"

test_non = pd.read_csv(test_non_path)
test_tryp = pd.read_csv(test_tryp_path)
test_tryp = test_tryp[test_tryp["PRECURSOR_CHARGE"] > 1]
df_test = pd.concat([test_non, test_tryp], axis=0)

hcd = pd.read_csv(HCD_path)
tims = pd.read_csv(TIMS_path)

hcd["label"] = "HCD Prosit 2020"
tims["label"] = "HCD TIMS Prosit 2023"

df_prediction = pd.concat([hcd, tims], axis=0)
new_columns = {col: col.replace(" ", "_").upper() for col in df_prediction.columns}
df_prediction.rename(columns=new_columns, inplace=True)
df_prediction_map = pd.merge(df_prediction, df_test, on="SCAN_NUMBER", how="left")

df_prediction_filtered = df_prediction_map.loc[df_prediction_map.apply(lambda row: row["RAW_FILE"] in row["RAWFILE"], axis=1)]
df_prediction_filtered["type"] = df_prediction_filtered["pool_name"].apply(lambda x: x[4:8])
df_prediction_filtered["type"] = df_prediction_filtered["type"].replace({"firs":"Tryptic", "HLA_":"MHC-I", "HLA2":"MHC-II", "lysn":"LysN", "aspn":"AspN"})

# Data for mirror plots
charge = "1"
peptide = "GVDAANSAAQQY"

mgf1_filename="/Users/adams/Projects/300K/2022-library-run/Annotation/mapped-summed-mgf/TUM_HLA_133-3.mgf"
mgf2_filename="/Users/adams/Projects/300K/PXD038782-comparison/mgf/HLA_133_hcd_pred.mgf"
mgf3_filename="/Users/adams/Projects/300K/PXD038782-comparison/mgf/HLA_133_tims_pred.mgf"

msms_path="/Users/adams/Projects/300K/2022-library-run/msms-txt/TUM_HLA_133.txt"
mgf1_id = "3378"
mgf2_id = "14565"

mgf1_spectrum = None
for spec in read_mgf(mgf1_filename):
    if spec.identifier == mgf1_id:
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

mgf3_spectrum = None
for spec in read_mgf(mgf3_filename):
    if spec.identifier == mgf2_id:
        mgf3_spectrum = spec
        mgf3_spectrum.precursor_charge = mgf1_spectrum.precursor_charge
        break

if mgf3_spectrum is None:
    raise ValueError('Could not find the specified query spectrum in mgf 2')

mgf1_spectrum.annotate_proforma(peptide, 10, "ppm", ion_types="by")
mgf2_spectrum.annotate_proforma(peptide, 10, "ppm", ion_types="by")
mgf3_spectrum.annotate_proforma(peptide, 10, "ppm", ion_types="by")

# ---------------------- Figure 2 ----------------------
sns.set_context("paper")
sns.set_style("whitegrid")

cm = 1/2.54  # centimeters in inches
width = 18*cm
height = 14*cm

fig = plt.figure(constrained_layout=True, figsize=(width, height))
# top, bottom = fig.subfigures(nrows=2, ncols=1, height_ratios=[1,2])
axes = fig.subplot_mosaic([['a', 'a'], ['b', 'c'], ['b', 'd']])
# axes = fig.subplot_mosaic([['a', 'a', 'a', 'a', '', '', '', '', '', ''],['b', 'b','b', 'b', 'b', 'c', 'c', 'c', 'c', 'c'],['b', 'b','b', 'b','b', 'd', 'd', 'd', 'd', 'd']])
# axes = fig.subplot_mosaic([['a', 'a','a', '', '', ''],['b', 'b','b','c', 'c', 'c'],['b', 'b','b','d', 'd', 'd']])

# Bar Plots
ax_top = axes['a'].subplots(1, 2, sharey=True, width_ratios=[2, 1])

sns.barplot(data=tryptic_melt, x='Condition', y='Peptide', hue='Charge',
                palette=["#0E1C36", "#7A8DB3"], ax = ax_top[0])
ax_top[0].set_title("Tryptic")
ax_top[0].set_ylabel("Number of PSMs")

sns.barplot(data=nontryptic_melt, x='Condition', y='Peptide', hue='Charge',
                palette=["#CDEAC0", "#0E1C36", "#7A8DB3"], ax = ax_top[1])
ax_top[1].set_title("Non-Tryptic")
ax_top[1].set_ylabel('')
axes_list = [ax_top[0], ax_top[1]]
ax_top[1].legend(frameon=False)

for ax in axes_list:
    ax.set_ylim(0, 125000)
    ax.set_xlabel('\n\n')
    # ax.set_xticklabels(["Training", "Validation", "Test"], rotation=45)
    ylabels = ["{:,.0f}".format(x) + "K" for x in ax.get_yticks()/1000]
    ax.set_yticklabels(ylabels)
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.yaxis.grid(True, linewidth=0.5, which="major", color="lightgrey", alpha=0.5)
    ax.tick_params(axis="x", colors="#464646")
    ax.tick_params(axis="y", colors="#464646")

ax_top[0].get_legend().remove()
ax_top[1].yaxis.set_ticklabels([])

# Violin Plot
order = ["MHC-I", "MHC-II", "LysN", "AspN", "Tryptic"]

sns.violinplot(x="type",y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                    # edgecolor=None, linewidth=0, palette=["#0E1C36", "#CDEAC0"],
                    edgecolor=None, linewidth=0, palette=["#022873", "#7a8db3"],
                    hue_order=["HCD TIMS Prosit 2023", "HCD Prosit 2020"],
                    order=order,
                    data=df_prediction_filtered, split=True,
                    ax=axes['b'])

axes['b'].legend(
    loc='upper center', 
    bbox_to_anchor=(0.5, 1.1),
    ncol=2, frameon=False
)
axes['b'].set_ylabel("Spectral angle", color="black")
axes['b'].set(xlabel=None)
axes['b'].spines["right"].set_color("none")
axes['b'].spines["top"].set_color("none")
axes['b'].spines["left"].set_color("none")
axes['b'].spines["bottom"].set_color("none")
axes['b'].yaxis.grid(True, linewidth=0.5, which="major", color="lightgrey", alpha=0.5)
axes['b'].tick_params(axis="x", colors="#464646")
axes['b'].tick_params(axis="y", colors="#464646")
ax_x = plt.gca().xaxis
ax_x.set_tick_params(pad=-5)

# Mirror Plot
sup.mirror(mgf1_spectrum, mgf2_spectrum, ax = axes['c'])
axes['c'].grid(False, which="both")

sup.mirror(mgf1_spectrum, mgf3_spectrum, ax = axes['d'])
axes['d'].grid(False, which="both")

sns.despine()
plt.tight_layout()

for label, ax in axes.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='x-large', va='bottom', weight='bold')

plot_path = "/Users/adams/Projects/300K/Results/Figures/paper-fig-2.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")