# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
import matplotlib.transforms as mtransforms


sns.set_context("paper")
sns.set_style("whitegrid")

cm = 1/2.54  # centimeters in inches
width = 18*cm
height = 14*cm

fig = plt.figure(constrained_layout=True, figsize=(width, height))
axes = fig.subplot_mosaic([['a', 'a', '', '', ''],['b', 'b','b', 'b', 'b']])

# axes = fig.subplot_mosaic([['0', '0'],['1', '1']],
#                           gridspec_kw={'width_ratios':[1, 2]})

# Bar Plots
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

sns.barplot(data=tryptic_melt, x='Condition', y='Peptide', hue='Charge',
                palette=["#0E1C36", "#7A8DB3"], ax = axes['a'])
axes['a'].set_title("Tryptic")
axes['a'].set_ylabel("Number of PSMs")

sns.barplot(data=nontryptic_melt, x='Condition', y='Peptide', hue='Charge',
                palette=["#CDEAC0", "#0E1C36", "#7A8DB3"], ax = axes[''])
axes[''].set_title("Non-Tryptic")
axes[''].set_ylabel('')

axes_list = [axes['a'], axes['']]

for ax in axes_list:
    ax.set_ylim(0, 125000)
    ax.set_xlabel('')
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

axes['a'].get_legend().remove()

# Violin Plot
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
    ncol=2, 
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

# for label, ax in axes.items():
#     # label physical distance in and down:
#     # trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
#     # ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
#     #         fontsize='medium', verticalalignment='top', fontfamily='serif',
#     #         bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
#     ax.set_title(label, loc='left', fontsize='large', weight='bold')

for label, ax in axes.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='x-large', va='bottom', weight='bold')

# for i, (ax, c) in enumerate(zip(axes, 'abcd')):
#     ax.annotate(c, xy=(-0.17, 1.1), xycoords='axes fraction',
#                 fontsize='xx-large', weight='bold')

sns.despine()
fig.tight_layout()

# plt.xticks(rotation = 45) # Rotates X-Axis Ticks by 45-degrees
plot_path = "/Users/adams/Projects/300K/Results/Figures/paper-fig-2.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")