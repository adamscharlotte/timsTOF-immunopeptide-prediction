# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
import statistics
import pandas as pd
import numpy as np
import seaborn as sns
import glob

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/"  # nolint
file_1_path = base_path + "20230315104417_prediction.csv"
file_2_path = base_path + "20230315111004_prediction.csv"
file_3_path = base_path + "20230316105504_prediction.csv"

base = "/Users/adams/Projects/300K/2022-library-run/"
test_non_path = base + "Annotation/total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-test-non-tryptic.csv"
train_non_path = base + "Annotation/total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-train-non-tryptic.csv"
test_tryp_path = base + "Annotation/total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-test-tryptic.csv"
train_tryp_path = base + "Annotation/total-scan-consensus/calibrated-linear-40-ppm/calibrated-40-ppm-train-tryptic.csv"

test_non = pd.read_csv(test_non_path)
train_non = pd.read_csv(train_non_path)
test_tryp = pd.read_csv(test_tryp_path)
test_tryp = test_tryp[test_tryp["PRECURSOR_CHARGE"] > 1]
train_tryp = pd.read_csv(train_tryp_path)
train_tryp = train_tryp[train_tryp["PRECURSOR_CHARGE"] > 1]

df_test = pd.concat([test_non, test_tryp], axis=0)
df_train = pd.concat([train_non, train_tryp], axis=0)

# 2020 VS 2023 Model
df1 = pd.read_csv(file_1_path)
df2 = pd.read_csv(file_2_path)

df1["label"] = "HCD Prosit 2020"
df2["label"] = "HCD TIMS Prosit 2023"

df_prediction = pd.concat([df1, df2], axis=0)
new_columns = {col: col.replace(" ", "_").upper() for col in df_prediction.columns}
df_prediction.rename(columns=new_columns, inplace=True)
df_prediction_map = pd.merge(df_prediction, df_test, on="SCAN_NUMBER", how="left")

df_prediction_filtered = df_prediction_map.loc[df_prediction_map.apply(lambda row: row["RAW_FILE"] in row["RAWFILE"], axis=1)]
df_prediction_filtered["type"] = df_prediction_filtered["pool_name"].apply(lambda x: x[4:8])
df_prediction_filtered["type"] = df_prediction_filtered["type"].replace({"firs":"Tryptic", "HLA_":"MHC-I", "HLA2":"MHC-II", "lysn":"LysN", "aspn":"AspN"})

order = ["MHC-I", "MHC-II", "LysN", "AspN", "Tryptic"]

# # Calculate the medians
# path = "/Users/adams/Projects/300K/2022-library-run/Annotation/spectral-angle/pairwise/test-set/*.csv"
# all_files = glob.glob(path)
# dfs_sa = []

# for filename in all_files:
#     df = pd.read_csv(filename, header=0)
#     new_file_name = filename[88:-6]
#     df['file_name'] = filename  # add a new column with the file name
#     df['new_file_name'] = new_file_name
#     dfs_sa.append(df)

# sa_df = pd.concat(dfs_sa, ignore_index=True)
# median_df = sa_df.groupby('new_file_name')['SPECTRAL_ANGLE'].median().round(3)

# Plot
plt.rcParams["figure.figsize"] = [8, 5]
# median_sa1 = 0.678
# median_sa2 = 0.726
# median_sa3 = 0.625
# median_sa4 = 0.744
# median_sa5 = 0.676
fig = plt.figure()
sns.set_context("talk")

ax = fig.gca()
ax = sns.violinplot(x="type",y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                    # edgecolor=None, linewidth=0, palette=["#0E1C36", "#CDEAC0"],
                    edgecolor=None, linewidth=0, palette=["#022873", "#7a8db3"],
                    hue_order=["HCD TIMS Prosit 2023", "HCD Prosit 2020"],
                    order=order,
                    data=df_prediction_filtered, split=True)
# ax.axhline(median_sa1, color="black", linestyle="--", linewidth=1, xmax=0.2)
# ax.axhline(median_sa2, color="black", linestyle="--", linewidth=1, xmin=0.2, xmax= 0.4)
# ax.axhline(median_sa3, color="black", linestyle="--", linewidth=1, xmin=0.4, xmax= 0.6)
# ax.axhline(median_sa4, color="black", linestyle="--", linewidth=1, xmin=0.6, xmax=0.8)
# ax.axhline(median_sa5, color="black", linestyle="--", linewidth=1, xmin=0.8)
# plt.legend(loc="lower right")
ax.legend(loc="lower left", bbox_to_anchor=(1.0, 0.1), frameon=False)
ax.set_ylabel("Spectral angle", color="black")
ax.set(xlabel=None)
ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")
ax.spines["left"].set_color("none")
ax.spines["bottom"].set_color("none")
sns.set_style("whitegrid")
ax.yaxis.grid(True, linewidth=0.5, which="major", color="lightgrey", alpha=0.5)
plt.yticks([i/10 for i in range(0, 11)])
ax.tick_params(axis="x", colors="#464646")
ax.tick_params(axis="y", colors="#464646")
ax_x = plt.gca().xaxis
ax_x.set_tick_params(pad=-5)
plot_path = "/Users/adams/Projects/300K/Results/Figures/model-performance/violin-py/2020-vs-2023/type-ppt.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
# plt.show()


def violin_plot(df, x_label, label_name, size, path, h_order, p_color, x_order):
    plt.rcParams["figure.figsize"] = size
    fig = plt.figure()
    ax = fig.gca()
    ax = sns.violinplot(x=x_label, y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                        edgecolor=None, linewidth=0, palette=p_color,
                        hue_order=h_order, order=x_order,
                        data=df, split=True)
    ax.legend(loc="lower left", bbox_to_anchor=(1.0, 0.1), frameon=False)
    ax.set_ylabel("Spectral angle", color="black")
    ax.set_xlabel(x_label, color="black")
    print(ax.set_xlabel(label_name, color="black"))
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["bottom"].set_color("none")
    sns.set_style("whitegrid")
    ax.yaxis.grid(True, linewidth=0.5, which="major", color="lightgrey", alpha=0.5)
    plt.yticks([i/10 for i in range(0, 11)])
    ax.tick_params(axis="x", colors="#464646")
    ax.tick_params(axis="y", colors="#464646")
    ax_x = plt.gca().xaxis
    ax_x.set_tick_params(pad=-5)
    plot_path = "/Users/adams/Projects/300K/Results/Figures/model-performance/violin-py/" + path + ".png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")

df_prediction_filtered.columns
df_prediction_filtered["n_term"] = df_prediction_filtered["OBS_SEQUENCE"].apply(lambda x: x[0:1])
df_prediction_filtered["c_term"] = df_prediction_filtered["OBS_SEQUENCE"].str[-1]
df_prediction_filtered["length"] = df_prediction_filtered["OBS_SEQUENCE"].str.len()

n_term = df_prediction_filtered["n_term"].unique().tolist()
sort_n_term = sorted(n_term)

c_term = df_prediction_filtered["c_term"].unique().tolist()
sort_c_term = sorted(c_term)
order = ["MHC-I", "MHC-II", "LysN", "AspN", "Tryptic"]

violin_plot(df_prediction_filtered, "type", "", [8, 5],
    "2020-vs-2023/type_1", ["HCD TOF Prosit 2023", "HCD Prosit 2020"],
    ["#0E1C36", "#CDEAC0"], order)
violin_plot(df_prediction_filtered, "n_term", "N-terminus", [14, 5],
    "2020-vs-2023/n_term", ["HCD TOF Prosit 2023", "HCD Prosit 2020"],
    ["#0E1C36", "#CDEAC0"], n_term)
violin_plot(df_prediction_filtered, "c_term", "C-terminus", [14, 5],
    "2020-vs-2023/c_term", ["HCD TOF Prosit 2023", "HCD Prosit 2020"],
    ["#0E1C36", "#CDEAC0"])
violin_plot(df_prediction_filtered, "length", "Peptide length", [14, 5],
    "2020-vs-2023/length", ["HCD TOF Prosit 2023", "HCD Prosit 2020"],
    ["#0E1C36", "#CDEAC0"])
violin_plot(df_prediction_filtered, "COLLISION_ENERGY_ALIGNED_NORMED", "Normalized Aligned Collision Energy", [14, 5],
    "2020-vs-2023/COLLISION_ENERGY_ALIGNED_NORMED", ["HCD TOF Prosit 2023", "HCD Prosit 2020"],
    ["#0E1C36", "#CDEAC0"])
violin_plot(df_prediction_filtered, "aligned_collision_energy", "Aligned Collision Energy", [14, 5],
    "2020-vs-2023/aligned_collision_energy", ["HCD TOF Prosit 2023", "HCD Prosit 2020"],
    ["#0E1C36", "#CDEAC0"])

# -----------------------------------------------------------------------------
df2 = pd.read_csv(file_2_path)
df3 = pd.read_csv(file_3_path)

df2["label"] = "Test set"
df3["label"] = "Training set"

new_columns = {col: col.replace(" ", "_").upper() for col in df2.columns}
df2.rename(columns=new_columns, inplace=True)
new_columns = {col: col.replace(" ", "_").upper() for col in df3.columns}
df3.rename(columns=new_columns, inplace=True)

df_test_map = pd.merge(df2, df_test, on="SCAN_NUMBER", how="left")
df_train_map = pd.merge(df3, df_train, on="SCAN_NUMBER", how="left")
df_test_train_map = pd.concat([df_test_map, df_train_map], axis=0)

df_test_train_filtered = df_test_train_map.loc[df_test_train_map.apply(lambda row: row["RAW_FILE"] in row["RAWFILE"], axis=1)]
df_test_train_filtered["type"] = df_test_train_filtered["pool_name"].apply(lambda x: x[4:8])
df_test_train_filtered["type"] = df_test_train_filtered["type"].replace({"firs":"Tryptic", "HLA_":"MHC-I", "HLA2":"MHC-II", "lysn":"LysN", "aspn":"AspN"})
df_test_train_filtered["n_term"] = df_test_train_filtered["OBS_SEQUENCE"].apply(lambda x: x[0:1])
df_test_train_filtered["c_term"] = df_test_train_filtered["OBS_SEQUENCE"].str[-1]
df_test_train_filtered["length"] = df_test_train_filtered["OBS_SEQUENCE"].str.len()

n_term = df_test_train_filtered["n_term"].unique().tolist()
sort_n_term = sorted(n_term)

c_term = df_test_train_filtered["c_term"].unique().tolist()
sort_c_term = sorted(c_term)

df_test_train_map["RAW_FILE"]
df_test_train_filtered["length"]

df_test_train_filtered = df_test_train_filtered[df_test_train_filtered["length"] < 30]
df_test_train_filtered["aligned_collision_energy"] = df_test_train_filtered["aligned_collision_energy"].astype(int)
max(df_test_train_filtered["length"])
min(df_test_train_filtered["aligned_collision_energy"])

sns.set_context("paper")
sns.set_style("whitegrid")

cm = 1/2.54  # centimeters in inches
width = 18*cm
height = 12*cm
# constrained_layout=True, 
fig = plt.figure(figsize=(width, height))
axes = fig.subplot_mosaic([
     ['a', 'b', 'b', 'b'],
     ['c', 'c', 'c', 'c']])

sns.violinplot(x="PRECURSOR_CHARGE", y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                    edgecolor=None, linewidth=0,
                    palette=["#022873", "#7a8db3"],
                    hue_order=["Training set", "Test set"], order=[1,2,3],
                    data=df_test_train_filtered, split=True, ax=axes['a'])
axes['a'].set_xlabel("Precursor charge\n")

sns.violinplot(x="length", y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                    edgecolor=None, linewidth=0,
                    palette=["#022873", "#7a8db3"],
                    hue_order=["Training set", "Test set"], #order=[1,2,3],
                    data=df_test_train_filtered, split=True, ax=axes['b'])
axes['b'].set_xlabel("Peptide length")

sns.violinplot(x="aligned_collision_energy", y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                    edgecolor=None, linewidth=0,
                    palette=["#022873", "#7a8db3"],
                    hue_order=["Training set", "Test set"], #order=[1,2,3],
                    data=df_test_train_filtered, split=True, ax=axes['c'])
axes['c'].set_xlabel("Aligned collision energy")

ax_list = [axes['a'], axes['b'], axes['c']]
for ax in ax_list:
    ax.set_ylabel("Spectral angle")
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.yaxis.grid(True, linewidth=0.5, which="major", color="lightgrey", alpha=0.5)
    ax.get_legend().remove()
    ax.set_ylim(0, 1)

sns.despine()
fig.tight_layout()

handles, labels = axes["b"].get_legend_handles_labels()
fig.legend(handles=handles, loc='upper center', ncol=2,
                # handlelength=1, handleheight=1, handletextpad=0.4,
                bbox_to_anchor=(0.5, 1.07),
                # fontsize="7",
                frameon=False, columnspacing = 0.7)

for label, ax in axes.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    ax.text(0.0, 1, label, transform=ax.transAxes + trans,
            fontsize='x-large', va='bottom', weight='bold')

plot_path = "/Users/adams/Projects/300K/Results/Figures/supplementary-violin-1.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")



sns.set_context("paper")
sns.set_style("whitegrid")

cm = 1/2.54  # centimeters in inches
width = 18*cm
height = 12*cm
# constrained_layout=True, 
fig = plt.figure(figsize=(width, height))
axes = fig.subplot_mosaic([['a'], ['b']])

sns.violinplot(x="n_term", y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                    edgecolor=None, linewidth=0,
                    palette=["#022873", "#7a8db3"],
                    hue_order=["Training set", "Test set"], order=sort_n_term,
                    data=df_test_train_filtered, split=True, ax=axes['a'])
axes['a'].set_xlabel("N-terminus\n")

sns.violinplot(x="c_term", y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                    edgecolor=None, linewidth=0,
                    palette=["#022873", "#7a8db3"],
                    hue_order=["Training set", "Test set"], order=sort_c_term,
                    data=df_test_train_filtered, split=True, ax=axes['b'])
axes['b'].set_xlabel("C-terminus")


ax_list = [axes['a'], axes['b']]
for ax in ax_list:
    ax.set_ylabel("Spectral angle")
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.yaxis.grid(True, linewidth=0.5, which="major", color="lightgrey", alpha=0.5)
    ax.get_legend().remove()
    ax.set_ylim(0, 1)

sns.despine()
fig.tight_layout()

handles, labels = axes["b"].get_legend_handles_labels()
fig.legend(handles=handles, loc='upper center', ncol=2,
                # handlelength=1, handleheight=1, handletextpad=0.4,
                bbox_to_anchor=(0.5, 1.07),
                # fontsize="7",
                frameon=False, columnspacing = 0.7)

for label, ax in axes.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    ax.text(0.0, 1, label, transform=ax.transAxes + trans,
            fontsize='x-large', va='bottom', weight='bold')

plot_path = "/Users/adams/Projects/300K/Results/Figures/supplementary-violin-2.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")

