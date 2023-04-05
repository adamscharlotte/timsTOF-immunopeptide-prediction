# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import matplotlib.pyplot as plt
import statistics
import pandas as pd
import numpy as np
import seaborn as sns

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/"  # nolint
file_1_path = base_path + "20230315104417_prediction.csv"
file_2_path = base_path + "20230315111004_prediction.csv"
file_3_path = base_path + "20230316105504_prediction.csv"

base = "/Users/adams/Projects/300K/2022-library-run/"
test_non_path = base + "Annotation/total-scan-consensus/split/non-tryptic/annotated-40-ppm-test.csv"
train_non_path = base + "Annotation/total-scan-consensus/split/non-tryptic/annotated-40-ppm-train.csv"
test_tryp_path = base + "Annotation/total-scan-consensus/split/tryptic/annotated-40-ppm-test.csv"
train_tryp_path = base + "Annotation/total-scan-consensus/split/tryptic/annotated-40-ppm-train.csv"

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
df2["label"] = "HCD TOF Prosit 2023"

df_prediction = pd.concat([df1, df2], axis=0)
new_columns = {col: col.replace(" ", "_").upper() for col in df_prediction.columns}
df_prediction.rename(columns=new_columns, inplace=True)
df_prediction_map = pd.merge(df_prediction, df_test, on="SCAN_NUMBER", how="left")

df_prediction_filtered = df_prediction_map.loc[df_prediction_map.apply(lambda row: row["RAW_FILE"] in row["RAWFILE"], axis=1)]
df_prediction_filtered["type"] = df_prediction_filtered["pool_name"].apply(lambda x: x[4:8])
df_prediction_filtered["type"] = df_prediction_filtered["type"].replace({"firs":"Tryptic", "HLA_":"MHC-I", "HLA2":"MHC-II", "lysn":"LysN", "aspn":"AspN"})

df_prediction_filtered.columns
df_prediction_map.columns

order = ["MHC-I", "MHC-II", "LysN", "AspN", "Tryptic"]

plt.rcParams["figure.figsize"] = [8, 5]
median_sa1 = 0.96
median_sa2 = 0.36
median_sa3 = 0.99
median_pcorr1 = 0.56
median_pcorr2 = 0.87
fig = plt.figure()
ax = fig.gca()
ax = sns.violinplot(x="type",y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                    edgecolor=None, linewidth=0, palette=["#0E1C36", "#CDEAC0"],
                    hue_order=["HCD TOF Prosit 2023", "HCD Prosit 2020"],
                    order=order,
                    data=df_prediction_filtered, split=True)
ax.axhline(median_sa1, color="black", linestyle="--", linewidth=1, xmax=0.2)
ax.axhline(median_sa2, color="black", linestyle="--", linewidth=1, xmin=0.2, xmax= 0.4)
ax.axhline(median_sa3, color="black", linestyle="--", linewidth=1, xmin=0.4, xmax= 0.6)
ax.axhline(median_pcorr1, color="black", linestyle="--", linewidth=1, xmin=0.6, xmax=0.8)
ax.axhline(median_pcorr2, color="black", linestyle="--", linewidth=1, xmin=0.8)
# plt.legend(loc="lower right")
ax.legend(loc="lower left", bbox_to_anchor=(1.0, 0.1), frameon=False)
#sns.move_legend(ax, "lower right")
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
plot_path = "/Users/adams/Projects/300K/Results/Figures/model-performance/violin-py/test-vs-train/type.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
# plt.show()


def violin_plot(df, x_label, label_name, size):
    plt.rcParams["figure.figsize"] = size
    fig = plt.figure()
    ax = fig.gca()
    ax = sns.violinplot(x=x_label, y="SPECTRAL_ANGLE", hue="LABEL", inner=None,
                        edgecolor=None, linewidth=0, palette=["#0E1C36", "#CDEAC0"],
                        hue_order=["HCD TOF Prosit 2023", "HCD Prosit 2020"],
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
    plot_path = "/Users/adams/Projects/300K/Results/Figures/model-performance/violin-py/test-vs-train/" + x_label + ".png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")

violin_plot(df_prediction_filtered, "PRECURSOR_CHARGE", "Precursor charge", [6, 5])

df_prediction_filtered.columns
df_prediction_filtered["n_term"] = df_prediction_filtered["OBS_SEQUENCE"].apply(lambda x: x[0:1])
df_prediction_filtered["c_term"] = df_prediction_filtered["OBS_SEQUENCE"].str[-1]
df_prediction_filtered["length"] = df_prediction_filtered["OBS_SEQUENCE"].str.len()

violin_plot(df_prediction_filtered, "n_term", "N-terminus", [14, 5])
violin_plot(df_prediction_filtered, "c_term", "C-terminus", [14, 5])
violin_plot(df_prediction_filtered, "length", "Peptide length", [14, 5])

