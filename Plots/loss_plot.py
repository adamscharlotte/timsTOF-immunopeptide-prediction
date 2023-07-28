# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
import seaborn as sns
import os

# ---------------------- Import data ----------------------
file_path = "/Users/adams/Projects/300K/2022-library-run/total-calibrated-tryp.csv"
loss_df = pd.read_csv(file_path, sep=",")
loss_df.columns

sns.set_context("paper")
sns.set_style("whitegrid", {'axes.grid' : False})

cm = 1/2.54  # centimeters in inches
width = 18*cm
height = 8*cm
fig = plt.figure(figsize=(width, height))
axes = fig.subplot_mosaic([['a']])

loss_df.plot(x="Step", y=["total-calibrated-tryp - loss", "total-calibrated-tryp - val_loss"], ax=axes['a'], color=["#022873", "#7a8db3"],)
axes['a'].set_xlabel("Epoch")
axes['a'].set_ylabel("Loss")
axes['a'].set_xlim(0, 33.5)
axes['a'].set_ylim(0, 0.30)
axes['a'].legend(['Train', 'Validation'])
axes['a'].tick_params(axis='both', which='major', pad=-0.5)


sns.despine()
fig.tight_layout(pad=0.5)
plot_path = "/Users/adams/Projects/300K/Results/Figures/supplementary-fig-4.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
