# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logomaker

sns.set_context("paper")
sns.set_style("whitegrid")

cm = 1/2.54  # centimeters in inches
width = 9*cm
height = 6*cm

fig = plt.figure(constrained_layout=True, figsize=(width, height))
axes = fig.subplot_mosaic([['a'],['a']])

# load ww information matrix
ww_df = logomaker.get_example_matrix('ww_information_matrix',
                                     print_description=False)

# create Logo object
logomaker.Logo(ww_df,
                color_scheme='chemistry',
                vpad=.1,
                width=.8,
                ax = axes['a'])

sns.despine()
fig.tight_layout()

# plt.xticks(rotation = 45) # Rotates X-Axis Ticks by 45-degrees
plot_path = "/Users/adams/Projects/300K/Results/Figures/logo-test.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")

# style using Logo methods
axes['a'].style_xticks(anchor=0, spacing=5, rotation=45)
axes['a'].highlight_position(p=4, color='gold', alpha=.5)
axes['a'].highlight_position(p=26, color='gold', alpha=.5)

# style using Axes methods
ww_logo.ax.set_ylabel('information (bits)')
ww_logo.ax.set_xlim([-1, len(ww_df)])

