# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random

data = {'Category': ['Tryptic', 'Tryptic', 'Non-Tryptic', 'Non-Tryptic', 'Non-Tryptic'],
        'Charge': [2, 3, 1, 2, 3],
        'Training': [125075, 28734, 20652, 45263, 11662],
        'Validation': [13316, 3167, 1745, 4420, 1613],
        'Test': [10794, 3468, 1943, 4623, 1306]}

df = pd.DataFrame(data)
df_melted = pd.melt(df, id_vars=['Category', 'Charge'], var_name='Condition', value_name='Peptide')

sns.set(style="whitegrid")

g = sns.catplot(x='Condition', y='Peptide', hue='Charge', col='Category', data=df_melted, kind='bar', 
                height=3.5, aspect=1, palette=["#CDEAC0", "#0E1C36", "#7A8DB3"], legend=False)

g.despine(left=True)
g.set_ylabels("Number of PSMs")
g.set_xlabels("")
g.set_titles("{col_name}")
plt.legend(loc="lower left", bbox_to_anchor=(1.0, 0.1), frameon=False, title='Charge')
plt.tight_layout()
plot_path = "/Users/adams/Projects/300K/Results/Figures/bar/data_info_final.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
