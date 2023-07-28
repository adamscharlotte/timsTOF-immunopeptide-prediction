# /Users/adams/opt/miniconda3/envs/spectrum_utils/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import logomaker
from logomaker.src.matrix import transform_matrix

def transform_psfm(matrix):
    t_matrix = transform_matrix(matrix,
                 from_type = "probability",
                 to_type = "information")
    return t_matrix

# Load data

psfm_g1 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_g1.txt", sep=" ")
psfm_g2 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_g2.txt", sep=" ")
psfm_g3 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_g3.txt", sep=" ")
psfm_g4 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_g4.txt", sep=" ")
psfm_s1 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_s1.txt", sep=" ")
psfm_s2 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_s2.txt", sep=" ")
psfm_s3 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_s3.txt", sep=" ")
psfm_s4 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_s4.txt", sep=" ")
psfm_l1 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_l1.txt", sep=" ")
psfm_l2 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_l2.txt", sep=" ")
psfm_l3 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_l3.txt", sep=" ")
psfm_l4 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_l4.txt", sep=" ")
psfm_0101 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_0101.txt", sep=" ")
psfm_0202 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_0202.txt", sep=" ")
psfm_4403 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_4403.txt", sep=" ")
psfm_5701 = pd.read_csv("/Users/adams/Projects/300K/MSV000091456-SCP/PSFM_5701.txt", sep=" ")

# Process data
inf_0101 = transform_psfm(psfm_0101)
inf_0202 = transform_psfm(psfm_0202)
inf_4403 = transform_psfm(psfm_4403)
inf_5701 = transform_psfm(psfm_5701)

inf_g1 = transform_psfm(psfm_g1)
inf_g2 = transform_psfm(psfm_g2)
inf_g3 = transform_psfm(psfm_g3)
inf_g4 = transform_psfm(psfm_g4)
inf_s1 = transform_psfm(psfm_s1)
inf_s2 = transform_psfm(psfm_s2)
inf_s3 = transform_psfm(psfm_s3)
inf_s4 = transform_psfm(psfm_s4)
inf_l1 = transform_psfm(psfm_l1)
inf_l2 = transform_psfm(psfm_l2)
inf_l3 = transform_psfm(psfm_l3)
inf_l4 = transform_psfm(psfm_l4)

# Figure out which plot where

def kl(p, q):
    """Kullback-Leibler divergence D(P || Q) for discrete distributions
    Parameters
    ----------
    p, q : array-like, dtype=float, shape=n
    Discrete probability distributions.
    """
    p = np.asarray(p, dtype=np.float)
    q = np.asarray(q, dtype=np.float)
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def get_hla_match(psfm_list, psfm_hla):
    kl_dict = {}
    kl_dict["1"]=kl(psfm_list[0], psfm_hla)
    kl_dict["2"]=kl(psfm_list[1], psfm_hla)
    kl_dict["3"]=kl(psfm_list[2], psfm_hla)
    kl_dict["4"]=kl(psfm_list[3], psfm_hla)
    return min(kl_dict.items(), key=lambda x:x[1])

get_hla_match([psfm_g1, psfm_g2, psfm_g3, psfm_g4], psfm_0101)
get_hla_match([psfm_g1, psfm_g2, psfm_g3, psfm_g4], psfm_0202)
get_hla_match([psfm_g1, psfm_g2, psfm_g3, psfm_g4], psfm_4403)
get_hla_match([psfm_g1, psfm_g2, psfm_g3, psfm_g4], psfm_5701)

get_hla_match([psfm_s1, psfm_s2, psfm_s3, psfm_s4], psfm_0101)
get_hla_match([psfm_s1, psfm_s2, psfm_s3, psfm_s4], psfm_0202)
get_hla_match([psfm_s1, psfm_s2, psfm_s3, psfm_s4], psfm_4403)
get_hla_match([psfm_s1, psfm_s2, psfm_s3, psfm_s4], psfm_5701)

get_hla_match([psfm_l1, psfm_l2, psfm_l3, psfm_l4], psfm_5701)
get_hla_match([psfm_l1, psfm_l2, psfm_l3, psfm_l4], psfm_0101)
get_hla_match([psfm_l1, psfm_l2, psfm_l3, psfm_l4], psfm_0202)
get_hla_match([psfm_l1, psfm_l2, psfm_l3, psfm_l4], psfm_4403)


# Plot

color_dict =  {
        'A': '#0e1c36',
        'C': '#cdeac0',
        'D': '#f33b16',
        'E': '#f33b16',
        'F': '#0e1c36',
        'G': '#cdeac0',
        'H': '#7a8db3',
        'I': '#0e1c36',
        'K': '#7a8db3',
        'L': '#0e1c36',
        'M': '#0e1c36',
        'N': '#b8b3e9',
        'P': '#0e1c36',
        'Q': '#b8b3e9',
        'R': '#7a8db3',
        'S': '#cdeac0',
        'T': '#cdeac0',
        'V': '#0e1c36',
        'W': '#0e1c36',
        'Y': '#cdeac0'
    }

# sns.set_context("paper")
sns.set_context("talk")
sns.set_style("whitegrid", {'axes.grid' : False})

cm = 1/2.54  # centimeters in inches
width = 11*cm
height = 8*cm
fig = plt.figure(constrained_layout=True, figsize=(width, height))
axes = fig.subplot_mosaic([['hla']])

logomaker.Logo(inf_l1,
            color_scheme=color_dict,
            center_values=True,
            vpad=.1,
            width=.8,
            ax = axes['hla'])

# axes['hla'].set_title("HLA-A*01:01", fontweight = "bold", y=1.02, pad=-1)
axes['hla'].set_title("Lost", fontweight = "bold", y=1.02, pad=-1)
axes['hla'].set(xticks=np.arange(1, 10))
axes['hla'].set_ylabel("Bits", labelpad=-0.5)
axes['hla'].tick_params(axis='both', which='major', pad=-1)
axes['hla'].set_xlabel("Position")

plot_path = "/Users/adams/Projects/300K/Results/Figures/ppt-logos-l1.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")


cm = 1/2.54  # centimeters in inches
# width = 18*cm
# height = 12*cm
width = 38*cm
height = 24*cm

fig = plt.figure(constrained_layout=True, figsize=(width, height))
axes = fig.subplot_mosaic([
    ['a1', 'b1', 'c1', 'd1'],
    ['a2', 'b2', 'c2', 'd2'],
    ['a3', 'b3', 'c3', 'd3'],
    ['a4', 'b4', 'c4', 'd4']])

motif_axes = [axes['a1'], axes['a2'], axes['a3'], axes['a4'],
              axes['b1'], axes['b2'], axes['b3'], axes['b4'],
              axes['c1'], axes['c2'], axes['c3'], axes['c4'],
              axes['d1'], axes['d2'], axes['d3'], axes['d4']]

matrix_list = [inf_0101, inf_0202, inf_4403, inf_5701,
    inf_g1, inf_g4, inf_g2, inf_g3,
    inf_s1, inf_s4, inf_s2, inf_s3,
    inf_l3, inf_l4, inf_l2, inf_l1]

for ax_logo, matrix in zip(motif_axes, matrix_list):
    logomaker.Logo(matrix,
                color_scheme=color_dict,
                center_values=True,
                vpad=.1,
                width=.8,
                ax = ax_logo)

for ax in motif_axes:
    ax.set(xticks=np.arange(1, 10))
    ax.set_ylabel("Bits", labelpad=-0.5)
    ax.tick_params(axis='both', which='major', pad=-1)

axes['a1'].set_title("HLA-A*01:01", fontweight = "bold", y=1.02, pad=-1)
axes['a2'].set_title("HLA-B*44:03", fontweight = "bold", y=1.02, pad=-1)
axes['a3'].set_title("HLA-B*57:01", fontweight = "bold", y=1.02, pad=-1)
axes['a4'].set_title("HLA-A*02:02", fontweight = "bold", y=1.02, pad=-1)
axes['b1'].set_title("Gained", fontweight = "bold", y=1.02, pad=-1)
axes['c1'].set_title("Shared", fontweight = "bold", y=1.02, pad=-1)
axes['d1'].set_title("Lost", fontweight = "bold", y=1.02, pad=-1)
axes['a4'].set_xlabel("Position")
axes['b4'].set_xlabel("Position")
axes['c4'].set_xlabel("Position")
axes['d4'].set_xlabel("Position")

list_color = ["#f33b16", "#7a8db3", "#0e1c36", "#b8b3e9", "#cdeac0"]
list_lab = ['Acidic','Basic','Hydrophobic', 'Neutral', 'Polar']

handles, labels = axes["a1"].get_legend_handles_labels()
for col, lab in zip(list_color, list_lab):
    patch = mpatches.Patch(color=col, label=lab)
    handles.append(patch) 

# plot the legend
fig.legend(handles=handles, loc='upper center', ncol=5,
                handlelength=1, handleheight=1, handletextpad=0.4,
                bbox_to_anchor=(0.5, 1.04),
                # fontsize="7",
                frameon=False, columnspacing = 0.7)

sns.despine()
fig.tight_layout()

# plt.xticks(rotation = 45) # Rotates X-Axis Ticks by 45-degrees
plot_path = "/Users/adams/Projects/300K/Results/Figures/ppt-logos.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
