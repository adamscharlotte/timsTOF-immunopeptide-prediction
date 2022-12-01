# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

from math import comb
import pandas as pd
import numpy as np
import os
from itertools import combinations

import argparse, pathlib

def get_spectral_angle(intensities):
    pred= np.array(intensities[0])
    true = np.array(intensities[1])
    epsilon = 1e-7
    list_1 = np.argwhere(true>0)
    list_2 = np.argwhere(pred>0)
    indices = np.union1d(list_1,list_2)
    pred_masked = pred[indices]
    true_masked = true[indices]
    true_masked += epsilon
    pred_masked += epsilon
    true_norm = true_masked*(1/np.sqrt(np.sum(np.square(true_masked), axis=0)))
    pred_norm = pred_masked*(1/np.sqrt(np.sum(np.square(pred_masked), axis=0)))
    product = np.sum(true_norm*pred_norm, axis=0)
    arccos = np.arccos(product)
    return 1-2*arccos/np.pi

parser = argparse.ArgumentParser()
parser.add_argument('pool', type=str)					# Filename
args = parser.parse_args()

pool = args.pool
# pool = "TUM_aspn_36"

base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
# unsummed_path = base_path + "full-truncated-qc/annotated-40-ppm/" + pool + ".csv"
sum_path = base_path + "scan-consensus/annotated-20-ppm/" + pool + ".csv"
sa_path = base_path + "spectral-angle/"

# unsummed_df = pd.read_csv(unsummed_path)
sum_df = pd.read_csv(sum_path)

# unsummed_df.INTENSITIES = unsummed_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
# unsummed_df.MZ = unsummed_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
# unsummed_df.SEQUENCE_INT = unsummed_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

sum_df.INTENSITIES = sum_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
sum_df.MZ = sum_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
sum_df.SEQUENCE_INT = sum_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

# # ----------------------------------------- UN-SUMMED -------------------------------------------

# col_filter = ['SCORE']
# unsummed_df[col_filter] = unsummed_df[unsummed_df[col_filter] >= 70][col_filter]
# filtered_unsummed_df = unsummed_df.dropna()

# filtered_unsummed_df['PRECURSOR_CHARGE'].value_counts()
# filtered_unsummed_df['COLLISION_ENERGY'].min()
# filtered_unsummed_df.columns
# grouped_charge_df = filtered_unsummed_df.groupby('PRECURSOR_CHARGE')

# for charge, df_charge in grouped_charge_df:
#     df_grouped = df_charge.groupby('PRECURSOR')
#     df_pairs = pd.DataFrame({'combination' : [], 'A' : '', 'B': '', 'INTENSITIES_A' : [], 'INTENSITIES_B' : [], 'PRECURSOR': ''})
#     for group, df_group in df_grouped:
#         if len(df_group) == 1:
#             continue
#         else:
#             comb_list = list(combinations(df_group.FRAME, 2))
#             combdf = pd.DataFrame({'combination':comb_list})
#             combdf['combination'] = combdf['combination'].astype(str)
#             combdf.combination = combdf.combination.str.strip(")(").str.split(", ").apply(lambda s: [int(x) for x in s])
#             combdf[['A','B']] = pd.DataFrame(combdf.combination.tolist(), index= combdf.index)
#             mergdfA = pd.merge(combdf,df_group[['FRAME','INTENSITIES']],left_on='A', right_on='FRAME', how='left')
#             mergdfA = mergdfA.drop(columns=['FRAME'])
#             mergdfA = mergdfA.rename(columns={'INTENSITIES': 'INTENSITIES_A'})
#             mergdfAB = pd.merge(mergdfA,df_group[['FRAME','INTENSITIES', 'PRECURSOR']],left_on='B', right_on='FRAME', how='left')
#             mergdfAB = mergdfAB.drop(columns=['FRAME'])
#             mergdfAB = mergdfAB.rename(columns={'INTENSITIES': 'INTENSITIES_B'})
#             df_pairs = pd.concat([df_pairs, mergdfAB], ignore_index=True)
#     df_pairs["SPECTRAL_ANGLE"] = df_pairs[['INTENSITIES_A','INTENSITIES_B']].apply(lambda x : get_spectral_angle(x), axis=1)
#     df_pairs["SPECTRAL_ANGLE"].fillna(0, inplace=True)
#     df_pairs.to_csv(sa_path + "pairwise/frames/all-40-ppm/" + pool + '_' + str(charge) + '.csv')

# ------------------------------------------- SUMMED --------------------------------------------

# col_filter = ['SCORE']
# sum_df[col_filter] = sum_df[sum_df[col_filter] >= 70][col_filter]
# filtered_sum_df = sum_df.dropna()

# filtered_sum_df.columns

# filtered_sum_df['PRECURSOR_CHARGE'].value_counts()
# filtered_sum_df['SCORE'].unique()
# filtered_sum_df['COLLISION_ENERGY'].min()

# grouped_charge_df = filtered_sum_df.groupby('PRECURSOR_CHARGE')

# for charge, df_charge in grouped_charge_df:
#     df_grouped = df_charge.groupby('SEQUENCE')
#     df_pairs = pd.DataFrame({'combination' : [], 'A' : '', 'B': '', 'INTENSITIES_A' : [], 'INTENSITIES_B' : [], 'SEQUENCE': ''})
#     for group, df_group in df_grouped:
#         if len(df_group) == 1:
#             continue
#         else:
#             comb_list = list(combinations(df_group.PRECURSOR, 2))
#             combdf = pd.DataFrame({'combination':comb_list})
#             combdf['combination'] = combdf['combination'].astype(str)
#             combdf.combination = combdf.combination.str.strip(")(").str.split(", ").apply(lambda s: [int(x) for x in s])
#             combdf[['A','B']] = pd.DataFrame(combdf.combination.tolist(), index= combdf.index)
#             mergdfA = pd.merge(combdf,df_charge[['PRECURSOR','INTENSITIES']],left_on='A', right_on='PRECURSOR', how='left')
#             mergdfA = mergdfA.drop(columns=['PRECURSOR'])
#             mergdfA = mergdfA.rename(columns={'INTENSITIES': 'INTENSITIES_A'})
#             mergdfAB = pd.merge(mergdfA,df_charge[['PRECURSOR','INTENSITIES', 'SEQUENCE']],left_on='B', right_on='PRECURSOR', how='left')
#             mergdfAB = mergdfAB.drop(columns=['PRECURSOR'])
#             mergdfAB = mergdfAB.rename(columns={'INTENSITIES': 'INTENSITIES_B'})
#             df_pairs = pd.concat([df_pairs, mergdfAB], ignore_index=True)
#     df_pairs["SPECTRAL_ANGLE"] = df_pairs[['INTENSITIES_A','INTENSITIES_B']].apply(lambda x : get_spectral_angle(x), axis=1)
#     df_pairs["SPECTRAL_ANGLE"].fillna(0, inplace=True)
#     df_pairs.to_csv(sa_path + "pairwise/scans/20-ppm/" + pool + '_' + str(charge) + '.csv')

# ------------------------------------------- SCANS --------------------------------------------

col_filter = ['SCORE']
sum_df[col_filter] = sum_df[sum_df[col_filter] >= 70][col_filter]
filtered_sum_df = sum_df.dropna()

filtered_sum_df.columns

filtered_sum_df['PRECURSOR_CHARGE'].value_counts()
filtered_sum_df['SCORE'].unique()
# filtered_sum_df['COLLISION_ENERGY'].min()
filtered_sum_df['mean_CE'].min()

grouped_charge_df = filtered_sum_df.groupby('PRECURSOR_CHARGE')

for charge, df_charge in grouped_charge_df:
    df_grouped = df_charge.groupby('SEQUENCE')
    df_pairs = pd.DataFrame({'combination' : [], 'A' : '', 'B': '', 'INTENSITIES_A' : [], 'INTENSITIES_B' : [], 'SEQUENCE': ''})
    for group, df_group in df_grouped:
        if len(df_group) == 1:
            continue
        else:
            comb_list = list(combinations(df_group.SCAN_NUMBER, 2))
            combdf = pd.DataFrame({'combination':comb_list})
            combdf['combination'] = combdf['combination'].astype(str)
            combdf.combination = combdf.combination.str.strip(")(").str.split(", ").apply(lambda s: [int(x) for x in s])
            combdf[['A','B']] = pd.DataFrame(combdf.combination.tolist(), index= combdf.index)
            mergdfA = pd.merge(combdf,df_charge[['SCAN_NUMBER','INTENSITIES']],left_on='A', right_on='SCAN_NUMBER', how='left')
            mergdfA = mergdfA.drop(columns=['SCAN_NUMBER'])
            mergdfA = mergdfA.rename(columns={'INTENSITIES': 'INTENSITIES_A'})
            mergdfAB = pd.merge(mergdfA,df_charge[['SCAN_NUMBER','INTENSITIES', 'SEQUENCE']],left_on='B', right_on='SCAN_NUMBER', how='left')
            mergdfAB = mergdfAB.drop(columns=['SCAN_NUMBER'])
            mergdfAB = mergdfAB.rename(columns={'INTENSITIES': 'INTENSITIES_B'})
            df_pairs = pd.concat([df_pairs, mergdfAB], ignore_index=True)
    df_pairs["SPECTRAL_ANGLE"] = df_pairs[['INTENSITIES_A','INTENSITIES_B']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_pairs["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    df_pairs.to_csv(sa_path + "pairwise/scans/20-ppm/" + pool + '_' + str(charge) + '.csv')
