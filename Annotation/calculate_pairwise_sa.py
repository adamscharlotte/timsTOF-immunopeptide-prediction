# ssh cadams@10.152.135.57
# /home/cadams/anaconda3/envs/prosit-annotate/bin/python3
# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

from math import comb
import pandas as pd
import numpy as np
import os
from itertools import combinations

# from sqlalchemy import true

# from fundamentals import constants
# from fundamentals.fragments import initialize_peaks
# from fundamentals.annotation.annotation import annotate_spectra
# from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
# import h5py
import argparse, pathlib
# from fundamentals.mod_string import parse_modstrings, maxquant_to_internal
# from fundamentals.constants import ALPHABET
# import prosit_grpc
# from prosit_grpc.predictPROSIT import PROSITpredictor

def get_spectral_angle(intensities):
    pred= np.array(intensities[0])
    true = np.array(intensities[1])
    epsilon = 1e-7
    list_1 = np.argwhere(true>0)
    list_2 = np.argwhere(pred>0)
    indices = np.union1d(list_1,list_1)
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
calibrated_unsummed_path = base_path + "full-truncated-qc/calibrated-20ppm/" + pool + ".csv"
calibrated_sum_path = base_path + "precursor-consensus/calibrated-20ppm/" + pool + ".csv"
sa_path = base_path + "spectral-angle/"

unsummed_df = pd.read_csv(calibrated_unsummed_path)
sum_df = pd.read_csv(calibrated_sum_path)

unsummed_df.INTENSITIES = unsummed_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
unsummed_df.MZ = unsummed_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
unsummed_df.SEQUENCE_INT = unsummed_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

sum_df.INTENSITIES = sum_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
sum_df.MZ = sum_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
sum_df.SEQUENCE_INT = sum_df.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

# ----------------------------------------- UN-SUMMED -------------------------------------------

col_filter = ['SCORE']
unsummed_df[col_filter] = unsummed_df[unsummed_df[col_filter] >= 70][col_filter]
filtered_unsummed_df = unsummed_df.dropna()

filtered_unsummed_df['PRECURSOR_CHARGE'].value_counts()
filtered_unsummed_df['ORIG_COLLISION_ENERGY'].min()
filtered_unsummed_df['aligned_collision_energy'].min()
filtered_unsummed_df.columns
grouped_charge_df = filtered_unsummed_df.groupby('PRECURSOR_CHARGE')

filtered_unsummed_df[['PRECURSOR', 'aligned_collision_energy', 'ORIG_COLLISION_ENERGY']].sort_values(by = ['aligned_collision_energy'])

for charge, df_charge in grouped_charge_df:
    random_10_precursors = df_charge['PRECURSOR'].drop_duplicates().sample(n = 10, random_state = 43)
    df_random_10 = df_charge[df_charge['PRECURSOR'].isin(random_10_precursors)]
    grouped_df_random_10 = df_random_10.groupby('PRECURSOR')
    df_random_10_pairs = pd.DataFrame({'combination' : [], 'A' : '', 'B': '', 'INTENSITIES_A' : [], 'INTENSITIES_B' : [], 'PRECURSOR': ''})
    for group, df_group in grouped_df_random_10:
        if len(df_group) == 1:
            continue
        else:
            comb_list = list(combinations(df_group.FRAME, 2))
            combdf = pd.DataFrame({'combination':comb_list})
            combdf['combination'] = combdf['combination'].astype(str)
            combdf.combination = combdf.combination.str.strip(")(").str.split(", ").apply(lambda s: [int(x) for x in s])
            combdf[['A','B']] = pd.DataFrame(combdf.combination.tolist(), index= combdf.index)
            mergdfA = pd.merge(combdf,df_random_10[['FRAME','INTENSITIES']],left_on='A', right_on='FRAME', how='left')
            mergdfA = mergdfA.drop(columns=['FRAME'])
            mergdfA = mergdfA.rename(columns={'INTENSITIES': 'INTENSITIES_A'})
            mergdfAB = pd.merge(mergdfA,df_random_10[['FRAME','INTENSITIES', 'PRECURSOR']],left_on='B', right_on='FRAME', how='left')
            mergdfAB = mergdfAB.drop(columns=['FRAME'])
            mergdfAB = mergdfAB.rename(columns={'INTENSITIES': 'INTENSITIES_B'})
            df_random_10_pairs = pd.concat([df_random_10_pairs, mergdfAB], ignore_index=True)
    df_random_10_pairs["SPECTRAL_ANGLE"] = df_random_10_pairs[['INTENSITIES_A','INTENSITIES_B']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10_pairs["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    df_random_10_pairs.to_csv(sa_path + "pairwise/frames/20-ppm/" + pool + '_' + str(charge) + '.csv')

df_random_10_pairs["SPECTRAL_ANGLE"] = df_random_10_pairs[['INTENSITIES_A','INTENSITIES_B']].apply(lambda x : get_spectral_angle(x), axis=1)

# ------------------------------------------- SUMMED --------------------------------------------

col_filter = ['SCORE']
sum_df[col_filter] = sum_df[sum_df[col_filter] >= 70][col_filter]
filtered_sum_df = sum_df.dropna()

filtered_sum_df.columns

filtered_sum_df['PRECURSOR_CHARGE'].value_counts()
filtered_sum_df['SCORE'].unique()
filtered_sum_df['ORIG_COLLISION_ENERGY'].min()
filtered_sum_df['aligned_collision_energy'].min()

grouped_charge_df = filtered_sum_df.groupby('PRECURSOR_CHARGE')

for charge, df_charge in grouped_charge_df:
    random_10_precursors = df_charge['SEQUENCE'].drop_duplicates().sample(n = 10, random_state = 43)
    df_random_10 = df_charge[df_charge['SEQUENCE'].isin(random_10_precursors)]
    grouped_df_random_10 = df_random_10.groupby('SEQUENCE')
    df_random_10_pairs = pd.DataFrame({'combination' : [], 'A' : '', 'B': '', 'INTENSITIES_A' : [], 'INTENSITIES_B' : [], 'SEQUENCE': ''})
    for group, df_group in grouped_df_random_10:
        if len(df_group) == 1:
            continue
        else:
            comb_list = list(combinations(df_group.PRECURSOR, 2))
            combdf = pd.DataFrame({'combination':comb_list})
            combdf['combination'] = combdf['combination'].astype(str)
            combdf.combination = combdf.combination.str.strip(")(").str.split(", ").apply(lambda s: [int(x) for x in s])
            combdf[['A','B']] = pd.DataFrame(combdf.combination.tolist(), index= combdf.index)
            mergdfA = pd.merge(combdf,df_random_10[['PRECURSOR','INTENSITIES']],left_on='A', right_on='PRECURSOR', how='left')
            mergdfA = mergdfA.drop(columns=['PRECURSOR'])
            mergdfA = mergdfA.rename(columns={'INTENSITIES': 'INTENSITIES_A'})
            mergdfAB = pd.merge(mergdfA,df_random_10[['PRECURSOR','INTENSITIES', 'SEQUENCE']],left_on='B', right_on='PRECURSOR', how='left')
            mergdfAB = mergdfAB.drop(columns=['PRECURSOR'])
            mergdfAB = mergdfAB.rename(columns={'INTENSITIES': 'INTENSITIES_B'})
            df_random_10_pairs = pd.concat([df_random_10_pairs, mergdfAB], ignore_index=True)
    df_random_10_pairs["SPECTRAL_ANGLE"] = df_random_10_pairs[['INTENSITIES_A','INTENSITIES_B']].apply(lambda x : get_spectral_angle(x), axis=1)
    df_random_10_pairs["SPECTRAL_ANGLE"].fillna(0, inplace=True)
    df_random_10_pairs.to_csv(sa_path + "pairwise/precursors/20-ppm/" + pool + '_' + str(charge) + '.csv')


# ------------------------------------- FIGURE OUT HOW TO --------------------------------------

# df = pd.DataFrame({'year': ['2001', '2004', '2005', '2006', '2007', '2008', '2009',
#                              '2003', '2004', '2005', '2006', '2007', '2008', '2009',
#                             '2003', '2004', '2005', '2006', '2007', '2008', '2009'],
#                    'id': ['1', '1', '1', '1', '1', '1', '1', 
#                           '2', '2', '2', '2', '2', '2', '2',
#                          '5', '5', '5','5', '5', '5', '5'],
#                    'money': ['15', '15', '15', '21', '21', '21', '21', 
#                              '17', '17', '17', '20', '17', '17', '17',
#                             '25', '30', '22', '25', '8', '7', '12']}).astype(int)

# grouped_df = df.groupby('id')
# for id, df_id in grouped_df:
#     comb_list = list(combinations(df_id.year, 2))
#     combdf = pd.DataFrame({'combination':comb_list})
#     combdf['combination'] = combdf['combination'].astype(str)
#     combdf.combination = combdf.combination.str.strip(")(").str.split(", ").apply(lambda s: [int(x) for x in s])
#     combdf[['A','B']] = pd.DataFrame(combdf.combination.tolist(), index= combdf.index)
#     df_id
#     mergdfA = pd.merge(combdf,df_id[['year','money']],left_on='A', right_on='year', how='left')
#     mergdfA = mergdfA.drop(columns=['year'])
#     mergdfA = mergdfA.rename(columns={'money': 'money_A'})
#     mergdfAB = pd.merge(mergdfA,df_id[['year','money']],left_on='B', right_on='year', how='left')
#     mergdfAB = mergdfAB.drop(columns=['year'])
#     mergdfAB = mergdfAB.rename(columns={'money': 'money_B'})
