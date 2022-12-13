# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

import numpy as np
import pandas as pd
import os
# from itertools import combinations

def get_spectral_angle(intensities):
    pred= np.array(intensities[0])
    true = np.array(intensities[1])
    epsilon = 1e-7
    list_1 = np.argwhere(true>0)
    list_2 = np.argwhere(pred>0)
    indices = np.intersect1d(list_1,list_2)
    pred_masked = pred[indices]
    true_masked = true[indices]
    true_masked += epsilon
    pred_masked += epsilon
    true_norm = true_masked*(1/np.sqrt(np.sum(np.square(true_masked), axis=0)))
    pred_norm = pred_masked*(1/np.sqrt(np.sum(np.square(pred_masked), axis=0)))
    product = np.sum(true_norm*pred_norm, axis=0)
    arccos = np.arccos(product)
    return 1-2*arccos/np.pi

pool = "TUM_HLA2_7"
base_path = "/Users/adams/Projects/300K/2022-library-run/Annotation/"
# base_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/" # nolint
calibrated_path = base_path + "scan-consensus/calibrated-40-ppm/" + pool + ".csv"
map_path = base_path + "full-length-map/" + pool + ".csv"

map_df = pd.read_csv(map_path)
map_df = map_df.rename({"OBS_SEQUENCE": "OBS_SEQUENCE_map"}, axis='columns')

calibrated_annot_df = pd.read_csv(calibrated_path)
calibrated_annot_df = calibrated_annot_df.rename({"SEQUENCE": "OBS_SEQUENCE", "MODIFIED_SEQUENCE": "OBS_MODIFIED_SEQUENCE", "SEQUENCE_INT": "OBS_SEQUENCE_INT"}, axis='columns')
calibrated_annot_df.INTENSITIES = calibrated_annot_df.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
calibrated_annot_df.MZ = calibrated_annot_df.MZ.str.split(";").apply(lambda s: [float(x) for x in s])

calibrated_annot_df_merged = pd.merge(left=calibrated_annot_df, right=map_df, how='left', left_on='SCAN_NUMBER', right_on='SCAN_NUMBER')

calibrated_annot_df_merged.columns

len(calibrated_annot_df_merged[calibrated_annot_df_merged["OBS_SEQUENCE"] == calibrated_annot_df_merged["OBS_SEQUENCE_map"]])
len(calibrated_annot_df_merged.dropna())

df_merged = calibrated_annot_df_merged.dropna()

df_merged.groupby('SEQUENCE')['OBS_SEQUENCE'].transform('nunique')
df_merged['count'] = df_merged.groupby('SEQUENCE')['SEQUENCE'].transform('count')
df_merged.sort_values('count', ascending=False)
df_merged[df_merged['count'] > 3]

df_filter = df_merged[df_merged["SEQUENCE"].isin(["HTAYSDFLSDK", "GDGTFQKWAAVVVPSGEEQ"])]
# df_filter = df_merged[df_merged["SEQUENCE"] == "HTAYSDFLSDK"]
df_filter[['SEQUENCE','OBS_SEQUENCE', 'PRECURSOR_CHARGE']]

df_filter.columns

grouped_charge_df = df_filter[['SEQUENCE','OBS_SEQUENCE','PRECURSOR_CHARGE','MZ','INTENSITIES']].groupby(['PRECURSOR_CHARGE','SEQUENCE'])

for charge, df_charge in grouped_charge_df:
    df_fullength = df_charge[df_charge["SEQUENCE"] == df_charge["OBS_SEQUENCE"]]
    df_trunc = df_charge[(df_charge["SEQUENCE"].isin(df_fullength["SEQUENCE"])) &
        (df_charge["SEQUENCE"] != df_charge["OBS_SEQUENCE"])].add_prefix('trunc_')
    df_merge = pd.merge(df_trunc,df_fullength[['MZ','INTENSITIES','SEQUENCE']],left_on='trunc_SEQUENCE', right_on='SEQUENCE', how='left')
    df_merge["SPECTRAL_ANGLE"] = df_merge[['INTENSITIES','trunc_INTENSITIES']].apply(lambda x : get_spectral_angle(x), axis=1)


    df_random_10["SPECTRAL_ANGLE_CAL"].fillna(0, inplace=True)
    df_random_10.to_csv(sa_path + "frames/20-ppm/" + pool + '_' + str(charge) + '.csv')



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
