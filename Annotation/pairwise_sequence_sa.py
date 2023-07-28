# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3

from math import comb
import pandas as pd
import numpy as np
import os
from itertools import combinations

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

base = "/Users/adams/Projects/300K/2022-library-run/"
sa_path = base + "Annotation/spectral-angle/"

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

df_test["type"] = df_test["pool_name"].apply(lambda x: x[4:8])

df_test.INTENSITIES = df_test.INTENSITIES.str.split(";").apply(lambda s: [float(x) for x in s])
df_test.MZ = df_test.MZ.str.split(";").apply(lambda s: [float(x) for x in s])
df_test.SEQUENCE_INT = df_test.SEQUENCE_INT.str.strip("][").str.split(", ").apply(lambda s: [int(x) for x in s])

df_test['PRECURSOR_CHARGE'].value_counts()
df_test['SCORE'].unique()
# filtered_sum_df['COLLISION_ENERGY'].min()
df_test['median_CE'].min()

df_test_filtered = df_test.loc[~((df_test['type'] == 'lysn') & (df_test['PRECURSOR_CHARGE'] == 1))]

grouped_type_df = df_test_filtered.groupby('type')

for types, df_type in grouped_type_df:
    # grouped_charge_df = df_type.groupby('PRECURSOR_CHARGE')
    # for charge, df_charge in grouped_charge_df:
    #     counts = df_charge['OBS_SEQUENCE'].value_counts()
    #     counting = counts > 1
    #     print(types)
    #     print(counting.value_counts())
    grouped_charge_df = df_type.groupby('PRECURSOR_CHARGE')
    for charge, df_charge in grouped_charge_df:
        df_grouped = df_charge.groupby('OBS_SEQUENCE')
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
        df_pairs.to_csv(sa_path + "pairwise/test-set/test_" + str(types) + '_' + str(charge) + '.csv')
