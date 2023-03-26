# /home/cadams/anaconda3/envs/prosit-annotate/bin/python3

import sqlite3
import pandas as pd
import argparse, pathlib
parser = argparse.ArgumentParser()
parser.add_argument("pool", type=str)
args = parser.parse_args()
pool = args.pool

# pool = "TUM_first_pool_124"
# raw_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/Analysis/" + pool + "_01_01/TIMS-30min-R1-TIMS_unspecific/230719_f1-" + pool + "_01_01-TIMS-30min-R1.d"
# raw_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/Analysis/" + pool + "_01_01/TIMS-30min-R1-TIMS_semitryptic/230719_f1-" + pool + "_01_01-TIMS-30min-R1.d"
raw_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/d-folder/" + pool + ".d"
tdf_path = raw_path + "/analysis.tdf"
con = sqlite3.connect(tdf_path)
cur = con.cursor()
pasef_df = pd.read_sql_query("SELECT * from PasefFrameMsMsInfo", con)
con.close()

# output_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/extract-pasef/" + pool + ".csv"
output_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/extract-pasef/" + pool + ".csv"
pasef_df.to_csv(output_path)