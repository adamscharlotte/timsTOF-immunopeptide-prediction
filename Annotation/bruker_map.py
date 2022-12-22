# /home/cadams/anaconda3/envs/prosit-annotate/bin/python3

import sqlite3
import pandas as pd
import argparse, pathlib
parser = argparse.ArgumentParser()
parser.add_argument("pool", type=str)
args = parser.parse_args()
pool = args.pool

pool = "TUM_first_pool_126"
raw_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/Analysis/" + pool + "_01_01/TIMS-30min-R1-TIMS_unspecific/230719_f1-" + pool + "_01_01-TIMS-30min-R1.d"
# raw_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/Analysis/" + pool + "_01_01/TIMS-30min-R1-TIMS_semitryptic/230719_f1-" + pool + "_01_01-TIMS-30min-R1.d"
tdf_path = raw_path + "/analysis.tdf"
con = sqlite3.connect(tdf_path)
cur = con.cursor()
pasef_df = pd.read_sql_query("SELECT * from PasefFrameMsMsInfo", con)
con.close()

min_frame = pasef_df.Frame.min()
max_frame = pasef_df.Frame.max()

f = open("/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/pool-path-frames.txt", "a")
f.write(f"{pool};{raw_path};{min_frame};{max_frame}" + '\n')
f.close()