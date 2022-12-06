# ssh vsc20709@login-leibniz.hpc.uantwerpen.be

import sqlite3
import pandas as pd

# import argparse, pathlib
# parser = argparse.ArgumentParser()
# parser.add_argument("folder_path", type=str)
# parser.add_argument("mgf_path", type=str)
# args = parser.parse_args()
# folder_path = args.folder_path
# mgf_path = args.mgf_path
folder_path = "/Users/adams/Projects/300K/2022-library-run/raw-folders/HLAI_1_96_p2-A3_S2-A3_1_6928.d"
mgf_path = "/Users/adams/Projects/300K/2022-library-run/precursor-mgf/HLAI_1_96_p2-A3_S2-A3_1_6928.mgf"

# Get the pasef information
tdf_path = folder_path + "/analysis.tdf"
con = sqlite3.connect(tdf_path)
cur = con.cursor()
pasef_df = pd.read_sql_query("SELECT * from PasefFrameMsMsInfo", con)
con.close()

# Get all frames
frames = pasef_df.Frame.to_numpy()

