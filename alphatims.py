# pip install alphatims
import alphatims
import os, sys
import bruker as timsbr

cwd = os.getcwd()
print(cwd)

alphatims.bruker.open_bruker_d_folder("/Users/adams/Projects/300K/2022-library-run/MaxQuant/LysN/AspNLysN-A5_S1-A5_1_6546.d")

