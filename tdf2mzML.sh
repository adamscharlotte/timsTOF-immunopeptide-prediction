conda create -n tdf2mzml python=3.8
conda activate tdf2mzml

# Changed the requirements.txt file

docker build -t tdf2mzml .

PWD=/Users/adams/Projects/300K/20220613-library-test
docker run --rm -it -v $PWD:/data tdf2mzml tdf2mzml.py -i /data/E1_50fmol_bruker_S1-E2_1_6339.d 