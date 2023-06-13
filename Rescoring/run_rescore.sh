ssh cadams@10.152.135.57

scp Names/d-hlai-rescore.txt cadams@10.152.135.57:/home/cadams
scp Names/scp-rescore.txt cadams@10.152.135.57:/home/cadams
scp Annotation/d_master_spec.py cadams@10.152.135.57:/home/cadams

scp -r oktoberfest-0.1.0-tims cadams@10.152.135.57:/home/cadams
cd /Users/adams/Code/oktoberfest-0.1.0-tims/oktoberfest-0.1.0-tims-1/
scp -r oktoberfest cadams@10.152.135.57:/home/cadams/oktoberfest-0.1.0-tims/
scp oktoberfest/ce_calibration.py cadams@10.152.135.57:/home/cadams/oktoberfest-0.1.0-tims/oktoberfest
scp -r cadams@10.152.135.57:/home/cadams/spectrum_io .

# conda create --name oktoberfest-tims python=3.8
# conda activate oktoberfest-tims
conda activate oktoberfest-0_1
cd /home/cadams/oktoberfest-0.1.0-tims
pip install .


screen -r rescore-hcd
conda activate oktoberfest-0_1

name_file=d-hlai-rescore.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison
    mkdir $searchdir/reresults/d-hcd/${array[0]}
    scp $searchdir/reresults/d-tims/${array[0]}/${array[0]}.pkl $searchdir/reresults/d-hcd/${array[0]}
    scp $searchdir/reresults/d-tims/${array[0]}/msms.txt $searchdir/reresults/d-hcd/${array[0]}
    scp $searchdir/reresults/d-hcd/test.json $searchdir/reresults/d-hcd/${array[0]}
    time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/d-hcd/${array[0]} --config_path $searchdir/reresults/d-hcd/${array[0]}/test.json
done

name_file=d-hlai-rescore.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison
    mkdir $searchdir/reresults/d-cid/${array[0]}
    scp $searchdir/reresults/d-tims/${array[0]}/${array[0]}.pkl $searchdir/reresults/d-cid/${array[0]}
    scp $searchdir/reresults/d-tims/${array[0]}/msms.txt $searchdir/reresults/d-cid/${array[0]}
    scp $searchdir/reresults/d-cid/test.json $searchdir/reresults/d-cid/${array[0]}
    time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/d-cid/${array[0]} --config_path $searchdir/reresults/d-cid/${array[0]}/test.json
done


screen -r rescore-scp

name_file=scp-rescore.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/A375
    mkdir $searchdir/reresults/A375-low-input-HLAI/${array[0]}
    scp $searchdir/reresults/test.json $searchdir/reresults/A375-low-input-HLAI/${array[0]}
    scp $searchdir/Searches/A375-low-input-HLAI/${array[0]}/combined/txt/msms.txt $searchdir/reresults/A375-low-input-HLAI/${array[0]}
    scp $searchdir/Searches/A375-low-input-HLAI/${array[0]}/combined/txt/pasefMsmsScans.txt $searchdir/reresults/A375-low-input-HLAI/${array[0]}
    scp $searchdir/Searches/A375-low-input-HLAI/${array[0]}/combined/txt/accumulatedMsmsScans.txt $searchdir/reresults/A375-low-input-HLAI/${array[0]}
    time Rscript d_extract.R "$searchdir/A375-low-input-HLAI/${array[0]}.d" "$searchdir/reresults/A375-low-input-HLAI/${array[0]}" "/media/c2m/5_user_files/cadams/libtimsdata.so"
    time python d_master_spec_2.py "$searchdir/reresults/A375-low-input-HLAI/${array[0]}/${array[0]}.csv"
    time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/A375-low-input-HLAI/${array[0]} --config_path $searchdir/reresults/A375-low-input-HLAI/${array[0]}/test.json
done



python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir 220404_NHG_benign_UDN10_Liver_W6-32_17%_DDA_Rep1-test --config_path 220404_NHG_benign_UDN10_Liver_W6-32_17%_DDA_Rep1-test/test.json


csv_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/test-d/211113_SS_malignant_HNSCC_Tue39L243_20%_DDA_Rep1.csv" # nolint
base_path = os.path.dirname(csv_path)
csv_df = pd.read_csv(csv_path)

screen -r rescore-loop-tims

name_file=d-tims-rescore.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison
    mkdir $searchdir/reresults/d-tims/${array[0]}
    scp $searchdir/reresults/d-tims/test.json $searchdir/reresults/d-tims/${array[0]}
    scp $searchdir/reresults/d/${array[0]}/msms.txt $searchdir/reresults/d-tims/${array[0]}
    scp $searchdir/reresults/d/${array[0]}/pasefMsmsScans.txt $searchdir/reresults/d-tims/${array[0]}
    scp $searchdir/reresults/d/${array[0]}/accumulatedMsmsScans.txt $searchdir/reresults/d-tims/${array[0]}
    scp $searchdir/reresults/d/${array[0]}/${array[0]}.pkl $searchdir/reresults/d-tims/${array[0]}
    time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/d-tims/${array[0]} --config_path $searchdir/reresults/d-tims/${array[0]}/test.json
done

name_file=d-tims-rescore.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison
    time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/d-tims/${array[0]} --config_path $searchdir/reresults/d-tims/${array[0]}/test.json
done


searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison
file=211113_SS_malignant_HNSCC_Tue39L243_20%_DDA_Rep1
time python d_master_spec_3.py "$searchdir/reresults/d/$file/$file.csv"


conda activate oktoberfest-0_1

searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/Jurkat-A549
file=IPX72_HLAI_ElA_S3-H2_1_11987
file=IPX_A549_B_S1-A8_1_11150
file=IPX_A549_A_S1-A2_1_11144
mkdir $searchdir/reresults/$file
scp $searchdir/Searches/A549/$file/combined/txt/msms.txt $searchdir/reresults/$file
scp $searchdir/Searches/A549/$file/combined/txt/pasefMsmsScans.txt $searchdir/reresults/$file
scp $searchdir/Searches/A549/$file/combined/txt/accumulatedMsmsScans.txt $searchdir/reresults/$file
scp $searchdir/Searches/A549/$file/combined/txt/accumulatedMsmsScans.txt $searchdir/reresults/$file
scp $searchdir/reresults/test.json $searchdir/reresults/$file

time Rscript d_extract.R "$searchdir/d-folder/$file.d" "$searchdir/reresults/$file" "/media/c2m/5_user_files/cadams/libtimsdata.so"
time python d_master_spec_3.py "$searchdir/reresults/$file/$file.csv"
time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/$file --config_path $searchdir/reresults/$file/test.json

time python d_master_spec_3.py "$searchdir/reresults/$file/$file.csv"


searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/Jurkat-A549
# file=IPX69_HLA1_B_S4-A6_1_10893
# file=IPX_A549_B_S1-A8_1_11150
# file=IPX_A549_A_S1-A2_1_11144
file=IPX67_HLAI_A_S1-A2_1_10830

scp cadams@10.152.135.57:$searchdir/reresults/$file/msms.txt /Users/adams/Projects/300K/Jurkat-A549/reresults/msms/$file.txt
scp cadams@10.152.135.57:$searchdir/reresults/$file/results-tims/percolator/rescore_target.psms /Users/adams/Projects/300K/Jurkat-A549/reresults/rescore-tims/$file.psms
scp cadams@10.152.135.57:$searchdir/reresults/$file/results-hcd/percolator/rescore_target.psms /Users/adams/Projects/300K/Jurkat-A549/reresults/rescore-hcd/$file.psms
scp cadams@10.152.135.57:$searchdir/reresults/$file/results-cid/percolator/rescore_target.psms /Users/adams/Projects/300K/Jurkat-A549/reresults/rescore-cid/$file.psms
scp cadams@10.152.135.57:$searchdir/reresults/$file/results-tims/percolator/original_target.psms /Users/adams/Projects/300K/Jurkat-A549/reresults/rescore-tims/$file-andromeda.psms

220706_NHG_benign_UDN10_Spleen_W6-32_17%_DDA_Rep1 IIv
220404_NHG_malignant_CLL_01_W6-32_17%_DDA_Rep1 IIv
220404_NHG_benign_UDN31_PBMC_W6-32_17%_Rep3 II
220331_NHG_malignant_CLL_02_W6-32_17%_DDA_Rep1 II
220331_NHG_malignant_CLL_02_Tue39L243_17%_DDA_Rep1 I

211113_SS_malignant_HNSCC_W6-32_20%_DDA_Rep1

file="211113_SS_malignant_HNSCC_W6-32_20%_DDA_Rep1"
searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison
mkdir $searchdir/reresults/d-tims/$file
scp $searchdir/reresults/d-tims/test.json $searchdir/reresults/d-tims/$file
scp $searchdir/Searches/d-hla-i/$file/combined/txt/msms.txt $searchdir/reresults/d-tims/$file
scp $searchdir/Searches/d-hla-i/$file/combined/txt/pasefMsmsScans.txt $searchdir/reresults/d-tims/$file
scp $searchdir/Searches/d-hla-i/$file/combined/txt/accumulatedMsmsScans.txt $searchdir/reresults/d-tims/$file
# scp $searchdir/Searches/d-hla-i/$file/$file.pkl $searchdir/reresults/d-tims/$file
time Rscript d_extract.R "$searchdir/d-folder/$file.d" "$searchdir/reresults/d-tims/$file" "/media/c2m/5_user_files/cadams/libtimsdata.so"
time python d_master_spec_3.py "$searchdir/reresults/$file/$file.csv"
time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/d-tims/$file --config_path $searchdir/reresults/d-tims/$file/test.json



# 220706_NHG_benign_UDN12_Thymus_Tue39L243_17%_DDA_Rep3
220404_NHG_benign_UDN10_Liver_Tue39L243_17%_DDA_Rep2 xx
220331_NHG_malignant_CLL_02_Tue39L243_17%_DDA_Rep1 xx
211113_SS_malignant_HNSCC_Tue39L243_20%_DDA_Rep1 xx

file="220404_NHG_benign_UDN10_Liver_Tue39L243_17%_DDA_Rep2"
searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison
mkdir $searchdir/reresults/d-tims-ii/$file
scp $searchdir/reresults/d-tims/test.json $searchdir/reresults/d-tims-ii/$file
scp $searchdir/reresults/d/$file/msms.txt $searchdir/reresults/d-tims-ii/$file
scp $searchdir/reresults/d/$file/pasefMsmsScans.txt $searchdir/reresults/d-tims-ii/$file
scp $searchdir/reresults/d/$file/accumulatedMsmsScans.txt $searchdir/reresults/d-tims-ii/$file
scp $searchdir/reresults/d/$file/$file.pkl $searchdir/reresults/d-tims-ii/$file
time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/d-tims-ii/$file --config_path $searchdir/reresults/d-tims-ii/$file/test.json

# file="220404_NHG_malignant_CLL_01_W6-32_17%_DDA_Rep1"
file="211113_SS_malignant_HNSCC_W6-32_20%_DDA_Rep1"
searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison
scp $searchdir/reresults/d/$file/$file.pkl $searchdir/reresults/d-tims/$file
time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/d-tims/$file --config_path $searchdir/reresults/d-tims/$file/test.json


scp -r /Users/adams/Projects/300K/MSV000091456-SCP/A375_lowInput_IP/A375-low-input-HLAI cadams@10.152.135.57:/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/A375/

rsync -au /Users/adams/Projects/300K/MSV000091456-SCP/A375_lowInput_IP/A375-low-input-HLAI cadams@10.152.135.57:/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/A375/A375-low-input-HLAI



E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep1_Slot1-06_1_3498
E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep3_Slot2-1_1_3508
E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep4_Slot2-1_1_3509
E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep1_Slot2-3_1_3512

E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep2_Slot2-5_1_3520
E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep4_Slot2-5_1_3522

file="E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496"
searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/A375
mkdir $searchdir/reresults/A375-low-input-HLAI/$file
scp $searchdir/reresults/test.json $searchdir/reresults/A375-low-input-HLAI/$file
scp $searchdir/Searches/A375-low-input-HLAI/$file/combined/txt/msms.txt $searchdir/reresults/A375-low-input-HLAI/$file
scp $searchdir/Searches/A375-low-input-HLAI/$file/combined/txt/pasefMsmsScans.txt $searchdir/reresults/A375-low-input-HLAI/$file
scp $searchdir/Searches/A375-low-input-HLAI/$file/combined/txt/accumulatedMsmsScans.txt $searchdir/reresults/A375-low-input-HLAI/$file
# scp $searchdir/Searches/A375-low-input-HLAI/$file/$file.pkl $searchdir/reresults/A375-low-input-HLAI/$file
# time Rscript d_extract.R "$searchdir/A375-low-input-HLAI/$file.d" "$searchdir/reresults/A375-low-input-HLAI/$file" "/media/c2m/5_user_files/cadams/libtimsdata.so"
# time python d_master_spec_3.py "$searchdir/reresults/$file/$file.csv"
time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/A375-low-input-HLAI/$file --config_path $searchdir/reresults/A375-low-input-HLAI/$file/test.json

rm $searchdir/reresults/A375-low-input-HLAI/$file/$file.csv

E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep1_Slot1-06_1_3498
E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep3_Slot2-1_1_3508
E_20221201_NO30_400nL_HLAc1_1e6_directIP_titration_rep4_Slot2-1_1_3509
E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep1_Slot2-3_1_3512
E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep2_Slot2-3_1_3513
E_20221201_NO30_400nL_HLAc1_1e7_directIP_titration_rep3_Slot2-3_1_3514
# E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496
E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep2_Slot2-5_1_3520
E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep3_Slot2-5_1_3521
E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep4_Slot2-5_1_3522

file="E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep1_Slot1-04_1_3496"
searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/A375

scp cadams@10.152.135.57:$searchdir/reresults/A375-low-input-HLAI/$file/msms.txt /Users/adams/Projects/300K/MSV000091456-SCP/reresults/msms/$file.txt
scp cadams@10.152.135.57:$searchdir/reresults/A375-low-input-HLAI/$file/results/percolator/rescore_target.psms /Users/adams/Projects/300K/MSV000091456-SCP/reresults/rescore-tims/$file.psms
scp cadams@10.152.135.57:$searchdir/reresults/A375-low-input-HLAI/$file/results/percolator/original_target.psms /Users/adams/Projects/300K/MSV000091456-SCP/reresults/andromeda/$file.psms

Proteomics22!!

rm -r msms/ mzML/ proc/ results/ msms.prosit 