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

# Rescore SCP data
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

file="E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep1_Slot1-37_1_3502"
file="E_20221201_NO30_400nL_HLAc1_5e6_directIP_titration_rep2_Slot2-2_1_3510"
file="E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep2_Slot2-5_1_3520"
searchdir=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/A375
mkdir $searchdir/reresults/A375-low-input-HLAI/$file
scp $searchdir/reresults/test.json $searchdir/reresults/A375-low-input-HLAI/$file
scp $searchdir/Searches/A375-low-input-HLAI/$file/combined/txt/msms.txt $searchdir/reresults/A375-low-input-HLAI/$file
scp $searchdir/Searches/A375-low-input-HLAI/$file/combined/txt/pasefMsmsScans.txt $searchdir/reresults/A375-low-input-HLAI/$file
scp $searchdir/Searches/A375-low-input-HLAI/$file/combined/txt/accumulatedMsmsScans.txt $searchdir/reresults/A375-low-input-HLAI/$file
# scp $searchdir/Searches/A375-low-input-HLAI/$file/$file.pkl $searchdir/reresults/A375-low-input-HLAI/$file
time Rscript d_extract.R "$searchdir/A375-low-input-HLAI/$file.d" "$searchdir/reresults/A375-low-input-HLAI/$file" "/media/c2m/5_user_files/cadams/libtimsdata.so"
time python d_master_spec_2.py "$searchdir/reresults/A375-low-input-HLAI/$file/$file.csv"
time python /home/cadams/oktoberfest-0.1.0/oktoberfest/run_oktoberfest.py --search_dir $searchdir/reresults/A375-low-input-HLAI/$file --config_path $searchdir/reresults/A375-low-input-HLAI/$file/test.json


rm $searchdir/reresults/A375-low-input-HLAI/$file/$file.csv
