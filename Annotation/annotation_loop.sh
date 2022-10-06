ssh cadams@10.152.135.57
scp Annotation/create_annotation_input.R cadams@10.152.135.57:/home/cadams
scp Annotation/filter_annotation_input.R cadams@10.152.135.57:/home/cadams
scp Annotation/filter_and_sum.R cadams@10.152.135.57:/home/cadams

# cd /Users/adams/Projects/300K/2022-library-run/metadata/poolnames
# scp all-pool-names.txt cadams@10.152.135.57:/home/cadams
# scp mapped_done.txt cadams@10.152.135.57:/home/cadams
# scp Annotation/filter_and_add_ce.R cadams@10.152.135.57:/home/cadams

grep -vxFf mapped_done.txt all-pool-names.txt > map_queue.txt

name_file=map_queue.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript create_annotation_input.R "${array[0]}"
done

grep -vxFf done_done.txt done_pool.txt > filter-queue.txt

name_file=filter-queue.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript filter_annotation_input.R "${array[0]}"
done

Rscript filter_annotation_input.R "TUM_lysn_25"

filter_annotation_input.R "${array[0]}"
name_file=pool-queue.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time Rscript filter_and_add_ce.R "${array[0]}" "${array[1]}"
done

name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time Rscript filter_and_sum.R "${array[0]}"
done

scp Plots/prepare_ce_df.R cadams@10.152.135.57:/home/cadams
# scp Plots/prepare_qc_in_trap.R cadams@10.152.135.57:/home/cadams
# scp done_pool.txt cadams@10.152.135.57:/home/cadams

# cp done_pool.txt ce_queue.txt
grep -vxFf ce_done.txt done_pool.txt > ce-queue.txt

name_file=ce_queue.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript prepare_ce_df.R "${array[0]}"
done

name_file=done_pool.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript prepare_ce_mz.R "${array[0]}"
done

name_file=done_pool.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript prepare_qc_in_trap.R "${array[0]}"
done


# cd /Users/adams/Projects/300K

name_file=/Users/adams/Code/timsTOF-immunopeptide-prediction/done_pool.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 /Users/adams/Code/timsTOF-immunopeptide-prediction/Annotation/prosit_annotate.py "${array[0]}"
done

# /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 /Users/adams/Code/timsTOF-immunopeptide-prediction/Annotation/shuffle_split.py "test"

# --------------------------------------- MAP BRUKER ----------------------------------------
scp Annotation/bruker_map.R cadams@10.152.135.57:/home/cadams

Rscript bruker_map.R "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/PXD030334/190926_TIMSiDE_LCMS02_sample-3_90min_R2_Slot1-30_01_3622.d" "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/PXD030334/combined/txt" "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/PXD030334/precursor-mapped/S3_R2.csv"




# ---------------------------------------- ANNOTATE -----------------------------------------

name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/prosit_annotate.py "${array[0]}"
done

grep -vxFf ann_done.txt all-pool-names.txt > ann_queue.txt

name_file=sample-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/frame_annotate.py "${array[0]}"
done

name_file=all-pool-names.txt
# name_file=sample-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/sum_annotate.py "${array[0]}"
done

# ---------------------------------------- CALIBRATE ----------------------------------------

scp Annotation/frame_calibrate.py cadams@10.152.135.57:/home/cadams
scp Annotation/sum_calibrate.py cadams@10.152.135.57:/home/cadams

scp ann_queue.txt cadams@10.152.135.57:/home/cadams
# scp Annotation/prosit_annotate_calibrate.py cadams@10.152.135.57:/home/cadams

name_file=done_pool.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 prosit_annotate_calibrate.py "${array[0]}"
    printf '%d %s' "$i" "${array[i]}"
done

name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 calibrate_ce.py "${array[0]}"
    printf '%d %s' "$i" "${array[i]}"
done



scp /Users/adams/Projects/300K/fundamentals/fundamentals/fragments.py cadams@10.152.135.57:/home/cadams/fundamentals/fundamentals

grep -vxFf cal_done.txt  all-pool-names.txt > cal_queue.txt

name_file=sample-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 frame_calibrate.py "${array[0]}"
    printf '%d %s' "$i" "${array[i]}"
done

name_file=all-pool-names.txt
# name_file=sample-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 sum_calibrate.py "${array[0]}"
    printf '%d %s' "$i" "${array[i]}"
done

# --------------------------------------- SUM SPECTRA ---------------------------------------

# scp cadams@10.152.135.57:/home/cadams/all-pool-names.txt .
name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time Rscript Annotation/filter_and_sum.R "${array[0]}"
done

grep -vxFf un-annotated_done.txt all-pool-names.txt > annotation-queue.txt

# name_file=all-pool-names.txt
# lines=`tail -n+1 $name_file`
# for line in $lines
# do
#     IFS=';' read -r -a array <<< "$line"
#     time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/prosit_sum_annotate.py "${array[0]}"
# done

name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/master_spectrum_precursor.py "${array[0]}"
done

name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/annotate.py "${array[0]}"
done

name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/filter_generate_hsf5.py "${array[0]}"
done

# --------------------------------------- GENERATE HDF5 ---------------------------------------

name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/calibrated_hdf5.py "${array[0]}"
done


# ---------------------------------------- CALCULATE SA ----------------------------------------
scp Annotation/calculate_spectral_angle.py cadams@10.152.135.57:/home/cadams
# scp cadams@10.152.135.57:/home/cadams/sample-pool-names.txt .

# shuf -n 10 all-pool-names.txt > sample-pool-names.txt
name_file=sample-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 calculate_spectral_angle.py "${array[0]}"
done

name_file=sample-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/calculate_pairwise_sa.py "${array[0]}"
done

# ------------------------------------------- PLOT SA -------------------------------------------

name_file=sample_charge.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript Plots/sa_boxplot.R "${array[0]}"
done

name_file=sample_charge.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript Plots/sa_boxplot_pairwise.R "${array[0]}"
done