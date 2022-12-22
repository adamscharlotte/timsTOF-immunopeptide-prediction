
# --------------------------------------- MAP BRUKER ----------------------------------------
ssh cadams@10.152.135.57

grep -vxFf mapped_done.txt all-pool-names.txt > map_queue.txt
grep -vxFf cali_done.txt all-pool-names.txt > annot_queue.txt
scp Annotation/bruker_map.py cadams@10.152.135.57:/home/cadams
scp Names/tryptic_pool_names.txt cadams@10.152.135.57:/home/cadams

# screen -S annotate
screen -r annotate

# name_file=map_queue.txt
name_file=test.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # Rscript create_annotation_input.R "${array[0]}"
    Rscript tryptic_annotation_input.R "${array[0]}"
done

screen -r extract_d

name_file=tryptic_pool_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 bruker_map.py "${array[0]}"
done

# name_file=/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/pool-path-frames.txt
name_file=test-1.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript bruker_map.R ${array[1]} "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/extract-d/${array[0]}.csv" ${array[2]} ${array[3]}
done

# -------------------------------------- FILTER DATA ----------------------------------------

name_file=Names/all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # Rscript filter_annotation_input.R "${array[0]}"
    # time Rscript filter_and_add_ce.R "${array[0]}" "${array[1]}"
    # time Rscript filter_and_sum.R "${array[0]}"
    time Rscript Annotation/scan_filter_sum.R "${array[0]}"
done

# -------------------------------- CREATE CONSENSUS SPECTRA ---------------------------------

name_file=Names/all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/master_spectrum_precursor.py "${array[0]}"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/scan_master_spectrum.py "${array[0]}"
done

# ---------------------------------------- ANNOTATE -----------------------------------------

name_file=Names/all-pool-names.txt
# name_file=sample-pool-names.txt
# name_file=annot_queue
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/prosit_annotate.py "${array[0]}"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/frame_annotate.py "${array[0]}"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/sum_annotate.py "${array[0]}"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/scan_annotate.py "${array[0]}"
done

/Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/sum_annotate.py "TUM_HLA2_154"

# ---------------------------------------- CALIBRATE ----------------------------------------
ssh cadams@10.152.135.57

scp Annotation/frame_calibrate.py cadams@10.152.135.57:/home/cadams
scp annot_queue.txt cadams@10.152.135.57:/home/cadams

# name_file=all-pool-names.txt
name_file=annot_queue.txt
# name_file=sample-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 frame_calibrate.py "${array[0]}"
    # time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 sum_calibrate.py "${array[0]}"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 scan_calibrate.py "${array[0]}"
    printf '%d %s' "$i" "${array[i]}"
done

# --------------------------------- MAP FULL-LENGTH SEQUENCES --------------------------------

name_file=Names/all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time Rscript Annotation/full_length_mapper.R "${array[0]}"
done

# --------------------------------------- GENERATE HDF5 ---------------------------------------

name_file=Names/all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/calibrated_df_precursor_scans.py "${array[0]}"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/scan_df_precursor_scans.py "${array[0]}"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/trunc_peptides.py "${array[0]}"
done

# --------------------------------------- SHUFFLE SPLIT ----------------------------------------
ssh cadams@10.152.135.57
cd /home/cadams/ozapft/scripts

# Adjust the name of the df_precursor_scans file
view shuffle.sh

# Run script
bash shuffle.sh

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

name_file=sample_charge.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript Plots/sa_boxplot_compare.R "${array[0]}"
done

# ------------------------------------ CALIBRATION CE PLOTS ---------------------------------
# scp Annotation/scan_calibrate.py cadams@10.152.135.57:/home/cadams

name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 scan_calibrate.py "${array[0]}"
done

# ----------------------------------------- CE PLOTS ----------------------------------------
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

# ------------------------------------ GENERAL BRUKER MAP -----------------------------------

scp Annotation/bruker_map.R cadams@10.152.135.57:/home/cadams

Rscript bruker_map.R "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/PXD030334/190926_TIMSiDE_LCMS02_sample-3_90min_R2_Slot1-30_01_3622.d" "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/PXD030334/combined/txt" "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/PXD030334/precursor-mapped/S3_R2.csv"

# -------------------------------------- GENERAL STATS --------------------------------------

name_file=Names/all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time Rscript Analysis-code/identification_stats.R "${array[0]}"
done

# ------------------------------------- GENERATE FASTA --------------------------------------

name_file=Names/tryptic_pool_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time Rscript Fasta/create_fasta.R "${array[0]}"
done
