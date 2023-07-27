
# --------------------------------------- MAP BRUKER ----------------------------------------
name_file=map_queue.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # Rscript create_annotation_input.R "${array[0]}"
    Rscript tryptic_annotation_input.R "${array[0]}"
done

name_file=tryptic_pool_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 extract_pasef.py "${array[0]}"
done

name_file=d-missing.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript extract_d.R ${array[1]} "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/extract-d/${array[0]}.csv" "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/extract-pasef/${array[0]}.csv"
done

name_file=psm-missing.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript extract_psm.R "${array[0]}"
done

# -------------------------------------- FILTER DATA ----------------------------------------
name_file=psm-missing.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # Rscript filter_annotation_input.R "${array[0]}"
    # time Rscript filter_and_add_ce.R "${array[0]}" "${array[1]}"
    # time Rscript filter_and_sum.R "${array[0]}"
    time Rscript scan_filter_sum.R "${array[0]}"
done

# -------------------------------- CREATE CONSENSUS SPECTRA ---------------------------------
name_file=anno_queu.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 scan_master_spectrum.py "${array[0]}"
done

name_file=sa_queue.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 scan_prediction_sa.py "${array[0]}"
done

# ---------------------------------------- ANNOTATE -----------------------------------------
name_file=anno_queu.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/prosit_annotate.py "${array[0]}"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/frame_annotate.py "${array[0]}"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/sum_annotate.py "${array[0]}"
    time /home/cadams/anaconda3/envs/prosit-annotate/bin/python3 scan_annotate.py "${array[0]}"
done

# ---------------------------------------- CALIBRATE ----------------------------------------
name_file=annot_queue.txt
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
name_file=Names/train_set.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/calibrated_df_precursor_scans.py "${array[0]}"
    # time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/scan_df_precursor_scans.py "${array[0]}"
    time /Users/adams/opt/miniconda3/envs/prosit-annotate/bin/python3 Annotation/scan_hdf5_nocal.py "${array[0]}"
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

# ------------------------------------- GENERATE FASTA --------------------------------------
name_file=Names/tryptic_pool_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time Rscript Fasta/create_fasta.R "${array[0]}"
done