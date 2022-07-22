## How to use this Queue

1. Create a file `queue_tbd.csv` based on `queue_template.csv` with your search parameters
2. Append lines to `RawFilesAnnotation.csv` based on `RawFilesAnnotation_template.csv`
3. To start processing the searches in `queue_tbd.csv`, copy the lines into the `queue.csv` file.

### queue.csv

The comma separated file should start with a header line (this header is not parsed though). The following columns are mandatory, do not leave columns empty in the middle of a row as this will break the queue (you can however leave columns empty at the end of a row):

1. `_name`: Specifies subfolder of MQ processing folder (results will be written to `<prdb1_drive>\<_project>\<_name>`)
2. `_MQ_version`: MaxQuant version, currently only 1.5.3.8 is supported. 1.6.12.0 is in progress
3. `_pre_payload`: Bash function to run before the MQ search (see `additionalPayloads.sh`)
4. `_post_payload`: Bash function to run after the MQ search (see `additionalPayloads.sh`)
5. `_threads`: Maximum number of threads MQ uses. Note that for certain steps MQ will not be able to use more than a couple of cores. In some cases it's therefore better to process multiple datasets in parallel rather than allocating many cores per dataset.
6. `_project`: Specifies main folder of MQ processing folder (results will be written to `<prdb1_drive>\<_project>\<_name>`)
7. `_fasta_file`: Fasta file, should be present in `<prdb1_drive>\MaxQuant`
8. `_raw_folder`: Location of the raw files specified in `RawFilesAnnotation.csv` (drives: `w`=`prdb1`, `y`=`raw_files`, `z`=`kusterlab`)
9. `_phospho`: Set to `1` to search for phospho (STY) and `0` to not search phospho.
10. `_protease`: Set to protease string. Typically `Trypin` or `Trypsin/P`, for other proteases check the `Configuration` tab of the MQ GUI.
11. `_mqpar_base`: Specify which template to use to generate the `mqpar.xml`, e.g. `mqpar_base_<mq_version>.xml`. If left empty, it will use `mqpar_base_<mq_version>.xml` (but do not leave empty if you have a value in the `_tmt_corr_factors` column!)
12. `_tmt_corr_factors`: TMT correction factor batch number (including `.txt` extension). The correction factors hould be placed as a `.txt` file in `<prdb1_drive>\MaxQuant\TMT_correction_factors`. Leave empty in case of LFQ or if no correction factors are available.
13. `_peptide_fdr`: Peptide FDR threshold. Default value: 0.01 (=1%)
14. `_protein_fdr`: Protein FDR threshold. Default value: 1 (=100%)


### RawFilesAnnotation.csv

The comma separated file should start with a header line (not parsed). The following columns are mandatory:

1. `File`: Raw file name without full path (the path is specified in `queue.csv`) and without `.raw` extension
2. `Output_folder`: Should correspond to `_name` field of the corresponding row in `queue.csv`. This (partly) specifies the MQ processing folder.
3. `Experiment`: `Experiment` field in `mqpar.xml`. Fractionated samples should have the same experiment name
4. `MS_Plate`: Should be the same for the RAW files searched together in one search, e.g. can typically be set to a fixed string for an entire project or, to be more safe, taken as the same as `Output_folder`. This column can be used to reuse the same RAW file in multiple searches, e.g. use the same DMSO for different `Output_folder`s.
5. `Fraction`: `Fraction` field in `mqpar.xml`
6. `Parameter_group`: `Parameter group` field in `mqpar.xml`
7. `PTM`: `PTM` field in `mqpar.xml`. If set to `True` this run will not be used for protein quantification.
8. `Reference_channel`: `Reference channel` field in `mqpar.xml`. Used to specify the bridge channel for isobaric match between runs. Leave empty if you don't want to use the isobaric match between runs normalization.
9. `File_order`: Not used
10. `Compound_order`: Not used
11. `Project`: Should correspond to `_project` field of the corresponding row in `queue.csv`. This (partly) specifies the MQ processing folder and acts as a filter for which Raw files to add to the `mqpar.xml`.
12. `Comments`: Comments

Parameter groups in `mqpar_base_1.6.12.0.xml`: 

0. TMT Full proteome (variable mods: M(ox) and (ac)N-term)
1. TMT Phospho search (extra variable mods: Phospho (STY))
2. TMT Acetylation search (extra variable mods: Acetyl (K))

Parameter groups in `mqpar_base_ubi_1.6.12.0.xml`: 
0. LFQ Ubiquitination search (extra variable mods: GlyGly (K))

Parameter groups in `mqpar_base_var_1.6.12.0.xml`: 
0. LFQ-like search with TMT as variable modification, used for QC (labeling efficiency, etc.)

## How the Queue works

The `runK.sh` script runs an endless loop, reading the `queue.csv` file every `X` seconds (e.g. every 2 minutes). It processes the `queue.csv` file line by line and checks if there are enough cores available (`NCORES` in `runK-<server>.cfg`) to start a job. If this is the case, it writes a line to `queue_done.csv` with the status `processing` and runs the specified `_pre_payload` function (e.g. copying data and creating an `mqpar.xml` file), runs MaxQuant and finally runs the specified `_post_payload` function.

To restart a failed job, just remove the line from `queue_done.csv` and the Queue will try to re-run this job the next time the `queue.csv` file is checked.

You can make changes to the configuration file (`runK-<server>.cfg`) without having to restart the Queue, the new configuration will automatically be loaded the next time `queue.csv` is checked. Changes to `additionalPayloads.sh` and the `payload` folder also don't require restarting of the Queue. Only if the `runK.sh` script is changed, the Queue needs to be restarted (see below).

### Set-up the Queue in a new folder

```
git clone https://gitlab.lrz.de/proteomics/infrastructure/queue.git
```

1. Create a file `queue_tbd.csv` based on `queue_template.csv` with your searches
2. Create `RawFilesAnnotation.csv` based on `RawFilesAnnotation_template.csv`
3. Rename `queue_done_template.csv` to `queue_done.csv`
4. Change the configuration in `runK-<server>.cfg`
5. If `LOCAL_PROCESSING` is turned on, remember to create a project folder on the scratch drive `D:/PrDB/<project>`.
6. Start the Queue by running `./runK.sh start` in cygwin`
7. To start processing the searches in `queue_tbd.csv`, either rename it to `queue.csv`, or copy the lines into a `queue.csv` file.

To track this folder in git with the custom configuration file:

1. Create a new branch
2. Merge in changes from the master branch with `git pull; git merge master` (make sure this does not overwrite changes you made to your local configuration file).
3. Merge changes into the master branch with `git cherry-pick <my_commit>`

### How to start this Queue (bioinformaticians only!)

1. Make sure that the proper drives are mapped (`W`=`prdb1`, `Y`=`raw_files`, `Z`=`kusterlab`)
2. Open the `cygwin` terminal.
3. Navigate to this repository `cd /cygwin/w/Queue`
4. Start the queue `./runK.sh start`
5. You can safely close the terminal

### How to kill the Queue  (bioinformaticians only!)

1. Find the process id in `.runK.sh-<server>.pid`
2. Open `cygwin` (if the queue was started under another user, run `cygwin` as administrator)
3. Kill the process `kill <pid>`

## How to troubleshoot the Queue

You can find log files for each search in the `log-<server>` folder.
