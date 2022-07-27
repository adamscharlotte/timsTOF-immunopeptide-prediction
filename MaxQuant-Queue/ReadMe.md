## How to use this Queue

1. Create a file `queue.csv` using create_queue.R
2. Define the path to the queue in the cfg file
3. Go to the remote desktop and open the cygwin terminal
4. Copy the code to a folder you can access there and go to it
5. Run ./runK.sh start to start the queue
6. In the folder a .pid file will be created
7. To end the queue run kill `number in the pid`

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
11. `_mqpar_base`: Specify which template to use to generate the `mqpar.xml`, e.g. `mqpar_base_<mq_version>.xml`. If left empty, it will use `mqpar_base_<mq_version>.xml` (but do not leave empty if you have a value in the `_tmt_corr_factors` column!)

## How the Queue works

The `runK.sh` script runs an endless loop, reading the `queue.csv` file every `X` seconds (e.g. every 2 minutes). It processes the `queue.csv` file line by line and checks if there are enough cores available (`NCORES` in `runK-<server>.cfg`) to start a job. If this is the case, it writes a line to `queue_done.csv` with the status `processing` and runs the specified `_pre_payload` function (e.g. copying data and creating an `mqpar.xml` file), runs MaxQuant and finally runs the specified `_post_payload` function.

To restart a failed job, just remove the line from `queue_done.csv` and the Queue will try to re-run this job the next time the `queue.csv` file is checked.

You can make changes to the configuration file (`runK-<server>.cfg`) without having to restart the Queue, the new configuration will automatically be loaded the next time `queue.csv` is checked. Changes to `additionalPayloads.sh` and the `payload` folder also don't require restarting of the Queue. Only if the `runK.sh` script is changed, the Queue needs to be restarted (see below).

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
