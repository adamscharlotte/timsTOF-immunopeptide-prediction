#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#
#SBATCH --mem-per-cpu=60
#SBATCH --time=06:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=charlotte.adams@uantwerpen.be
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j

module purge
module load calcua/2020a
module load calcua/supported
module load calcua/x86_64
module load calcua/system
module load monitor
module use /data/antwerpen/200/vsc20003/easybuild/modules/vaughan/2021b
module use /data/antwerpen/203/vsc20366/EasyBuild/modules/centos8/rome/2021b
module load dotNET-SDK/3.1.300-linux-x64
module load MaxQuant/2.3.1.0-GCCcore-11.2.0

monitor -l mem-raw.log -d 5 -- mono MaxQuant_2.0.3.1/bin/MaxQuantCmd.exe $VSC_SCRATCH/PXD038782-comp/mqpar-raw-worked.xml
