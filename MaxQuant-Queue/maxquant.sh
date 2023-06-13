#   CalcUA
# ssh -i /Users/adams/.ssh/id_rsa_ua vsc20709@login-leibniz.hpc.uantwerpen.be
ssh -i /Users/adams/.ssh/id_rsa_ua vsc20709@login-vaughan.hpc.uantwerpen.be

scp -i /Users/adams/.ssh/id_rsa_ua -r HLAI_p2_97_178_p2-D6_S1-D6_1_6871.d vsc20709@login-leibniz.hpc.uantwerpen.be:/scratch/antwerpen/207/vsc20709/Kuster_test
scp -i /Users/adams/.ssh/id_rsa_ua -r MaxQuant-Queue/runMaxQuant.sh vsc20709@login-leibniz.hpc.uantwerpen.be:/user/antwerpen/207/vsc20709/
scp -i /Users/adams/.ssh/id_rsa_ua MaxQuant-Queue/payloads/maxdel/MaxDEL.sh vsc20709@login-leibniz.hpc.uantwerpen.be:/user/antwerpen/207/vsc20709/

scp -i /Users/adams/.ssh/id_rsa_ua -r UP000005640_9606.fasta vsc20709@login-leibniz.hpc.uantwerpen.be:/scratch/antwerpen/207/vsc20709/PXD038782-comp
scp -i /Users/adams/.ssh/id_rsa_ua -r 220331_NHG_malignant_CLL_02_W6-32_17%_DDA_Rep3.d vsc20709@login-leibniz.hpc.uantwerpen.be:/scratch/antwerpen/207/vsc20709/PXD038782-comp
scp -i /Users/adams/.ssh/id_rsa_ua -r 220329_NHG_benign_UDN31_PBMC_Tue39L243_17%_orbitrap_DDA_Rep1.raw vsc20709@login-leibniz.hpc.uantwerpen.be:/scratch/antwerpen/207/vsc20709/PXD038782-comp
scp -i /Users/adams/.ssh/id_rsa_ua -r MaxQuant_2.0.3.1 vsc20709@login-leibniz.hpc.uantwerpen.be:/user/antwerpen/207/vsc20709
scp -i /Users/adams/.ssh/id_rsa_ua mqpar-raw-worked.xml vsc20709@login-leibniz.hpc.uantwerpen.be:/scratch/antwerpen/207/vsc20709/PXD038782-comp

scp -i /Users/adams/.ssh/id_rsa_ua -r vsc20709@login-leibniz.hpc.uantwerpen.be:/user/antwerpen/207/vsc20709/mem.log .
scp -i /Users/adams/.ssh/id_rsa_ua -r vsc20709@login-leibniz.hpc.uantwerpen.be:/scratch/antwerpen/207/vsc20709/Kuster_test/mqpar-20.xml .

cd /Users/adams/Projects/300K/Jurkat/Jurkat_A549

scp -r Jurkat cadams@10.152.135.57:/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/

cd /scratch/antwerpen/grp/abiodatamining_students/
cd $VSC_SCRATCH

sbatch runMaxQuant.sh
squeue
scancel 675328

srun --pty bash
ssh r1c06cn4

module purge
module load monitor # Mogelijks moet purge weg voor deze
module use /data/antwerpen/200/vsc20003/easybuild/modules/vaughan/2021b
module use /data/antwerpen/203/vsc20366/EasyBuild/modules/centos8/rome/2021b
module load dotNET-SDK/3.1.300-linux-x64
module load MaxQuant/2.3.1.0-GCCcore-11.2.0
# module load MaxQuant/2.0.3.1-GCCcore-11.2.0

mono MaxQuant_2.0.3.1/bin/MaxQuantCmd.exe $VSC_SCRATCH/PXD038782-comp/mqpar-raw-worked.xml

monitor -l mem.log -d 5 -- mono $EBROOTMAXQUANT/bin/MaxQuantCmd.exe $VSC_SCRATCH/Kuster_test/mqpar-1.xml
dotnet $EBROOTMAXQUANT/bin/MaxQuantCmd.exe $VSC_SCRATCH/Kuster_test/mqpar-1.xml

#   TUM linux-ml
ssh cadams@10.152.135.57

scp -r /Users/adams/Projects/300K/PXD038782-comparison cadams@10.152.135.57:/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker

cd /media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison

for f in *.tar; do tar -m --no-overwrite-dir -xvf "$f"; done

cd /media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/d-folder
mkdir tar
mv *.d.tar tar/
mv mnt/f/Pride_Upload/umgenannt/*/*.d .

# Clean files

MaxDEL.sh -c "simple" -f "/scratch/antwerpen/207/vsc20709/PXD038782-comp/Searches"