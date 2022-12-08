ssh -i /Users/adams/.ssh/id_rsa_ua vsc20709@login-leibniz.hpc.uantwerpen.be

# ---------------------------------------- COPY FILES ---------------------------------------
scp -i /Users/adams/.ssh/id_rsa_ua MSFragger/TUM_HLA_3.params vsc20709@login-leibniz.hpc.uantwerpen.be:/user/antwerpen/207/vsc20709/msfragger
scp -i /Users/adams/.ssh/id_rsa_ua /Users/adams/Projects/300K/2022-library-run/fasta-decoy/TUM_HLA_3.fasta vsc20709@login-leibniz.hpc.uantwerpen.be:/data/antwerpen/207/vsc20709/fasta
scp -i /Users/adams/.ssh/id_rsa_ua /Users/adams/Projects/300K/2022-library-run/precursor-mgf/test/TUM_HLA_3.mgf vsc20709@login-leibniz.hpc.uantwerpen.be:/data/antwerpen/207/vsc20709/MGF/
scp -i /Users/adams/.ssh/id_rsa_ua -r /Users/adams/Projects/300K/2022-library-run/raw-folders/HLAI_1_96_p2-A3_S2-A3_1_6928.d vsc20709@login-leibniz.hpc.uantwerpen.be:/data/antwerpen/207/vsc20709/bruker/

# --------------------------------------- RUN MSFRAGGER -------------------------------------

module load Java

time java -jar -Xmx112G MSFragger-3.2/MSFragger-3.2.jar msfragger/TUM_HLA_3.params /data/antwerpen/207/vsc20709/MGF/TUM_HLA_3.mgf
time java -jar -Xmx112G MSFragger-3.2/MSFragger-3.2.jar msfragger/TUM_HLA_3.params /data/antwerpen/207/vsc20709/bruker/HLAI_1_96_p2-A3_S2-A3_1_6928.d

cd /data/antwerpen/207/vsc20709/MGF/

# ---------------------------------------- GET OUTPUT ---------------------------------------

scp -i /Users/adams/.ssh/id_rsa_ua vsc20709@login-leibniz.hpc.uantwerpen.be:/data/antwerpen/207/vsc20709/MGF/TUM_HLA_3*.tsv /Users/adams/Projects/300K/2022-library-run/msfragger-results/test/
scp -i /Users/adams/.ssh/id_rsa_ua vsc20709@login-leibniz.hpc.uantwerpen.be:/data/antwerpen/207/vsc20709/bruker/HLAI_1_96_p2-A3_S2-A3_1_6928.tsv /Users/adams/Projects/300K/2022-library-run/msfragger-results/test/d
