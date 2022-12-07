ssh -i /Users/adams/.ssh/id_rsa_ua vsc20709@login-leibniz.hpc.uantwerpen.be

module load Java

time java -jar -Xmx112G MSFragger-3.2/MSFragger-3.2.jar msfragger/CRC_closed.params /data/antwerpen/207/vsc20709/interim/mzML/I-CRC-05/*
