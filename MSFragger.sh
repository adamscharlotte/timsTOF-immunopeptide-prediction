# uninstall java
# sudo rm -fr /Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin
# sudo rm -fr /Library/PreferencesPanes/JavaControlPanel.prefPane
# sudo rm -fr ~/Library/Application\ Support/Java

#!/bin/bash
set -xe

dataDirPath="/Users/adams/Projects/300K/20220613-library-test/raw/"
fastaPath="/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_1.fasta"
msfraggerPath="/Users/adams/Downloads/MSFragger-3.5/MSFragger-3.5.jar" # download from http://msfragger-upgrader.nesvilab.org/upgrader/
fraggerParamsPath="/Users/adams/Projects/300K/20220613-library-test/params/msfragger_nonspecific.params"
philosopherPath="philosopher" # download from https://github.com/Nesvilab/philosopher/releases/latest
ionquantPath="IonQuant.jar" # download from https://github.com/Nesvilab/IonQuant/releases/latest
decoyPrefix="rev_"

java -Xmx64G -jar $msfraggerPath $fraggerParamsPath $dataDirPath/E1_50fmol_bruker_S1-E2_1_6339.d

java -jar MSFragger-3.5/MSFragger-3.5.jar

java -jar MSfragger-GUI.jar
