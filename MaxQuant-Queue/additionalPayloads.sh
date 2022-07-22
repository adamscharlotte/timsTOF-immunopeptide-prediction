#!/usr/bin/env bash

function cleanPostSearchAndQCTopasFP() {

	NAME=$1
	BASE=$2
	THREADS=$3
	
	cleanPostSearch $1 $2 $3
	runQCTopas $1 $2 $3 FP
	
}

function cleanPostSearchAndQCTopasPP() {

	NAME=$1
	BASE=$2
	THREADS=$3
	
	cleanPostSearch $1 $2 $3
	runQCTopas $1 $2 $3 PP
	
}

function cleanPostSearchAndQCTopasVar() {

	NAME=$1
	BASE=$2
	THREADS=$3
	
	cleanPostSearch $1 $2 $3
	runQCTopas $1 $2 $3 var
	
}

function cleanPostSearch() {

	NAME=$1
	BASE=$2
	THREADS=$3

	DIR="$DESTINATION/$BASE/$NAME"
	"$PAYLOAD_FOLDER/maxdel/MaxDEL.sh" -c "simpleraw" -f "$DIR" # keeps andromeda search results, removes raw files
	
}

function cleanSearch() {

	NAME=$1
	BASE=$2
	THREADS=$3

	DIR="$DESTINATION/$BASE/$NAME"
	
}

function changeDirectoryLocal() {

	NAME=$1
	BASE=$2
	THREADS=$3
	
	DIR_BASE="$DESTINATION_LOCAL/$BASE"

	if [ ! -d "$DIR_BASE" ]; then
		echo -e "\t\tCreating folder for $BASE"
		mkdir $DIR_BASE
	fi
	
	DIR="$DIR_BASE/$NAME"
	LOC="$DIR/mqpar.xml"

	if [ ! -d "$DIR" ]; then
		echo -e "\t\tCreating folder for $NAME"
		mkdir $DIR
	fi

	cd "$DIR"

}

function changeDirectory() {
	
	NAME=$1
	BASE=$2
	THREADS=$3
	
	DIR_BASE="$DESTINATION/$BASE"

	if [ ! -d "$DIR_BASE" ]; then
		echo -e "\t\tCreating folder for $BASE"
		mkdir $DIR_BASE
	fi
	
	DIR="$DIR_BASE/$NAME"
	LOC="$DIR/mqpar.xml"

	if [ ! -d "$DIR" ]; then
		echo -e "\t\tCreating folder for $NAME"
		mkdir $DIR
	fi

	cd "$DIR"
	
}

function moveToLocal() {

	NAME=$1
	BASE=$2
	THREADS=$3

	LOCAL_DIR_BASE="$DESTINATION_LOCAL/$BASE"

	if [ ! -d "$LOCAL_DIR_BASE" ]; then
		echo -e "\t\tCreating folder for $BASE"
		mkdir $LOCAL_DIR_BASE
	fi
	
	REMOTE_DIR="$DESTINATION/$BASE/$NAME"
	LOCAL_DIR="$LOCAL_DIR_BASE/$NAME"

	if [ ! -d "$LOCAL_DIR" ]; then
		echo -e "\t\tCreating folder for $NAME"
		mkdir $LOCAL_DIR
	fi
	
	echo "cp -Rup $REMOTE_DIR/* $LOCAL_DIR/"
	cp -Rup $REMOTE_DIR/* $LOCAL_DIR/
	
	echo "cd $LOCAL_DIR"
	cd $LOCAL_DIR
	
}

function moveToRemote() {

	NAME=$1
	BASE=$2
	THREADS=$3

	REMOTE_DIR="$DESTINATION/$BASE/$NAME"
	LOCAL_DIR="$DESTINATION_LOCAL/$BASE/$NAME"

	echo "cp -Rup $LOCAL_DIR/* $REMOTE_DIR/"
	cp -Rup $LOCAL_DIR/* $REMOTE_DIR/
	
	echo "rm -rf $LOCAL_DIR/*"
	rm -rf $LOCAL_DIR/*
	
	echo "cd $REMOTE_DIR"
	cd $REMOTE_DIR
	
}

# function testGenerateMQpar() {
# 	RAWFILE_POSITION="1"
# 	FOLDER_POSITION="2"
# 	CONDITION_POSITION="3"
# 	MSPLATE_POSITION="4"
# 	FRACTION_POSITION="5"
# 	PARAMETERGROUP_POSITION="6"
# 	PTM_QUANT_POSITION="7" # controls if PTMs are used for protein quantification
# 	REFCHANNEL_POSITION="8" # used for isobaric match between runs as bridge channel
	
# 	LOCAL_PROCESSING=0
# 	PATH_PREFIX="test\\\\"
# 	FASTA_FILE_ROOT="fasta_test\\\\"
# 	MAPPING="RawFilesAnnotation_template.csv"
# 	BASE_MQPAR="MaxQuant/mqpar_base_1.6.12.0.xml"
# 	generateMQpar Coleoptile RICE 12 Rice_UniProt.fasta 1 "Trypsin/P" "Maxquant/mqpar_base_1.6.12.0.xml" Maxquant/TMT_correction_factors/XX000000.txt 0.01 1 raw
# }


function generateMQpar() {

	NAME=$1
	BASE=$2
	THREADS=$3
	FASTA_FILE=$4
	RAW_FOLDER=$5
	MAX_LENGTH=$6
	
	# for backwards compatibility, in newer versions this is controlled by the parameter group column in RawFilesAnnotation
	# if [ "$PHOSPHO" = "1" ]; then
	# 	ADDMOD="|            <string>Phospho (STY)</string>"
	# else
	# 	ADDMOD=""
	# fi
	
	
	# if [ "$LOCAL_PROCESSING" = "1" ]; then
	# 	PATH_PREFIX_L=$(echo "$PATH_PREFIX_LOCAL$BASE")
	# else
	# 	PATH_PREFIX_L=$(echo "$PATH_PREFIX$BASE")
	# fi

	LOC="mqpar.xml"

	echo -e "\tGenerating $1"

	if [ -e $LOC ]; then
		echo -e "\t\tFound old file, saving ..."
		mv $LOC $LOC`date +"_%y%m%d-%H%M%S"`
	fi

	MSPLATE=$(grep ",$NAME," $MAPPING | awk -F "," -v x=$MSPLATE_POSITION '{print $x}' | uniq)
	LMAPPING="localMapping.txt"

	#add real files
	grep $MSPLATE $MAPPING | grep $NAME | awk -F "," '{print $1".d,"$2","$3","$4""}' > $LMAPPING
	# grep $MSPLATE $MAPPING | grep ",$NAME," | awk -v EXT="$RAW_FILE_EXT" -F "," '{print $1"."d","$2","$3","$4","$5","$6","$7","$8}' > $LMAPPING
	# grep $MSPLATE $MAPPING | grep ",DMSO," | awk -v EXT="$RAW_FILE_EXT" -F "," '{print $1"."EXT","$2","$2","$4","$5","$6","$7","$8}' >> $LMAPPING

	#add ctrl files
	#$DMSO_FOLDER/*.raw | awk -F "," -v x=$MSPLATE '{print $1","substr($1, 1, length($1)-4)","substr($1, 1, length($1)-4)","x",1,1"}' >> $LMAPPING
	
	# TMT_FACTORS_STRING=""
	# if [ ! -z $TMT_CORR_FACTORS_FILE ]; then
	# 	TMT_FACTORS_STRING=$(generateMQparTMTCorrectionFactors $TMTFACTORS_FOLDER/$TMT_CORR_FACTORS_FILE)
	# fi
	
	ACTIVE_BASE_MQPAR="$ROOT_WD/Queue/Parameters/$ALT_BASE_MQPAR"
	
	# temporarily replace newlines with "^^" (the record separator, which 
	# normally should not occur in text files) for sed to work
	sed -e "s#%%RAW_FOLDER%%#`cat $MAPPING-$NAME |
		awk -F "," -v x=$RAW_FOLDER_POSITION '{print $x}' |
		sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
		tr "\n" "|"`#g" |
		sed -e "s#%%FASTA_FILE%%#`cat $MAPPING-$NAME |
			awk -F "," -v x=$FASTA_FILE_POSITION '{print $x}' |
			sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
			tr "\n" "|"`#g" |
		sed -e "s#%%MAX_LENGT%%#`cat $MAPPING-$NAME |
			awk -F "," -v x=$MAX_LENGT_POSITION '{print $x}' |
			sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
			tr "\n" "|"`#g" |
		sed -e "s#%%OUTPUT_FOLDER%%#`cat $MAPPING-$NAME |
			awk -F "," -v x=$OUTPUT_FOLDER_POSITION '{print $x}' |
			sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
			tr "\n" "|"`#g" |
		sed -e "s#\!#/#g" > $LOC

	# sed -e "s#%%RAW_FILES%%#`cat $LMAPPING |
	# 		awk -F "," -v x=$RAWFILE_POSITION '{print $x}' |
	# 		sed -e "s/^/      \<string\>$PATH_PREFIX_L\!$NAME\!/g" |
	# 		sed -e "s/$/\<\/string\>/g" | tr "\n" "^^"`#g" $ACTIVE_BASE_MQPAR |
	# 	sed -e "s#%%EXPERIMENT_NAMES%%#`cat $LMAPPING |
	# 		awk -F "," -v x=$CONDITION_POSITION '{print $x}' |
	# 		sed -e "s/^/      \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
	# 		tr "\n" "^^"`#g" |
	# 	sed -e "s#%%IS_PTM_SEARCH%%#`cat $LMAPPING |
	# 		awk -F "," -v x=$PTM_QUANT_POSITION '{print $x}' |
	# 		sed -e "s/^/      \<boolean\>/g" | sed -e "s/$/\<\/boolean\>/g" |
	# 		tr "\n" "^^"`#g" |
	# 	sed -e "s#%%FRACTION_NUMBER%%#`cat $LMAPPING |
	# 		awk -F "," -v x=$FRACTION_POSITION '{print $x}' |
	# 		sed -e "s/^/      \<short\>/g" | sed -e "s/$/\<\/short\>/g" |
	# 		tr "\n" "^^"`#g" |
	# 	sed -e "s#%%PARAMETER_GROUP%%#`cat $LMAPPING |
	# 		awk -F "," -v x=$PARAMETERGROUP_POSITION '{print $x}' |
	# 		sed -e "s/^/      \<int\>/g" | sed -e "s/$/\<\/int\>/g" |
	# 		tr "\n" "^^"`#g" |
	# 	sed -e "s#%%REF_CHANNELS%%#`cat $LMAPPING |
	# 		awk -F "," -v x=$REFCHANNEL_POSITION '{print $x}' |
	# 		sed -e "s/^/      \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
	# 		tr "\n" "^^"`#g" |
	# 	sed -e "s#%%TMT_CORR_FACTORS%%#`echo "$TMT_FACTORS_STRING" | tr "\n" "^^"`#g" |
	# 	sed -e "s#%%FASTA_FILE%%#$FASTA_FILE_ROOT$FASTA#" |
	# 	sed -e "s#%%PEPTIDE_FDR%%#$PEPTIDE_FDR#" |
	# 	sed -e "s#%%PROTEIN_FDR%%#$PROTEIN_FDR#" |
	# 	sed -e "s#%%PHOSPHO%%#$ADDMOD#" |
	# 	sed -e "s#%%PROTEASE%%#$PROTEASE_IN#" |
	# 	sed -e "s#%%THREADS%%#$THREADS#" |
	# 	tr "^^" "\n" |
	# 	sed -e "s#\!#\\\\#g" > $LOC

	unix2dos -q $LOC

	diff $ACTIVE_BASE_MQPAR $LOC | sed -e "s/^/\t\t/g"
	
	# if [ "$RAW_FILE_EXT" != "raw" ]; then
	# 	sed -i -e "s#$RAW_FILE_EXT#raw#" $LMAPPING
	# fi

}

# function testGenerateMQparTMTCorrectionFactors() {
# 	TMT_VAR=$(generateMQparTMTCorrectionFactors Maxquant/TMT_correction_factors/XX000000.txt)
# 	echo "$TMT_VAR"
# }

# function generateMQparTMTCorrectionFactors() {
	
# 	TMT_FACTORS_FILE=$1
	
# 	TMT_FACTORS_TEMPLATE="
#             <IsobaricLabelInfo>
#                <internalLabel>%%INTERNAL_LABEL%%</internalLabel>
#                <terminalLabel>%%TERMINAL_LABEL%%</terminalLabel>
#                <correctionFactorM2>%%FACTOR_M2%%</correctionFactorM2>
#                <correctionFactorM1>%%FACTOR_M1%%</correctionFactorM1>
#                <correctionFactorP1>%%FACTOR_P1%%</correctionFactorP1>
#                <correctionFactorP2>%%FACTOR_P2%%</correctionFactorP2>
#                <tmtLike>%%TMT_LIKE%%</tmtLike>
#             </IsobaricLabelInfo>"
	
# 	TMT_FACTORS_STRING=""
# 	FIRST="0"
# 	while IFS='' read -r line || [[ -n "$line" ]]; do
# 		if [ "$FIRST" -eq "0" ]; then
# 			FIRST="1"
# 			continue
# 		fi
# 		INTERNAL_LABEL=$(echo $line | awk '{print $1}')
# 		TERMINAL_LABEL=$(echo $line | awk '{print $2}')
# 		FACTOR_M2=$(echo $line | awk '{print $3}')
# 		FACTOR_M1=$(echo $line | awk '{print $4}')
# 		FACTOR_P1=$(echo $line | awk '{print $5}')
# 		FACTOR_P2=$(echo $line | awk '{print $6}')
# 		TMT_LIKE=$(echo $line | awk '{print tolower($7)}' | tr -d " \r\n")
		
# 		TMT_FACTORS_STRING+=$(echo "$TMT_FACTORS_TEMPLATE" | 
# 			sed -e "s#%%INTERNAL_LABEL%%#$INTERNAL_LABEL#g" |
# 			sed -e "s#%%TERMINAL_LABEL%%#$TERMINAL_LABEL#g" |
# 			sed -e "s#%%FACTOR_M2%%#$FACTOR_M2#g" |
# 			sed -e "s#%%FACTOR_M1%%#$FACTOR_M1#g" |
# 			sed -e "s#%%FACTOR_P1%%#$FACTOR_P1#g" |
# 			sed -e "s#%%FACTOR_P2%%#$FACTOR_P2#g" |
# 			sed -e "s#%%TMT_LIKE%%#${TMT_LIKE^}#g")
# 	done < $TMT_FACTORS_FILE
# 	echo "$TMT_FACTORS_STRING"
# }

# function cpCTRLs() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3
	
# 	LMAPPING="localMapping.txt"

# 	echo "copying controls"
# 	for file in `ls $DMSO_FOLDER/*.raw`; do
# 		if [ ! -e `basename $file` ]; then
# 			cp -v $file .
# 		else
# 			echo -E "$file already located in here"
# 		fi
# 	done
# 	echo "done copying DMSOs"
# }

# function cpRaws() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3
# 	FASTA_FILE=$4
# 	RAW_FOLDER=$5
# 	MAX_LENGTH=$6

# 	LMAPPING="localMapping.txt"

# 	for file in `cat $LMAPPING | grep ",$NAME," | awk -F ',' -v x=$RAW_FOLDER_POSITION '{print $x}'`; do
# 		if [ ! -e `basename $file` ]; then
# 			# cp -v $RAW_LOCATION/$file .
# 			cp -v $RAW_LOCATION$RAW_FILE_FOLDER/$file .
# 		else
# 			echo -E "$file already located in here"
# 		fi
# 	done
# }

# function convertRaws() {	
# 	"$MSCONVERT" --32 --mzXML *.raw
# 	MSCONVERT_RT_CODE=$?
# 	echo -E "msconvert returned $MSCONVERT_RT_CODE"
# 	if [ $MSCONVERT_RT_CODE -ne 0 ]; then
# 		RET_CODE="problem converting RAW files to mzXML"
# 	else
# 		rm -f *.raw
# 	fi
# }

# function cpConvertAndGenerateMQpar() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3
# 	FASTA_IN=$4
# 	RAW_FILE_FOLDER=$5
# 	PHOSPHO=$6
# 	PROTEASE=$7
# 	ALT_BASE_MQPAR=$8
# 	TMT_CORR_FACTORS_FILE=$9
# 	PEPTIDE_FDR=${10}
# 	PROTEIN_FDR=${11}
	
# 	generateMQpar $1 $2 $3 $4 $6 $7 $8 $9 ${10} ${11} "mzXML"
# 	cpRaws $1 $2 $3 $5
# 	convertRaws
# 	#cpCTRLs $1 $2 $3
	
# }

# function cpAndGenerateMQpar() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3
# 	FASTA_IN=$4
# 	RAW_FILE_FOLDER=$5
# 	PHOSPHO=$6
# 	PROTEASE=$7
# 	ALT_BASE_MQPAR=$8
# 	TMT_CORR_FACTORS_FILE=$9
# 	PEPTIDE_FDR=${10}
# 	PROTEIN_FDR=${11}
	
# 	generateMQpar $1 $2 $3 $4 $6 $7 $8 $9 ${10} ${11} "raw"
# 	cpRaws $1 $2 $3 $5
# 	#cpCTRLs $1 $2 $3
	
# }

# function fit() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3

# 	echo "fit data"
# 	#if [ -e "$RFOLDER/post_EC50_R/EC_fit_$BASE-auto_annotation.R" ]; then
# 	if [ "$BASE" = "TwoDose" ]; then
# 		"$RSCRIPT" --arch 64 --vanilla `cygpath -w "$RFOLDER/post_EC50_R/EC_fit_$BASE-auto_annotation.R"`
# 	else
# 		"$RSCRIPT" --arch 64 --vanilla `cygpath -w "$RFOLDER/post_EC50_R/EC_fit_$BASE.R"`
# 	fi
	
# }

# function prepareEC50() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3

# 	echo "prepare dose data"
# 	"$RSCRIPT" --arch 64 --vanilla `cygpath -w "$RFOLDER/post_EC50_R/EC_prepare.R"`

# }

# function fitEC50Curves() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3

# 	prepareEC50 $1 $2 $3
# 	fit $1 $2 $3
	
# }

# function finalizeFit() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3

# 	echo "finalizing output"
# 	fit $1 $2 $3
	
# }

# function finalizeSupp() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3

# 	echo "finalizing output"
# 	"$RSCRIPT" --vanilla "/media/msdata10/Maxquant3Queue/payloads/post_EC50_R/EC_convert_suppfile.R"
# 	RET_CODE=$?
# 	echo "done with Supp"
# }

# # Topas QC
# # Do NOT use this function as a post payload as it requires an extra positional argument (MODE)
# # Modes: FP = Full proteome, PP = Phospho, var = variable mod search
# function runQCTopas() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3
# 	MODE=$4
	
# 	REMOTE_DIR=$(echo "$PATH_PREFIX$BASE!$NAME" | sed -e "s#\!#\\\\#g" )
# 	RFILE=$MODE"_QC_pipeline.R"
	
# 	echo "finalizing output"
# 	# "$RSCRIPT" --vanilla "$RFOLDER/post_TOPAS_QC_R/$RFILE" $REMOTE_DIR
# 	"$RSCRIPT" --arch 64 `cygpath -w "$RFOLDER/post_TOPAS_QC_R/$RFILE"` $REMOTE_DIR
# 	RT_CODE=$?
# 	if [ $RT_CODE -ne 0 ]; then
# 	  RET_CODE+="; QC script failed with error code $RT_CODE"
# 	else
# 	  RET_CODE+="; QC script executed without error"
# 	fi
# 	echo "done with QC"
# }
