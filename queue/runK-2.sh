#!/bin/bash

CONFIG_FILE="runK-$HOSTNAME.cfg"
source $CONFIG_FILE

if [ ! -d "$LOG" ]; then
	mkdir "$LOG/"
fi

function echoLog() {
	if [ "$VERBOSE" != "" ]; then
		if [ $VERBOSE -le $1 ]; then
			ECHOTIME=$(date +"%x %X")
			echo -E "$ECHOTIME: $2"
		fi
	fi
}

function run() {

	FULL=$(echo $1)
	NAME=$(echo $1 | awk -F "," '{print $1}')
	VERSION=$(echo $1 | awk -F "," '{print $2}' | tr -d " ")
	PRE_PAYLOAD=$(echo $1 | awk -F "," '{print $3}' | tr -d " ")
	POST_PAYLOAD=$(echo $1 | awk -F "," '{print $4}' | tr -d " ")
	THREADS=$(echo $1 | awk -F "," '{print $5}' | tr -d " ")
	BASE=$(echo $1 | awk -F "," '{print $6}' | tr -d " ")
	FASTA_FILE=$(echo $1 | awk -F "," '{print $7}' | tr -d " ")
	RAW_FOLDER=$(echo $1 | awk -F "," '{print $8}' | tr -d " ")
	MAX_LENGTH=$(echo $1 | awk -F "," '{print $9}' | tr -d " ")
	# PHOSPHO=$(echo $1 | awk -F "," '{print $9}' | tr -d " ")
	# PROTEASE=$(echo $1 | awk -F "," '{print $10}' | tr -d " \r\n")
	# ALT_BASE_MQPAR=$(echo $1 | awk -F "," '{print $11}' | tr -d " \r\n")
	# TMT_CORR_FACTORS_FILE=$(echo $1 | awk -F "," '{print $12}' | tr -d " \r\n")
	# PEPTIDE_FDR=$(echo $1 | awk -F "," '{print $13}' | tr -d " \r\n")
	# PROTEIN_FDR=$(echo $1 | awk -F "," '{print $14}' | tr -d " \r\n")
	TIME=$2

# PRE_PAYLOAD="cpAndGenerateMQpar"
# NAME="TUM_lysn_10"
# BASE="20220621_AspNlysN"
# THREADS="1"
# FASTA_FILE="/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/fasta/TUM_lysn_10.fasta"
# RAW_FOLDER="/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/raw-library-run/20220621_AspNlysN/AspNLysN-A10_S1-A10_1_6551.d"
# MAX_LENGTH="16"

#   if [ -z $PHOSPHO ]; then
#     PHOSPHO=0
#   fi
  
#   if [ -z $PROTEASE ]; then
#     PROTEASE="Trypsin"
#   fi
  
#   if [ -z $ALT_BASE_MQPAR ]; then
#     ALT_BASE_MQPAR=$MQPAR
#   fi
  
#   if [ -z $TMT_CORR_FACTORS_FILE ]; then
#     TMT_CORR_FACTORS_FILE="no_factors.txt"
#   fi
  
#   if [ -z $PEPTIDE_FDR ]; then
#     PEPTIDE_FDR=0.01
#   fi
  
#   if [ -z $PROTEIN_FDR ]; then
#     PROTEIN_FDR=1
#   fi
  
	CWD=$(pwd)
	
	echo -E "$NAME-$BASE,$THREADS,$TIME-$NAME.log,processing" >> "$QUEUE_DONE_FILE"
	echo $THREADS > $LOG/$TIME-$NAME.th

	source $CONFIG_FILE
	source $ADDITIONAL_PAYLOADS

	RT_CODE=0
	OTHREADS=$(echo $THREADS)
	
	changeDirectory $NAME $BASE $THREADS
	
	if [ "$PRE_PAYLOAD" != "" ]; then
		# echo -E "starting pre-payload $PRE_PAYLOAD $NAME $BASE $THREADS $FASTA_FILE $RAW_FOLDER $PHOSPHO $PROTEASE $ALT_BASE_MQPAR $TMT_CORR_FACTORS_FILE $PEPTIDE_FDR $PROTEIN_FDR"
		# $PRE_PAYLOAD $NAME $BASE $THREADS $FASTA_FILE $RAW_FOLDER $PHOSPHO $PROTEASE $ALT_BASE_MQPAR $TMT_CORR_FACTORS_FILE $PEPTIDE_FDR $PROTEIN_FDR
		echo -E "starting pre-payload $PRE_PAYLOAD $NAME $BASE $THREADS $FASTA_FILE $RAW_FOLDER $MAX_LENGTH"
		$PRE_PAYLOAD $NAME $BASE $THREADS $FASTA_FILE $RAW_FOLDER $MAX_LENGTH
		RT_CODE=$?
		echo -E "pre payload returned $RT_CODE"
		if [ $RT_CODE -eq 1 ]; then
			RET_CODE="empty fasta file"
		fi
	fi

	if [ $RT_CODE -eq 0 ]; then
	
		UUPATH=$(pwd)
		echo -E "Attempting to start $NAME with $THREADS threads in $MQPAR with version $VERSION ($PRE_PAYLOAD, $POST_PAYLOAD)"
		echo -E "Sitting in $UUPATH"
		
		# case $VERSION in
		# 	"1.5.3.8"|"1.6.12.0")
		# 		if [ "$LOCAL_PROCESSING" = "1" ]; then
		# 			moveToLocal $NAME $BASE $THREADS
		# 		fi

		# 		cleanSearch $NAME $BASE $THREADS
		# 		echo "`cygpath -u "$MAXQUANT_PREFIX$VERSION$MAXQUANT_SUFFIX"` mqpar.xml"
		# 		if [ "$DRY" == "0" ]; then
		# 			`cygpath -u "$MAXQUANT_PREFIX$VERSION$MAXQUANT_SUFFIX"` mqpar.xml 2>&1
		# 		fi

		# 		if [ "$LOCAL_PROCESSING" = "1" ]; then
		# 			moveToRemote $NAME $BASE $THREADS
		# 		fi

		# 		cleanPostSearch $NAME $BASE $THREADS
				
		# 		;;
		# 	*)
		# 		echo "Do not know maxquant version $VERSION. MQ not running."
		# 		;;
		# esac

		if [ -d "./combined" ]; then
			if [ -d "./combined/txt" ]; then
				if [ -s "./combined/txt/evidence.txt" ]; then
					RET_CODE="looks fine"
				else
					RET_CODE="evidence does not exist"
				fi
			else
				RET_CODE="txt does not exist"
			fi
		else
			RET_CODE="combined does not exist"
		fi
		
		if [ "$POST_PAYLOAD" != "" ]; then
			echo -E "starting post-payload $POST_PAYLOAD"
			
			source $ADDITIONAL_PAYLOADS
			
			$POST_PAYLOAD $NAME $BASE $THREADS
			RT_CODE=$?
			echo -E "post payload returned $RT_CODE"
		fi
		
	fi

	cd "$CWD" > /dev/null
	echo -E "Sitting in $CWD"
	
	rm $LOG/$TIME-$NAME.th
	echoLog 1 "End $FULL with return code $RET_CODE"

	echo "sed -i -e "s/$NAME-$BASE,$THREADS,$TIME-$NAME.log,processing/$NAME-$BASE,$THREADS,$TIME-$NAME.log,$RET_CODE/g" "$QUEUE_DONE_FILE""
	sed -i -e "s/$NAME-$BASE,$THREADS,$TIME-$NAME.log,processing/$NAME-$BASE,$THREADS,$TIME-$NAME.log,$RET_CODE/g" "$QUEUE_DONE_FILE"

}

function process() {

	NAME=$(echo $1 | awk -F "," '{print $1}')
	TIME=$(date +"%Y%m%d_%H%M%S")

	run "$1" "$TIME"  &> "$LOG/$TIME-$NAME.log" &
	
}

function runNext() {
	FIRST="0"
	#cat $QUEUE_FILE | tr "\\" "|" > $QUEUE_FILE-
	while IFS='' read -r line || [[ -n "$line" ]]; do
		CCORES=$(runningInstances)
		LEFTCORES=$(($NCORES - $CCORES))
		if [ "$FIRST" -eq "0" ]; then
			echoLog 0 "skipping first line in $QUEUE_FILE"
			FIRST="1"
			continue
		fi
		if [ "$line" == "" ]; then
			#skipping empty lines
			continue
		fi
		if [ -s $QUEUE_DONE_FILE ]; then
			NAME=$(echo $line | awk -F "," '{print $1}')
			BASE=$(echo $line | awk -F "," '{print $6}' | tr -d " \r\n")
			FOUND=$(grep "^$NAME-$BASE" $QUEUE_DONE_FILE | wc -l)
			if [ $FOUND -eq 0 ]; then
				TC=$(echo $line | tr -d " " | awk -F "," '{print $5}')
				if [ "$TC" == "" ] || [ "$1" == "" ]; then
					echoLog 1 "TC: $TC or 1: $1 is empty. Skipping line: $line"
					continue
				fi
				echoLog 0 "line: $line"
				echoLog 0 "TC: $TC"
				echoLog 0 "1: $1"
				if [ $TC -le $LEFTCORES ]; then
					echoLog 1 "$NAME-$BASE with [$line] not processed yet. Starting."
					echoLog 1 "Available now $1-$TC"
					process "$line"
					sleep $SLEEP_TIME
					#return 0
				else
					if [ "$INORDER" == "1" ]; then
						return 2
					fi
					echoLog 0 "Not enough resources for $NAME-$BASE with $TC threads (available $1). Skipping for now."
				fi
			fi
		else
			echoLog 1 "Cannot fine $QUEUE_DONE_FILE"
		fi
	done < $QUEUE_FILE
	return 1
}

function runningInstances() {
	FM=$(ls $LOG/ | grep th$ | wc -l)
	if [ $FM -gt 0 ]; then
		cat `ls $LOG/*.th` | awk '{ sum += $1 } END { print sum }'
	else
		echo 0
	fi
}

function payload() {
while [ true ]; do
	checkforterm
		
	source $CONFIG_FILE
	
	# check how many are running, stop if NCORES
	CCORES=$(runningInstances)
	if [ "$CCORES" -ge "$NCORES" ]; then
		echoLog 0 "$CCORES/$NCORES occupied. Waiting ..."
	else
		# run NCORS-CCORES procs
		while [ $CCORES -lt $NCORES ]; do
			echoLog 0 "Cores $CCORES/$NCORES"
			runNext $(($NCORES - $CCORES))
			RET_CODE=$?
			CCORES=$(runningInstances)
			echoLog 0 "$CCORES/$NCORES occupied."
			if [ $RET_CODE -eq 1 ]; then
				echoLog 0 "Nothing to do. Waiting ..."
				break
			fi
			if [ $RET_CODE -eq 2 ]; then
				echoLog 1 "Not enough resources for next tasks in line. Waiting ..."
				break
			fi
			
			source $CONFIG_FILE
			
		done
	fi

	sleep $CYCLE_TIME
done
}

source $DAEMONIZER
