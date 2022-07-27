#!/usr/bin/env bash

function cleanPostSearch() {

	NAME=$1
	BASE=$2
	THREADS=$3

	DIR="$DESTINATION/$BASE/$NAME"
	"$PAYLOAD_FOLDER/maxdel/MaxDEL.sh" -c "simple" -f "$DIR" # keeps andromeda search results, keeps raw files
	
}

function generateMQpar() {

	NAME=$1
	BASE=$2
	THREADS=$3
	FASTA_FILE=$4
	RAW_FOLDER=$5
	MAX_LENGTH=$6
	BASE_MQPAR=$7
	OUTPUT_FOLDER=$8

	ROOT="/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K"
	DESTINATION="$ROOT/Searches"
	ACTIVE_BASE_MQPAR="$ROOT/Queue/Parameters/$BASE_MQPAR"

	DIR="$DESTINATION/$BASE/$NAME"
	LOC="$DIR/mqpar.xml"

	if [ ! -d "$DIR" ]; then
		echo -e "\t\tCreating folder for $NAME"
		mkdir $DIR
	fi

	echo -e "\tGenerating $1"

	if [ -e $LOC ]; then
		echo -e "\t\tFound old file, saving ..."
		mv $LOC $LOC`date +"_%y%m%d-%H%M%S"`
	fi

	# MSPLATE=$(grep $NAME $MAPPING | awk -F "," -v x=$MSPLATE_POSITION '{print $x}' | uniq)
	# LMAPPING="localMapping.txt"

	# #add real files
	# grep $MSPLATE $MAPPING | grep $NAME | awk -F "," '{print $1".d,"$2","$3","$4""}' > $LMAPPING
	# grep $MSPLATE $MAPPING | grep ",$NAME," | awk -v EXT="$RAW_FILE_EXT" -F "," '{print $1"."d","$2","$3","$4","$5","$6","$7","$8}' > $LMAPPING
	# grep $MSPLATE $MAPPING | grep ",DMSO," | awk -v EXT="$RAW_FILE_EXT" -F "," '{print $1"."EXT","$2","$2","$4","$5","$6","$7","$8}' >> $LMAPPING
	
	# temporarily replace newlines with "^^" (the record separator, which 
	# normally should not occur in text files) for sed to work
	sed -e "s#%%RAW_FOLDER%%#$RAW_FOLDER#" $ACTIVE_BASE_MQPAR |
		sed -e "s#%%FASTA_FILE%%#$FASTA_FILE#" |
		sed -e "s#%%MAX_LENGTH%%#$MAX_LENGTH#" |
		sed -e "s#%%OUTPUT_FOLDER%%#$OUTPUT_FOLDER#" |
		sed -e "s#\!#\\\\#g" > $LOC

	# sed -e "s#%%RAW_FOLDER%%#`cat $LMAPPING |
	# 	awk -F "," -v x=$RAW_FOLDER_POSITION '{print $x}' |
	# 	sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
	# 	tr "\n" "|"`#g" |
	# 	sed -e "s#%%FASTA_FILE%%#`cat $LMAPPING |
	# 		awk -F "," -v x=$FASTA_FILE_POSITION '{print $x}' |
	# 		sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
	# 		tr "\n" "|"`#g" |
	# 	sed -e "s#%%MAX_LENGT%%#`cat $LMAPPING |
	# 		awk -F "," -v x=$MAX_LENGT_POSITION '{print $x}' |
	# 		sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
	# 		tr "\n" "|"`#g" |
	# 	sed -e "s#%%OUTPUT_FOLDER%%#`cat $LMAPPING |
	# 		awk -F "," -v x=$OUTPUT_FOLDER_POSITION '{print $x}' |
	# 		sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
	# 		tr "\n" "|"`#g" |

	unix2dos -q $LOC

	diff $ACTIVE_BASE_MQPAR $LOC | sed -e "s/^/\t\t/g"
}


# function cpRaws() {

# 	NAME=$1
# 	BASE=$2
# 	RAW_FOLDER=$3

# 	LMAPPING="localMapping.txt"

# 	for file in `cat $LMAPPING | grep ",$NAME," | awk -F ',' -v x=$RAW_FOLDER_POSITION '{print $x}'`; do
# 		if [ ! -e `basename $file` ]; then
# 			# cp -v $RAW_LOCATION/$file .
# 			cp -v $RAW_FOLDER/$file .
# 		else
# 			echo -E "$file already located in here"
# 		fi
# 	done
# }


# function cpAndGenerateMQpar() {

# 	NAME=$1
# 	BASE=$2
# 	THREADS=$3
# 	FASTA_FILE=$4
# 	RAW_FOLDER=$5
# 	MAX_LENGTH=$6
	
# 	generateMQpar $1 $2 $3 $4 $5 $6
# 	cpRaws $1 $2 $5
	
# }

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
