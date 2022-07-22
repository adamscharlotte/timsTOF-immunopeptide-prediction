#!/usr/bin/env bash

source vars.sh

NAME=$1
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

MSPLATE=$(grep $NAME $MAPPING | awk -F "," -v x=$MSPLATE_POSITION '{print $x}' | uniq)

#add real files
grep $MSPLATE $MAPPING | grep $NAME | awk -F "," '{print $1".d,"$2","$3","$4""}' > $MAPPING-$NAME
# grep $MSPLATE $MAPPING | grep DMSO | awk -F "," '{print $1".raw,"$2","$2","$4",1,1"}' >> $MAPPING-$NAME

# #add ctrl files
# ls $DMSO_FOLDER | awk -F "," -v x=$MSPLATE '{print $1","substr($1, 1, length($1)-4)","substr($1, 1, length($1)-4)","x",1,0"}' >> $MAPPING-$NAME

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


# sed -e "s#%%RAW_FILES%%#`cat $MAPPING-$NAME |
# 		awk -F "," -v x=$RAWFILE_POSITION '{print $x}' |
# 		sed -e "s/^/      \<string\>$PATH_PREFIX_TWODOSE\!$NAME\!/g" |
# 		sed -e "s/$/\<\/string\>/g" | tr "\n" "|"`#g" $BASE_MQPAR |
# 	sed -e "s#%%EXPERIMENT_NAMES%%#`cat $MAPPING-$NAME |
# 	# 	awk -F "," -v x=$CONDITION_POSITION '{print $x}' |
# 	# 	sed -e "s/^/      \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
# 	# 	tr "\n" "|"`#g" |
# 	# sed -e "s#%%FRACTION_NUMBER%%#`cat $MAPPING-$NAME |
# 	# 	awk -F "," -v x=$FRACTION_POSITION '{print $x}' |
# 	# 	sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
# 	# 	tr "\n" "|"`#g" |
# 	# sed -e "s#%%PARAMETER_GROUP%%#`cat $MAPPING-$NAME |
# 	# 	awk -F "," -v x=$PARAMETERGROUP_POSITION '{print $x}' |
# 	# 	sed -e "s/^/     \<string\>/g" | sed -e "s/$/\<\/string\>/g" |
# 	# 	tr "\n" "|"`#g" |
# 	sed -e "s#%%FASTA_FILE%%#$FASTA_FILE_ROOT#" |
# 	tr "|" "\n" |
# 	sed -e "s#\!#/#g" > $LOC

unix2dos -q $LOC

diff $BASE_MQPAR $LOC | sed -e "s/^/\t\t/g"

#rm $MAPPING-$NAME
