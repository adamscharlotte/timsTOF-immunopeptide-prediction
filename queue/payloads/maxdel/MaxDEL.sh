#!/usr/bin/env bash

echoLog() {
	echo $1 $2
}

run() {
	if [ $verbose == 1 ]; then
		echoLog -e "\t\t\t$1 $2"
	fi
	if [ $apply == 1 ]; then
		$1 "$2"
	fi
}

bzip2_file() {
	tar -cjf "$1.tar.bz2" "$1"
	retc=$(echo $?)
	if [ $retc == 0 ]; then
		run rm "$1"
	else
		if [ $verbose == 1 ]; then
			echoLog -e "\t\t\t\tCould not bzip2 '$1'. Return code $retc. Will not delete original file."
		fi
	fi
}

iterate_do() {
        objects=$1
        dtype="$2"
        action="$3"
	depth="$4"

	if [ "$objects" != "" ]; then
		if [ $verbose == 1 ]; then
	       	 	echoLog -e "\t\t$objects $dtype $action $depth"
		fi

	        IFS=',' read -a array <<< "$objects"
	        for element in "${array[@]}"
	        do
			if [ "$dtype" == "f" ]; then
				find . -maxdepth $4 -type $dtype -iname "*.$element" > $element-$dtype-$action.maxdellist
			else
				find . -maxdepth $4 -type $dtype -iname "$element" > $element-$dtype-$action.maxdellist
			fi
	                while read object
	                do
				if [ "$object" == "" ];then
					continue
				fi
				size=$(du -c "$object" | grep "total$" | awk -F " " '{print $1}')
	                        case "$action" in
					"bz2")	
						if [ $verbose == 1 ];then
							echoLog -e "\t\t\tbzip2 $object"
						else
							echoLog -n "+"
						fi
						while [ 1 ]; do
							rp=$(ps -aef | grep bz2 | grep tar | grep -v grep | wc -l)
							if [ $rp -lt $ncores ]; then
								AR_SIZE=$(($AR_SIZE + $size))
								if [ $apply == 1 ];then
									bzip2_file "$object" &
								fi
								break
							fi
							sleep 1
						done
						;;
					"rm")   
						RM_SIZE=$(($RM_SIZE + $size))
						run "rm -r" "$object"
						if [ $verbose == 0 ]; then
							echoLog -n "-"
						fi
	                                        ;;
	                                '?')    echoLog "do not know what to do"
	                                        ;;
	                        esac 
	                done < $element-$dtype-$action.maxdellist
			rm $element-$dtype-$action.maxdellist &> /dev/null
	        done
	fi
}

#size deleted
RM_SIZE=0
AR_SIZE=0

# Initialize our own variables:
verbose=0

apply=1
ncores=10

# complete
c_del_folder="proc,search,ps,ser"
c_del_file="ser,apl,res,apar,raw,Data"
c_zip_folder=""
c_zip_file=""
c_keep_file="txt"

# best practice
b_del_folder="proc,search,ps,ser"
b_del_file="ser,Data,apl,res,apar"
b_zip_folder=""
b_zip_file=""
b_keep_file="txt"

# best practice + remove raw
br_del_folder="proc,search,ps,ser"
br_del_file="ser,Data,apl,res,apar,raw"
br_zip_folder=""
br_zip_file=""
br_keep_file="txt"

# simple
#s_del_folder="proc,search,ps,ser"
s_del_folder="search,ps,ser"
s_del_file="ser,Data"
s_zip_folder=""
s_zip_file=""
s_keep_file="res,apl,apar,zip,bz2,txt"

# simple + remove raw
#s_del_folder="proc,search,ps,ser"
sr_del_folder="search,ps,ser"
sr_del_file="ser,Data,raw"
sr_zip_folder=""
sr_zip_file=""
sr_keep_file="res,apl,apar,zip,bz2,txt"

# archive
a_del_folder="proc,search,ps,ser"
a_del_file="ser,Data"
a_zip_folder=""
a_zip_file="res,apl,apar,txt"
a_keep_file="txt"

del_folder=$b_del_folder
del_file=$b_del_file
zip_folder=$b_zip_folder
zip_file=$b_zip_file
keep_file=$b_keep_file


show_help() {
cat << EOF
Usage: ${0##*/} [-hvl] [-c clean_level] [-d del_files] [-D del_folder] [-k keep_file] [-z zip_file] [-Z zip_folder] -f directory
Clean up script for Maxquant folders. Scripts scans directory for "combined" and "txt"
folders and deletes and zips specified files and folders.
    
    -h          display this help and exit
    -v          verbose mode.
    -f FOLDER   directory (folder) to scan
    -l          only list actions, do not run them
    -d PATTERN  delete files containing PATTERN (default: $del_file)
    -D PATTERN  delete directories containing PATTERN (default: $del_folder)
    -k PATTERN  keep files containing PATTERN (default: $keep_file)
    -z PATTERN  bzip2 files containing PATTERN (default: $zip_files)
    -Z PATTERN  bzip2 directories containing PATTERN (default: $zip_folder)
    -c level    clean folders according to levels:
                  "complete":
		    delete files: $c_del_file
		    delete directories: $c_del_folder
		    zip files: $c_zip_file
		    zip directories: $c_zip_folder
		    keep files: $c_keep_file
                  "bestpracticeraw":
		    delete files: $br_del_file
		    delete directories: $br_del_folder
		    zip files: $br_zip_file
		    zip directories: $br_zip_folder
		    keep files: $br_keep_file
                  "bestpractice":
		    delete files: $b_del_file
		    delete directories: $b_del_folder
		    zip files: $b_zip_file
		    zip directories: $b_zip_folder
		    keep files: $b_keep_file
                  "simple":
		    delete files: $s_del_file
		    delete directories: $s_del_folder
		    zip files: $s_zip_file
		    zip directories: $s_zip_folder
		    keep files: $s_keep_file
		              "simpleraw":
		    delete files: $sr_del_file
		    delete directories: $sr_del_folder
		    zip files: $sr_zip_file
		    zip directories: $sr_zip_folder
		    keep files: $sr_keep_file
                  "archive":
		    delete files: $a_del_file
		    delete directories: $a_del_folder
		    zip files: $a_zip_file
		    zip directories: $a_zip_folder
		    keep files: $a_keep_file
EOF
}


OPTIND=1 # Reset is necessary if getopts was used previously in the script.  It is a good idea to make this local in a function.
while getopts "hvlf:k:d:D:z:Z:c:" opt; do
    case "$opt" in
        h)
        	show_help
        	exit 0
        	;;
        v)	
		verbose=1
        	;;
        f)	
		input_dir=$OPTARG
        	;;
	l)	
		apply=0
		verbose=1
		;;
	d)	
		del_file=$OPTARG
		;;
	D)	
		del_folder=$OPTARG
		;;
	k)	
		keep_file=$OPTARG
		;;
	K)
		keep_folder=$OPTARG
		;;
	z)
		zip_file=$OPTARG
		;;
	Z)
		zip_folder=$OPTARG
		;;
	c)
		level=$OPTARG
		case "$level" in
			"bestpractice")	
				del_folder=$b_del_folder
				del_file=$b_del_file

				zip_folder=$b_zip_folder
				zip_file=$b_zip_file

				keep_file=$b_keep_file
				;;
			"bestpracticeraw")	
				del_folder=$br_del_folder
				del_file=$br_del_file

				zip_folder=$br_zip_folder
				zip_file=$br_zip_file

				keep_file=$br_keep_file
				;;
			"complete")
				del_folder=$c_del_folder
				del_file=$c_del_file

				zip_folder=$c_zip_folder
				zip_file=$c_zip_file

				keep_file=$c_keep_file
				;;
			"simple")
				del_folder=$s_del_folder
				del_file=$s_del_file

				zip_folder=$s_zip_folder
				zip_file=$s_zip_file

				keep_file=$s_keep_file
				;;
		  "simpleraw")
				del_folder=$sr_del_folder
				del_file=$sr_del_file

				zip_folder=$sr_zip_folder
				zip_file=$sr_zip_file

				keep_file=$sr_keep_file
				;;
			"archive")
				del_folder=$a_del_folder
				del_file=$a_del_file

				zip_folder=$a_zip_folder
				zip_file=$a_zip_file

				keep_file=$a_keep_file
				;;
			"?")
				show_help >&2
				exit 2
				;;
		esac
		;;
        '?')
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

if [ "$input_dir" == "" ]; then
	show_help >&2
	exit 0
fi

if [ ! -d "$input_dir" ]; then
	echo "${0##*/}: cannot find '$input_dir': not specified or does not exist or is not a folder" >&2
	exit -1
fi

old_pwd=$(pwd)
cd "$input_dir"

if [ $apply == 0 ]; then
	echoLog "list actions only"
fi
echoLog "scanning ..."
echoLog "-e" "\t$input_dir"

if [ $verbose == 1 ]; then
	stdbuf -o 0 find . -type d -name "combined" | tee combined_1.list
	stdbuf -o 0 find . -type d -name "txt" | sed -e "s/\/txt//g" | tee combined_2.list
else
	find . -type d -name "combined" > combined_1.list
	find . -type d -name "txt" | sed -e "s/\/txt//g" > combined_2.list
fi

cat combined_1.list combined_2.list | sort -u > combined.list
rm combined_1.list
rm combined_2.list
combc=$(cat combined.list | wc -l)

echoLog "found $combc folders to clean"

echoLog "start processing ..."

while read sub_dir;
do
        cd "$sub_dir"
        cd ..

	echo -en "\tcleaning $sub_dir"
	if [ $verbose == 1 ]; then
		echo ""
		echo "-e" "\t\trm raw folders"
	fi

	find . -maxdepth 2 -type d -name n0 | awk -F "/" '{print ( $(NF-1) )}' | uniq > to_del.list
	find . -maxdepth 1 -type f -name "*.index" >> to_del.list
	while read folder
	do
		if [ $verbose == 0 ]; then
			echo -n "-"
		fi
		run "rm -rf" "$folder"
	done < to_del.list
	rm to_del.list

	#delete files and folders in root of combined
	iterate_do "$del_file" "f" "rm" "1"

        cd - > /dev/null

	#delete files and folders
	iterate_do "$del_folder" "d" "rm" "1"
	iterate_do "$del_file" "f" "rm" "2"

	#archive files and folders
	iterate_do "$zip_folder" "d" "bz2" "1"
	iterate_do "$zip_file" "f" "bz2" "2"

	#waiting for zipping to finish
	while [ 1 ]; do
		rp=$(ps -aef | grep bz2 | grep tar | grep -v grep | wc -l)
		if [ $rp -eq 0 ]; then
			break;
		fi
	done

	#clean other files in combined
	if [ $verbose == 1 ]; then
		echoLog "-e" "\t\trm other files in combined"
	fi
	grep_string=$(echo $keep_file | sed -e "s/,/\\\\|/g")
	find . -maxdepth 2 -type f ! -name "*.*" | grep -v $grep_string > to_delete.list
	while read file; do
		if [ $verbose == 0 ]; then
			echo -n "-"
		fi
		size=$(du -c "$file" | grep "total$" | awk -F " " '{print $1}')
		RM_SIZE=$(($RM_SIZE + $size))
		run "rm" "$file"
	done < to_delete.list
	rm to_delete.list

	if [ $verbose == 1 ]; then
		echoLog -e "\t\tcleanup"
	fi
	iterate_do "maxdellist" "f" "rm" "2"
	rm maxdellist-f-rm.maxdellist &> /dev/null

	cd "$input_dir" 2> /dev/null

	echoLog

done < combined.list

cd "$input_dir"
rm combined.list

echoLog "... done"
if [ $apply == 1 ]; then
	echoLog "deleted $(($RM_SIZE/1024/1024))GB ($RM_SIZE kB)"
	echoLog "archive $(($AR_SIZE/1024/1024))GB ($AR_SIZE kB)"
else
	echoLog "would delete $(($RM_SIZE/1024/1024))GB ($RM_SIZE kB)"
	echoLog "would archive $(($AR_SIZE/1024/1024))GB ($AR_SIZE kB)"
fi

cd "$old_pwd"

exit 0

