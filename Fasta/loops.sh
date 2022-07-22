name_file=/Users/adams/Projects/300K/2022-library-run/metadata/poolnames/All.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript Fasta/create_fasta.R "${array[0]}"
done

name_file=/Users/adams/Projects/300K/2022-library-run/metadata/poolnames/AspN_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript Fasta/create_fasta.R "${array[0]}"
done

name_file=/Users/adams/Projects/300K/2022-library-run/metadata/poolnames/HLA-I_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript Fasta/create_fasta.R "${array[0]}"
done

name_file=/Users/adams/Projects/300K/2022-library-run/metadata/poolnames/LysN_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript Fasta/create_fasta.R "${array[0]}"
done
