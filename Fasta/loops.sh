name_file=/Users/adams/Projects/300K/2022-library-run/metadata/poolnames/HLA-II_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript create_fasta.R "${array[0]}" "HLA-II"
done

name_file=/Users/adams/Projects/300K/2022-library-run/metadata/poolnames/AspN_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript create_fasta.R "${array[0]}" "AspN"
done

name_file=/Users/adams/Projects/300K/2022-library-run/metadata/poolnames/HLA-I_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript create_fasta.R "${array[0]}" "HLA-I"
done

name_file=/Users/adams/Projects/300K/2022-library-run/metadata/poolnames/LysN_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    Rscript create_fasta.R "${array[0]}" "LysN"
done
