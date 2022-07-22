cd /Users/adams/Projects/300K/PXD021013/fasta

while read line; do echo $line; done < /Users/adams/Projects/300K/PXD021013/fasta/TUM_HLA_138_sequences.txt

name_file=/Users/adams/Projects/300K/PXD021013/fasta/TUM_HLA_138_sequences.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
  echo "$line"
done


name_file=/Users/adams/Projects/300K/PXD021013/fasta/TUM_HLA_138_sequences.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    echo ">sp|${array[1]}|${array[0]}" >> TUM_HLA_138.fasta
    echo "${array[1]}" >> TUM_HLA_138.fasta
done

