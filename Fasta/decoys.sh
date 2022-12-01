
name_file=all-pool-names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    name_file=/Users/adams/Projects/300K/2022-library-run/fasta/"${array[0]}".fasta
    lines=`tail -n+1 $name_file`
    for line in $lines
        do echo $line >> /Users/adams/Projects/300K/2022-library-run/fasta-decoy/"${array[0]}".fasta
    done
    for line in $lines
    do
        if [[ "$line" == ">"* ]]; then
            head="${line:1}"
            revhead=">rev_$head"
            echo $revhead >> /Users/adams/Projects/300K/2022-library-run/fasta-decoy/"${array[0]}".fasta
        else
            echo $line | rev >> /Users/adams/Projects/300K/2022-library-run/fasta-decoy/"${array[0]}".fasta
        fi
    done
done
