ssh cadams@10.152.135.57

scp cadams@10.152.135.57:/home/cadams
cd /media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/TUM-Searches

mkdir accumulatedMsmsScans

name_file=pool_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=',' read -r -a array <<< "$line"
    cp ${array[0]}/combined/txt/accumulatedMsmsScans.txt accumulatedMsmsScans/"${array[0]}".txt
done

scp cadams@10.152.135.57:/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Searches/accumulatedMsmsScans/*.txt .

# -----------------------------------------------------------------------------
cd /media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Searches

name_file=pool_plate.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=',' read -r -a array <<< "$line"
    cp "${array[1]}"/"${array[0]}"/combined/txt/pasefMsmsScans.txt pasefMsmsScans/"${array[0]}".txt
done
