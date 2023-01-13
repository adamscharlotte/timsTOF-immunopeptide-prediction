ssh cadams@10.152.135.57

scp cadams@10.152.135.57:/home/cadams
cd /media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/TUM-Searches

mkdir accumulatedMsmsScans-unsp

name_file=/home/cadams/tryptic_pool_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=',' read -r -a array <<< "$line"
    cp first_pool-unsp/${array[0]}/combined/txt/accumulatedMsmsScans.txt accumulatedMsmsScans-unsp/"${array[0]}".txt
done

cd /Users/adams/Projects/300K/2022-library-run/
scp -r cadams@10.152.135.57:/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/TUM-Searches/accumulatedMsmsScans-unsp .
# -----------------------------------------------------------------------------
# All msms.txt files
cd /media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K
mkdir msms-txt

name_file=/home/cadams/tryptic_pool_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=',' read -r -a array <<< "$line"
    cp TUM-Searches/first_pool-unsp/${array[0]}/combined/txt/msms.txt msms-txt/"${array[0]}".txt
done

name_file=/home/cadams/tryptic_pool_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=',' read -r -a array <<< "$line"
    cp TUM-Searches/first_pool-unsp/${array[0]}/combined/txt/msms.txt msms-txt/"${array[0]}".txt
done

scp Searches/pool_plate.txt /home/cadams/
name_file=/home/cadams/pool_plate.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=',' read -r -a array <<< "$line"
    cp Searches/"${array[1]}"/"${array[0]}"/combined/txt/msms.txt msms-txt/"${array[0]}".txt
done

cd /Users/adams/Projects/300K/2022-library-run
scp -r cadams@10.152.135.57:/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/msms-txt .

# -----------------------------------------------------------------------------
cd /media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Searches

name_file=pool_plate.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=',' read -r -a array <<< "$line"
    cp "${array[1]}"/"${array[0]}"/combined/txt/pasefMsmsScans.txt pasefMsmsScans/"${array[0]}".txt
done
