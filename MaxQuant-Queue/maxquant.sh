ssh cadams@10.152.135.57

scp -r /Users/adams/Projects/300K/PXD038782-comparison cadams@10.152.135.57:/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker

cd /media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison

for f in *.tar; do tar -m --no-overwrite-dir -xvf "$f"; done

cd /media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/d-folder
mkdir tar
mv *.d.tar tar/
mv mnt/f/Pride_Upload/umgenannt/*/*.d .
