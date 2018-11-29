
for f in /home/soniamai/Bureau/2018---2019-partage/new_data/HOMSTRAD_extended/*/*.map;
do
b=$(basename $f)
./runpsipred_single /home/soniamai/Bureau/2018---2019-partage/new_data/HOMSTRAD_extended/*/${b%.*}.fasta
rm *.horiz
rm *.ss
mv *.ss2 /home/soniamai/Bureau/psipred/
cd  /home/soniamai/Bureau/2018---2019-partage/new_data/HOMSTRAD_extended/${b%.*}/
python3 /home/soniamai/Bureau/2018---2019-Equipe-2/src/template.py ${b%.*}.map /home/soniamai/Bureau/psipred/${b%.*}.ss2 /home/soniamai/Bureau/pssm_templates_new/${b%.*}.aamtx
cd /home/soniamai/Bureau/2018---2019-Equipe-2/src
echo $b
done;
