
CURRDIR = $(pwd)

for f in $CURRDIR/HOMSTRAD/*/*.map;
do
b=$(basename $f)
./template.py $f /home/sdv/m2bi/stieo/Bureau/pssm_templates/${b%.*}.aatmx
echo $b

done; 
