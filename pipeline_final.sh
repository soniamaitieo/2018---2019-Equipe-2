#!/bin/bash
#formatdb -i /home/sdv/m2bi/dde_murat/MEET_U/BDD/uniref50.fasta -p T -o T

##### A changer en fonction de l'emplacement de uniref50.fasta ou autre base de donnees #######
DIR_BD=/home/madeleine/BD_90/uniref90.fasta
DIR_CURRENT=$(pwd)
####################################################

FASTA=$1
b=$(basename $FASTA)
name_fasta=${b%.*}
OUTPUT_file=$2
mkdir -p $OUTPUT_file

legacy_blast blastpgp -d $DIR_BD -i $FASTA -j 3 -C ff.chd.ckp -o $name_fasta.blast_out;
$DIR_CURRENT/src/parser_psiblast_output.py $name_fasta.blast_out $name_fasta.parseout
cat $FASTA $name_fasta.parseout > $name_fasta"2.parseout"
$DIR_CURRENT/bin/mafft-linux64/mafft.bat $name_fasta"2.parseout" > $name_fasta.aln
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' > $name_fasta.mfasta $name_fasta.aln

$DIR_CURRENT/src/runpsipred_single $1

$DIR_CURRENT/src/query_ss2.py $name_fasta.mfasta $name_fasta.ss2 $name_fasta.aamtx
python3 $DIR_CURRENT/src/align_pssm.py $name_fasta.aamtx
cat $name_fasta.foldrec1 $name_fasta.foldrec2 >> $name_fasta.foldrec
rm $name_fasta.foldrec1 $name_fasta.foldrec2
mv $name_fasta* $OUTPUT_file
mv output.ss2 $OUTPUT_file
mv ff.chd.ckp $OUTPUT_file
mkdir -p $OUTPUT_file/aln
mv $OUTPUT_file/$name_fasta.aln $OUTPUT_file/aln
