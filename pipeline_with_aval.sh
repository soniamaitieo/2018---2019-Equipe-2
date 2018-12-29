#!/bin/bash
FASTA=$1
OUTPUT_file=$2

b=$(basename $FASTA)
name_fasta=${b%.*}

echo "LAUNCH UPSTREAM CODE"
./pipeline_final.sh $FASTA $OUTPUT_file

echo "LAUNCH DOWNSTREAM CODE"
cd code_aval
./script/run_ccmpred.py --folder ../$OUTPUT_file/aln
mv data/aln/clustal/* ../$OUTPUT_file
mv data/ccmpred/* ../$OUTPUT_file
./fold_u ../$OUTPUT_file/$name_fasta.foldrec ../$OUTPUT_file/$name_fasta.clustal \
../$OUTPUT_file/$name_fasta.mat -o ../$OUTPUT_file
