### parser_psiblast_output.py

INPUT : OUTPUT PSI-BLAST of query  
OUTPUT : multi-fasta file with all match sequences with query  

> **Command-line PSI-BLAST:**

> - Install Legacy Blast : PSIBLAST 2.2.31+  
`sudo apt-get install ncbi-blast+-legacy`
OR
` conda install -c bioconda blast-legacy`


> - Formate database Uniref50  
`formatdb -i uniref50.fasta -p T -o T`

> - PSIBLAST (3 rounds)  
`legacy_blast blastpgp -d uniref50.fasta -i Cohesin.fasta -o Cohesin.out -j 3 -C ff.chd.ckp`
