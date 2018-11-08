# MEET-U Amont: team n°2

One Paragraph of project description goes here


### Prerequisites

All the software and database you need :

 **Downloading database (ex: Uniref50) :**
[Dowload here](https://www.uniprot.org/downloads)


**Installing PSI-BLAST:**

 - Install Legacy Blast : PSIBLAST 2.2.31+  
`sudo apt-get install ncbi-blast+-legacy`
or ` conda install -c bioconda blast-legacy`

 - Formate database Uniref50  
`formatdb -i uniref50.fasta -p T -o T`

**MAFFT :**

- [Dowload here](https://mafft.cbrc.jp/alignment/software/linux.html)  
OR
- Use the program in */prog* .


### How to use

**Clone the repository :**

``` ruby
git clone https://github.com/meetU-MasterStudents/2018---2019-Equipe-2.git
```


**/!\ CHANGE DIRECTORY FOR DATABASE: in *pipeline_final.sh***

``` ruby
##### A changer en fonction de l'emplacement de uniref50.fasta ou autre base de donnees #######
DIR_BD=/home/sdv/m2bi/dde_murat/MEET_U/BDD/uniref50.fasta
```

**Lauch from directory : *2018---2019-Equipe-2***
``` ruby
./pipeline_final.sh Query.fasta
```

## Authors

* **Daniel DE MURAT** - *M2BI - Université Paris Diderot* - [DanielBI](https://github.com/DanielBI)

* **Madeleine DE SOUSA-VIOLANTE** - *M2BI - Université Paris Diderot* - [madeleinevlt](https://github.com/madeleinevlt)

* **Meishiue KUO** - *M2BI - Université Paris Diderot* - [meishiue](https://github.com/meishiue)

* **Sonia TIEO** - *M2BI - Université Paris Diderot* - [soniamaitieo](https://github.com/soniamaitieo)
