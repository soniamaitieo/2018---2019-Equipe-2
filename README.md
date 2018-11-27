# MEET-U Upstream: team n°2


[![Foo](https://i.imgur.com/b4GC6hM.png)](http://google.com.au/)

### Our method

We are team n°2 from Paris Diderot. Here, we decide to make a profile-profile alignment in order to find a template which matches with a query. Our profile is composed by 20 amino-acids, gap and secondary structure informations. The alignment is a semi-global alignment and we made a dot-product score and a Pearson's correlation score. We will use those differents scores for our benchmark.


If you have any questions or problems, feel free to create an issue :)

### Prerequisites

All the software and database you need :

 **Downloading database (ex: Uniref50) :**

[Dowload database here](https://www.uniprot.org/downloads)


**Installing PSI-BLAST:**

 - Install Legacy Blast : PSIBLAST 2.2.31+  
`sudo apt-get install ncbi-blast+-legacy`
or ` conda install -c bioconda blast-legacy`


Sometimes legacy blast need a perl library. Check if it needs it with running legacy blast



 - Formate database Uniref50  
`formatdb -i uniref50.fasta -p T -o T`

**MAFFT :**

- [Dowload MAFFT here](https://mafft.cbrc.jp/alignment/software/linux.html)  
OR
- Use the program in */prog* .

**Python 3 :**

- [Dowload Python here](https://www.python.org/downloads/)  

**Python library :**

``` ruby
pip install numpy
pip install pandas
pip install sys
pip install os
pip install glob
```

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
