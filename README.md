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

**PSIPRED :**

We give you the binary file of psipred, so you don't need to dowload it. But you need to download csh

`sudo apt-get install csh`

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

# Team 4 dependancies :
pip install docopt
pip install numpy
pip install biopython
pip install pandas
pip install schema
pip install tqdm
pip install matplotlib
pip install m2r

```


**Team 4 dependancies :**


`sudo apt-get install dssp`

`conda install -c salilab modeller`

### How to use

**Clone the repository :**

``` ruby
git clone https://github.com/meetU-MasterStudents/2018---2019-Equipe-2.git
```


**/!\ CHANGE DIRECTORY FOR DATABASE: in *pipeline_final.sh***

``` ruby
DIR_BD=/home/sdv/m2bi/dde_murat/MEET_U/BDD/uniref50.fasta
```

**/!\ CHANGE DIRECTORY FOR DATABASE: in *src/align_pssm.py***

``` ruby
/home/madeleine/Documents/2018---2019-partage/Data/HOMSTRAD/
```

**Lauch from directory : *2018---2019-Equipe-2***
``` ruby
./pipeline_final.sh Query.fasta OUTPUT_FILE
```

**Full pipeline UPSTREAM and DOWNSTREAM (launch from directory : 2018---2019-Equipe-2)**
``` ruby
./pipeline_with_aval.sh Query.fasta OUTPUT_FILE
```

**Example**

``` ruby
./pipeline_with_aval.sh example/Cohesin.fasta Cohesin_output
```

## Authors

* **Daniel DE MURAT** - *M2BI - Université Paris Diderot* - [DanielBI](https://github.com/DanielBI)

* **Madeleine DE SOUSA-VIOLANTE** - *M2BI - Université Paris Diderot* - [madeleinevlt](https://github.com/madeleinevlt)

* **Meishiue KUO** - *M2BI - Université Paris Diderot* - [meishiue](https://github.com/meishiue)

* **Sonia TIEO** - *M2BI - Université Paris Diderot* - [soniamaitieo](https://github.com/soniamaitieo)
