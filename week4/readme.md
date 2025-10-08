# Solution to Assignment 4


## 1. Download Data
Identify the accession numbers and write the commands to download the data.


Accession number: GCF_000882815.3
```
micromamba activate bioinfo
datasets download genome accession GCF_000882815.3 --include genome,gff3,gtf
unzip ncbi_dataset.zip
cd ncbi_dataset.zip/data/GCF_000882815.3
```


## 2. Visualize the genome and annotations
How big is the genome?
How many features does the GFF file contain?
What is the longest gene?
What is its name and function?


## 4. Pick another gene and describe its name and function.
Are the features closely packed; is there a lot of intragenomic space?
Estimate how much of the genome is covered by coding sequences.


## 5. Find alternative genome builds that could be used to answer a different question.
Include their accession numbers.
What other questions could be answered using a different genome build?