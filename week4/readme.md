# Solution to Assignment 4

## 1. Download Data
Identify the accession numbers and write the commands to download the data.

Accession number: GCF_000882815.3

To download the genome data, I used the following commands:
```
micromamba activate bioinfo
datasets download genome accession GCF_000882815.3 --include genome,gff3,gtf
unzip ncbi_dataset.zip
cd ncbi_dataset/data/GCF_000882815.3/

# To ensure the files were downloaded, I checked the contents of the directory
ls
```

## 2. Visualize the genome and annotations.
A. How big is the genome?
The genome is 10,794 nucleotides long. 

B. How many features does the GFF file contain?
I used the following commands to give me the features and how many there are: 
```
cat genomic.gff | grep -v '#' > genomic.gff3
cat genomic.gff3 | cut -f 3 | sort-uniq-count-rank
```
Which gave me the following output:
```
14      mature_protein_region_of_CDS
1       CDS
1       five_prime_UTR
1       gene
1       region
1       three_prime_UTR
```

C. What is the longest gene?
The longest gene is id-YP_002790881.1:2517..3419.

D. What is its name and function?
This gene codes for RNA-dependent RNA polymerase (NS5) which is the protein that connects amino acids to form a polypeptide chain and thus form protein.

E. Pick another gene and describe its name and function.
Another gene is id-YP_002790881.1:291..790, which codes for envelope protein (E). According to NCBI, the envelope protein "facilitates assembly and release of the virus." It also has ion channel activity. 

## 3. Look at the genomic features.
A. Are the features closely packed; is there a lot of intragenomic space?
The features are very closely packed. It was hard to tell where one gene ends and another starts when looking at the whole genome because of how close the genes are. 

B. Estimate how much of the genome is covered by coding sequences.
The length of both untranslated regions is 437 bp long, and the whole genome is 10,794 bp long, so the percentage of the genome covered by coding sequences is about 96%.

## 4. Find alternative genome builds that could be used to answer a different question.
Include their accession numbers.

Two other Zika virus genome builds are Genome assembly clc4 contig length 11155 (GCF_004788155.1) and Genome assembly ViralProj411812 (GCF_002366285.1).

What other questions could be answered using a different genome build?
Both of the genomes builds mentioned above are newer genome assemblies than the one used in the Zika virus paper, so differences in the genomes could give insight into how the virus has evolved over time. The first genome I suggested is also a bit longer than the original Zika virus genome, so it could contain new genes that weren't present or found in other Zika virus genomes.