# Solution to Assignment 13

## 1. Write a Makefile that aligns the reads to the genome and creates BAM and BigWig files.
I decided to reuse the data and code from the Biostar Workflows: RNA-Seq with Hisat2 chapter in the website. This project contained two datasets: the Universal Human Reference (UHR), which contains RNA from 10 cancer cell lines, and the Human Brain Reference (HBR), which contains RNA from the brains of twenty-three individuals of European descent (both male and female of varying ages, but mainly sixty to eighty years old).

Download and unpack the data using the following command:
```
# Download the data
wget -nc  http://data.biostarhandbook.com/data/uhr-hbr.tar.gz

# Unpack the data
tar xzvf uhr-hbr.tar.gz
```
This commnand downloads the data and creates the `reads` and `refs` directories. For this assignment, I'm only using human chromosome twenty-two as a reference.

I then investigated the files present in the `refs` directory, a FASTA genome file, GTF file, and FASTA transcriptome file, using the following commands:
```
seqkit stats refs/*.fa
```
Which gave me the following output:
```
processed files:  2 / 2 [====================================] ETA: 0s. done
file                       format  type  num_seqs     sum_len     min_len     avg_len     max_len
refs/chr22.genome.fa       FASTA   DNA          1  50,818,468  50,818,468  50,818,468  50,818,468
refs/chr22.transcripts.fa  FASTA   DNA      4,506   7,079,970          33     1,571.2      84,332
```

I then investigated the FASTQ files for the samples using the following commands:
```
seqkit stats reads/*.fq
```
Which gave me the following output:
```
processed files:  6 / 6 [====================================] ETA: 0s. done
file               format  type  num_seqs     sum_len  min_len  avg_len  max_len
reads/HBR_1_R1.fq  FASTQ   DNA    118,571  11,857,100      100      100      100
reads/HBR_2_R1.fq  FASTQ   DNA    144,826  14,482,600      100      100      100
reads/HBR_3_R1.fq  FASTQ   DNA    129,786  12,978,600      100      100      100
reads/UHR_1_R1.fq  FASTQ   DNA    227,392  22,739,200      100      100      100
reads/UHR_2_R1.fq  FASTQ   DNA    162,373  16,237,300      100      100      100
reads/UHR_3_R1.fq  FASTQ   DNA    185,442  18,544,200      100      100      100
```

I used the following command to make a design.csv file:
```
make design.csv
```
Which generated a design.csv file with the six samples (three control samples, UHR, and three treatment samples, HBR).

Next, I indexed the reference genome using the following command:
```
make index
```

Then I aligned the samples using the following command:
```
make align
```
This generated BAM, BAI, BedGraph, and BigWig files for each of the samples. I visualized the BigWig and GTF files in IGV.

## 2. Run a feature counter to create a count matrix for your data. 
- The final result of your code should be a count matrix that summarizes read counts for each dataset.
- Include IGV screenshots that demonstrate your data is RNA-Seq data.

First, I made the counts.txt file using the following command:
```
make res/counts-hisat.txt
```
Which generated a counts-histat.txt file and a summary file.

I then made a design.csv file with the information from the counts-histat.txt file using the following command:
```
make res/counts-hisat.csv
```
Which gave me the following output:
```
micromamba run -n stats Rscript src/r/format_featurecounts.r -c res/counts-hisat.txt -o res/counts-hisat.csv
# Reformating featurecounts.
# Input: res/counts-hisat.txt
# Output: res/counts-hisat.csv
```
This generated the file counts-histat.csv in the `res` directory.

I then used the following commands to map the transcripts to the gene names:
```
micromamba activate stats
Rscript src/r/create_tx2gene.r -s > names.txt
Rscript src/r/create_tx2gene.r -d hsapiens_gene_ensembl
```
This created a names.txt file containing gene information that was mapped to the genome for this experiment. This was the output:
```
# Create tx2gene mapping
# Dataset: hsapiens_gene_ensembl
# Connecting to ensembl.
# Submitting the query.
# Output: tx2gene.csv
```
I used the following command to add the Ensembl gene IDs to the counts-histats.csv file:
```
Rscript src/r/format_featurecounts.r -c res/counts-hisat.txt -t tx2gene.csv -o res/counts-hisat.csv
```
Which gave me the following output:
```
# Reformating featurecounts.
# Input: res/counts-hisat.txt
# Tx2gene: tx2gene.csv
# Output: res/counts-hisat.csv
```

## 3. Discuss a few lines of the resulting count matrix. 
- Visually identify rows where the counts show consistent gene expression levels. 
- Try to visually verify that the counts in the count matrix are consistent with the numbers you can observe in the alignment tracks.

Below is a selection of counts from the counts-hisat.csv file:
```
name,gene,HBR_1,HBR_2,HBR_3,UHR_1,UHR_2,UHR_3
ENSG00000225255.6,PSLNR,0,0,0,2,1,0
ENSG00000235992.1,GRAMD4P2,0,0,0,0,1,0
ENSG00000206195.10,DUXAP8,10,4,7,184,156,168
ENSG00000271127.1,ENSG00000271127,0,1,0,16,9,19
ENSG00000232775.6,BMS1P22,0,1,3,21,21,22
ENSG00000272872.1,ENSG00000272872,3,2,0,96,36,49
ENSG00000271672.1,DUXAP8,0,0,0,0,0,0
ENSG00000223875.2,NBEAP3,0,0,0,4,3,0
ENSG00000215270.3,TOMM40P2,0,0,0,11,6,6
ENSG00000229286.1,ENSG00000229286,0,0,0,0,0,0
ENSG00000233866.1,ENSG00000233866,0,0,0,0,2,0
```
This selection is from lines 42-52, and demonstrates varying distribution of some genes in the samples. For example, the genes DUXAP8 and ENSG00000272872 are highly expressed in the UHR samples, but not in the HBR samples. This trend is present in several of the genes shown above. Otherwise, the genes are not present/expressed in any of the samples.