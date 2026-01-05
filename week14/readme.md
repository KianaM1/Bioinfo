# Solution ot Assignment 14
Perform a differential expression analysis on an RNA-Seq count matrix. This assignment builds upon your previous work, extending the analysis pipeline from read mapping and quantification to differential expression and functional enrichment.

First, I simulated the read counts using the following command in the stats environment:
```
bio code
Rscript src/r/simulate_counts.r
```
Which gave me the following output:
```
Initializing  PROPER ... done
# PROspective Power Evaluation for RNAseq
# Error level: 1 (bottomly)
# All genes: 20000
# Genes with data: 4755
# Genes that changed: 1000
# Changes we can detect: 218
# Replicates: 3
# Design: design.csv
# Counts: counts.csv
```
This created two design.csv files, one with the genes and the changes and expression of each gene (counts.csv), and one with the sample names (design.csv).

Next, I used the following command to determine the number of genes that pass the false discovery rate, or likely had genes changed:
```
Rscript src/r/edger.r -c counts.csv -d design.csv
```
Which gave me the following output:
```
Initializing edgeR tibble dplyr tools ... done
# Tool: edgeR
# Design: design.csv
# Counts: counts.csv
# Sample column: sample
# Factor column: group
# Factors: A B
# Group A has 3 samples.
# Group B has 3 samples.
# Method: glm
# Input: 20000 rows
# Removed: 15224 rows
# Fitted: 4776 rows
# Significant PVal:  388 ( 8.10 %)
# Significant FDRs:  141 ( 3.00 %)
# Results: edger.csv
```
This showed that the program was able to detect 141 significant gene changes. (Maybe mention the number of genes that changed compared to the number of changed genes detected initially)

Next, I wanted to compare the number of significant gene changes detected with EDGER to DESEQ. I used the following commands to determine this:
```
Rscript src/r/deseq2.r -c counts.csv -d design.csv
```
Which gave me the following output:
```
Running DESeq2
# Design: design.csv
# Counts: counts.csv
# Sample column: sample
# Factor column: group
# Group A has 3 samples.
# Group B has 3 samples.
# Initializing  DESeq2 tibble dplyr tools ... done
converting counts to integer mode
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
# Input: 20000 rows
# Removed: 13620 rows
# Fitted: 6380 rows
# Significant PVal: 377 ( 5.9 %)
# Significant FDRs: 168 ( 2.6 %)
# Results: deseq2.csv
```
This tool gave a slightly higher number for genes detected (FDR), than EDGER.

## 1. Generate PCA and heatmap visualizations of your data


## 2. Identify a set of differentially expressed genes or transcripts

## 3. Perform functional enrichment analysis on your differentially expressed genes