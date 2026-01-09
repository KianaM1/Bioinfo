# Solution to Assignment 14
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

Next, I used the following command to determine the number of genes that passed the false discovery rate, or likely had genes changed:
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
This tool gave a slightly higher number for genes detected (FDR) than EDGER.

## 1. Generate PCA and heatmap visualizations of your data
To generate the PCA plot, I used the following command:
```
src/r/plot_pca.r -c edger.csv -d design.csv -o pca.pdf
```
Which gave me the following output:
```
Generating PCA plot
# Design: design.csv
# Counts: edger.csv
# Sample column: sample
# Factor column: group
# Group A has 3 samples.
# Group B has 3 samples.
# Initializing  DESeq2 tibble dplyr ... done
using ntop=500 top features by variance
Warning message:
`aes_string()` was deprecated in ggplot2 3.0.0.
i Please use tidy evaluation idioms with `aes()`.
i See also `vignette("ggplot2-in-packages")` for more information.
i The deprecated feature was likely used in the DESeq2 package.
  Please report the issue to the authors.
# PCA plot: pca.pdf
```
Here is the PCA plot that was generated:
[pca.pdf](https://github.com/user-attachments/files/24532375/pca.pdf)

To generate the heatmap, I used the following commands:
```
src/r/plot_heatmap.r -c edger.csv -d design.csv -o heatmap.pdf
```
Which gave me the following output:
```
Initializing  gplots tibble dplyr tools ... done
# Tool: Create heatmap
# Design: design.csv
# Counts: edger.csv
# Sample column: sample
# Factor column: group
# Group A has 3 samples.
# Group B has 3 samples.
null device
          1
# Output: heatmap.pdf
```
Here is the generated heatmap:
[heatmap.pdf](https://github.com/user-attachments/files/24532378/heatmap.pdf)

## 2. Identify a set of differentially expressed genes or transcripts
To determine the differentially expressed genes, I used the following command:
```
cat edger.csv | cut -f 1 -d ,  | head -10
```
Which gave me the following output:
```
name
GENE-13258
GENE-15610
GENE-5747
GENE-5929
GENE-6396
GENE-10250
GENE-17463
GENE-1171
GENE-423
```

## 3. Perform functional enrichment analysis on your differentially expressed genes
To do the functional enrichment analysis, I used the following command:
```
# Download the data file 
wget http://data.biostarhandbook.com/books/rnaseq/data/edger.csv

# Run g:Profiler to find the functional enrichment of the genes
bio gprofiler -c edger.csv -d hsapiens
```
Note: I had to delete the preexisting edger.csv file. I could've just run the g:Profiler since I already had the data in the directory.

I got the following output:
```
Running g:Profiler
# Counts: edger.csv
# Organism: hsapiens
# Name column: gene
# Pval column: FDR < 0.05
# Gene count: 279
# Genes: IGLC2,SEPTIN3,SYNGR1,MIAT,SEZ6L,[...]
# Submitting to gProfiler
# Found 377 functions
# Output: gprofiler.csv
```
This created the file "gprofiler.csv."

Next, I ran Enrichr to determine the functions of the genes in the file. I used the following command:
```
bio enrichr -c edger.csv
```
Which gave me the following output:
```
Running Enrichr
# Counts: edger.csv
# Organism: mmusculus
# Name column: gene
# Pval column: FDR < 0.05
# Gene count: 279
# Genes: IGLC2,SEPTIN3,SYNGR1,MIAT,SEZ6L,[...]
# Submitting to Enrichr
# User list id: 114734475
# Entries: 95
# Output: enrichr.csv
```
This created the file "enrichr.csv," which lists different functional processes and the genes associated with them.
