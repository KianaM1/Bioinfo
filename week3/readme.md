# Solution to Assignment 3

## 1. Use IGV to visualize your genome and the annotations relative to the genome.
How big is the genome?

I selected a bacteria from the "bacteria 5 collection," acinetobacter_baumannii_atcc_19606_cip_70_34_jcm_6841_gca_000162295. I used the following commands to determine the length of the genome:

```
conda activate bioinfo
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/fasta/bacteria_5_collection/acinetobacter_baumannii_atcc_19606_cip_70_34_jcm_6841_gca_000162295/dna/Acinetobacter_baumannii_atcc_19606_cip_70_34_jcm_6841_gca_000162295.ASM16229v1_.dna.toplevel.fa.gz
gunzip Acinetobacter_baumannii_atcc_19606_cip_70_34_jcm_6841_gca_000162295.ASM16229v1_.dna.toplevel.fa.gz
mv Acinetobacter_baumannii_atcc_19606_cip_70_34_jcm_6841_gca_000162295.ASM16229v1_.dna.toplevel.fa abaum.fa
seqkit stats abaum.fa
```
The genome has a length of 3,971,516.

How many features of each type does the GFF file contain?

I used the following commands to give me the features and how many there are:

```
cat abaum.gff | grep -v '#' > abaum.gff3
cat abaum.gff3 | cut -f 3 | sort-uniq-count-rank
```
Which gave me the following output:

```
3761    exon
3722    CDS
3722    gene
3722    mRNA
39      ncRNA
39      ncRNA_gene
22      region
```

## 2. From your GFF file, separate the intervals of type "gene" or "transcript" into a different file. Show the commands.
I used the following script to put intervals with "gene" or "transcript" into a separate file.

```
cat abaum.gff | grep "gene" | wc -l > genes.gff
cat abaum.gff | grep "transcript" | wc -l > transcripts.gff
```

I used the "ls" command to double check that the new files were created.

## 3. Visualize the simplified GFF in IGV as a separate track. Compare the visualization of the original GFF with the simplified GFF.
Here are the collapsed and expanded views of the sequence for "supercont1.1" in the genome that I chose.


The expanded view presents the transcripts and genes present in the sequence in a way that is much easier to see and understand a quick glance. You can see where the different genes start and stop, whereas with the collapsed view, if two genes are next to each other they appear as one continous box. 

## 4. Zoom in to see the sequences, expand the view to show the translation table in IGV.
Note how the translation table needs to be displayed in the correct orientation for it to make sense.


## 5. Visually verify that the first coding sequence of a gene starts with a start codon and that the last coding sequence of a gene ends with a stop codon.