# Solution to Assignment 7

## 1. Add code to last week's Makefile to create BigWig coverage tracks.
In your README.md, demonstrate the use of your Makefile to generate a BAM file for:
- The original data
- The second sequencing data obtained with a different instrument

## 2. Visualize the GFF annotations and both wiggle and BAM files in IGV.


## 3. Answer the following questions:
For both readsets, I used the following command to generate statistics for the respective BAM files:
```
samtools flagstat (SRR).bam
```
I substituted SRR with the respective SRR number for each BAM file.

For the first readset, SRR3194431.bam, the following output was generated:
```
1500 + 0 in total (QC-passed reads + QC-failed reads)
1500 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
6 + 0 mapped (0.40% : N/A)
6 + 0 primary mapped (0.40% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

For the second readset, SRR26779677.bam, the following output was generated:
```
1660 + 0 in total (QC-passed reads + QC-failed reads)
1500 + 0 primary
0 + 0 secondary
160 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1407 + 0 mapped (84.76% : N/A)
1247 + 0 primary mapped (83.13% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

- Briefly describe the differences between the alignment in both files.
The SRR26779677 readset covers a much larger portion of the genome than the initial readset.

Briefly compare the statistics for the two BAM files.
- How many primary alignments does each of your BAM files contain?
I used the following commands to determine the number of primary alignments for the first readset:
```
samtools view -c -F 0x100 SRR3194431.bam
```
I used the same commands to get the number of primary alignments for the second readset as well but with The first readset contains 1500 reads, while the second readset contains 1660. 

- What coordinate has the largest observed coverage? (hint: samtools depth)
I used the following command to determine the observed coverage for the first readset:
```
samtools depth SRR3194431.bam | sort -k3,3nr | head -n 1
```
Which gave me the following output:
```
AY632535.2      10179   2
```

For the second readset, I used the following commands:
```
samtools depth SRR26779677.bam | sort -k3,3nr | head -n 1
```
Which gave me the following output:
```
AY632535.2      3106    118
```

- Select a gene of interest. How many alignments on the forward strand cover the gene?
Going back to the Week 4 assignment, there was only one gene listed in the .gff file, YP_002790881.1. This genic region does encode for different proteins that aid in the virus's function spanning the entire genome.