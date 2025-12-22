# Solution to Assignment 6

## 1. Transform last week's bash script into a makefile.
Include the following:
- Rules for obtaining the genome and downloading reads from SRA
- Targets to index the genome and generate a sorted and indexed BAM file using align

## 2. Explain the use of the Makefile in your project. 
For this assignment, the goal is to download reads for my assigned virus and align them to create a BAM file using a Makefile. The first step was to download the reference genome and verify the genome was downloaded correctly. The next step was to create a reads directory and download the reads, which were then indexed and aligned. Alignment statistics were generated to confirm the action was completed successfully.

Makefiles allow you to utilize a script but you can choose which steps in a pipeline you wish to run and change variables without having to rewrite the script. You can also print commands with a Makefile, which shows you the commands that would run without having to run them and if files are changed, only steps that depend on the file are rerun. (come back and change this)

## 3. Visualize the BAM files for both simulated reads and reads downloaded from SRA.
For the reads I downloaded, I didn't see any alignment with the BAM file, as shown below.


When I compared the BAM file to the whole SRA viral genome, I saw a few sites of overlap, which can be seen below.

## 4. Generate alignment stats for the BAM file.
I received the following stats were generated for the BAM file"
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
- What percentage of reads aligned to the genome?
0% (none) of the reads aligned to the genome.

- What was the expected average coverage?
I downloaded 1,500 reads, which would have given me roughly 10x coverage.

- What is the observed average coverage?
The observed average coverage is 0x.

- How much does the coverage vary across the genome? (Provide a visual estimate.)
There isn't much coverage observed across the genome. There are some small regions, but not a significant amount of coverage.