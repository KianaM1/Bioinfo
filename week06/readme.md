# Solution to Assignment 6

## 1. Transform last week's bash script into a makefile.
- Include rules for obtaining the genome and downloading reads from SRA
- Index the genome
- Generate a sorted and indexed BAM file using align

## 2. Explain the use of the Makefile.
Makefiles allow you to utilize a script but you can choose which steps in a pipeline you wish to run and change variables without having to rewrite the script. You can also print commands with a Makefile, which shows you the commands that would run without having to run them and if files are changed, only steps that depend on the file are rerun.

## 3. Visualize the BAM file for both simulated reads and reads downloaded from SRA.

## 4. Generate alignment stats for the BAM file.
- What percentage of reads aligned to the genome?
- What was the expected average coverage?
- What is the observed average coverage?
- How much does the coverage vary across the genome? (Provide a visual estimate.)