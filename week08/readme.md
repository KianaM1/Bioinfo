# Solution to Assignment 8

This assignment requires the presence of a Makefile, a README.md markdown file, and a design.csv file.

## 1. Reuse the last Makefile. Take the paper you were assigned to reproduce and create a design.csv file that connects the SRR numbers to the sample names.

1. Identify the sample names that connect the SRR numbers to the samples.
2. Create a design.csv file that connects the SRR numbers to the sample names.
3. Create a Makefile that can produce multiple BAM alignment files (you can reuse the one from the previous assignment) where from an SRR you can produce a BAM alignment file named after the sample name.
4. Using GNU parallel, run the Makefile on all (or at least 10) samples.
5. Create a README.md file that explains how to run the Makefile