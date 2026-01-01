# Solution to Assignment 11

## 1. Generate a VCF file.
I used the following commands to generate a VCF file for the first Zika virus sample, SRR3194431:
```
# Download files from the Bioinformatics Toolkit
bio code

# Run the Makefile commands to generate a VCF file
make run
```
I then viewed the genome, BAM, and VCF files in IGV, as shown below:

## 2. Evaluate the effects of at least three variants in your VCF file.
I decided to visually inspect the variants using IGV. Most of the variants present in the sample were homozygous variants, meaning both alleles in the sample for a particular gene are variants from the reference genome. Two homozygous variants are shown below:

I did find one other variant, a heterozygous variant, meaning one allele matches the reference genome and the other does not. It is shown below:
