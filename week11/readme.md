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
<img width="1536" height="808" alt="BMMB852_9" src="https://github.com/user-attachments/assets/90594504-682d-4812-b98f-bba60b9a6795" />

## 2. Evaluate the effects of at least three variants in your VCF file.
I decided to visually inspect the variants using IGV. Most of the variants present in the sample were homozygous, meaning that both alleles in the sample for a particular gene are variants from the reference genome. Two homozygous variants are shown below:
<img width="1536" height="808" alt="BMMB852_10" src="https://github.com/user-attachments/assets/ae3d90dc-8fc8-41c0-b87e-e709959108e9" />

I did find one other variant, a heterozygous variant, meaning one allele matches the reference genome and the other does not. It is shown below:
<img width="1536" height="808" alt="BMMB852_11" src="https://github.com/user-attachments/assets/08a50592-3e82-40a9-9cd5-382eb37d3fba" />
