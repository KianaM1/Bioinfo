# Solution to Assignment 10

## 1. Extend your existing Makefile to call variants on a single sample of your choice.
- Use your existing BAM file as input
- Generate a VCF file for this sample
- Follow best practices for variant calling
- Visualize the resulting VCF file data alongside the BAM file

## 2. Call variants for all samples.

Run the variant calling workflow for all samples using your design.csv file.

## 3. Create a multisample VCF.
- Merge all individual sample VCF files into a single multisample VCF file (hint: bcftools merge)
- Visualize the multisample VCF in the context of the GFF annotation file.

If any samples show poor alignment or no variants, identify and replace them with better samples. Ensure you have sufficient genome coverage across all samples