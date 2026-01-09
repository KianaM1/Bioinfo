# Solution to Assignment 12
I decided to go with option 2, produce and evaluate variant calls, which involves generating variant calls from the Cancer Genome in a Bottle data and evaluating their quality.

(Expand more on the data/study)

## 1. Call variants for normal and tumor samples in a region of interest
First, I downloaded the reference genome using the following command: 
```
# Enter the bioinfo environment
micromamba activate bioinfo

# Download code and data needed
bio code

# Download reference genome
make genome
```

Then, I extracted the reads from the region of interest needed for the alignment. I used the following command:
```
make bam

# Get statistics on the BAM file.
samtools flagstat bam/KRAS-N-P.bam
```
Which gave me the following output:
```
1127355 + 0 in total (QC-passed reads + QC-failed reads)
1126129 + 0 primary
0 + 0 secondary
1226 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1125885 + 0 mapped (99.87% : N/A)
1124659 + 0 primary mapped (99.87% : N/A)
1126129 + 0 paired in sequencing
562661 + 0 read1
563468 + 0 read2
1116978 + 0 properly paired (99.19% : N/A)
1123189 + 0 with itself and mate mapped
1470 + 0 singletons (0.13% : N/A)
5525 + 0 with mate mapped to a different chr
5077 + 0 with mate mapped to a different chr (mapQ>=5)
```

Next, I called variants for the region of interest:
```
make vcf
```
Which generated a vcf.gz file for the reference genome. 

I then used the following command to make the design.csv file for the experiment:
```
make design
```
With the design file, I called variants for both samples:
```
make batch
```
This created BAM and VCF files for the tumor genome. 

## 2. Compare variant calls between samples and identify tumor-specific variants
To evaluate the variant calls between the samples, both VCF files were merged into one VCF file. This was done using the following command:
```
# Run the variant evaluation recipe
make -f src/recipes/cancer-study-snpeval.mk merge
```
This generated control, tumor, and merged vcf.gz files.

Next, shared and unique variants were found between the two samples:
```
# Intersect the variant calls
bcftools isec -p vcf/isec \
    vcf/KRAS-N-P.vcf.gz \
    vcf/KRAS-T.vcf.gz
```
This resulted in the creation of several files in the vcf/isec directory, including:
- 0000.vcf: contains variants unique to the reference (control) sample
- 0001.vcf: contains variants unique to the tumor sample
- 0002.vcf: contains variants present in both samples

The three VCF files are seen in IGV here:
<img width="1536" height="806" alt="BMMB852_14" src="https://github.com/user-attachments/assets/5df040df-078f-429f-b8e2-aeece6e690f8" />

The tumor VCF file appears to have many more variants than the reference genome. 

## 3. Compare your results to the gold standard DeepVariant calls
I downloaded the DeepVariant VCF file for the gene region using the following command:
```
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz
```

I then loaded the file in IGV, which can be seen here:
<img width="1536" height="806" alt="BMMB852_15" src="https://github.com/user-attachments/assets/4ec7c30b-9f0f-41a2-9a50-0315ec625e9e" />

I believe there are fewer variants present in the gold standard genome compared to both the reference and tumor samples. Interestingly, every variant present in the gold standard is present in the tumor sample.
