# Solution to Assignment 5

## 1. Identify the BioProject and SRR accession numbers from last week's paper.
- BioProject Number: PRJNA313294
- SRR Number: SRR3194431

## 2. Write a bash shell script 
- Add commands to download a sequencing dataset using the SRR number. 
- Download a subset of the data that would provide 10x genome coverage
- Explain how you estimated the amount of data needed for 10x coverage.

### Coverage Calculation
Coverage = total number of sequenced bases/genome size

Goal: Want 10x genome coverage

Zika virus genome size = 10,794 bp
```
10 = x/10,794
10(10,794) = x(10,794)/10,794
107,940 = x 
Need 107,940 sequenced bases
```

Average read length = 75 bp (found by looking at the metadata for 1,000 reads)
```
107,940/75 = x
1,439.2 = x = number of reads needed for 10x genome coverage
```

Can download 1,500 reads to get 10x genome coverage.

## 3. Perform a quality assessment.
- Generate basic statistics on the downloaded reads (e.g., number of reads, total bases, average read length).
- Run FASTQC on the downloaded data to generate a quality report. Evaluate the report and summarize the findings.
- Perform any necessary quality control steps (e.g., trimming, filtering) and briefly describe your process.

### FastQC Report
- The per base sequence quality of the reads was generally high. No sequences were flagged as being of poor quality
- The reads had 47% GC content
- The per sequence GC content was quite varied, but somewhat followed a normal distribution
- The per base N content was consistently zero except for in position six in the reads
- Sequences varied in length from 35 to 76 bp. Most sequences were greater than 72 bp
- Two sequences appeared twice in the read set
- There was no adapter content present in the reads
FastQC report: [FastQC_Report1.pdf](https://github.com/user-attachments/files/23667183/FastQC_Report1.pdf)

### Quality Control Steps
I decided to trim and filter the downloaded reads using fastp. The reads are single-end, so I used the appropriate commands to trim and filter single-end reads.

## 4. Compare sequencing platforms.
Search the SRA for another dataset for the same genome, but generated using a different sequencing platform (e.g., if the original data was Illumina, select PacBio or Oxford Nanopore).
Briefly compare the quality or characteristics of the datasets from the two platforms.

### Zika virus genome sequenced using Oxford Nanopore
- BioProject number: PRJNA1035959
- SRR number: SRR26779677

This genome read set has over eight times as many sequences as the previous read set, despite my downloading 1,500 reads, indicating that I should probably have calculated how many reads would have given me 10x genome coverage. The average read length was 537.4 bp, which is much higher than the 35 bp average read length that my initial read set had.

The FastQC Report:
- The percentage GC content was 49%, which is slightly higher than my initial set.
- The majority of the sequences had a low per base sequence quality, and their quality scores were low
- There was no N content present in the read set
- The majority of sequences were about 500 bp in length
- There were a lot of sequences that appeared more than once in the read set
- Very little adapter content was detected, similar to the previous read set.
FastQC report: [FastQC_Report2.pdf](https://github.com/user-attachments/files/23667184/FastQC_Report2.pdf)

## Bash Script
```
# Obtain run metadata based on the SRR number
bio search SRR3194431

# Download the first 1500 reads based on SRR number
mkdir -p reads
cd reads/
fastq-dump -X 1500 -F --split-files SRR3194431

# View file stats
seqkit stats SRR3194431_1.fastq

# Generate a FastQC report
fastqc SRR3194431_1.fastq

# Single end quality control with fastp
fastp --cut_tail -i SRR3194431_1.fastq -o SRR3194431_1.trim.fq
```
