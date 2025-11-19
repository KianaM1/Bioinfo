# Solution to Assignment 5

## 1. Identify the BioProject and SRR accession numbers from last week's paper.
BioProject Number: PRJNA313294
SRR Number: SRR3194431

## 2. Write a bash shell script 
- Add commands to download a sequencing dataset using the SRR number. 
- Download a subset of the data that would provide 10x genome coverage
- Explain how you estimated the amount of data needed for 10x coverage.
```
# Obtain run metadata based on the SRR number
bio search SRR3194431

# Download the first 1000 reads based on SRR number
mkdir -p reads
fastq-dump -X 1000 -F --outdir reads --split-files SRR3194431

# View file stats
cd reads/
seqkit stats SRR3194431_1.fastq
```

## 3. Perform a quality assesment.
- Generate basic statistics on the downloaded reads (e.g., number of reads, total bases, average read length).
- Run FASTQC on the downloaded data to generate a quality report. Evaluate the report and summarize the findings.
- Perform any necessary quality control steps (e.g., trimming, filtering) and briefly describe your process.

## 4. Compare sequencing platforms. 
Search the SRA for another dataset for the same genome, but generated using a different sequencing platform (e.g., if original data was Illumina select PacBio or Oxford Nanopore).
Briefly compare the quality or characteristics of the datasets from the two platforms.