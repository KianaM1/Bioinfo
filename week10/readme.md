# Solution to Assignment 10

## 1. Extend your existing Makefile to call variants on a single sample of your choice.
- Use your existing BAM file as input
- Generate a VCF file for this sample
- Follow best practices for variant calling
- Visualize the resulting VCF file data alongside the BAM file

To download the Variant Calling Makefile from the Bioinformatics Toolbox, I used the following command:
```
bio code
```
I then went to the "src," then the "recipes" directory to copy and paste the "variant-calling" Makefile into my Makefile.

Using the "variant-calling" makefile template, I added the information to the Makefile with the Zika sample I have been working with this semester.

I then used the following command to generate the VCF file:
```
make run
```
I visualized the Zika SRR3194431 genome, BAM, and VCF files in IGV, as shown below.


## 2. Call variants for all samples.
Run the variant calling workflow for all samples using your design.csv file.

I used the following command to test run the variant calling workflow for all samples:
```
cat design.csv | \
    parallel --dry-run --colsep , --header : --eta --lb -j 4 \
             make \
             SRR={Run} \
             SAMPLE={Sample} \
             GROUP={Group} \
             CONDITION={Condition} \
             run
```
Which gave me the following output:
```
parallel: Warning: The first three values end in CR-NL. Consider using -d "\r\n"

Computers / CPU threads / Max jobs to run
1:local / 16 / 4

Computer:jobs running/jobs completed/%of started jobs/Average seconds to complete
make SRR=SRR3191542 SAMPLE=Sample1 GROUP=Control CONDITION=Mock1_1 run
ETA: 0s Left: 7 AVG: 0.00s  local:2/1/100%/1.0s make SRR=SRR3194428 SAMPLE=Sample2 GROUP=Control CONDITION=Mock1_2 run
ETA: 0s Left: 6 AVG: 0.00s  local:1/2/100%/0.5s make SRR=SRR3191543 SAMPLE=Sample3 GROUP=Control CONDITION=Mock2_1 run
ETA: 0s Left: 5 AVG: 0.00s  local:1/3/100%/0.3s make SRR=SRR3194429 SAMPLE=Sample4 GROUP=Control CONDITION=Mock2_2 run
ETA: 0s Left: 4 AVG: 0.00s  local:1/4/100%/0.2s make SRR=SRR3191544 SAMPLE=Sample5 GROUP=Virus CONDITION=Zika1_1 run
ETA: 0s Left: 3 AVG: 0.00s  local:1/5/100%/0.2s make SRR=SRR3194430 SAMPLE=Sample6 GROUP=Virus CONDITION=Zika1_2 run
ETA: 0s Left: 2 AVG: 0.00s  local:1/6/100%/0.2s make SRR=SRR3191545 SAMPLE=Sample7 GROUP=Virus CONDITION=Zika2_1 run
ETA: 0s Left: 1 AVG: 0.00s  local:1/7/100%/0.1s make SRR=SRR3194431 SAMPLE=Sample8 GROUP=Virus CONDITION=Zika2_2 run
ETA: 0s Left: 0 AVG: 0.00s  local:0/8/100%/0.1s
```

I then used the following command to run the workflow and generate the VCF files for all samples:
```
cat design.csv | \
    parallel --colsep , --header : --eta --lb -j 4 \
             make \
             SRR={Run} \
             SAMPLE={Sample} \
             GROUP={Group} \
             CONDITION={Condition} \
             run
 ```
This command generated VCF files for samples 5 through 8, which I believe is due to the first four samples being mock, or control, genomes that wouldn't have any variants. 

## 3. Create a multisample VCF.
- Merge all individual sample VCF files into a single multisample VCF file (hint: bcftools merge)
- Visualize the multisample VCF in the context of the GFF annotation file.

I used the following command to merge the VCF files:
```
make merge
```
Which generated the merged VCF file. I visualized it in IGV, which can be seen below:
