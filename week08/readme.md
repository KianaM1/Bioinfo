# Solution to Assignment 8

## 1. Identify the sample names that connect the SRR numbers to the samples.
To do this, I went to the Zika paper, found the accession number, and went to the list of SRA experiments that were performed. I then listed the SRR numbers according to their condition, as shown in the design.csv file.

## 2. Create a design.csv file that connects the SRR numbers to the sample names.
Following the example provided in the handbook, I made a design.csv file for the eight SRR numbers from the paper.

## 3. Create a Makefile that can produce multiple BAM alignment files (you can reuse the one from the previous assignment), where from an SRR you can produce a BAM alignment file named after the sample name.
I reused the Makefile from last week's assignment, but I added the commands to make a BigWig file to the "all" command. This ensures the BigWig files are made when running the GNU parallel.

## 4. Using GNU parallel, run the Makefile on all samples.
I used the following command to test the GNU parallel before running it:
```
cat design.csv | \
    parallel --dry-run --colsep , --header : --eta --lb -j 4 \
             make \
             SRR={Run} \
             SAMPLE={Sample} \
             GROUP={Group} \
             CONDITION={Condition} \
             all
```
Which gave me the following output:
```
parallel: Warning: The first three values end in CR-NL. Consider using -d "\r\n"

Computers / CPU threads / Max jobs to run
1:local / 16 / 4

Computer:jobs running/jobs completed/%of started jobs/Average seconds to complete
make SRR=SRR3191542 SAMPLE=Sample1 GROUP=Control CONDITION=Mock1_1 all
ETA: 0s Left: 7 AVG: 0.00s  local:1/1/100%/1.0s make SRR=SRR3194428 SAMPLE=Sample2 GROUP=Control CONDITION=Mock1_2 all
ETA: 0s Left: 6 AVG: 0.00s  local:1/2/100%/0.5s make SRR=SRR3191543 SAMPLE=Sample3 GROUP=Control CONDITION=Mock2_1 all
ETA: 0s Left: 5 AVG: 0.00s  local:1/3/100%/0.3s make SRR=SRR3194429 SAMPLE=Sample4 GROUP=Control CONDITION=Mock2_2 all
ETA: 0s Left: 4 AVG: 0.00s  local:1/4/100%/0.2s make SRR=SRR3191544 SAMPLE=Sample5 GROUP=Virus CONDITION=Zika1_1 all
ETA: 0s Left: 3 AVG: 0.00s  local:1/5/100%/0.2s make SRR=SRR3194430 SAMPLE=Sample6 GROUP=Virus CONDITION=Zika1_2 all
ETA: 0s Left: 2 AVG: 0.00s  local:1/6/100%/0.2s make SRR=SRR3191545 SAMPLE=Sample7 GROUP=Virus CONDITION=Zika2_1 all
ETA: 0s Left: 1 AVG: 0.00s  local:1/7/100%/0.1s make SRR=SRR3194431 SAMPLE=Sample8 GROUP=Virus CONDITION=Zika2_2 all
ETA: 0s Left: 0 AVG: 0.00s  local:0/8/100%/0.1s
```

I did see one warning (shown above), but I didn't know what it meant, so I ran the parallel just to see if it would cause any issues. I used the following command to run all of the functions:
```
cat design.csv |     parallel --colsep , --header : --eta --lb -j 4
       make              SRR={Run}              SAMPLE={Sample}
 GROUP={Group}              CONDITION={Condition}              all
```

I didn't see any errors, and the only potential issue I saw was that two samples, Sample2 and Sample4, didn't have a BAI, BEDGRAPH, or BW file.
