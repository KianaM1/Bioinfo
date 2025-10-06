# Solution to Assignment 2


## 1. My Organism: The Channel Bull Blenny
I chose Cottoperca gobio, which is commonly known as the channel bull blenny. This is an ocean fish found below South America, off the coasts of Chile and Argentina.  


## 2. What is the number of sequence regions?
Before I answer this question, I have to download and view the GFF3 file:
```
 micromamba activate bioinfo
 wget https://ftp.ensembl.org/pub/current_gff3/cottoperca_gobio/Cottoperca_gobio.fCotGob3.1.115.gff3.gz
 gunzip Cottoperca_gobio.fCotGob3.1.115.gff3.gz
 ```
 Then I view the file columns:
 ```
 cat Cottoperca_gobio.fCotGob3.1.115.gff3 | more
 ```
 Using this feature, I can determine that the channel bull blenny has 24 chromosomes, as seen below.
```
##sequence-region   1 1 27055867
##sequence-region   10 1 27438269
##sequence-region   11 1 22193419
##sequence-region   12 1 22902521
##sequence-region   13 1 27742916
##sequence-region   14 1 25704503
##sequence-region   15 1 24963210
##sequence-region   16 1 26581398
##sequence-region   17 1 25156145
##sequence-region   18 1 14930383
##sequence-region   19 1 21060126
##sequence-region   2 1 12923129
##sequence-region   20 1 17599799
##sequence-region   21 1 24104793
##sequence-region   22 1 22609112
##sequence-region   23 1 15932363
##sequence-region   24 1 22441820
##sequence-region   3 1 30028995
##sequence-region   4 1 28949045
##sequence-region   5 1 30479438
##sequence-region   6 1 27680345
##sequence-region   7 1 23069680
##sequence-region   8 1 23432536
##sequence-region   9 1 30072547
```
However, using the following code, I observe that the file contains 67,657 individual sequence regions.
```
cat Cottoperca_gobio.fCotGob3.1.115.gff3 | grep CAAAFJ | wc -l
```


## 3. How many features does the file contain?
To view the different features in the file, I used the following scripts:
```
cat Cottoperca_gobio.fCotGob3.1.115.gff3 | grep -v '#' > cotto.gff3
cat cotto.gff3 | cut -f 3 | sort-uniq-count-rank
```
This gave me a list of the file's different features, ranked by number:
```
773175  exon
753229  CDS
57074   mRNA
22564   five_prime_UTR
21662   gene
15877   three_prime_UTR
14099   biological_region
2824    ncRNA_gene
2334    lnc_RNA
563     rRNA
322     region
278     pseudogene
278     pseudogenic_transcript
251     snoRNA
140     snRNA
90      miRNA
31      transcript
26      V_gene_segment
14      J_gene_segment
9       scRNA
1       Y_RNA
```


## 4. How many genes are listed for this organism?
This organism contains 21,662 genes.


## 5. Is there a feature type you haven't heard about before?
Yes, I have never heard of snoRNAs before. From what I have gathered online, snoRNA stands for small nucleolar RNA, which is a group of small RNAs that guide chemical modification of other RNAs.


## 6. What are the top 10 most annotated feature types across the genome?
To answer this question, I could look at the first 10 lines of the last script, or I could use the following script:
```
cat cotto.gff3 | cut -f 3 | sort-uniq-count-rank | head
```
Which gave me the following output:
```
773175  exon
753229  CDS
57074   mRNA
22564   five_prime_UTR
21662   gene
15877   three_prime_UTR
14099   biological_region
2824    ncRNA_gene
2334    lnc_RNA
563     rRNA
```


## 7. Does this seem like a complete and well-annotated organism?
I'm not really sure. There seems to be a lot of chromosome fragments, but I don't know what that means for how well-annotated the organism is. Different features of the organism are labeled, so in that regard, I do think it is well annotated.


## 8. Share any other insights you have.
Most of the chromosome appears to be either exons or coding sequences, which makes sense since most DNA is noncoding. I did notice that the ratio of exons to coding sequences is pretty close, so I wonder if only about half of the DNA is noncoding.