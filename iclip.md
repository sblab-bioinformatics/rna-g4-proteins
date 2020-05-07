
## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Trim illumina adapters and quality trimming](#trim-illumina-adapters-and-quality-trimming)
- [Remove duplicates and reads with stretches of several identical nucleotides](#remove-duplicates-and-reads-with-stretches-of-several-identical-nucleotides)
- [Alignment](#alignment)


## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- Standard Unix tools: cat, awk, sort, uniq, paste, grep, pigz
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [sambamba v0.6.5](https://academic.oup.com/bioinformatics/article/31/12/2032/214758)
- [igvtools v2.3.91](https://software.broadinstitute.org/software/igv/igvtools)
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)
- [macs2 v2.1.1.20160309](https://github.com/taoliu/MACS)
- [fastaRegexFinder.py v0.1.1](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder)
- [gat-run.py](http://gat.readthedocs.io/en/latest/contents.html)
- [python v2.7.12](https://www.python.org/). Libraries:
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [GenomicFeatures v1.26.4](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [edgeR v3.16.5](https://bioconductor.org/packages/release/bioc/html/edgeR.html)


## Libraries

Protein | Genotype | Replicate | File name | GEO GSE | GEO GSM
:------:|:--------:|:---------:|:---------:|:-------:|:-------:
DHX36 | WT | 1 | DHX36_WT1.fastq.gz | - | -
DHX36 | WT | 2 | DHX36_WT2.fastq.gz | - | -
DHX36 | WT | 3 | DHX36_WT3.fastq.gz | - | -
DHX36 | EA | 1 | DHX36_EA1.fastq.gz | - | -
DHX36 | EA | 2 | DHX36_EA2.fastq.gz | - | -
DHX36 | EA | 3 | DHX36_EA3.fastq.gz | - | -
GRSF1 | WT | 1 | GRSF1_WT1.fastq.gz | - | -
GRSF1 | WT | 2 | GRSF1_WT2.fastq.gz | - | -
GRSF1 | WT | 3 | GRSF1_WT3.fastq.gz | - | -
GRSF1 | WT | 3 | GRSF1_WT3.fastq.gz | - | -
DDX3X | WT | 1 | DDX3X_WT1.fastq.gz | GSE106476 | GSM2838585
DDX3X | WT | 2 | DDX3X_WT2.fastq.gz | GSE106476 | GSM2838586
DDX3X | WT | 3 | DDX3X_WT3.fastq.gz | GSE106476 | GSM2838587


## Trim illumina adapters and quality trimming

```bash
cd fastq

mkdir ../fastq_trimmed

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4G --wrap "cutadapt -a AGATCGGAAGAGC -m 10 -q 10 -O 5 -o ../fastq_trimmed/$fq $fq > ../fastq_trimmed/$bname.txt"
done
```


## Remove duplicates and reads with stretches of several identical nucleotides

```bash
cd fastq_trimmed

mkdir ../fastq_trimmed_nodup/

for fq in *.fastq.gz
do
  nohup zcat $fq | \
  paste - - - - | \
  awk 'BEGIN {FS = "\t|_"}; {print $1 "\t" $3 "_" $2 "\t" $4 "\t" $5}' | \
  sort -k 2,2 | \
  awk -F "\t" '!_[$2]++' | \
  awk  'BEGIN {FS="\t|_"}; {print $1"_"$3 "\n" $2 "\n" $4 "\n" $5}' | \
  paste -d "\t" - - - - | \
  grep -P -v 'T{10,}' | \
  grep -P -v 'A{10,}' | \
  grep -P -v 'C{10,}' | \
  grep -P -v 'G{10,}' | \
  grep -P -v '\tT{8,}' | \
  grep -P -v '\tA{8,}' | \
  grep -P -v '\tG{8,}' | \
  grep -P -v '\tC{8,}' | \
  awk 'BEGIN{ FS = "\t" }; {print $1 "\n" $2 "\n" $3 "\n" $4}' | \
  pigz > ../fastq_trimmed_nodup/$fq &
done
```


## Alignment

### Index reference genome

Already performed before



### Align, sort, index and flagstat

Following the [iCLIP CTK tutorial bwa options](https://zhanglab.c2b2.columbia.edu/index.php/ICLIP_data_analysis_using_CTK#Read_mapping_.26_parsing).

```bash
cd /scratchb/sblab/martin03/repository/20180117_Ummi/data/20190605/fastq_trimmed_nodup

mkdir ../bam
mkdir ../flagstat

ref=/scratcha/sblab/martin03/reference_data/genomes/gencode/Gencode_human/release_28/fasta/GRCh38.p12.genome.clean.fa

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
#  echo $fq, $bname
  sbatch -J $bname -o ../bam/$bname.log --mem 32G --wrap "bwa aln -t 20 -n 0.06 -q 20 $ref $fq | \
  bwa samse $ref - $fq | \
  samtools view -@ 20 -bS - | \
  samtools sort -@ 20 -T /scratchb/sblab/martin03/tmp/$bname -o ../bam/$bname.bam - && \
  samtools index ../bam/$bname.bam && \
  samtools flagstat ../bam/$bname.bam > ../flagstat/$bname.txt"
done

cd ../bam
tail *.log

for bam in *.bam
do
  bname=${bam%.bam}
  sbatch -J $bname --mem 32G --wrap "samtools flagstat $bam > ../flagstat/$bname.txt"
done

cd ../flagstat
grep "mapped (" *.txt
# DDX17_WT1.txt:4149539 + 0 mapped (84.23% : N/A)
# DDX17_WT2.txt:8138904 + 0 mapped (88.88% : N/A)
# DDX17_WT3.txt:12240199 + 0 mapped (90.79% : N/A)
# DDX3X_mRG1.txt:3907263 + 0 mapped (73.70% : N/A)
# DDX3X_mRG2.txt:5819734 + 0 mapped (75.91% : N/A)
# DDX3X_mRG3.txt:5425604 + 0 mapped (70.66% : N/A)
# DDX3X_WT1.txt:6035717 + 0 mapped (75.02% : N/A)
# DDX3X_WT2.txt:2900577 + 0 mapped (67.66% : N/A)
# DDX3X_WT3.txt:4575647 + 0 mapped (66.32% : N/A)
# DHX36_EA1.txt:16084667 + 0 mapped (85.92% : N/A)
# DHX36_EA2.txt:4257481 + 0 mapped (73.84% : N/A)
# DHX36_EA3.txt:11294049 + 0 mapped (85.25% : N/A)
# DHX36_WT1.txt:3682542 + 0 mapped (80.74% : N/A)
# DHX36_WT2.txt:3560509 + 0 mapped (78.89% : N/A)
# DHX36_WT3.txt:5006384 + 0 mapped (78.87% : N/A)
# GRSF1_WT1.txt:4171362 + 0 mapped (75.34% : N/A)
# GRSF1_WT2.txt:6256568 + 0 mapped (79.45% : N/A)
# GRSF1_WT3.txt:5661404 + 0 mapped (76.11% : N/A)
# GRSF1_WT4.txt:9718016 + 0 mapped (78.48% : N/A)
```



### Filter duplicates, sort and index

```bash
srun --mem 16G --pty /usr/bin/bash

cd /scratchb/sblab/martin03/repository/20180117_Ummi/data/20190605/bam

ref=/scratcha/sblab/martin03/reference_data/genomes/gencode/Gencode_human/release_28/fasta/GRCh38.p12.genome.clean.fa.fai

for bam in *.bam
do
  bname=${bam%.bam}
  nohup samtools view -@ 20 -F 4 $bam | \
  awk 'BEGIN {FS="\t|_"}; {print $4 "_" $5 "_" $2 "_" $3 "\t" $0}' | \
  sort -k 1,1 | \
  awk -F "\t" '!_[$1]++' | \
  cut -f 1 --complement | \
  samtools view -@ 20 -bS -t $ref - | \
  samtools sort -@ 20 -T /scratchb/sblab/martin03/tmp/$bname -o $bname.clean.bam - &
done

for bam in *.clean.bam
do
  nohup samtools index $bam &
done

rm nohup.out

for bam in *.clean.bam
do
  n=`samtools view -@ 20 $bam -c`
  echo $bam, $n
done
# DDX17_WT1.clean.bam, 2772030
# DDX17_WT2.clean.bam, 6161796
# DDX17_WT3.clean.bam, 10516220
# DDX3X_mRG1.clean.bam, 1417348
# DDX3X_mRG2.clean.bam, 2537033
# DDX3X_mRG3.clean.bam, 2039541
# DDX3X_WT1.clean.bam, 2928061
# DDX3X_WT2.clean.bam, 837188
# DDX3X_WT3.clean.bam, 1516372
# DHX36_EA1.clean.bam, 13669267
# DHX36_EA2.clean.bam, 2693858
# DHX36_EA3.clean.bam, 9182754
# DHX36_WT1.clean.bam, 2489367
# DHX36_WT2.clean.bam, 2324378
# DHX36_WT3.clean.bam, 2964051
# GRSF1_WT1.clean.bam, 2340663
# GRSF1_WT2.clean.bam, 3437582
# GRSF1_WT3.clean.bam, 2779663
# GRSF1_WT4.clean.bam, 4743784

exit
```

Protein | Sample | Raw | Trimmed | Clean | Aligned | Filtered | %
:------:|:------:|:----:|:------:|:-----:|:-------:|:--------:|:-:
DHX36 | WT1 | 50,880,139 | 50,348,978 | 4,560,786 | 3,682,542 | 2,489,367 | 4.9
DHX36 | WT2 | 53,371,798 | 52,234,294 | 4,513,415 | 3,560,509 | 2,324,378 | 4.4
DHX36 | WT3 | 68,693,632 | 68,399,312 | 6,347,369 | 5,006,384 | 2,964,051 | 4.3
DHX36 | EA1 | 72,180,436 | 72,060,200 | 18,721,182 | 16,084,667 | 13,669,267 | 19
DHX36 | EA2 | 61,775,848 | 61,575,095 | 5,765,857 | 4,257,481 | 2,693,858 | 4.4
DHX36 | EA3 | 68,026,536 | 67,945,342 | 13,248,802 | 11,294,049 | 9,182,754 | 13.5
GRSF1 | WT1 | 98,449,831 | 97,536104 | 5,536,812 | 4,171,362 | 2,340,663 | 2.4
GRSF1 | WT2 | 116,295,996 | 115,008684 | 7,875,038 | 6,256,568 | 3,437,582 | 3.0
GRSF1 | WT3 | 96,728,879 | 96,218238 | 7,438,519 | 5,661,404 | 2,779,663 | 2.9
GRSF1 | WT4 | 122,375,817 | 121,825246 | 12,382,333 | 9,718,016 | 4,743,784 | 3.9
DDX3X | WT1 | 74,249,222 | 74,121,070 | 8,045,617 | 6,035,717 | 2,928,061 | 4.0
DDX3X | WT2 | 64,455,519 | 64,093,750 | 4,286,683 | 2,900,577 | 837,188 | 1.3
DDX3X | WT3 | 68,927,680 | 68,729,583 | 6,899,370 | 4,575,647 | 1,516,372 | 2.2
DDX3X | mRG1 | 67,402,219 | 67,151,641 | 5,301,572 | 3,907,263 | 1,417,348 | 2.1
DDX3X | mRG2 | 70,529,397 | 70,445,995 | 7,666,258 | 5,819,734 | 2,537,033 | 3.6
DDX3X | mRG3 | 87,728,012 | 87,604,141 | 7,678,671 | 5,425,604 | 2,039,541 | 2.3
DDX17 | WT1 | 29,693,100 | 29,669,719 | 4,926,485 | 4,149,539 | 2,772,030 | 9
DDX17 | WT2 | 41,499,120 | 41,477,279 | 9,157,529 | 8,138,904 | 6,161,796 | 15
DDX17 | WT3 | 43,146,889 | 43,077,633 | 13,481,282 | 12,240,199 | 10,516,220 | 24
