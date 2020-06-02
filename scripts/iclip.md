
## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Trim illumina adapters and quality trimming](#trim-illumina-adapters-and-quality-trimming)
- [Remove duplicates and reads with stretches of several identical nucleotides](#remove-duplicates-and-reads-with-stretches-of-several-identical-nucleotides)
- [Alignment](#alignment)
- [Genomic signal correlation](#genomic-signal-correlation)


## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- Standard Unix tools: cat, awk, sort, uniq, paste, grep, pigz, cut, sbatch, nohup
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [deeptools v3.3.0](https://deeptools.readthedocs.io/en/develop/)


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

### Prepare reference genome

Human reference genome and annotations downloaded from [GENCODE](https://www.gencodegenes.org/)

```bash
cd ~/reference
awk '{print $1}' GRCh38.p12.genome.fa > GRCh38.p12.genome.clean.fa
samtools faidx GRCh38.p12.genome.clean.fa
```

### Align, sort, index and flagstat

```bash
cd fastq_trimmed_nodup

mkdir ../bam
mkdir ../flagstat

ref=~/reference/GRCh38.p12.genome.clean.fa

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../bam/$bname.log --mem 32G --wrap "bwa aln -t 20 -n 0.06 -q 20 $ref $fq | \
  bwa samse $ref - $fq | \
  samtools view -@ 20 -bS - | \
  samtools sort -@ 20 -T ~/tmp/$bname -o ../bam/$bname.bam - && \
  samtools index ../bam/$bname.bam && \
  samtools flagstat ../bam/$bname.bam > ../flagstat/$bname.txt"
done
```

### Filter duplicates, sort and index

```bash
cd bam

ref=~/reference/GRCh38.p12.genome.clean.fa.fai

for bam in *.bam
do
  bname=${bam%.bam}
  nohup samtools view -@ 20 -F 4 $bam | \
  awk 'BEGIN {FS="\t|_"}; {print $4 "_" $5 "_" $2 "_" $3 "\t" $0}' | \
  sort -k 1,1 | \
  awk -F "\t" '!_[$1]++' | \
  cut -f 1 --complement | \
  samtools view -@ 20 -bS -t $ref - | \
  samtools sort -@ 20 -T ~/tmp/$bname -o $bname.clean.bam - &
done

for bam in *.clean.bam
do
  nohup samtools index $bam &
done
```


## Genomic signal correlation

### deeptools

multiBamSummary and plotCorrelation:

```bash
cd bam

mkdir ../deeptools

nohup multiBamSummary bins -b DHX36_WT1.clean.bam \
DHX36_WT2.clean.bam \
DHX36_WT3.clean.bam \
GRSF1_WT1.clean.bam \
GRSF1_WT2.clean.bam \
GRSF1_WT3.clean.bam \
GRSF1_WT4.clean.bam \
DDX3X_WT1.clean.bam \
DDX3X_WT2.clean.bam \
DDX3X_WT3.clean.bam \
-out ../deeptools/DXH36_WT.GRSF1_WT.DDX3X_WT.genome.bins.npz \
-l DHX36_WT_rep1 DHX36_WT_rep2 DHX36_WT_rep3 GRSF1_WT_rep1 GRSF1_WT_rep2 GRSF1_WT_rep3 GRSF1_WT_rep4 DDX3X_WT_rep1 DDX3X_WT_rep2 DDX3X_WT_rep3 \
-p "max" > ../deeptools/DXH36_WT.GRSF1_WT.DDX3X_WT.genome.bins.log &

cd ../deeptools/

nohup plotCorrelation -in DXH36_WT.GRSF1_WT.DDX3X_WT.genome.bins.npz \
-o DXH36_WT.GRSF1_WT.DDX3X_WT.genome.bins_pearson_heatmap.png \
-c pearson \
-p heatmap \
-l "DHX36_WT_rep1" "DHX36_WT_rep2" "DHX36_WT_rep3" "GRSF1_WT_rep1" "GRSF1_WT_rep2" "GRSF1_WT_rep3" "GRSF1_WT_rep4" "DDX3X_WT_rep1" "DDX3X_WT_rep2" "DDX3X_WT_rep3" \
--removeOutliers \
--colorMap Reds \
--plotNumbers > DXH36_WT.GRSF1_WT.DDX3X_WT.genome.bins_pearson_heatmap.log &
```



