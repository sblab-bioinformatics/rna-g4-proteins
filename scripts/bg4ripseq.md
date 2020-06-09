
## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Trim illumina adapters and quality trimming](#trim-illumina-adapters-and-quality-trimming)
- [Alignment](#alignment)
- [Genomic signal normalization](#genomic-signal-normalization)
- [Peak calling](#peak-calling)


## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- Standard Unix tools: cat, awk, sort, uniq, paste, grep, pigz, cut, sbatch, nohup, sed
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [sambamba v0.6.5](https://academic.oup.com/bioinformatics/article/31/12/2032/214758)
- [deeptools v3.3.0](https://deeptools.readthedocs.io/en/develop/)
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [edgeR v3.16.5](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html)


## Libraries

Sample | Replicate | File name | GEO GSE | GEO GSM
:--------:|:---------:|:---------:|:-------:|:-------:
BG4 | 1 | BG4_1.fastq.gz | - | -
BG4 | 2 | BG4_2.fastq.gz | - | -
BG4 | 3 | BG4_3.fastq.gz | - | -
A9 | 1 | A9_1.fastq.gz | - | -
A9 | 2 | A9_2.fastq.gz | - | -
A9 | 3 | A9_3.fastq.gz | - | -
INPUT | 1 | INPUT_1.fastq.gz | - | -
INPUT | 2 | INPUT_2.fastq.gz | - | -
INPUT | 3 | INPUT_3.fastq.gz | - | -


## Trim illumina adapters and quality trimming

```bash
cd ~/fastq

mkdir ../fastq_trimmed

for fq in *.fq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4G --wrap "cutadapt -a AGATCGGAAGAGC -m 20 --max-n=20 -q 20 -o ../fastq_trimmed/${bname}.fq.gz $fq > ../fastq_trimmed/$bname.txt"
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

### Align and sort

```bash
cd ~/fastq_trimmed

mkdir ../bam

ref=~/reference/GRCh38.p12.genome.clean.fa

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../bam/$bname.log --mem 32G --wrap "bwa mem -t 20 -M $ref $fq | \
  samtools view -@ 20 -b - | \
  samtools sort -@ 20 -T ~/tmp/$bname -o ../bam/$bname.tmp.bam -"
done
```

### Mark duplicates, index and flagstat

```bash
cd ~/bam

mkdir ../flagstat

for bam in *.tmp.bam
do
  id=${bam%.tmp.bam}
  sbatch -J $id -o $id.markdup.log --mem 16G --wrap "sambamba markdup -t 20 $bam ${id}.bam 2> ${id}.markdup.txt && \
  samtools flagstat ${id}.bam > ../flagstat/${id}.txt"
done
```

### Filter and index

```bash
cd ~/bam

for bam in `ls *.bam | grep -v "tmp"`
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.clean.log --mem 16G --wrap "samtools idxstats $bam | cut -f1 | grep '^chr' | xargs samtools view -@ 20 -b -F 516 -q 1 $bam > $bname.clean.bam && samtools index $bname.clean.bam"
done
```


## Genomic signal normalization

```bash
cd ~/bam

mkdir ../bw

for bam in *.clean.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../bw/$bname.log --mem 32G --wrap "bamCoverage -b $bam -o ../bw/$bname.bw -of bigwig --binSize 1 -p 20 --normalizeUsing CPM"
done

cd ~/bam

mkdir ../bedgraph

for bam in *.clean.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../bedgraph/$bname.log --mem 32G --wrap "bamCoverage -b $bam -o ../bedgraph/$bname.bedgraph -of bedgraph --binSize 1 -p 20 --normalizeUsing CPM"
done
```


## Peak calling

### define target regions

with enough CPM and width:

```bash
cd ~/bedgraph

for bdg in BG4_*.bedgraph A9_*.bedgraph
do
  bname=${bdg%.bedgraph}
  nohup awk '{if($4 >= 1){print $0}}' $bdg | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{if(($3-$2)>=50){print $0}}' > $bname.bed &
done

bedtools multiinter -i A9_*.bed | awk '$4>=2' | bedtools sort -i - | bedtools merge -i - > A9.clean.bed
bedtools multiinter -i BG4_*.bed | awk '$4>=2' | bedtools sort -i - | bedtools merge -i - > BG4.clean.bed

bedtools multiinter -i A9.clean.bed BG4.clean.bed | awk '$4>=1' | bedtools sort -i - | bedtools merge -i - > A9.BG4.clean.bed
```

### count

```bash
cd ~/bam

sbatch -J counting -o counting.log --mem 64G --wrap "bedtools multicov -bams BG4_1.clean.bam \
BG4_2.clean.bam \
BG4_3.clean.bam \
A9_1.clean.bam \
A9_2.clean.bam \
A9_3.clean.bam \
INPUT_1.clean.bam \
INPUT_2.clean.bam \
INPUT_3.clean.bam \
-bed ../bedgraph/A9.BG4.clean.bed > ../bedgraph/A9.BG4.clean.counts.bed"
```

### edgeR

```R
library(data.table)
library(edgeR)
library(ggplot2)

# Enlarge the view width when printing tables
options(width = 250)

# Load data
data <- fread("~/bedgraph_i2/A9.BG4.clean.counts.bed")
setnames(data, c("chr", "start", "end", "bg4_1", "bg4_2", "bg4_3", "a9_1", "a9_2", "a9_3", "input_1", "input_2", "input_3"))

# Define group
group <- factor(c('bg4', 'bg4', 'bg4', 'a9', 'a9', 'a9', 'input', 'input', 'input'))

# Define DGEList object
y <- DGEList(counts = data[,-c(1,2,3)], group = group, genes = data[,c(1,2,3)])

# Define design matrix
des <- model.matrix(~ 0 + group, data = y$samples)
colnames(des) <- levels(factor(y$samples$group))

# Calculate normalization factors
y <- calcNormFactors(y, method = "TMM")

# Estimate dispersion
y <- estimateDisp(y, des)

# Fit linear model
fit <- glmFit(y, des)

# Define matrix of contrasts
my.contrasts <- makeContrasts(bg4vsinput = bg4 - input, a9vsinput = a9 - input, levels=des)

# Obtain likelihoods
lrt_bg4vsinput <- glmLRT(fit, contrast=my.contrasts[,"bg4vsinput"])
lrt_a9vsinput <- glmLRT(fit, contrast=my.contrasts[,"a9vsinput"])

# Tables
table_bg4vsinput <- data.table(data.frame(topTags(lrt_bg4vsinput, Inf)))[logFC > log2(1.5) & FDR < 0.05][order(-logFC)]
write.table(table_bg4vsinput[, c("chr", "start", "end", "logFC", "PValue", "FDR")], file = "~/bedgraph/table_bg4vsinput.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

table_a9vsinput <- data.table(data.frame(topTags(lrt_a9vsinput, Inf)))[logFC > log2(1.5) & FDR < 0.05][order(-logFC)]
write.table(table_a9vsinput[, c("chr", "start", "end", "logFC", "PValue", "FDR")], file = "~/bedgraph/table_a9vsinput.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```

Further filtering and keep BG4 - A9 only:

```bash
cd ~/bedgraph

bedtools intersect \
-a <(awk '$4 > 0.8 {print $0}' table_bg4vsinput.bed | bedtools sort -i) \
-b <(awk '$4 > 0.8 {print $0}' table_a9vsinput.bed | bedtools sort -i) \
-v > bg4_rip-seq.bed
```
