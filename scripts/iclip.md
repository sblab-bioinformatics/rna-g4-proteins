
## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Trim illumina adapters and quality trimming](#trim-illumina-adapters-and-quality-trimming)
- [Remove duplicates and reads with stretches of several identical nucleotides](#remove-duplicates-and-reads-with-stretches-of-several-identical-nucleotides)
- [Alignment](#alignment)
- [Genomic signal correlation](#genomic-signal-correlation)
- [Genomic signal normalization](#genomic-signal-normalization)
- [Peak calling](#peak-calling)
- [Venn diagram](#venn-diagram)
- [GAT](#gat)
- [Overlap with genomic features and PQS](#overlap-with-genomic-features-and-pqs)
- [Binding profiles around PQS](#binding-profiles-around-pqs)
- [Transcripts analysis](#transcripts-analysis)
- [Motifs](#motifs)
- [G4 structural subclasses](#g4-structural-subclasses)
- [Hierarchical overlap of peaks with genomic features](#hierarchical-overlap-of-peaks-with-genomic-features)


## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- Standard Unix tools: cat, awk, sort, uniq, paste, grep, pigz, cut, sbatch, nohup, sed
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [deeptools v3.3.0](https://deeptools.readthedocs.io/en/develop/)
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [piranha v1.2.1](https://github.com/smithlabcode/piranha)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)
- [fastaRegexFinder.py v0.1.1](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder)
- [gat-run.py](http://gat.readthedocs.io/en/latest/contents.html)
- [meme suite v4.11.2](http://meme-suite.org/)
- [Python v2.7.12](https://www.python.org/)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [VennDiagram v1.6.20](https://cran.r-project.org/web/packages/VennDiagram/index.html)
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [GenomicFeatures v1.26.4](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)


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
cd ~/fastq

mkdir ../fastq_trimmed

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4G --wrap "cutadapt -a AGATCGGAAGAGC -m 10 -q 10 -O 5 -o ../fastq_trimmed/$fq $fq > ../fastq_trimmed/$bname.txt"
done
```


## Remove duplicates and reads with stretches of several identical nucleotides

```bash
cd ~/fastq_trimmed

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
cd ~/fastq_trimmed_nodup

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
cd ~/bam

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
cd ~/bam

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


## Genomic signal normalization

```bash
cd ~/bam

mkdir ../bw

for bam in *.clean.bam
do
  bname=${bam%.bam}
  sbatch -J $bname --mem 32G --wrap "bamCoverage -b $bam -o ../bw/$bname.bw -of bigwig --binSize 10 -p 20 --normalizeUsing CPM"
done
```


## Peak calling

```bash
cd ~/bam

mkdir ../piranha

for bam in *.clean.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../piranha/$bname.log --mem 8G --wrap "bedtools bamtobed -i $bam > ../piranha/${bname}.bed && \
  Piranha ../piranha/${bname}.bed -o ../piranha/${bname}.peaks.bed -s -p 0.0001 -b 200 -u 0 -d ZeroTruncatedPoisson && \
  rm ../piranha/${bname}.bed"
done
```

Consensus peaks:

```bash
cd ~/piranha

# DHX36 WT - 3 in 3
tableCat.py -i DHX36_{WT1,WT2,WT3}.clean.peaks.bed | \
grep "^chr*" | \
sort -k1,1 -k2,2n | \
awk -v OFS="\t" 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1, $2, $3, $4, $5, $6, abs($7), $8}' | \
bedtools merge -s -d -1 -c 6,7,8,8 -o distinct,min,distinct,count_distinct -i - | \
awk -v OFS="\t" '$7 > 2 {print $1, $2, $3, $5, $7, $4}' > DHX36_WT.clean.peaks.consensus.bed

# DHX36 EA - 3 in 3
tableCat.py -i DHX36_{EA1,EA2,EA3}.clean.peaks.bed | \
grep "^chr*" | \
sort -k1,1 -k2,2n | \
awk -v OFS="\t" 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1, $2, $3, $4, $5, $6, abs($7), $8}' | \
bedtools merge -s -d -1 -c 6,7,8,8 -o distinct,min,distinct,count_distinct -i - | \
awk -v OFS="\t" '$7 > 2 {print $1, $2, $3, $5, $7, $4}' > DHX36_EA.clean.peaks.consensus.bed

# GRSF1 WT - 3 in 4
tableCat.py -i GRSF1_{WT1,WT2,WT3,WT4}.clean.peaks.bed | \
grep "^chr*" | \
sort -k1,1 -k2,2n | \
awk -v OFS="\t" 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1, $2, $3, $4, $5, $6, abs($7), $8}' | \
bedtools merge -s -d -1 -c 6,7,8,8 -o distinct,min,distinct,count_distinct -i - | \
awk -v OFS="\t" '$7 > 2 {print $1, $2, $3, $5, $7, $4}' > GRSF1_WT.clean.peaks.consensus.bed

# DDX3X WT - 3 in 3
tableCat.py -i DDX3X_{WT1,WT2,WT3}.clean.peaks.bed | \
grep "^chr*" | \
sort -k1,1 -k2,2n | \
awk -v OFS="\t" 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1, $2, $3, $4, $5, $6, abs($7), $8}' | \
bedtools merge -s -d -1 -c 6,7,8,8 -o distinct,min,distinct,count_distinct -i - | \
awk -v OFS="\t" '$7 > 2 {print $1, $2, $3, $5, $7, $4}' > DDX3X_WT.clean.peaks.consensus.bed
```

Union of DHX36_WT, GRSF1_WT and DDX3X_WT peaks:

```bash
cd ~/piranha

cat DHX36_{WT1,WT2,WT3}.clean.peaks.bed GRSF1_{WT1,WT2,WT3,WT4}.clean.peaks.bed DDX3X_{WT1,WT2,WT3}.clean.peaks.bed | \
grep "^chr*" | \
sort -k1,1 -k2,2n | \
bedtools merge -s -d -1 -c 6 -o distinct -i - > DHX36_WT_GRSF1_WT_DDX3X_WT_union.clean.peaks.bed
```


## Venn diagram

Intersections:

```bash
cd ~/piranha

# DHX36 EA and WT intersection
bedtools intersect \
-a DHX36_EA.clean.peaks.consensus.bed \
-b DHX36_WT.clean.peaks.consensus.bed \
-wa -u > DHX36_EAiWT.clean.peaks.consensus.bed

# DHX36 WT and EA intersection
bedtools intersect \
-a DHX36_WT.clean.peaks.consensus.bed \
-b DHX36_EA.clean.peaks.consensus.bed \
-wa -u > DHX36_WTiEA.clean.peaks.consensus.bed

# DHX36 EA unique
bedtools intersect \
-a DHX36_EA.clean.peaks.consensus.bed \
-b DHX36_WT.clean.peaks.consensus.bed \
-v > DHX36_EAu.clean.peaks.consensus.bed

# DHX36 WT unique
bedtools intersect \
-a DHX36_WT.clean.peaks.consensus.bed \
-b DHX36_EA.clean.peaks.consensus.bed \
-v > DHX36_WTu.clean.peaks.consensus.bed
```

Venn diagrams:

```r
library(VennDiagram)

# EAu, EAiWT and WTu - all peaks
venn.plot <- draw.pairwise.venn(
  area1 = 24290,
  area2 = 10891,
  cross.area = 9075,
  category = c("EA\n (24290)", "WT\n (10891)"),
  ext.percent = 0.1,
  fill = c("darkgoldenrod1", "darkgoldenrod"),
  cex = 2,
  fontfamily = "sans",
  cat.pos = c(-50, 30),
  cat.dist = c(0.08, 0.08),
  cat.cex = 2,
  cat.fontfamily = "sans",
  ext.pos = 90,
  ext.dist = -0.05,
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  margin = 0.075)

pdf("~/figures/EAu_EAiWT_WTu_all.pdf")
g <- grid.draw(venn.plot)
dev.off()
```


## GAT

Create PQS files:

```bash
cd ~/annotation

ref=~/reference/GRCh38.p12.genome.clean.fa

nohup fastaRegexFinder.py -f $ref -q | \
bedtools sort -i | \
sed 's/$/&\tG3L7/' > GRCh38.p12.genome.clean.g3l7.bed &

nohup fastaRegexFinder.py -f $ref -r '([gG]{3,}\w{1,12}){3,}[gG]{3,}' -q | \
bedtools sort -i | \
sed 's/$/&\tG3L12/' > GRCh38.p12.genome.clean.g3l12.bed &

nohup fastaRegexFinder.py -f $ref -r '([gG]{2,}\w{1,7}){3,}[gG]{2,}' -q | \
bedtools sort -i | \
sed 's/$/&\tG2L7/' > GRCh38.p12.genome.clean.g2l7.bed &

nohup fastaRegexFinder.py -f $ref -r '([gG]{2,}\w{1,12}){3,}[gG]{2,}' -q | \
bedtools sort -i | \
sed 's/$/&\tG2L12/' > GRCh38.p12.genome.clean.g2l12.bed &
```

Create genomics feature files:

```r
library(data.table)
library(GenomicFeatures)

# change width
options(width = 250)

# prepare coordinates table
txdb <- makeTxDbFromGFF("~/reference/gencode.v28.annotation.sorted.gtf", format="gtf")

# genes
genes <- data.table(data.frame(genes(txdb)))[, c("seqnames", "start", "end", "gene_id", "strand")][order(seqnames, start)]
genes <- genes[, dummy := 0][, c("seqnames", "start", "end", "gene_id", "dummy", "strand")]
write.table(genes, "~/reference/gencode.v28.annotation.sorted.genes.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# transcripts
transcripts <- data.table(data.frame(transcripts(txdb, columns=c("tx_name", "gene_id"))))[,.(seqnames, start, end, tx_name, gene_id, strand)]
write.table(transcripts, "~/reference/gencode.v28.annotation.sorted.transcripts.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# gene promoters
gene_promoters <- data.table(data.frame(promoters(genes(txdb), upstream=1000, downstream=0)))[, c("seqnames", "start", "end", "gene_id", "strand")][order(seqnames, start)]
gene_promoters <- gene_promoters[, dummy := 0][, c("seqnames", "start", "end", "gene_id", "dummy", "strand")]
write.table(gene_promoters, "~/reference/gencode.v28.annotation.sorted.gene_promoters.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# transcript 5'UTR
transcript_utr5 <- data.table(data.frame(fiveUTRsByTranscript(txdb, use.names = TRUE)))[, c("seqnames", "start", "end", "group_name", "strand")][order(seqnames, start)]
transcript_utr5 <- transcript_utr5[, dummy := 0][, c("seqnames", "start", "end", "group_name", "dummy", "strand")]
write.table(transcript_utr5, "~/reference/gencode.v28.annotation.sorted.transcript_utr5.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# transcript 3'UTR
transcript_utr3 <- data.table(data.frame(threeUTRsByTranscript(txdb, use.names = TRUE)))[, c("seqnames", "start", "end", "group_name", "strand")][order(seqnames, start)]
transcript_utr3 <- transcript_utr3[, dummy := 0][, c("seqnames", "start", "end", "group_name", "dummy", "strand")]
write.table(transcript_utr3, "~/reference/gencode.v28.annotation.sorted.transcript_utr3.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# transcript exons
transcript_exons <- data.table(data.frame(exonsBy(txdb, by = "tx", use.names = TRUE)))[, c("seqnames", "start", "end", "group_name", "strand")][order(seqnames, start)]
transcript_exons <- transcript_exons[, dummy := 0][, c("seqnames", "start", "end", "group_name", "dummy", "strand")]
write.table(transcript_exons, "~/reference/gencode.v28.annotation.sorted.transcript_exons.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# transcript introns
transcript_introns <- data.table(data.frame(intronsByTranscript(txdb, use.names = TRUE)))[, c("seqnames", "start", "end", "group_name", "strand")][order(seqnames, start)]
transcript_introns <- transcript_introns[, dummy := 0][, c("seqnames", "start", "end", "group_name", "dummy", "strand")]
write.table(transcript_introns, "~/reference/gencode.v28.annotation.sorted.transcript_introns.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```

Combine annotation files:

```bash
cd ~/annotation

cat \
<(grep -h "^chr" GRCh38.p12.genome.clean.g*l*.bed | cut -f1-3,8) \
<(cut -f1-3 gencode.v28.annotation.sorted.gene_promoters.bed | sed 's/$/&\tpromoter/') \
<(cut -f1-3 gencode.v28.annotation.sorted.transcript_utr5.bed | sed 's/$/&\t5UTR/') \
<(cut -f1-3 gencode.v28.annotation.sorted.transcript_utr3.bed | sed 's/$/&\t3UTR/') \
<(cut -f1-3 gencode.v28.annotation.sorted.transcript_exons.bed | sed 's/$/&\texon/') \
<(cut -f1-3 gencode.v28.annotation.sorted.transcript_introns.bed | sed 's/$/&\tintron/') | \
awk -v OFS="\t" '$2 > -1' | \
bedtools sort -i - > GRCh38.p12.genome.clean.gencode.v28.annotation.sorted.bed
```

Run gat-run.py:

```bash
cd ~/piranha

mkdir ../gat/

annotations=~/annotation/GRCh38.p12.genome.clean.gencode.v28.annotation.sorted.bed
mappable_g=~/annotation/gencode.v28.annotation.sorted.genes.bed

for bed in DHX36_EA.clean.peaks.consensus.bed \
DHX36_WT.clean.peaks.consensus.bed \
GRSF1_WT.clean.peaks.consensus.bed \
DDX3X_WT.clean.peaks.consensus.bed
do
  bname=`basename ${bed%.bed}`
  nohup gat-run.py -a $annotations -s $bed -w <(cut -f1-3 $mappable_t) --ignore-segment-tracks -n 10000 -L ../gat/$bname.transcripts.log > ../gat/$bname.transcripts.txt &
  nohup gat-run.py -a $annotations -s $bed -w <(cut -f1-3 $mappable_g) --ignore-segment-tracks -n 10000 -L ../gat/$bname.genes.log > ../gat/$bname.genes.txt &
done
```

Plotting:

```r
library(data.table)
library(ggplot2)

# Enlarge the view width when printing tables
options(width = 400)

# Load data
data <- fread("tableCat.py -i ~/gat/*.transcripts.txt -S 1 -r .clean.peaks.consensus.transcripts.txt | cut -f2-6,8-12,15,18,21-25")
setnames(data, c("annotation", "observed", "expected", "CI95low", "CI95high", "fold", "l2fold", "pvalue", "qvalue", "track_nsegments", "annotation_nsegments", "overlap_nsegments", "percent_overlap_nsegments_track", "percent_overlap_size_track", "percent_overlap_nsegments_annotation", "percent_overlap_size_annotation", "library"))


# Boxplot all
data_selected <- data[annotation %in% c("5UTR", "exon", "intron", "3UTR", "G2L7", "G3L7")]
data_selected$annotation <- factor(data_selected$annotation, levels = c("5UTR", "exon", "intron", "3UTR", "G2L7", "G3L7"))
data_selected$library <- factor(data_selected$library, levels = c("DHX36_EA", "DHX36_WT", "GRSF1_WT", "DDX3X_WT"))

gg <- ggplot(data_selected, aes(x=annotation, y=fold, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab("Fold Change) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("darkgoldenrod1", "darkgoldenrod", "deepskyblue3", "seagreen3", "seagreen4"), labels = c("DHX36 EA", "DHX36 WT", "GRSF1 WT", "DDX3X_WT")) +
geom_errorbar(aes(ymin=observed/CI95low,ymax=observed/CI95high),position=position_dodge(0.9),width=0.2)
ggsave("~/figures/dhx36_grsf1_ddx3x_gat.pdf", width = 20, height = 14, units= 'cm')


# Boxplot genomic features DHX36 WT and EA
data_genomic_dhx36 <- data[annotation %in% c("5UTR", "exon", "intron", "3UTR") & library %in% c("DHX36_EA", "DHX36_WT")]
data_genomic_dhx36[annotation == "5UTR", annotation := "5'UTR"]
data_genomic_dhx36[annotation == "3UTR", annotation := "3'UTR"]
data_genomic_dhx36$annotation <- factor(data_genomic_dhx36$annotation, levels = c("5'UTR", "exon", "intron", "3'UTR"))
data_genomic_dhx36$library <- factor(data_genomic_dhx36$library, levels = c("DHX36_EA", "DHX36_WT"))

gg <- ggplot(data_genomic_dhx36, aes(x=annotation, y=fold, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab("Fold Change") +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("darkgoldenrod1", "darkgoldenrod"), labels = c("EA", "WT")) +
geom_errorbar(aes(ymin=observed/CI95low,ymax=observed/CI95high),position=position_dodge(0.9),width=0.2)
ggsave("~/figures/dhx36_genomic_gat.pdf", width = 10, height = 10, units= 'cm')


# Boxplot PQS DHX36 WT and EA
data_pqs_dhx36 <- data[annotation %in% c("G2L7", "G3L7") & library %in% c("DHX36_EA", "DHX36_WT")]
data_pqs_dhx36$library <- factor(data_pqs_dhx36$library, levels = c("DHX36_EA", "DHX36_WT"))

gg <- ggplot(data_pqs_dhx36, aes(x=annotation, y=fold, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab("Fold Change") +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("darkgoldenrod1", "darkgoldenrod"), labels = c("EA", "WT")) +
geom_errorbar(aes(ymin=observed/CI95low,ymax=observed/CI95high),position=position_dodge(0.9),width=0.2)
ggsave("~/figures/dhx36_pqs_gat.pdf", width = 8, height = 10, units= 'cm')


# Boxplot genomic features GRSF1 WT
data_genomic_grsf1 <- data[annotation %in% c("5UTR", "exon", "intron", "3UTR") & library %in% c("GRSF1_WT")]
data_genomic_grsf1[annotation == "5UTR", annotation := "5'UTR"]
data_genomic_grsf1[annotation == "3UTR", annotation := "3'UTR"]
data_genomic_grsf1$annotation <- factor(data_genomic_grsf1$annotation, levels = c("5'UTR", "exon", "intron", "3'UTR"))

gg <- ggplot(data_genomic_grsf1, aes(x=annotation, y=fold, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5, show.legend = FALSE) +
theme_classic() +
ylab("Fold Change") +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("deepskyblue3")) +
geom_errorbar(aes(ymin=observed/CI95low,ymax=observed/CI95high),position=position_dodge(0.9),width=0.2)
ggsave("~/figures/grsf1_transcript_gat.pdf", width = 7, height = 10, units= 'cm')


# Boxplot PQS GRSF1 WT
data_pqs_grsf1 <- data[annotation %in% c("G2L7", "G3L7") & library %in% c("GRSF1_WT")]

gg <- ggplot(data_pqs_grsf1, aes(x=annotation, y=fold, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5, show.legend = FALSE) +
theme_classic() +
ylab("Fold Change") +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("deepskyblue3")) +
geom_errorbar(aes(ymin=observed/CI95low,ymax=observed/CI95high),position=position_dodge(0.9),width=0.2)
ggsave("~/figures/grsf1_pqs_gat.pdf", width = 5, height = 10, units= 'cm')


# Boxplot genomic features DDX3X WT
data_genomic_ddx3x <- data[annotation %in% c("5UTR", "exon", "intron", "3UTR") & library %in% c("DDX3X_WT")]
data_genomic_ddx3x[annotation == "5UTR", annotation := "5'UTR"]
data_genomic_ddx3x[annotation == "3UTR", annotation := "3'UTR"]
data_genomic_ddx3x$annotation <- factor(data_genomic_ddx3x$annotation, levels = c("5'UTR", "exon", "intron", "3'UTR"))

gg <- ggplot(data_genomic_ddx3x, aes(x=annotation, y=fold, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab("Fold Change") +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("seagreen3", "seagreen4"), labels = c("WT")) +
geom_errorbar(aes(ymin=observed/CI95low,ymax=observed/CI95high),position=position_dodge(0.9),width=0.2)
ggsave("~/figures/ddx3x_genomic_gat.pdf", width = 10, height = 10, units= 'cm')


# Boxplot PQS DDX3X WT
data_pqs_ddx3x <- data[annotation %in% c("G2L7", "G3L7") & library %in% c("DDX3X_WT")]
data_pqs_ddx3x$library <- factor(data_pqs_ddx3x$library, levels = c("DDX3X_WT"))

gg <- ggplot(data_pqs_ddx3x, aes(x=annotation, y=l2fold, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab(expression("log"[2]*"FC")) +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("seagreen3", "seagreen4"), labels = c("mRG", "WT")) +
geom_errorbar(aes(ymin=observed/CI95low,ymax=observed/CI95high),position=position_dodge(0.9),width=0.2)
ggsave("~/figures/ddx3x_pqs_gat.pdf", width = 8, height = 10, units= 'cm')
```


## Overlap with genomic features and PQS

```bash
cd ~/piranha

for peaks in DHX36_EA.clean.peaks.consensus.bed \
DHX36_WT.clean.peaks.consensus.bed \
GRSF1_WT.clean.peaks.consensus.bed \
DDX3X_WT.clean.peaks.consensus.bed
do
  t1=`cat $peaks | wc -l`
  echo -e "$peaks\t$t1"
  for pqs in ../annotation/GRCh38.p12.genome.clean.g*l*.bed
  do
    t2=`cat $pqs | wc -l`
    o1=`bedtools intersect \
    -a $peaks \
    -b <(grep -h "^chr" $pqs) \
    -wa -u -s | wc -l`
    pct1=`echo "scale=1; 100*$o1/$t1" | bc`
    o2=`bedtools intersect \
    -b $peaks \
    -a <(grep -h "^chr" $pqs) \
    -wa -u -s | wc -l`
    pct2=`echo "scale=1; 100*$o2/$t2" | bc`
    echo -e "`basename $pqs`\t$o1\t$pct1%\t$o2\t$pct2%"
  done
  for ann in ../annotation/gencode.v28.annotation.sorted.transcript_*.bed
  do
    t2=`cat $ann | wc -l`
    o1=`bedtools intersect \
    -a $peaks \
    -b $ann \
    -wa -u -s | wc -l`
    pct1=`echo "scale=1; 100*$o1/$t1" | bc`
    o2=`bedtools intersect \
    -b $peaks \
    -a $ann \
    -wa -u -s | wc -l`
    pct2=`echo "scale=1; 100*$o2/$t2" | bc`
    echo -e "`basename $ann`\t$o1\t$pct1%\t$o2\t$pct2%"
  done
  echo -e "----------"
done | column -t
```

Plotting:

```r
library(data.table)
library(ggplot2)

# Enlarge the view width when printing tables
options(width = 400)

# Barplot percentage overlap DHX36 WT and EA
# Load data
data_dhx36 <- data.table(annotation = c("G2L7", "G3L7", "G2L7", "G3L7"), library = c("DHX36_EA", "DHX36_EA", "DHX36_WT", "DHX36_WT"), pct_overlap = c(52.6, 8.7, 50.9, 8.2))

data_dhx36$library <- factor(data_dhx36$library, levels = c("DHX36_EA", "DHX36_WT"))

gg <- ggplot(data_dhx36, aes(x=annotation, y=pct_overlap, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab("% overlap") +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("darkgoldenrod1", "darkgoldenrod"), labels = c("EA", "WT")) +
coord_cartesian(ylim = c(0, 100))
ggsave("~/figures/dhx36_pqs_pct_overlap.pdf", width = 8, height = 10, units= 'cm')


# Barplot percentage overlap GRSF1
# Load data
data_grsf1 <- data.table(annotation = c("G2L7", "G3L7"), pct_overlap = c(63.1, 12.4))

gg <- ggplot(data_grsf1, aes(x=annotation, y=pct_overlap)) +
geom_bar(stat="identity", color="black", fill = "deepskyblue3", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab("% overlap") +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
coord_cartesian(ylim = c(0, 100))
ggsave("~/figures/grsf1_pqs_pct_overlap.pdf", width = 5, height = 10, units= 'cm')
```


## Binding profiles around PQS

### Merge bam files and normalise signal by CPM

```bash
cd ~/bam

mkdir ../bam_merge
mkdir ../bw_merge

for id in DDX3X_WT DHX36_EA DHX36_WT GRSF1_WT
do
  bams=`echo $id*.clean.bam`
  sbatch -J $id -o ../bam_merge/$id.log --mem 8G --wrap "samtools merge -@ 20 ../bam_merge/$id.clean.bam $bams && \
  samtools index ../bam_merge/$id.clean.bam && \
  bamCoverage -b ../bam_merge/$id.clean.bam -o ../bw_merge/$id.clean.bw -of bigwig --binSize 1 -p 20 --normalizeUsing CPM"
done
```

### Select PQSs with overlap

```bash
cd ~/piranha

for pqs in ../annotation/GRCh38.p12.genome.clean.g3l7.bed \
../annotation/GRCh38.p12.genome.clean.g3l12.bed \
../annotation/GRCh38.p12.genome.clean.g2l7.bed \
../annotation/GRCh38.p12.genome.clean.g2l12.bed
do
  bname_pqs=`basename ${pqs%.bed}`
  for peaks in *.clean.peaks.consensus.bed
  do
    bname_peaks=${peaks%.bed}
    nohup bedtools intersect \
    -a <(grep -h "^chr" $pqs) \
    -b $peaks \
    -wa -u -s > ../deeptools/$bname_pqs.$bname_peaks.bed &
  done
done
```

### deeptools

```bash
cd ~/bw_merge

for bw in *.clean.bw
do
  bname=${bw%.clean.bw}
  sbatch -J $bname -o ../deeptools/$bname.log --mem 8G --wrap "computeMatrix scale-regions -S $bw -R ../deeptools/GRCh38.p12.genome.clean.g3l7.$bname.clean.peaks.consensus.bed ../deeptools/GRCh38.p12.genome.clean.g3l12.$bname.clean.peaks.consensus.bed ../deeptools/GRCh38.p12.genome.clean.g2l7.$bname.clean.peaks.consensus.bed ../deeptools/GRCh38.p12.genome.clean.g2l12.$bname.clean.peaks.consensus.bed -out ../deeptools/$bname.mat.gz -m 50 --startLabel 's' --endLabel 'e' -b 500 -a 500 -bs 1 --sortRegions no --skipZeros -p 'max' && \
  plotProfile -m ../deeptools/$bname.mat.gz -o ../deeptools/20190704_$bname.png --dpi 300 --plotHeight 9 --plotWidth 10 --plotType 'se' --numPlotsPerRow 2 --startLabel 's' --endLabel 'e' --regionsLabel 'G3L7' 'G3L12' 'G2L7' 'G2L12' --samplesLabel '' --plotTitle $bname --colors 'black' 'black' 'black' 'black' -y 'CPM' --legendLocation none --perGroup"
done
```


## Transcripts analysis

### Parsing gff

```python
ifile = open("~/annotation/gencode.v28.annotation.sorted.gtf")
ilines = ifile.readlines()
ifile.close()

ofile = open("~/annotation/gencode.v28.annotation.sorted.transcripts.type.bed", "w")

for l in ilines:
  if not l.startswith("##"):
    fields = l.split("\t")
    if fields[2] == "transcript":
      chr = fields[0]
      start = fields[3]
      end = fields[4]
      strand = fields[6]
      annotations = fields[8].split(";")[:-1]
      gene_id = [a.split('"')[1] for a in annotations if a.startswith("gene_id")][0]
      transcript_id = [a.split('"')[1] for a in annotations if a.startswith(" transcript_id")][0]
      transcript_type = [a.split('"')[1] for a in annotations if a.startswith(" transcript_type")][0]
      ofile.write("\t".join([chr, start, end, transcript_id, gene_id, strand, transcript_type])+"\n")

ofile.close()
```

### Overlap with annotated transcripts

```bash
cd ~/piranha

for peaks in DHX36_EA.clean.peaks.consensus.bed \
DHX36_WT.clean.peaks.consensus.bed \
GRSF1_WT.clean.peaks.consensus.bed \
DDX3X_WT.clean.peaks.consensus.bed
do
  # calc
  n=`cat $peaks | wc -l`
  t=`bedtools intersect -a $peaks -b ../annotation/gencode.v28.annotation.sorted.transcripts.type.bed -loj -s | wc -l`
  c=`bedtools intersect -a $peaks -b ../annotation/gencode.v28.annotation.sorted.transcripts.type.bed -loj -s | awk '$13 == "protein_coding" || $13 == "retained_intron" || $13 == "nonsense_mediated_decay" || $13 == "processed_transcript"' | wc -l`
  c_pct=`echo "scale=1; 100*$c/$t" | bc`
  nc=`bedtools intersect -a $peaks -b ../annotation/gencode.v28.annotation.sorted.transcripts.type.bed -loj -s | awk '$7 != "." && $13 != "protein_coding" && $13 != "retained_intron" && $13 != "nonsense_mediated_decay" && $13 != "processed_transcript"' | wc -l`
  nc_pct=`echo "scale=1; 100*$nc/$t" | bc`
  o=`bedtools intersect -a $peaks -b ../annotation/gencode.v28.annotation.sorted.transcripts.type.bed -loj -s | awk '$7 == "."' | wc -l`
  o_pct=`echo "scale=1; 100*$o/$t" | bc`
  nc_types=`bedtools intersect -a $peaks -b ../annotation/gencode.v28.annotation.sorted.transcripts.type.bed -loj -s | awk '$7 != "." && $13 != "protein_coding" && $13 != "retained_intron" && $13 != "nonsense_mediated_decay" && $13 != "processed_transcript"' | cut -f13 | sort | uniq -c | sort -k1,1nr`
  # print
  echo -e "### $peaks ###\n"
  echo -e "total number of peaks:\t$n\n"
  echo -e "total number of binding sites:\t$t"
  echo -e "\t- coding:\t$c\t($c_pct%)"
  echo -e "\t- non-coding:\t$nc\t($nc_pct%)"
  echo -e "\t- other:\t$o\t($o_pct%)\n"
  echo -e "most common non-coding binding sites are:"
  echo -e "$nc_types"
  echo -e "-----------------------------------------\n"
done
```


## Motifs

### Overlap

```bash
cd ~/piranha

# 5'UTR
for bed in *.clean.peaks.consensus.bed
do
  bname=${bed%.bed}
  nohup bedtools intersect \
  -a $bed \
  -b ~/annotation/gencode.v28.annotation.sorted.transcript_utr5.bed \
  -wa -u > $bname.utr5.bed &  
done

# 3'UTR
for bed in *.clean.peaks.consensus.bed
do
  bname=${bed%.bed}
  nohup bedtools intersect \
  -a $bed \
  -b ~/annotation/gencode.v28.annotation.sorted.transcript_utr3.bed \
  -wa -u > $bname.utr3.bed &  
done
```

### dreme

Create fasta files:

```bash
cd ~/piranha

mkdir ../dreme

g=~/reference/GRCh38.p12.genome.clean.fa

# bed files
for bed in *.clean.peaks.consensus*.bed
do
  bname=${bed%.bed}
  echo $bname
  nohup bedtools sort -i $bed | \
  bedtools getfasta -fi $g -bed - -s > ../dreme/$bname.fasta &
  nohup bedtools sort -i $bed | \
  bedtools shuffle -i - -g <(grep '^chr' $g.fai) -seed 123 -incl ../annotation/gencode.v28.annotation.sorted.genes.bed | \
  bedtools getfasta -fi $g -bed - > ../dreme/$bname.shuffle.fasta &
done
```

Run dreme with negative sequences:

```bash
cd ~/dreme

mkdir run_neg

for k in 8 10 15 20 25 30
do
  mkdir $k
  for fasta in `ls *.fasta | grep -v "shuffle"`
  do
    bname=${fasta%.fasta}
    sbatch -J $k.${bname} -o run_neg/${k}/${bname}.log --mem 32G --wrap "dreme -oc run_neg/${k}/${bname} -p $fasta -n $bname.shuffle.fasta -rna -norc -s 123 -png -maxk $k"
  done
done
```

Run dreme with scrambled negative sequences:

```bash
cd ~/dreme

mkdir -p run_scr

for k in 8 10 15 20 25 30
do
  mkdir run_scr/$k
  for fasta in `ls *.fasta | grep -v "shuffle"`
  do
    bname=${fasta%.fasta}
    sbatch -J $k.${bname} -o run_scr/${k}/${bname}.log --mem 32G --wrap "dreme -oc run_scr/${k}/${bname} -p $fasta -rna -norc -s 123 -png -maxk $k"
  done
done
```

### meme-chip

```bash
cd ~/dreme

mkdir ../meme-chip

for fasta in `ls *.fasta | grep -v "shuffle"`
do
  bname=${fasta%.fasta}
  sbatch -J ${bname} -o ../meme-chip/${bname}.log --mem 32G --wrap "meme-chip -oc ../meme-chip/${bname} -neg $bname.shuffle.fasta -seed 123 -norc -spamo-skip -fimo-skip $fasta"
done
```


## G4 structural subclasses

```bash
cd ~/annotation

ref=~/reference/GRCh38.p12.genome.clean.fa

# Canonical G4s: G3L7
#GRCh38.p12.genome.clean.g3l7.bed

# Long loops: G3L12
#GRCh38.p12.genome.clean.g3l12.bed

# Bulges
nohup fastaRegexFinder.py -f $ref -r '([gG]{1}[aAcCtT]{2,7}[gG]{2,}|[gG]{2}[aAcCtT]{2,7}[gG]{1,})\w{1,7}[gG]{3,}\w{1,7}[gG]{3,}\w{1,7}[gG]{3,}|
[gG]{3,}\w{1,7}([gG]{1}[aAcCtT]{2,7}[gG]{2,}|[gG]{2}[aAcCtT]{2,7}[gG]{1,})\w{1,7}[gG]{3,}\w{1,7}[gG]{3,}|
[gG]{3,}\w{1,7}[gG]{3,}\w{1,7}([gG]{1}[aAcCtT]{2,7}[gG]{2,}|[gG]{2}[aAcCtT]{2,7}[gG]{1,})\w{1,7}[gG]{3,}|
[gG]{3,}\w{1,7}[gG]{3,}\w{1,7}[gG]{3,}\w{1,7}([gG]{1}[aAcCtT]{2,7}[gG]{2,}|[gG]{2}[aAcCtT]{2,7}[gG]{1,})|
([gG]{1}[aAcCtT]{1}[gG]{2,}|[gG]{2}[aAcCtT]{1,}[gG]{2,}|[gG]{3,})\w{1,7}([gG]{1}[aAcCtT]{1}[gG]{2,}|[gG]{2}[aAcCtT]{1,}[gG]{2,}|[gG]{3,})\w{1,7}([gG]{1}[aAcCtT]{1}[gG]{2,}|[gG]{2}[aAcCtT]{1,}[gG]{2,}|[gG]{3,})\w{1,7}([gG]{1}[aAcCtT]{1}[gG]{2,}|[gG]{2}[aAcCtT]{1,}[gG]{2,}|[gG]{3,})' -q | \
bedtools sort -i | \
sed 's/$/&\tbulges/' > GRCh38.p12.genome.clean.bulges.bed &

# 2-tetrads
#GRCh38.p12.genome.clean.g2l12.bed

## G3L12 - G3L7 ##
bedtools intersect \
-a <(bedtools sort -i GRCh38.p12.genome.clean.g3l12.bed | grep "^chr") \
-b <(bedtools sort -i GRCh38.p12.genome.clean.g3l7.bed | grep "^chr") \
-v -s > GRCh38.p12.genome.clean.g3l12minusg3l7.bed

## G3L7 and G3L12 ##
cat <(grep "^chr" GRCh38.p12.genome.clean.g3l7.bed) <(grep "^chr" GRCh38.p12.genome.clean.g3l12minusg3l7.bed) | bedtools sort -i > GRCh38.p12.genome.clean.g3l7.g3l12.bed

## G3L7 and G3L12 - bulges ##
bedtools intersect \
-a <(bedtools sort -i GRCh38.p12.genome.clean.bulges.bed | grep "^chr") \
-b <(bedtools sort -i GRCh38.p12.genome.clean.g3l7.g3l12.bed | grep "^chr") \
-v -s > GRCh38.p12.genome.clean.g3l7.g3l12minusbulges.bed

## G3L7 and G3L12 and bulges ##
cat <(grep "^chr" GRCh38.p12.genome.clean.g3l7.g3l12.bed) <(grep "^chr" GRCh38.p12.genome.clean.g3l7.g3l12minusbulges.bed) | bedtools sort -i > GRCh38.p12.genome.clean.g3l7.g3l12.bulges.bed

## G3L7 and G3L12 and bulges - G2L12 ##
bedtools intersect \
-a <(bedtools sort -i GRCh38.p12.genome.clean.g2l12.bed | grep "^chr") \
-b <(bedtools sort -i GRCh38.p12.genome.clean.g3l7.g3l12.bulges.bed | grep "^chr") \
-v -s > GRCh38.p12.genome.clean.g3l7.g3l12.bulgesminusg2l12.bed

## G3L7 and G3L12 and bulges and G2L12 ##
cat <(grep "^chr" GRCh38.p12.genome.clean.g3l7.g3l12.bulges.bed) <(grep "^chr" GRCh38.p12.genome.clean.g3l7.g3l12.bulgesminusg2l12.bed) | bedtools sort -i > GRCh38.p12.genome.clean.g3l7.g3l12.bulges.g2l12.bed
```

### Overlap with G4 structural subclasses

```bash
cd ~/piranha

for peaks in DHX36_EA.clean.peaks.consensus.bed \
DHX36_WT.clean.peaks.consensus.bed \
GRSF1_WT.clean.peaks.consensus.bed \
DDX3X_WT.clean.peaks.consensus.bed
do
  # calc
  n=`cat $peaks | wc -l`
  t=`bedtools intersect -a $peaks -b ../annotation/GRCh38.p12.genome.clean.g3l7.g3l12.bulges.g2l12.bed -loj -s | wc -l`
  c=`bedtools intersect -a $peaks -b ../annotation/GRCh38.p12.genome.clean.g3l7.g3l12.bulges.g2l12.bed -loj -s | cut -f14 | sort | uniq -c`
  # print
  echo -e "### $peaks ###\n"
  echo -e "total number of peaks:\t$n\n"
  echo -e "total number of overlaps:\t$t"
  echo -e "$c"
done
```

The legend is:

- G3L7 as canonical
- G3L12 as long loops
- bulges
- G2L12 as 2 quartets
- . as other (none of the above)


## Hierarchical overlap of peaks with genomic features

5'UTR > 3'UTR > exons > introns

```bash
cd ~/piranha

for peaks in DHX36_EA.clean.peaks.consensus.bed \
DHX36_WT.clean.peaks.consensus.bed \
GRSF1_WT.clean.peaks.consensus.bed \
DDX3X_WT.clean.peaks.consensus.bed
do
  t1=`cat $peaks | wc -l`
  echo -e "$peaks\t$t1"
  o_utr5=`bedtools intersect \
  -a $peaks \
  -b ../annotation/gencode.v28.annotation.sorted.transcript_utr5.bed \
  -wa -u -s | wc -l`
  o_utr3=`bedtools intersect \
  -a <(bedtools intersect -a $peaks -b ../annotation/gencode.v28.annotation.sorted.transcript_utr5.bed -v -s) \
  -b ../annotation/gencode.v28.annotation.sorted.transcript_utr3.bed \
  -wa -u -s | wc -l`
  o_exons=`bedtools intersect \
  -a <(bedtools intersect -a <(bedtools intersect -a $peaks -b ../annotation/gencode.v28.annotation.sorted.transcript_utr5.bed -v -s) -b ../annotation/gencode.v28.annotation.sorted.transcript_utr3.bed -v -s) \
  -b ../annotation/gencode.v28.annotation.sorted.transcript_exons.bed \
  -wa -u -s | wc -l`
  o_introns=`bedtools intersect \
  -a <(bedtools intersect -a <(bedtools intersect -a <(bedtools intersect -a $peaks -b ../annotation/gencode.v28.annotation.sorted.transcript_utr5.bed -v -s) -b ../annotation/gencode.v28.annotation.sorted.transcript_utr3.bed -v -s) -b ../annotation/gencode.v28.annotation.sorted.transcript_exons.bed -v -s) \
  -b ../annotation/gencode.v28.annotation.sorted.transcript_introns.bed \
  -wa -u -s | wc -l`
  o_other=`bedtools intersect \
  -a <(bedtools intersect -a <(bedtools intersect -a <(bedtools intersect -a $peaks -b ../annotation/gencode.v28.annotation.sorted.transcript_utr5.bed -v -s) -b ../annotation/gencode.v28.annotation.sorted.transcript_utr3.bed -v -s) -b ../annotation/gencode.v28.annotation.sorted.transcript_exons.bed -v -s) \
  -b ../annotation/gencode.v28.annotation.sorted.transcript_introns.bed \
  -v -s | wc -l`
  pct_utr5=`echo "scale=1; 100*$o_utr5/$t1" | bc`
  pct_utr3=`echo "scale=1; 100*$o_utr3/$t1" | bc`
  pct_exons=`echo "scale=1; 100*$o_exons/$t1" | bc`
  pct_introns=`echo "scale=1; 100*$o_introns/$t1" | bc`
  pct_other=`echo "scale=1; 100*$o_other/$t1" | bc`
  echo -e "gencode.v28.annotation.sorted.transcript_utr5.bed\t$o_utr5\t$pct_utr5%"
  echo -e "gencode.v28.annotation.sorted.transcript_utr3.bed\t$o_utr3\t$pct_utr3%"
  echo -e "gencode.v28.annotation.sorted.transcript_exons.bed\t$o_exons\t$pct_exons%"
  echo -e "gencode.v28.annotation.sorted.transcript_introns.bed\t$o_introns\t$pct_introns%"
  echo -e "other\t$o_other\t$pct_other%"
  echo -e "----------"
done | column -t
```

Plotting:

```r
library(data.table)
library(ggplot2)

# Enlarge the view width when printing tables
options(width = 400)


# Barplot percentage peak overlap
# Load data
data <- data.table(library = c("DHX36", "DHX36", "DHX36", "DHX36", "DHX36", "GRSF1", "GRSF1", "GRSF1", "GRSF1", "GRSF1", "DDX3X", "DDX3X", "DDX3X", "DDX3X", "DDX3X"), annotation = c("5'UTR", "3'UTR", "exons", "introns", "other", "5'UTR", "3'UTR", "exons", "introns", "other", "5'UTR", "3'UTR", "exons", "introns", "other"), pct_overlap = c(15.7, 42.6, 23.9, 5.1, 12.5, 7.5, 24.3, 15.4, 16.4, 36.1, 47.8, 6.8, 29.2, 4.6, 11.3))

data$library <- factor(data$library, levels = c("DHX36", "GRSF1", "DDX3X"))
data$annotation <- factor(data$annotation, levels = c("5'UTR", "3'UTR", "exons", "introns", "other"))

gg <- ggplot(data, aes(x=annotation, y=pct_overlap, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab("% peak overlap") +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("darkgoldenrod", "deepskyblue3", "seagreen4"), labels = c("DHX36", "GRSF1", "DDX3X")) +
coord_cartesian(ylim = c(0, 50))
ggsave("~/figures/pct_peak_overlap.pdf", width = 14, height = 14, units= 'cm')


# Barplot peak count overlap
# Load data
data <- data.table(library = c("DHX36", "DHX36", "DHX36", "DHX36", "DHX36", "GRSF1", "GRSF1", "GRSF1", "GRSF1", "GRSF1", "DDX3X", "DDX3X", "DDX3X", "DDX3X", "DDX3X"), annotation = c("5'UTR", "3'UTR", "exons", "introns", "other", "5'UTR", "3'UTR", "exons", "introns", "other", "5'UTR", "3'UTR", "exons", "introns", "other"), count_overlap = c(1719, 4640, 2610, 557, 1365, 1249, 4031, 2564, 2732, 5995, 2182, 312, 1333, 213, 517))

data$library <- factor(data$library, levels = c("DHX36", "GRSF1", "DDX3X"))
data$annotation <- factor(data$annotation, levels = c("5'UTR", "3'UTR", "exons", "introns", "other"))

gg <- ggplot(data, aes(x=annotation, y=count_overlap, fill=library)) +
geom_bar(stat="identity", color="black", position=position_dodge(), alpha = 0.5) +
theme_classic() +
ylab("peak count overlap") +
xlab("") +
theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16, color = "black"), axis.text.x = element_text(angle = 45, size = 16, color = "black", hjust = 1), legend.text = element_text(size = 16, color = "black")) +
scale_fill_manual(values = c("darkgoldenrod", "deepskyblue3", "seagreen4"), labels = c("DHX36", "GRSF1", "DDX3X")) +
coord_cartesian(ylim = c(0, 6000))
ggsave("~/figures/peak_count_overlap.pdf", width = 14, height = 14, units= 'cm')
```
