
## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Trim illumina adapters and quality trimming](#trim-illumina-adapters-and-quality-trimming)
- [Remove duplicates and reads with stretches of several identical nucleotides](#remove-duplicates-and-reads-with-stretches-of-several-identical-nucleotides)
- [Alignment](#alignment)
- [Genomic signal correlation](#genomic-signal-correlation)
- [Genomic signal normalization](#genomic-signal-normalization)
- [Peak calling](#peak-calling)
- [Overlap with PQS](#overlap-with-pqs)
- [Binding profiles around PQS](#binding-profiles-around-pqs)



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


- [igvtools v2.3.91](https://software.broadinstitute.org/software/igv/igvtools)
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


## Genomic signal normalization

```bash
cd bam

mkdir ../bw

for bam in *.clean.bam
do
  bname=${bam%.bam}
  sbatch -J $bname --mem 32G --wrap "bamCoverage -b $bam -o ../bw/$bname.bw -of bigwig --binSize 10 -p 20 --normalizeUsing CPM"
done
```


## Peak calling

```bash
cd bam

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
cd ../piranha

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
cd ../piranha

cat DHX36_{WT1,WT2,WT3}.clean.peaks.bed GRSF1_{WT1,WT2,WT3,WT4}.clean.peaks.bed DDX3X_{WT1,WT2,WT3}.clean.peaks.bed | \
grep "^chr*" | \
sort -k1,1 -k2,2n | \
bedtools merge -s -d -1 -c 6 -o distinct -i - > DHX36_WT_GRSF1_WT_DDX3X_WT_union.clean.peaks.bed
```


## Overlap with PQS

Create PQS files:

```bash
cd annotation

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

Overlap:

```bash
cd ../piranha

for peaks in *.consensus.bed
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
  echo -e "----------"
done | column -t
```


## Binding profiles around PQS

### Merge bam files and normalise signal by CPM

```bash
cd bam

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
cd ../piranha

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
cd ../bw_merge

for bw in *.clean.bw
do
  bname=${bw%.clean.bw}
  sbatch -J $bname -o ../deeptools/$bname.log --mem 8G --wrap "computeMatrix scale-regions -S $bw -R ../deeptools/GRCh38.p12.genome.clean.g3l7.$bname.clean.peaks.consensus.bed ../deeptools/GRCh38.p12.genome.clean.g3l12.$bname.clean.peaks.consensus.bed ../deeptools/GRCh38.p12.genome.clean.g2l7.$bname.clean.peaks.consensus.bed ../deeptools/GRCh38.p12.genome.clean.g2l12.$bname.clean.peaks.consensus.bed -out ../deeptools/$bname.mat.gz -m 50 --startLabel 's' --endLabel 'e' -b 500 -a 500 -bs 1 --sortRegions no --skipZeros -p 'max' && \
  plotProfile -m ../deeptools/$bname.mat.gz -o ../deeptools/20190704_$bname.png --dpi 300 --plotHeight 9 --plotWidth 10 --plotType 'se' --numPlotsPerRow 2 --startLabel 's' --endLabel 'e' --regionsLabel 'G3L7' 'G3L12' 'G2L7' 'G2L12' --samplesLabel '' --plotTitle $bname --colors 'black' 'black' 'black' 'black' -y 'CPM' --legendLocation none --perGroup"
done
```
