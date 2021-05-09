
## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Trim illumina adapters and quality trimming](#trim-illumina-adapters-and-quality-trimming)
- [Alignment](#alignment)
- [Deeptools](#deeptools)
- [Bigwigs](#bigwigs)


## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- Standard Unix tools: cat, awk, sort, uniq, paste, grep, pigz, cut, sbatch, nohup, sed
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [sambamba v0.6.5](https://lomereiter.github.io/sambamba/)
- [deeptools v3.3.0](https://deeptools.readthedocs.io/en/develop/)


## Libraries

Sample | Replicate | File name | GEO GSE | GEO GSM
:--------:|:---------:|:---------:|:-------:|:-------:
DHX36_WT | 1 | DHX36_WT_rep1.fastq.gz | - | -
DHX36_WT | 2 | DHX36_WT_rep2.fastq.gz | - | -
DHX36_WT | 3 | DHX36_WT_rep3.fastq.gz | - | -
DHX36_DAIH | 1 | DHX36_DAIH_rep1.fastq.gz | - | -
DHX36_DAIH | 2 | DHX36_DAIH_rep2.fastq.gz | - | -
DHX36_DAIH | 3 | DHX36_DAIH_rep3.fastq.gz | - | -
GRSF1 | 1 | GRSF1_rep1.fastq.gz | - | -
GRSF1 | 2 | GRSF1_rep2.fastq.gz | - | -
GRSF1 | 3 | GRSF1_rep3.fastq.gz | - | -


## Trim illumina adapters and quality trimming

```bash
cd ~/fastq

mkdir ../fastq_trimmed

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4G --wrap "cutadapt -a AGATCGGAAGAGC -m 20 --max-n=20 -q 20 -o ../fastq_trimmed/${bname}.fq.gz $fq > ../fastq_trimmed/$bname.txt"
done
```


## Alignment

### Align and sort

```bash
cd ~/fastq
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

Removing unmmaped reads and filter by alignment quality too but KEEPING DUPLICATES:

```bash
cd ~/bam

for bam in `ls *.bam | grep -v "tmp"`
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.clean.log --mem 16G --wrap "samtools idxstats $bam | cut -f1 | grep '^chr' | xargs samtools view -@ 20 -b -F 516 -q 1 $bam > $bname.clean.bam && samtools index $bname.clean.bam"
done
```


## Deeptools

Comparing RNAseq to iCLIP:

```bash
cd ~/bam

mkdir -p ../deeptools/{multiBamSummary,plotCorrelation}

nohup multiBamSummary bins -b pc5_rep*.clean.bam DHX36_WT_rep*.clean.bam DHX36_DAIH_rep*.clean.bam GRSF1_rep*.clean.bam \
-out ../deeptools/multiBamSummary/bins.npz \
-l pcDNA5_1 pcDNA5_2 pcDNA5_3 DHX36_WT_1 DHX36_WT_2 DHX36_WT_3 DHX36_DAIH_1 DHX36_DAIH_2 DHX36_DAIH_3 GRSF1_1 GRSF1_2 GRSF1_3 \
-p "max" &

cd ../deeptools/multiBamSummary

plotCorrelation -in bins.npz -o ../plotCorrelation/20180918_bins_pearson_heatmap.png -c pearson -p heatmap \
-l "pcDNA5_1" "pcDNA5_2" "pcDNA5_3" "DHX36_WT_1" "DHX36_WT_2" "DHX36_WT_3" "DHX36_DAIH_1" "DHX36_DAIH_2" \
"DHX36_DAIH_3" "GRSF1_1" "GRSF1_2" "GRSF1_3" \
--removeOutliers --colorMap Reds --plotNumbers &

sbatch --mem 32G -o DHX36_WT.multiBamSummary.out -e DHX36_WT.multBamSummary.err -J DHX36_WT --wrap "multiBamSummary BED-file \
--BED ~/reference/gencode.v28.annotation.sorted.transcripts.bed \
-b ~/bam/DHX36_WT*.clean.bam ~/bam_iclip/DHX36_WT_rep*.clean.bam \
-out DHX36_WT_transcripts.npz \
-l DHX36_WT_CLIP_rep1 DHX36_WT_CLIP_rep2 DHX36_WT_CLIP_rep3 DHX36_WT_RNAseq_rep1 DHX36_WT_RNASeq_rep2 DHX36_WT_RNASeq_rep3 -p \"max\""

plotCorrelation -in DHX36_WT_transcripts.npz -o ../plotCorrelation/DHX36_WT_transcripts.scatter.png -c pearson -p scatterplot \
-l "CLIP_rep1" "CLIP_rep2" "CLIP_rep3" "RNAseq_rep1" "RNASeq_rep2" "RNASeq_rep3" --removeOutliers --log1p &

sbatch --mem 32G -o DHX36_DAIH.multiBamSummary.out -e DHX36_DAIH.multBamSummary.err -J DHX36_DAIH --wrap "multiBamSummary BED-file \
--BED ~/reference/gencode.v28.annotation.sorted.transcripts.bed \
-b ~/bam/DHX36_EA*.clean.bam ~/bam_iclip/DHX36_DAIH_rep*.clean.bam \
-out DHX36_DAIH_transcripts.npz \
-l DHX36_DAIH_CLIP_rep1 DHX36_DAIH_CLIP_rep2 DHX36_DAIH_CLIP_rep3 DHX36_DAIH_RNAseq_rep1 DHX36_DAIH_RNASeq_rep2 DHX36_DAIH_RNASeq_rep3 \
-p \"max\""

plotCorrelation -in DHX36_DAIH_transcripts.npz -o ../plotCorrelation/DHX36_DAIH_transcripts.scatter.png -c pearson -p scatterplot \
-l "CLIP_rep1" "CLIP_rep2" "CLIP_rep3" "RNAseq_rep1" "RNASeq_rep2" "RNASeq_rep3" --removeOutliers --log1p &

sbatch --mem 32G -o GRSF1_WT.multiBamSummary.out -e GRSF1_WT.multBamSummary.err -J GRSF1_WT --wrap "multiBamSummary BED-file \
--BED ~/reference/gencode.v28.annotation.sorted.transcripts.bed \
-b ~/bam/GRSF1_WT*.clean.bam ~/bam_iclip/GRSF1_rep*.clean.bam \
-out GRSF1_WT_transcripts.npz \
-l GRSF1_WT_CLIP_rep1 GRSF1_WT_CLIP_rep2 GRSF1_WT_CLIP_rep3 GRSF1_WT_CLIP_rep4 GRSF1_WT_RNAseq_rep1 GRSF1_WT_RNASeq_rep2 \
GRSF1_WT_RNASeq_rep3 -p \"max\""

plotCorrelation -in GRSF1_WT_transcripts.npz -o ../plotCorrelation/GRSF1_WT_transcripts.scatter.png -c pearson -p scatterplot \
-l "CLIP_rep1" "CLIP_rep2" "CLIP_rep3" "CLIP_rep4" "RNAseq_rep1" "RNASeq_rep2" "RNASeq_rep3" --removeOutliers --log1p &
```

## Bigwigs

```bash
cd ~/bam

for bam in *.clean.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../bw/$bname.log --mem 32G --wrap "bamCoverage -b $bam -o ../bw/$bname.bw -of bigwig --binSize 1 -p 20 --normalizeUsing CPM"
done
```
