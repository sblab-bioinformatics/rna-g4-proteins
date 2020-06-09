## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Trim illumina adapters and quality trimming](#trim-illumina-adapters-and-quality-trimming)
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


