
## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Trim illumina adapters and quality trimming](#trim-illumina-adapters-and-quality-trimming)


## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [sambamba v0.6.5](https://academic.oup.com/bioinformatics/article/31/12/2032/214758)
- Standard Unix tools: awk, sort, uniq
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
