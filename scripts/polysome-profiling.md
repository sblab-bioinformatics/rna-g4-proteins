
## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Trim illumina adapters and quality trimming](#trim-illumina-adapters-and-quality-trimming)
- [Alignment](#alignment)
- [Count reads in genes](#count-reads-in-genes)
- [Downstream analysis](#downstream-analysis)


## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- Standard Unix tools: sbatch
- [tophat v2.1.1](https://ccb.jhu.edu/software/tophat/index.shtml)
- [htseq v0.7.2](https://htseq.readthedocs.io/en/master/)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [edgeR v3.16.5](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
  - [ggfortify v0.4.10](https://cran.r-project.org/web/packages/ggfortify/index.html)
  - [ggrepel v0.6.5](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)
  - [biomaRt v2.30.0](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
  - [gridExtra v2.2.1](https://cran.r-project.org/web/packages/gridExtra/index.html)
  - [limma v3.30.11](http://bioconductor.org/packages/release/bioc/html/limma.html)


## Libraries

Polysomal fraction | File name
:------:|:---------:
AML | AML.fastq.gz
AM11 | AM11.fastq.gz
A2L | A2L.fastq.gz
A211 | A211.fastq.gz
A10L | A10L.fastq.gz
A1011 | A1011.fastq.gz
BML | BML.fastq.gz
BM11 | BM11.fastq.gz
B2L | B2L.fastq.gz
B211 | B211.fastq.gz
B10L | B10L.fastq.gz
B1011 | B1011.fastq.gz
CML | CML.fastq.gz
CM11 | CM11.fastq.gz
C2L | C2L.fastq.gz
C211 | C211.fastq.gz
C10L | C10L.fastq.gz
C1011 | C1011.fastq.gz


## Trim illumina adapters and quality trimming

```bash
cd fastq

mkdir ../trimmed

for f in *.fastq.gz
do 
  sbatch --mem 8G -J $f -o ../trimmed/$f.out -e ../trimmed/$f.err --wrap "cutadapt -f fastq -m 10 -e 0.1 -q 20 -O 4 -a AGATCGGAAGAGC -o ../trimmed/${f%%.fastq.gz}.trimmed.fq.gz $f"
done
```


## Alignment

Using `hg19` reference genome:

```bash
cd ../trimmed

mkdir ../tophat_out

gtf='../reference/genes.gtf'
fa='../reference/Bowtie2Index/genome'

for f in *.trimmed.fq.gz
do
  sbatch --mem 16G -J $f -o ../tophat_out/$f.out -e ../tophat_out/$f.err --wrap "tophat -o ../tophat_out/${f}_hg19 -p 8 --library-type fr-unstranded -G $gtf $fa $f"
done
```


## Count reads in genes

```bash
cd ../tophat_out

mkdir merged_counts

gtf='../reference/genes.gtf'

for bam in *.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o merged_counts/$bname.log --mem 16G --wrap "htseq-count -f bam -r name -s no -t exon -i gene_id -m intersection-strict $bam $gtf > merged_counts/$bname.htseq"
done
```


## Visualization

```r
library(data.table)
library(ggplot2)
library(edgeR)
library(ggfortify)
library(ggrepel)
library(biomaRt)
library(gridExtra)
library(limma)

setwd("~/")

# gets gene symbol, transcript_id and go_id for all genes annotated with ribosomal Go terms
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl") # uses human ensembl annotations
gene.data <- getBM(attributes = c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'), filters = 'go', values = c('GO:0022625','GO:0022627'), mart = ensembl)
go_ribosome <- unique(gene.data$hgnc_symbol)
remove(ensembl, gene.data)

# read in table
cnt <- fread('./tableCat.py -i tophat_out/merged_counts/*.htseq')
setnames(cnt, names(cnt), c('gene_id', 'count', 'library_id'))
raw.data <- data.frame(dcast.data.table(data = cnt, gene_id ~ library_id, value.var = 'count'))
names(raw.data) <- gsub(".genes.htseq", "", names(raw.data))
names(raw.data) <- sapply(strsplit(names(raw.data), "_"), `[`, 3)
rownames(raw.data) <- raw.data[,1]
raw.data[,1] <- NULL

# edgeR and limma modelling
y <- DGEList(counts = raw.data)
y <- calcNormFactors(y)
raw.data <- data.frame(cpm(y))
remove(y, cnt)
raw.data <- raw.data[rowSums(raw.data) > 1,]

ratios <- data.frame(U_1_ratio=raw.data$AM11/raw.data$AML,U_2_ratio=raw.data$BM11/raw.data$BML,
                    U_3_ratio=raw.data$CM11/raw.data$CML,P2_1_ratio=raw.data$A211/raw.data$A2L,
                    P2_2_ratio=raw.data$B211/raw.data$B2L,P2_3_ratio=raw.data$C211/raw.data$C2L,
                    P10_1_ratio=raw.data$A1011/raw.data$A10L,P10_2_ratio=raw.data$B1011/raw.data$B10L,
                    P10_3_ratio=raw.data$C1011/raw.data$C10L,row.names=row.names(raw.data))

ratios <- na.omit(ratios)
ratios <- ratios[is.finite(rowSums(ratios)),]

group <- factor(c("DMSO","DMSO","DMSO","PDS2","PDS2","PDS2","PDS10","PDS10","PDS10"))
experiment <- factor(c(1,2,3,1,2,3,1,2,3))
design <- model.matrix(~0+group+experiment)
design
fit <- lmFit(ratios, design, group = group)
plotMDS(fitted(fit))
contrast.matrix <- makeContrasts(groupPDS2-groupDMSO, groupPDS10-groupDMSO, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
PDS2 <- topTable(fit2, coef = 1, adjust = "none", number = 1e6, p.value = 1)
PDS10 <- topTable(fit2, coef = 2, adjust = "none", number = 1e6, p.value = 1)

PDS2$gene_id <- rownames(PDS2)
rownames(PDS2) <- NULL
PDS10$gene_id <- rownames(PDS10)
rownames(PDS10) <- NULL
labs <- c(paste("down", nrow(PDS2[PDS2$logFC < -0.1375 & PDS2$adj.P.Val < 0.1,]), sep = " = "),
          paste("up", nrow(PDS2[PDS2$logFC > 0.1375 & PDS2$adj.P.Val < 0.1,]), sep = " = "),
          paste("down", nrow(PDS10[PDS10$logFC < -0.1375 & PDS10$adj.P.Val < 0.1,]), sep = " = "),
          paste("up", nrow(PDS10[PDS10$logFC > 0.1375 & PDS10$adj.P.Val < 0.1,]), sep = " = "))

plots <- list(go_ribosome)
names <- c("Ribosome")


# plotting
ggplot(PDS2) + geom_point(aes(x = logFC, y= -log10(adj.P.Val)),size=0.5, color="grey80") + 
  geom_point(data=PDS2[PDS2$gene_id %in% plots[[1]],], aes(x = logFC, y= -log10(adj.P.Val)),size = 1, color="red") + 
  theme_bw(base_family = "Arial", base_size = 12 ) + 
  scale_x_continuous(limits = c(-2.5, 2.5)) + scale_y_continuous(limits = c(0, 2)) +
  geom_text_repel(data=PDS2[PDS2$gene_id %in% plots[[1]] & PDS2$adj.P.Val < 0.1 ,],
                   aes(logFC,-log10(adj.P.Val), label=gene_id) , box.padding = 0.5, size = 4,force=2, 
                   point.padding = 0.1, segment.color = 'grey20') +
  xlab("") + ylab("") + theme(axis.text.x=element_text(colour="black"),rect=element_rect(color = "grey99",size = 1),
                              legend.position = "none")

ggsave("PDS2_ribosome.png",width=4.13,height=3.1,device="png",dpi=600)


ggplot(PDS10) + geom_point(aes(x = logFC, y= -log10(adj.P.Val)),size=0.5, color="grey80") + 
  geom_point(data=PDS10[PDS10$gene_id %in% plots[[1]],], aes(x = logFC, y= -log10(adj.P.Val)),size = 1, color="red") + 
  theme_bw(base_family = "Arial") + 
  scale_x_continuous(limits = c(-2.5, 2.5)) + scale_y_continuous(limits = c(0, 2)) +
  geom_text_repel(data=PDS10[PDS10$gene_id %in% plots[[1]] & PDS10$adj.P.Val < 0.1,],
                  aes(logFC,-log10(adj.P.Val), label=gene_id) , box.padding = 0.5, size = 4,force=2, 
                  point.padding = 0.1, segment.color = 'grey20') +
  xlab("") + ylab("") + theme(axis.text.x=element_text(colour="black"),rect=element_rect(color = "black",size = 1))

ggsave("PDS10_ribosome.png",width=4.13,height=3.1,device="png",dpi=600)


# output tables
PDS2_list <- PDS2[PDS2$adj.P.Val < 0.1 & PDS2$logFC < -0.1375,]
PDS10_list <- PDS2[PDS10$adj.P.Val < 0.1 & PDS10$logFC < -0.1375,]
write.table('PDS2_polysome_down_genes.txt', x = PDS2_list$gene_id, row.names = F, col.names = F, quote = F)
write.table('all_genes.txt', x = PDS2$gene_id, row.names = F, col.names = F, quote = F)
write.table('PDS10_polysome_down_genes.txt', x = PDS10_list$gene_id, row.names = F, col.names = F, quote = F)
```
