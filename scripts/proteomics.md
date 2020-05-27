
## Table of contents

- [Software requirements](#software-requirements)
- [Data](#data)
- [Analysis](#analysis)


## Software requirements

- [R v3.5.2](https://www.r-project.org/). Libraries:
  - [tidyverse v1.3.0](https://www.tidyverse.org/)
  - [pander v0.6.3](https://cran.r-project.org/web/packages/pander/index.html)
  - [qPLEXanalyzer v1.0.4](https://www.bioconductor.org/packages/release/bioc/html/qPLEXanalyzer.html)
  - [readxl v1.3.1](https://cran.r-project.org/web/packages/readxl/index.html)
  - [gridExtra v2.2.1](https://cran.r-project.org/web/packages/gridExtra/index.html)
  - [biomaRt v2.30.0](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [ggrepel v0.6.5](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)


## Data

See (DMSO_PDS_Peptideintensities.xlsx)[../data/DMSO_PDS_Peptideintensities.xlsx]


## Analysis

Script created by CRUK bioinformatics core:

```r
library(tidyverse)
library(pander)
library(qPLEXanalyzer)
library(readxl)
library(gridExtra)
library(biomaRt)
library(ggplot2)
library(ggrepel)

# format data set and create sample sheet
set1 <-  read_excel( '../data/DMSO_PDS_Peptideintensities.xlsx', sheet = 1) 

names(set1)[c(1:6)] <- c( 'Confidence', 'Sequence', 'Modifications', 'Protein_Groups', 'Proteins', 'Accessions' )

set1 <- set1 %>% 
  dplyr::select( -Confidence, -Protein_Groups, -Proteins) 

sample_names <- names(set1)[4:ncol(set1)] %>% 
  str_remove_all(.,'-') %>% 
  str_replace(.,'.Rep.', '_')

names(set1)[4:ncol(set1)] <- sample_names


SampleGroup <- str_remove(sample_names, '_[1-9].*')

BioRep <- str_remove( sample_names, '[A-Z0-9]*_')

s_sheet <- data.frame( SampleName=sample_names, SampleGroup=SampleGroup, BioRep=BioRep, TechRep=NA, Run=1, stringsAsFactors = F) %>% 
  arrange( SampleGroup)

write.csv( x=s_sheet, file='../data/samplesheet_v1.0.csv', row.names = F)


set1 <- set1 %>% 
  dplyr::select( Sequence, Modifications, Accessions, one_of(s_sheet$SampleName)) %>% 
  na.omit()


write.csv( x=set1, file='../data/DMSO_PDS_formatted_v1.0.csv', row.names = F)

## Reading in data and annotation
s_sheet <- read.csv( file='../data/samplesheet_v1.0.csv', stringsAsFactors = F) 

anno_tab <- read.delim( file='../data/protein_gene.human.txt', header = T, stringsAsFactors = F) %>% 
  dplyr::rename( Accessions=Protein, GeneSymbol=Gene_Symbol)

int_tab <- read.csv( file='../data/DMSO_PDS_formatted_v1.0.csv', stringsAsFactors = F) %>% 
  dplyr::select( -Modifications)

## Removing zeroes and normalising using qPLEXanalyzer

int_tab[,c(3:ncol(int_tab))] <- int_tab[,c(3:ncol(int_tab))] + 1

MSnset_int <- convertToMSnset( ExpObj = int_tab, metadata =s_sheet, Sequences =1, 
                               Accessions =2, indExpData = c(3:ncol(int_tab)), 
                               rmMissing=TRUE )

## Plotting mean intensities
p1 <- intensityPlot(MSnset_int, title = "Raw peptide intensity distribution")
MSnset_norm_int <- qPLEXanalyzer::normalizeScaling(MSnset_int,  median)

p2 <- intensityPlot(MSnset_norm_int, title = "Normalised Peptide intensity distribution")

grid.arrange(p1, p2,ncol = 2, nrow = 1)

## Plotting relative peptide intensity distribution box-plots
p4 <- rliPlot(MSnset_int, title = "Raw : Relative Peptide intensity")
p5 <- rliPlot(MSnset_norm_int, title = "Normalised : Relative Peptide intensity")
grid.arrange(p4, p5, ncol = 2, nrow = 1)

## Calculating summarised protein intensities as sum of all peptide intensities
MSnset_Summarised <- summarizeIntensities(MSnset_norm_int, sum, anno_tab)
p3 <- intensityPlot(MSnset_Summarised, title = "Summarised intensity distribution")
grid.arrange(p1, p2, p3, ncol = 3, nrow = 1)


## Hierarchical clustering based on peptide intensities before normalisation and after normalisation
p1 <- hierarchicalPlot(MSnset_int)
p2 <- hierarchicalPlot(MSnset_norm_int)
p3 <- hierarchicalPlot(MSnset_Summarised)
grid.arrange(p1, p2, p3, ncol = 1, nrow = 3)

## Differential analysis

out_path <- '../results'
de_groups <- unique(s_sheet$SampleGroup)
treat_group <- de_groups[2]
ref_group <- de_groups[1]
cot_name <- paste( treat_group, 'Vs', ref_group, sep='_')
contrasts <- paste(treat_group, '-', ref_group, sep=' ')
names(contrasts) <- cot_name

de_out_file <- paste(  treat_group,'Vs', ref_group, 'v1.0.csv', sep='_')

de_out_file_with_path <- paste( out_path, de_out_file, sep='/')

diffstats <- computeDiffStats(MSnset_Summarised, contrasts=contrasts, batchEffect=c( 'SampleGroup'))

diffexp <- getContrastResults(diffstats=diffstats, contrast=contrasts, writeFile= FALSE)
diffexp <- diffexp %>% 
  arrange( adj.P.Val )
write.csv( x=diffexp, file=de_out_file_with_path, row.names = F)

#gets gene symbol, transcript_id and go_id for all genes annotated with ribosomal Go terms
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = c('GO:0022625','GO:0022627'), mart = ensembl)
go_ribosome <- unique(gene.data$hgnc_symbol)
remove(ensembl,gene.data)

ggplot(diffexp) + geom_point(aes(x = log2FC, y= -log10(P.Value)),size=.5, color="grey") + 
  geom_point(data=diffexp[diffexp$GeneSymbol %in% go_ribosome,], aes(x = log2FC, y= -log10(P.Value), color = "red"),size = 1) +  
  theme_bw(base_family = "Arial", base_size = 12 ) + xlim(-1.5,1.5) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "grey29") + ylim(0,5) + 
  theme(axis.text.x=element_text(colour="black"),rect=element_rect(color = "grey99",size = 1),legend.position = "none") 
ggsave("../data/PDS_ribosome.png",width=4.13,height=3.1,device="png",dpi=600)

## Write rnk file for GSEA
path <- '../results'

files_list <- list.files( path = path)

for(each_file in files_list[grep('csv', files_list)]){
  cat( each_file, '\n')
  in_file <- paste( path, each_file, sep='/')
  out_file <- paste( str_replace(each_file, '\\.csv', '\\.rnk'), sep='/')
  
  tab <- read_csv( file =in_file ) %>% 
    dplyr::select(GeneSymbol, t ) %>% 
    filter( !duplicated(GeneSymbol)) %>% 
    arrange( desc(t))
  write.table( file=out_file, x=tab, col.names = F, row.names = F, sep='\t', quote = F)
}
```
