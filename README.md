# Final Project Outline
## Title
Differential Gene Expression in TCGA-PAAD within Stage II pancrease ductal and lobular neoplasms comparing drinkers and non-drinkers using DeSEQ2.
## Author
Henry Xu
## Overview of Project
I will identify differentially expressed genes for pancrease ductal and lobular neoplasms between drinkers and non-drinkers. This analysis will utilize the package DeSEQ2 and follow the specific vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html. For this analysis, I'll use the TCGA cohort and have identified 35 STAR-Counts files for tumors that fit within my cohort with 25 drinkers and 10 non-drinkers. Within the analysis, I will control for race (white), gender (male), ethnicity (not hispanic or latino).

Vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
## Data
I will use the data from https://portal.gdc.cancer.gov/repository. Examining clinical data, there are 35 tumor samples, and 25 are defined by me as drinkers, and 10 are identified as non-drinkers. The specific files are available are here 
https://github.com/Silverkey0/Final/tree/main/pancrease_sample
## Milestone 1
**Due Date**: Thursday November 22th
**Data fully loaded into vignette through DeSEQ2 steps.** I will complete an entire first draft of analysis analyzed through the vignette.
## Milestone 2
**Due Date**: Thursday November 29th
**An initial completion of vignette.** I will do further analysis through the vignette.Data loaded into vignette (through STAR), for seeking feedback. Not all sections in the writing will be completed, but will be final project.
## Deliverable
**Due Date**: December 3rd
A complete repository with clear documentation and description of your analysis and results.
## Method
**Regarding the method, I would like to give special thanks to Eun Sung. He helped me to solve many problems. I use some of his codes in my method, and I'll cite them.**

### Change File Names ###
I changed all the filename of gene_counts.txt to match with the case_id in the sample.txt. So, after I process all of them and merge them into a matrix, the first row of matrix matches the first column of sample.txt.

### Process gene_counts.tsv and Merge them into Matrix ###
Extract the column #4 which is unstranded counts of the sample

```{r}
awk '{print $4}' /Users/henry/Documents/GitHub/Final/Gene_Counts/feec3523-951e-4541-8664-2a71adfdecb5.rna_seq.augmented_star_gene_counts.tsv > /Users/henry/Documents/GitHub/Final/Gene_Counts_Merged/feec.txt
```

Remove the first line of the file without printing [Eun]

```{r}
tail -n +2 feec.txt > feec.tmp && mv feec.tmp feec.txt

Replace_header.sh to change the column of each .txt file to its file name [Eun]

./replace_header.sh
```

Extract the gene_id, the column #1, from one of .tsv file

```{r}
awk '{print $1}' feec.tsv > gene_id.txt
```

Remove the first line of the file without printing [Eun]

```{r}
tail -n +2 gene_id.txt > gene_id.tmp && mv gene_id.tmp gene_id.txt
```

Remove row #2 to #4 [Eun]

```{r}
awk '!/^N_*/' gene_id.txt > gene_id.tmp && mv gene_id.tmp gene_id.txt
```

Merged counts

```{r}
paste gene_id.txt 00bb.txt 0c7c.txt 1fbb.txt 002d.txt 47b1.txt 5d3a.txt 7d5d.txt 7f20.txt 7fc4.txt 9ed8.txt 13db.txt 033f.txt 65ac.txt 67be.txt 78b8.txt 90fe.txt 92aa.txt 875e.txt 1536.txt 3330.txt 6633.txt 6813.txt 6855.txt 9077.txt a5f5.txt a66b.txt afa6.txt c9b3.txt d5d8.txt d92e.txt e1f5.txt ea66.txt f98d.txt f658.txt feec.txt > Merged_Gene_Counts.txt
``` 

### Check Packages in Your Rstudio library ###
```{r}
library(tibble)
library(tidyverse)
library(apeglm)
library(ggplot2)
library(vsn)
library(pheatmap)
library(ReportingTools)
library(DESeq2)
```

I use the commands below to install packages I don't have.

```{r}
install.packages("BiocManager")
```

```{r}
library(BiocManager)
BiocManager::install("Package")
```

### Setting Working Directory and Wrangling the Raw Count Matrix and Sample Sheet ### [Eun]
```{r}
count_matrix <- read.delim("/Users/henry/Documents/GitHub/Final/Gene_Counts_Merged/Merged_Gene_Counts.txt", header=T, sep="\t")

# Change column #1 into row name and delete the column #1 after that
row.names(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[-c(1)]

# Read in the sample sheet
sampletable <- read_tsv('/Users/henry/Documents/GitHub/Final/pancrease_sample/sample.tsv')

# Change column #1 (sample_id) into row name
row.names(sampletable) <- sampletable$case_id

# Read in the sample sheet
sampletable$alcohol_history <- as.factor(sampletable$alcohol_history)
```

### Create DESeq2 object ###
```{r}
DES_dataset <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = sampletable,
                                         design = ~ alcohol_history)
```

### Filtering ###
```{r}
# Number of gene before filtering
nrow(DES_dataset)
```

```{r}
## [1] 60660
```

```{r}
# Filtering to keep only rows that have at least 10 reads total
DES_dataset <- DES_dataset[rowSums(counts(DES_dataset)) > 10, ]

# Number of gene after filtering
nrow(DES_dataset)
```

```{r}
## [1] 38297
```

### Performing standard differential expression analysis ###
```{r}
DES_dataset <- DESeq(DES_dataset)
```

```{r}
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
## -- replacing outliers and refitting for 1503 genes
## -- DESeq argument 'minReplicatesForReplace' = 7 
## -- original counts are preserved in counts(dds)
## estimating dispersions
## fitting model and testing
```

### Get result table ### [Eun]
```{r}
DES2Report <- HTMLReport(shortName = 'RNAseq_Analysis_with_DEseq2', title = 'Differential Expression Analysis in Pancrease', reportDirectory = "./reports")
publish(DES_dataset,DES2Report, pvalueCutoff=0.05, annotation.db="org.Mm.eg.db", factor = colData(DES_dataset)$alcohol_history, reportDir="./reports")
finish(DES2Report)
```

```{r}
## [1] "./reports/RNAseq_Analysis_with_DEseq2.html"
```

### Geberate Result table ###
```{r}
result_table <- results(DES_dataset)
result_table
```

```{r}
## log2 fold change (MLE): alcohol history Yes vs No 
## Wald test p-value: alcohol history Yes vs No 
## DataFrame with 38297 rows and 6 columns
##                      baseMean log2FoldChange     lfcSE
##                     <numeric>      <numeric> <numeric>
## ENSG00000000003.15 1847.93968      0.0101748  0.230122
## ENSG00000000005.6     2.74651     -0.6736128  0.785834
## ENSG00000000419.13 1520.79722     -0.2372484  0.187057
## ENSG00000000457.14  726.29566      0.0229574  0.172018
## ENSG00000000460.17  225.54164     -0.1487979  0.186415
## ...                       ...            ...       ...
## ENSG00000288660.1     1.34019       0.507832  0.760994
## ENSG00000288663.1    23.81928       0.236015  0.352715
## ENSG00000288670.1   212.39178      -0.110868  0.205451
## ENSG00000288674.1     5.82581       0.754388  0.405270
## ENSG00000288675.1    28.09440       0.504508  0.299792
##                         stat    pvalue      padj
##                    <numeric> <numeric> <numeric>
## ENSG00000000003.15  0.044215  0.964733  0.994355
## ENSG00000000005.6  -0.857194  0.391337        NA
## ENSG00000000419.13 -1.268319  0.204684  0.719505
## ENSG00000000457.14  0.133459  0.893830  0.978175
## ENSG00000000460.17 -0.798209  0.424749  0.851692
## ...                      ...       ...       ...
## ENSG00000288660.1   0.667328 0.5045629        NA
## ENSG00000288663.1   0.669137 0.5034080  0.881044
## ENSG00000288670.1  -0.539634 0.5894498  0.906777
## ENSG00000288674.1   1.861445 0.0626813  0.536226
## ENSG00000288675.1   1.682858 0.0924026  0.595163
```

### Use “lfcShrink” to Shrink the Effect Size with apeglm Method ### [Eun]
```{r}
resultLFC <- lfcShrink(DES_dataset, coef = "alcohol_history_Yes_vs_No", type = "apeglm")
```

```{r}
## using 'apeglm' for LFC shrinkage. If used in published research, please cite:
##     Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
##     sequence count data: removing the noise and preserving large differences.
##     Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
```

```{r}
resultLFC
```

```{r}
## log2 fold change (MAP): alcohol history Yes vs No 
## Wald test p-value: alcohol history Yes vs No 
## DataFrame with 38297 rows and 5 columns
##                      baseMean log2FoldChange      lfcSE
                    <numeric>      <numeric>  <numeric>
## ENSG00000000003.15 1847.93968   -3.75131e-04 0.00146853
## ENSG00000000005.6     2.74651   -4.73781e-06 0.00144270
## ENSG00000000419.13 1520.79722   -7.27350e-06 0.00144266
## ENSG00000000457.14  726.29566    2.47539e-06 0.00144265
## ENSG00000000460.17  225.54164   -2.34388e-05 0.00144275
## ...                       ...            ...        ...
## ENSG00000288660.1     1.34019    9.00561e-07 0.00144269
## ENSG00000288663.1    23.81928    1.91597e-06 0.00144268
## ENSG00000288670.1   212.39178    4.29716e-06 0.00144266
## ENSG00000288674.1     5.82581    4.37509e-06 0.00144269
## ENSG00000288675.1    28.09440    6.77654e-06 0.00144269
##                       pvalue      padj
##                    <numeric> <numeric>
## ENSG00000000003.15  0.964733  0.994355
## ENSG00000000005.6   0.391337        NA
## ENSG00000000419.13  0.204684  0.719505
## ENSG00000000457.14  0.893830  0.978175
## ENSG00000000460.17  0.424749  0.851692
## ...                      ...       ...
## ENSG00000288660.1  0.5045629        NA
## ENSG00000288663.1  0.5034080  0.881044
## ENSG00000288670.1  0.5894498  0.906777
## ENSG00000288674.1  0.0626813  0.536226
## ENSG00000288675.1  0.0924026  0.595163
## ```

### Exploring and Exporting Results ###
#### MA-plot ####
```{r}
plotMA(result_table, ylim=c(-2,2))
```
![Result_table](https://github.com/Silverkey0/Final/blob/main/Vigentte/result_table.png)