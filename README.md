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
awk '{print $1}' TCGA-FB-AAQ0.tsv > gene_id.txt
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
paste gene_id.txt TCGA-FB-AAQ0.txt TCGA-FB-AAPY.txt TCGA-3A-A9IU.txt TCGA-IB-7893.txt TCGA-HZ-7918.txt TCGA-FB-AAPQ.txt TCGA-IB-A6UF.txt TCGA-IB-7654.txt TCGA-F2-A8YN.txt TCGA-IB-7886.txt TCGA-IB-7646.txt TCGA-HZ-A49I.txt TCGA-3A-A9J0.txt TCGA-OE-A75W.txt TCGA-IB-A7LX.txt TCGA-XD-AAUL.txt TCGA-FB-AAPP.txt TCGA-FB-AAPZ.txt TCGA-IB-AAUP.txt TCGA-FB-AAQ1.txt TCGA-3E-AAAY.txt TCGA-IB-AAUR.txt TCGA-IB-A7M4.txt TCGA-S4-A8RM.txt TCGA-HZ-7925.txt TCGA-IB-A6UG.txt TCGA-IB-AAUM.txt TCGA-3A-A9I9.txt TCGA-IB-AAUU.txt TCGA-RB-AA9M.txt TCGA-YB-A89D.txt TCGA-IB-7647.txt TCGA-IB-AAUQ.txt TCGA-US-A77G.txt TCGA-HZ-A8P0.txt > Merged_Gene_Counts.txt
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
row.names(sampletable) <- sampletable$case_submitter_id

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
##                     baseMean log2FoldChange     lfcSE
##                    <numeric>      <numeric> <numeric>
## ENSG00000000003.15  1847.940      0.3432972  0.223171
## ENSG00000000005.6     17.548     -4.4119953  0.916820
## ENSG00000000419.13  1520.797     -0.0576332  0.191393
## ENSG00000000457.14   726.296     -0.0335386  0.172066
## ENSG00000000460.17   225.542      0.2176643  0.185136
## ...                      ...            ...       ...
## ENSG00000288660.1    1.34019     0.41225385  0.761056
## ENSG00000288663.1   23.81928     0.00502119  0.354229
## ENSG00000288670.1  212.39178    -0.08557983  0.205988
## ENSG00000288674.1    5.82581     0.19588790  0.406215
## ENSG00000288675.1   28.09440     0.46080399  0.302039
##                         stat      pvalue       padj
##                    <numeric>   <numeric>  <numeric>
## ENSG00000000003.15  1.538273 1.23982e-01 0.81971045
## ENSG00000000005.6  -4.812280 1.49218e-06 0.00403182
## ENSG00000000419.13 -0.301124 7.63320e-01 0.98175320
## ENSG00000000457.14 -0.194917 8.45458e-01 0.98772154
## ENSG00000000460.17  1.175701 2.39714e-01 0.89130663
...                      ...         ...        ...
## ENSG00000288660.1   0.541687    0.588034   0.963093
## ENSG00000288663.1   0.014175    0.988690   0.999113
## ENSG00000288670.1  -0.415460    0.677805   0.976402
## ENSG00000288674.1   0.482227    0.629645   0.971329
## ENSG00000288675.1   1.525642    0.127099   0.820672

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
##                     baseMean log2FoldChange      lfcSE
##                    <numeric>      <numeric>  <numeric>
## ENSG00000000003.15  1847.940   -9.66068e-06 0.00144268
## ENSG00000000005.6     17.548   -5.45537e-06 0.00144270
## ENSG00000000419.13  1520.797    3.21209e-05 0.00144283
## ENSG00000000457.14   726.296   -2.78548e-06 0.00144264
## ENSG00000000460.17   225.542   -4.81844e-07 0.00144265
## ...                      ...            ...        ...
## ENSG00000288660.1    1.34019    7.19334e-07 0.00144269
## ENSG00000288663.1   23.81928    3.16623e-08 0.00144268
## ENSG00000288670.1  212.39178   -7.45090e-06 0.00144267
## ENSG00000288674.1    5.82581    1.22911e-06 0.00144269
## ENSG00000288675.1   28.09440    2.50409e-06 0.00144268
##                         pvalue       padj
##                      <numeric>  <numeric>
## ENSG00000000003.15 1.23982e-01 0.81971045
## ENSG00000000005.6  1.49218e-06 0.00403182
## ENSG00000000419.13 7.63320e-01 0.98175320
## ENSG00000000457.14 8.45458e-01 0.98772154
## ENSG00000000460.17 2.39714e-01 0.89130663
## ...                        ...        ...
## ENSG00000288660.1     0.588034   0.963093
## ENSG00000288663.1     0.988690   0.999113
## ENSG00000288670.1     0.677805   0.976402
## ENSG00000288674.1     0.629645   0.971329
## ENSG00000288675.1     0.127099   0.820672
```

### Exploring and Exporting Results ###
#### MA-plot ####
```{r}
# Normal data
plotMA(result_table, ylim=c(-2,2))
```
![result_table](https://github.com/Silverkey0/Final/blob/main/Vigentte/result_table.png)

```{r}
# Shrink data
plotMA(resultLFC, ylim=c(-2,2))
```

![resultLFC](https://github.com/Silverkey0/Final/blob/main/Vigentte/resultLFC.png)

#### Set p-value and Print the Summary ####
```{r}
result05 <- results(DES_dataset, alpha = 0.05)
summary(result05)
```

```{r}
## out of 38283 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 23, 0.06%
## LFC < 0 (down)     : 11, 0.029%
## outliers [1]       : 0, 0%
## low counts [2]     : 13374, 35%
## (mean count < 4)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```{r}
sum(result05$padj < 0.05, na.rm = TRUE)
```

```{r}
## [1] 34
```

#### Use “plotCounts" to Make a Plot for the Read Counts of Single Gene Across the Groups ####
```{r}
plotCounts(DES_dataset, gene = which.min(result05$padj), intgroup = "alcohol_history")
```

![gene.png](https://github.com/Silverkey0/Final/blob/main/Vigentte/gene.png)

#### Extracting transformed values ####
```{r}
vsd <- vst(DES_dataset, blind = FALSE)
rld <- rlog(DES_dataset, blind = FALSE)
```

```{r}
## rlog() may take a few minutes with 30 or more samples,
## vst() is a much faster transformation
```

```{r}
head(assay(vsd), 3)
```

```{r}
## TCGA.FB.AAQ0 TCGA.FB.AAPY TCGA.3A.A9IU
## ENSG00000000003.15    10.319431     10.46900    11.681672
## ENSG00000000005.6      2.488669      3.51606     3.061314
## ENSG00000000419.13    10.694605      9.85416    10.514349
##                    TCGA.IB.7893 TCGA.HZ.7918 TCGA.FB.AAPQ
## ENSG00000000003.15    11.117784    11.473718    10.152505
## ENSG00000000005.6      2.488669     2.488669     2.488669
## ENSG00000000419.13    10.698456    11.018299    10.921165
##                    TCGA.IB.A6UF TCGA.IB.7654 TCGA.F2.A8YN
## ENSG00000000003.15    11.171759    10.208675    11.635816
## ENSG00000000005.6      2.488669     5.038297     2.488669
## ENSG00000000419.13    11.562186    10.190192    10.488326
                   TCGA.IB.7886 TCGA.IB.7646 TCGA.HZ.A49I
## ENSG00000000003.15    11.583447    11.977760    10.165152
## ENSG00000000005.6      3.023692     2.488669     3.559424
## ENSG00000000419.13    10.947554    11.474217    10.374246
##                    TCGA.3A.A9J0 TCGA.OE.A75W TCGA.IB.A7LX
## ENSG00000000003.15     9.663161    10.540009    11.528823
## ENSG00000000005.6      6.791195     3.539016     2.488669
## ENSG00000000419.13    10.721324    10.443468    10.590175
##                    TCGA.XD.AAUL TCGA.FB.AAPP TCGA.FB.AAPZ
## ENSG00000000003.15    10.504364    10.952590    11.027196
## ENSG00000000005.6      3.356005     3.570885     6.433774
## ENSG00000000419.13    10.460089    11.187699    10.788415
##                    TCGA.IB.AAUP TCGA.FB.AAQ1 TCGA.3E.AAAY
## ENSG00000000003.15    10.044796    10.983121    10.689385
## ENSG00000000005.6      3.548781     2.488669     4.106236
## ENSG00000000419.13    10.115541    10.438022    10.354839
##                    TCGA.IB.AAUR TCGA.IB.A7M4 TCGA.S4.A8RM
## ENSG00000000003.15    10.079768    11.089879    10.251064
## ENSG00000000005.6      3.982598     2.488669     3.414931
## ENSG00000000419.13    10.076477    11.331012    10.273562
##                    TCGA.HZ.7925 TCGA.IB.A6UG TCGA.IB.AAUM
## ENSG00000000003.15     10.47531    11.421067    10.895244
## ENSG00000000005.6       3.33652     2.488669     4.074597
## ENSG00000000419.13     10.50466    10.192412    10.115033
##                    TCGA.3A.A9I9 TCGA.IB.AAUU TCGA.RB.AA9M
## ENSG00000000003.15    10.828242    10.182738    11.276271
## ENSG00000000005.6      4.415881     2.488669     3.796845
## ENSG00000000419.13    10.502043    10.387122    10.237731
##                    TCGA.YB.A89D TCGA.IB.7647 TCGA.IB.AAUQ
## ENSG00000000003.15    10.712182     9.900150    10.316825
## ENSG00000000005.6      8.481342     3.114056     3.758927
## ENSG00000000419.13    10.339602    10.517144    10.155995
##                    TCGA.US.A77G TCGA.HZ.A8P0
## ENSG00000000003.15    10.119572    10.474424
## ENSG00000000005.6      3.305171     2.488669
## ENSG00000000419.13    10.241197     8.815162 
```

#### Effects of Transformations on the Variance ####

```{r}
ntd <- normTransform(DES_dataset)
meanSdPlot(assay(ntd))
```

![ntd.png](https://github.com/Silverkey0/Final/blob/main/Vigentte/ntd.png)

```{r}
meanSdPlot(assay(vsd))
```

![vsd.png](https://github.com/Silverkey0/Final/blob/main/Vigentte/vsd.png)

```{r}
meanSdPlot(assay(rld))
```

![rld.png](https://github.com/Silverkey0/Final/blob/main/Vigentte/rld.png)

#### Data Quality Assessment by Sample Clustering and Visualization ####

```{r}
select <- order(rowMeans(counts(DES_dataset,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(DES_dataset)[,c("case_submitter_id", "alcohol_history")])
df$case_submitter_id <- as.factor(df$case_submitter_id)
df <- df[-c(1)]
```

#### Heatmap of the Count Matrix ####
```{r}
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
```

![Heatmap(ntd).png](https://github.com/Silverkey0/Final/blob/main/Vigentte/Heatmap(ntd).png)

```{r}
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
```

![Heatmap(vsd).png](https://github.com/Silverkey0/Final/blob/main/Vigentte/Heatmap(vsd).png)

```{r}
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE,annotation_col=df)
```

![Heatmap(rld).png](https://github.com/Silverkey0/Final/blob/main/Vigentte/rld.png)

Sample-to-Sample distances
```{r}
sampleDists <- dist(t(assay(vsd)))
DistMatrix <- as.matrix(sampleDists)
pheatmap(DistMatrix)
```

![Heatmap(DistMatrix).png](https://github.com/Silverkey0/Final/blob/main/Vigentte/Heatmap(DistMatrix).png)