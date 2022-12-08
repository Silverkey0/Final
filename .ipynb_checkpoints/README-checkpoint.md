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
**Data fully loaded into vignette through DeSEQ2 steps.** I will complete the first draft of the processed data.
## Milestone 2
**Due Date**: Thursday November 29th
**An initial completion of vignette.** I will do further analysis through the processed data.load data into vignette (through STAR), for seeking feedback. Not all sections in the writing will be completed, but will be final project.
## Deliverable
**Due Date**: December 3rd
A complete repository with clear formatted vigneette and the data about the relationship between pancrease and drinking alcohol.
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
```

Replace_header.sh to change the column of each .txt file to its file name [Eun]

```{r}
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
setwd('/Users/henry/Documents/GitHub/Final')

# Read in the matrix
count_matrix <- read.delim("/Users/henry/Documents/GitHub/Final/DES_dataset/Merged_Gene_Counts.txt", header=T, sep="\t")

# Change column #1 into row name and delete the column #1 after that
row.names(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[-c(1)]

# Read in the sample sheet
sampletable <- read_tsv('/Users/henry/Documents/GitHub/Final/DES_dataset/sample.tsv')

# Change column #1 into row name
row.names(sampletable) <- sampletable$case_submitter_id

sampletable$alcohol_history <- as.factor(sampletable$alcohol_history)

sampletable$gender <- as.factor(sampletable$gender)
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
## [1] 41102
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
## DataFrame with 41102 rows and 6 columns
##                     baseMean log2FoldChange     lfcSE
##                     <numeric>      <numeric> <numeric>
## ENSG00000000003.15 1887.97781      0.0804583  0.142874
## ENSG00000000005.6     8.38228     -2.6045422  0.657433
## ENSG00000000419.13 1517.18200     -0.0159072  0.114282
## ENSG00000000457.14  708.13087     -0.0784189  0.112094
## ENSG00000000460.17  227.57440      0.1691936  0.123464
## ...                       ...            ...       ...
## ENSG00000288660.1     1.09663      0.9471235  0.606703
## ENSG00000288663.1    22.66617      0.0804520  0.246530
## ENSG00000288670.1   221.73127      0.0842019  0.137095
## ENSG00000288674.1     5.39217      0.2671591  0.284744
## ENSG00000288675.1    28.32291      0.3696539  0.193836
##                         stat      pvalue      padj
##                    <numeric>   <numeric> <numeric>
## ENSG00000000003.15  0.563141 5.73339e-01 0.8932726
## ENSG00000000005.6  -3.961684 7.44229e-05 0.0431611
## ENSG00000000419.13 -0.139192 8.89298e-01 0.9802883
## ENSG00000000457.14 -0.699584 4.84187e-01 0.8591191
## ENSG00000000460.17  1.370383 1.70567e-01 0.6892353
## ...                      ...         ...       ...
## ENSG00000288660.1   1.561100    0.118500        NA
## ENSG00000288663.1   0.326337    0.744169  0.948706
## ENSG00000288670.1   0.614186    0.539093  0.878788
## ENSG00000288674.1   0.938242    0.348120  0.799474
## ENSG00000288675.1   1.907043    0.056515  0.558374
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
## DataFrame with 41102 rows and 5 columns
##                      baseMean log2FoldChange      lfcSE
##                     <numeric>      <numeric>  <numeric>
## ENSG00000000003.15 1887.97781    2.94178e-04 0.00145817
## ENSG00000000005.6     8.38228   -7.79633e-06 0.00144270
## ENSG00000000419.13 1517.18200   -8.25095e-05 0.00144376
## ENSG00000000457.14  708.13087   -6.52219e-06 0.00144258
## ENSG00000000460.17  227.57440    1.35846e-05 0.00144263
## ...                       ...            ...        ...
## ENSG00000288660.1     1.09663    2.66383e-06 0.00144269
## ENSG00000288663.1    22.66617    2.52994e-06 0.00144267
## ENSG00000288670.1   221.73127   -4.76731e-05 0.00144301
## ENSG00000288674.1     5.39217    3.92397e-06 0.00144268
## ENSG00000288675.1    28.32291    9.94458e-06 0.00144267
##                         pvalue      padj
##                      <numeric> <numeric>
## ENSG00000000003.15 5.73339e-01 0.8932726
## ENSG00000000005.6  7.44229e-05 0.0431611
## ENSG00000000419.13 8.89298e-01 0.9802883
## ENSG00000000457.14 4.84187e-01 0.8591191
## ENSG00000000460.17 1.70567e-01 0.6892353
## ...                        ...       ...
## ENSG00000288660.1     0.118500        NA
## ENSG00000288663.1     0.744169  0.948706
## ENSG00000288670.1     0.539093  0.878788
## ENSG00000288674.1     0.348120  0.799474
## ENSG00000288675.1     0.056515  0.558374
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
## out of 41095 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 44, 0.11%
## LFC < 0 (down)     : 18, 0.044%
## outliers [1]       : 0, 0%
## low counts [2]     : 11161, 27%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```{r}
sum(result05$padj < 0.05, na.rm = TRUE)
```

```{r}
## [1] 62
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
##                    TCGA.IB.A5ST TCGA.IB.7889 TCGA.FB.AAQ3
## ENSG00000000003.15    10.600048    11.308229    11.501187
## ENSG00000000005.6      6.664012     3.594849     4.209485
## ENSG00000000419.13    10.394696    10.640558    10.983005
##                    TCGA.IB.AAUO TCGA.FB.AAPS TCGA.HZ.7919
## ENSG00000000003.15    10.097235    10.559839    10.520871
## ENSG00000000005.6      2.707905     6.345453     2.707905
## ENSG00000000419.13    10.491123    10.142043    10.627850
##                    TCGA.IB.7885 TCGA.3A.A9IC TCGA.IB.7651
## ENSG00000000003.15    10.703114    10.618378    10.596153
## ENSG00000000005.6      4.722643     2.707905     3.137335
## ENSG00000000419.13    10.342351    10.182673    10.224360
##                    TCGA.XD.AAUH TCGA.XD.AAUI TCGA.HZ.A77O
## ENSG00000000003.15    10.315429    10.893892    10.890884
## ENSG00000000005.6      4.611004     3.260566     3.318963
## ENSG00000000419.13    10.130414    10.488733    10.732432
##                    TCGA.US.A779 TCGA.LB.A9Q5 TCGA.IB.7887
## ENSG00000000003.15    10.735429    11.173607    10.893195
## ENSG00000000005.6      2.707905     3.479841     2.707905
## ENSG00000000419.13    10.275133    10.559516    10.696392
##                    TCGA.HZ.8317 TCGA.HZ.7922 TCGA.FB.AAQ2
## ENSG00000000003.15    11.257950    11.533112    11.163157
## ENSG00000000005.6      4.101546     3.246727     3.792674
## ENSG00000000419.13    10.768382    10.898406    10.863802
##                    TCGA.FB.A78T TCGA.IB.A5SS TCGA.HZ.A49H
## ENSG00000000003.15    10.285506    11.248745    10.975332
## ENSG00000000005.6      3.213246     3.236189     4.464828
## ENSG00000000419.13    10.079566    10.183576    10.625458
##                    TCGA.YY.A8LH TCGA.S4.A8RO TCGA.IB.7891
## ENSG00000000003.15    10.771852     9.967486    11.074258
## ENSG00000000005.6      2.707905     2.707905     3.451118
## ENSG00000000419.13    10.482645    11.067733    10.436659
##                    TCGA.3A.A9IB TCGA.FB.AAQ0 TCGA.FB.AAPY
## ENSG00000000003.15    11.239867    10.360883    10.510575
## ENSG00000000005.6      3.432783     2.707905     3.675724
## ENSG00000000419.13    10.324485    10.735677     9.896529
##                    TCGA.3A.A9IU TCGA.IB.7893 TCGA.HZ.7918
## ENSG00000000003.15    11.718699    11.132901    11.503508
## ENSG00000000005.6      3.245796     2.707905     2.707905
## ENSG00000000419.13    10.552200    10.713945    11.048385
##                    TCGA.FB.AAPQ TCGA.IB.A6UF TCGA.IB.7654
## ENSG00000000003.15    10.191084    11.202248    10.240826
## ENSG00000000005.6      2.707905     2.707905     5.125891
## ENSG00000000419.13    10.958960    11.592442    10.222368
##                    TCGA.F2.A8YN TCGA.IB.7886 TCGA.IB.7646
## ENSG00000000003.15    11.666558    11.613194    12.000182
## ENSG00000000005.6      2.707905     3.209143     2.707905
## ENSG00000000419.13    10.519929    10.977711    11.496885
##                    TCGA.HZ.A49I TCGA.3A.A9J0 TCGA.OE.A75W
## ENSG00000000003.15    10.201092     9.694395    10.572147
## ENSG00000000005.6      3.714737     6.839115     3.694336
## ENSG00000000419.13    10.409929    10.751108    10.475710
##                    TCGA.IB.A7LX TCGA.XD.AAUL TCGA.FB.AAPP
## ENSG00000000003.15    11.568193    10.532636    10.963881
## ENSG00000000005.6      2.707905     3.520731     3.717468
## ENSG00000000419.13    10.630210    10.488411    11.198800
##                    TCGA.FB.AAPZ TCGA.IB.AAUP TCGA.FB.AAQ1
## ENSG00000000003.15    11.062667    10.080222    11.012250
## ENSG00000000005.6      6.491947     3.704451     2.707905
## ENSG00000000419.13    10.824076    10.150868    10.467668
##                    TCGA.3E.AAAY TCGA.IB.AAUR TCGA.IB.A7M4
## ENSG00000000003.15    10.725913    10.118494    11.106114
## ENSG00000000005.6      4.234104     4.117073     2.707905
## ENSG00000000419.13    10.391711    10.115208    11.347074
##                    TCGA.S4.A8RM TCGA.HZ.7925 TCGA.IB.A6UG
## ENSG00000000003.15    10.299225    10.504500    11.455591
## ENSG00000000005.6      3.581969     3.502649     2.707905
## ENSG00000000419.13    10.321697    10.533819    10.228016
##                    TCGA.IB.AAUM TCGA.3A.A9I9 TCGA.IB.AAUU
## ENSG00000000003.15    10.934569    10.865270    10.216059
## ENSG00000000005.6      4.205384     4.530319     2.707905
## ENSG00000000419.13    10.155166    10.539375    10.420189
##                    TCGA.RB.AA9M TCGA.YB.A89D TCGA.IB.7647
## ENSG00000000003.15    11.300528    10.738449     9.929323
## ENSG00000000005.6      3.935167     8.512636     3.293500
## ENSG00000000419.13    10.262981    10.366276    10.545484
##                    TCGA.IB.AAUQ TCGA.US.A77G TCGA.HZ.A8P0
## ENSG00000000003.15    10.350067    10.154877    10.515422
## ENSG00000000005.6      3.902537     3.474657     2.707905
## ENSG00000000419.13    10.189443    10.276342     8.859356
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
df <- as.data.frame(colData(DES_dataset)[,c("case_submitter_id", "gender", "alcohol_history")])
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

![Heatmap(rld).png](https://github.com/Silverkey0/Final/blob/main/Vigentte/Heatmap(rld).png)

Sample-to-Sample distances

```{r}
sampleDists <- dist(t(assay(vsd)))
DistMatrix <- as.matrix(sampleDists)
pheatmap(DistMatrix)
```

![Heatmap(DistMatrix).png](https://github.com/Silverkey0/Final/blob/main/Vigentte/Heatmap(DistMatrix).png)

### Principal Component Analysis Plot ###

```{r}
plotPCA(vsd, intgroup="alcohol_history")
```

![principal_analysis.png](https://github.com/Silverkey0/Final/blob/main/Vigentte/principal_analysis.png.png)