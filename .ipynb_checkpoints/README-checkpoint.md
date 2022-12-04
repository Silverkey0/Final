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
## out of 38284 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 103, 0.27%
## LFC < 0 (down)     : 30, 0.078%
## outliers [1]       : 0, 0%
## low counts [2]     : 14115, 37%
## (mean count < 5)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```{r}
sum(result05$padj < 0.05, na.rm = TRUE)
```

```{r}
## [1] 133
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
##                    X0944ce65.a89a.4916.b90a.f674b334281e
## ENSG00000000003.15                              11.11784
## ENSG00000000005.6                                2.50032
## ENSG00000000419.13                              10.69853
##                    X953359da.e535.48e6.97dc.1b4ccad6a671
## ENSG00000000003.15                             10.469092
## ENSG00000000005.6                               3.523737
## ENSG00000000419.13                              9.854301
##                    X45e48080.0b27.45f4.8149.1ae7e36c6933
## ENSG00000000003.15                              11.42111
## ENSG00000000005.6                                2.50032
## ENSG00000000419.13                              10.19252
##                    X6a0a280d.4a30.4071.afe7.de6b82a6e34c
## ENSG00000000003.15                              11.08994
## ENSG00000000005.6                                2.50032
## ENSG00000000419.13                              11.33106
##                    ff8dac97.42ff.4f48.81c4.50315ad8de58
## ENSG00000000003.15                             10.98319
## ENSG00000000005.6                               2.50032
## ENSG00000000419.13                             10.43812
##                    X74aa6675.28b7.46dc.9015.f07a07c7b018
## ENSG00000000003.15                             10.119690
## ENSG00000000005.6                               3.313616
## ENSG00000000419.13                             10.241305
##                    f2313562.c894.486d.a9b6.ed362e12a32a
## ENSG00000000003.15                            11.027259
## ENSG00000000005.6                              6.435199
## ENSG00000000419.13                            10.788489
##                    a643d0b8.aa50.481c.9f12.f2e634d61d90
## ENSG00000000003.15                            10.689464
## ENSG00000000005.6                              4.111973
## ENSG00000000419.13                            10.354938
##                    f73dfd7f.e054.46ff.9485.e8d05afe2da6
## ENSG00000000003.15                            11.583489
## ENSG00000000005.6                              3.033212
## ENSG00000000419.13                            10.947620
##                    X9d5c5303.550f.454f.8fd4.8abcdd311215
## ENSG00000000003.15                             11.276324
## ENSG00000000005.6                               3.803556
## ENSG00000000419.13                             10.237839
##                    X1422bc35.4597.4e76.aa75.8e79613e67d5
## ENSG00000000003.15                             10.952656
## ENSG00000000005.6                               3.578367
## ENSG00000000419.13                             11.187756
##                   ecdd0e44.0add.4a08.a3f8.ab2f51df7afd
## ENSG00000000003.15                             10.15262
## ENSG00000000005.6                               2.50032
## ENSG00000000419.13                             10.92123
##                    X27286a38.1035.4a22.9ebc.874df8c49026
## ENSG00000000003.15                             10.165266
## ENSG00000000005.6                               3.566946
## ENSG00000000419.13                             10.374345
##                    X396ff766.a383.4f42.bce3.f2f32b7f1151
## ENSG00000000003.15                             10.828314
## ENSG00000000005.6                               4.420741
## ENSG00000000419.13                             10.502134
##                    afe89625.b355.454d.8c0b.b4161edd69f8
## ENSG00000000003.15                            10.316927
## ENSG00000000005.6                              3.765765
## ENSG00000000419.13                            10.156109
##                    X98f1d0eb.0977.4f53.a3b1.e6875a34c27b
## ENSG00000000003.15                              10.47452
## ENSG00000000005.6                                2.50032
## ENSG00000000419.13                               8.81545
##                    X5dbb4642.6636.4d67.9e73.25c7312aec44
## ENSG00000000003.15                              11.63586
## ENSG00000000005.6                                2.50032
## ENSG00000000419.13                              10.48842
##                    X56a1b8ff.a71c.4e5f.9c3c.d8162ff896aa
## ENSG00000000003.15                              11.17182
## ENSG00000000005.6                                2.50032
## ENSG00000000419.13                              11.56223
##                    a53c919a.4e08.46f1.af3f.30b16b597c33
## ENSG00000000003.15                             10.18285
## ENSG00000000005.6                               2.50032
## ENSG00000000419.13                             10.38722
##                    d1e0f479.77dd.42cd.89cb.c1eeac60c10a
## ENSG00000000003.15                             11.97779
## ENSG00000000005.6                               2.50032
## ENSG00000000419.13                             11.47426
##                    d14154a5.9c8a.4c55.8c2d.64f3cad7b7d5
## ENSG00000000003.15                             9.900287
## ENSG00000000005.6                              3.123225
## ENSG00000000419.13                            10.517233
##                    X620e0648.ec20.4a12.a6cb.5546fe829c77
## ENSG00000000003.15                              11.47376
## ENSG00000000005.6                                2.50032
## ENSG00000000419.13                              11.01836
##                    cfda26b9.d417.425c.9a72.fa76ca4b296c
## ENSG00000000003.15                            10.895312
## ENSG00000000005.6                              4.080428
## ENSG00000000419.13                            10.115151
##                    X81b6b844.d760.4146.9b01.1c4b2596eb77
## ENSG00000000003.15                              10.31953
## ENSG00000000005.6                                2.50032
## ENSG00000000419.13                              10.69468
##                    X65e5bf5f.f22a.4ec2.8839.aeda37212b98
## ENSG00000000003.15                             10.208785
## ENSG00000000005.6                               5.041704
## ENSG00000000419.13                             10.190304
                   cd9dc5c5.1161.43ea.a617.1356a3b23e16
## ENSG00000000003.15                            11.681712
## ENSG00000000005.6                              3.070687
## ENSG00000000419.13                            10.514438
##                    f91eb33d.76f0.418c.9962.11a5b31ca1dc
## ENSG00000000003.15                            10.475401
## ENSG00000000005.6                              3.344849
## ENSG00000000419.13                            10.504750
##                    b84b58c7.95b8.4162.8e61.414f8fe422c6
## ENSG00000000003.15                            10.504454
## ENSG00000000005.6                              3.364262
## ENSG00000000419.13                            10.460182
##                    X4a60e87e.6bc1.4ae0.9ccb.4cdc3e99c608
## ENSG00000000003.15                             10.540097
## ENSG00000000005.6                               3.546611
## ENSG00000000419.13                             10.443562
##                    c70621cc.c8e1.441d.93ec.bd6922cc7f61
## ENSG00000000003.15                            10.251171
## ENSG00000000005.6                              3.422972
## ENSG00000000419.13                            10.273668
##                   b84b58c7.95b8.4162.8e61.414f8fe422c6.1
## ENSG00000000003.15                              10.712260
## ENSG00000000005.6                                8.481704
## ENSG00000000419.13                              10.339703
                   a308a5ee.8ea1.47ee.823d.14736d74a925
## ENSG00000000003.15                             9.663322
## ENSG00000000005.6                              6.792323
## ENSG00000000419.13                            10.721401
##                    X6788e0e5.b67f.4fe6.aabb.6d0430ae2d06
## ENSG00000000003.15                             10.079888
## ENSG00000000005.6                               3.988713
## ENSG00000000419.13                             10.076598
##                    X0a3c8161.4186.4c95.8288.9928d6db0355
## ENSG00000000003.15                             10.044919
## ENSG00000000005.6                               3.556341
## ENSG00000000419.13                             10.115659
##                    X712c2f78.c736.42ce.b689.a954c5290987
## ENSG00000000003.15                              11.52887
## ENSG00000000005.6                                2.50032
## ENSG00000000419.13                              10.59026
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
df <- as.data.frame(colData(DES_dataset)[,c("case_id", "alcohol_history")])
row.names(df) <- df$case_id
df <- df[-c(1)]
```

#### Heatmap of the Count Matrix ####
```{r}
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
```

```{r}
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
```

```{r}
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE,annotation_col=df)
```

```{r}
sampleDists <- dist(t(assay(vsd)))
DistMatrix <- as.matrix(sampleDists)
pheatmap(DistMatrix)
```