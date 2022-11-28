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

Extract the column #4 which is unstranded counts of the sample

awk '{print $4}' /Users/henry/Documents/GitHub/Final/Gene_Counts/feec3523-951e-4541-8664-2a71adfdecb5.rna_seq.augmented_star_gene_counts.tsv > /Users/henry/Documents/GitHub/Final/Gene_Counts_Merged/feec.txt

Remove the first line of the file without printing [Eun]

tail -n +2 feec.txt > feec.tmp && mv feec.tmp feec.txt

Replace_header.sh to change the column of each .txt file to its file name [Eun]

./replace_header.sh

Extract the gene_id, the column #1, from one of .tsv file

awk '{print $1}' feec.tsv > gene_id.txt

Remove the first line of the file without printing [Eun]

tail -n +2 gene_id.txt > gene_id.tmp && mv gene_id.tmp gene_id.txt

Remove row #2 to #4 [Eun]

awk '!/^N_*/' gene_id.txt > gene_id.tmp && mv gene_id.tmp gene_id.txt

Merged counts

paste gene_id.txt 00bb.txt 0c7c.txt 1fbb.txt 002d.txt 47b1.txt 5d3a.txt 7d5d.txt 7f20.txt 7fc4.txt 9ed8.txt 13db.txt 033f.txt 65ac.txt 67be.txt 78b8.txt 90fe.txt 92aa.txt 875e.txt 1536.txt 3330.txt 6633.txt 6813.txt 6855.txt 9077.txt a5f5.txt a66b.txt afa6.txt c9b3.txt d5d8.txt d92e.txt e1f5.txt ea66.txt f98d.txt f658.txt feec.txt > Merged_Gene_Counts.txt