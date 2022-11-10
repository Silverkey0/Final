# Final Project Outline
## Title
Differential Gene Expression in TCGA-BRCA within Stage II breast ductal and lobular neoplasms comparing middle-aged and elderly using DeSEQ2.
## Author
Henry Xu
## Overview of Project
I will identify differentially expressed genes for breast ductal and lobular neoplasms between middle-aged and elderly. This analysis will utilize the package DeSEQ2 and follow the specific vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html. For this analysis, I'll use the TCGA cohort and have identified 284 STAR-Counts files for tumors that fit within my cohort with 77 middle-aged (20-49) and 207 elderly (50-90). Within the analysis, I will control for race (white), gender (female), ethnicity (not hispanic or latino).

Vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
## Data
I will use the data from https://portal.gdc.cancer.gov/repository. Examining clinical data, there are 284 tumor samples, and 77 are defined by me as middle-aged, and 207 are identified as elderly. The specific files are available are here https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22set_id%3AtaMsYIQBh-9x-zpXkd47%22%5D%7D%2C%22op%22%3A%22IN%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22STAR%20-%20Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D&searchTableTab=cases.
## Milestone 1
**Due Date**: Thursday November 22th
**Data fully loaded into vignette through htseq steps.** I will complete an entire first draft of analysis analyzed through the vignette.
## Milestone 2
**Due Date**: Thursday November 29th
**An initial completion of vignette.** I will complete an entire first draft of analysis analyzed through the vignette.Data loaded into vignette (through htseq), for seeking feedback. Not all sections in the writing will be completed, but will be final project.
## Deliverable
**Due Date**: December 3rd
A complete repository with clear documentation and description of your analysis and results.