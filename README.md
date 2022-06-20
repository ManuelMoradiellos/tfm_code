# Recopilation of code used for the UAM Master's Degree in Bioinformatics and Computational Biology's Dissertation

## _Systematic understanding of higher-order interactions between inherited variants and somatic mutations_

*Author*: Manuel Moradiellos Corpus

### Abstract

One the most well-known mechanisms behind cancer development is Knudson’s Two Hit Hypothesis which states that inactivations
of both alleles of tumor suppressor genes or only one for an oncogene are enough to initiate a tumorigenic event. 

Given recent evidence about how the genomic mutational background may alter this classic model into scenarios in which
the cancer evolutionary path is different and that commonly cancer genomic research is highly biased towards the study
of somatic alterations, more works needs to be done to understand the interactions between mutations being contingent
upon the context and how it may lead to genes switching from one-hit drivers to two hit-drivers. Also, studying these
higher-order genetic interactions may lead to the findings of new candidate genes to be further studied in their role
as possible cancer driver genes.

This study shows a comprehensive methodology on how to identify significant co-occurrences between two different genomic
alterations within the same gene and a third interaction in a different gene focusing on rare pathogenic germline
data of $\sim$10.000 patients collected and curated from the TCGA project, and studied using model-based statistical
methods in a pan-cancer cohort. 

Although no higher-order genetic interactions were found to be significant, this framework poses a large-scale computational
analysis method that can be further improved to study this rare combination of events, highlighting  possible new interactions
between alternative tumor development paths and new not-cancer-described genes to be further investigated upon.

### Description 
The code presented here was created in order to tackle the following objectives:

- Collect all possible relevant genomic alterations in TCGA including rare damaging germline variants, LOH and somatic mutations across more than 30 cancer types.
    
- Systematically address two-hit models between germline variants and LOH within the same gene (_cis_-interactions) to understand tissue-specific cancer fitness landscape and elucidate the possibility of identifying novel cancer predisposition genes using two-hit models.
    
- Evaluating the evidence of high-order genetic interactions between two genes (_trans_-interactions) using three genomic alterations (_geneA_: germline variants, LOH and _geneB_: somatic mutations).


This repo is divided in:

- LOH_Calling_Data/ : Code used to obtain binary matrix for loss of heterozygosity calls on TCGA patients obtained from [Nichols _et al (2020)](https://doi.org/10.1038/s41467-020-16399-y)



_Only the code used is show, not data_
