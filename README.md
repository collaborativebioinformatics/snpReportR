# expressed-variant-reporting

![](https://github.com/collaborativebioinformatics/expressed-variant-reporting/blob/main/snpReportR.png)

# Contributors


# Goal
Develop a tool to facilitate genetic reporting, which takes as input a vcf file and gene name or id (optional) and produces a report file that informs the user of all genetic variants overlapping the gene, including any non-coding regulatory elements affecting expression of the gene.

The tool generates two reports, one is aimed for clinical use and the second aimed for researchers, informing the interpretation of genetic variants pertaining to the gene provided by the user.

# Introduction
Producing sequencing data, whether it is whole genome sequencing (WGS), whole exome sequencing (WES) or targeted gene panels, is common practice for the study of genetic bases of biological processes. In biomedical research, NGS data are widely used to investigate the genetic causes of disease, allowing for the study of genomic variants including single nucleotide variations, small insertions or deletions of a few bases, as well as structural variants.
There are several practical challenges when processing NGS data. For example, 40x WGS data for one sample produced on the Illumina Hiseq 2000, one of the most popular sequencers, is about 400 gigabytes in its raw format (fastq format). Such big files are not easy to handle for the average non-specialised scientist or lab, since they require sophisticated tools, bioinformatics skills and high-performance computing for their analysis. Furthermore, with the increasing availability of next-generation sequencing data, non-specialists, including health care professionals and patients, are obtaining their genomic information without a corresponding ability to analyse and interpret it as the relevance of novel or existing variants in genes of interest is not always apparent. In this context, research and medical workers, with different scientific backgrounds and levels of bioinformatics skills, deal with NGS data constantly. Here we describe SNP-ReprteR, an extremely fast, accurate and computationally light bioinformatics pipeline for the analysis, annotation and visualisation of DNA NGS data. 



# Methods

# Implementation
The SNP-ReportR tool was developed to be implemented in a clinical setting to inform the interpretation of genetic variants. SNP-ReportR software is an open access, gene centric data browser for genetic analysis. SNP-ReportR is a package developed using R. After entering the gene name (HGNC, Ensembl gene (ENSG), or transcript SNP-ReportR will produce a comprehensive genetic report.  SNP-ReportR output includes a summary of expressed variants found in the gene, RIN score, allowing both clinicians and researchers to assess the accuracy of the report. Furthermore, SNP-ReportR provides visualisation of associated haplotype, gene location, and the reported chromosomal abnormalities.  
  
SNP-ReportR is available on GitHub(https://github.com/collaborativebioinformatics/expressed-variant-reporting). The repository provides detailed instructions for tool usage and installation. 

Inputs
https://docs.google.com/spreadsheets/d/1pcB_bI_83B__sJ_Qw3tYDUhAYTz7Bh9SBvxjMzd8L4U/edit?usp=sharing



# Operation

Flow Chart
![](https://github.com/collaborativebioinformatics/expressed-variant-reporting/blob/main/flowchart.v2.png)

Working flow chart link
https://docs.google.com/presentation/d/1GwtkY6gtXj9l-k23cuZHhgFGPFtRZ0RkTLCqfZa2p-o/edit?usp=sharing

SNPReprter creates two reports. The first report is aimed at patients, non-specialist clinicians. The second report is aimed for genetic researchers.


