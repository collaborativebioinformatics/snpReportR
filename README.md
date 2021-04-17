# expressed-variant-reporting

snpReportR <img src='man/figures/logo.png' align="right" height="244" />

# Contributors
Jenny Smith

Brandon Michael Blobner

Ahmad Al Khleifat


# Goal
Develop a tool to facilitate genetic reporting, which takes as input a vcf file and gene name or id (optional) and produces a report file that informs the user of all genetic variants overlapping the gene, including any non-coding regulatory elements affecting expression of the gene.

The tool generates two reports, one is aimed for clinical use and the second aimed for researchers, informing the interpretation of genetic variants pertaining to the gene provided by the user.

# Introduction
Producing sequencing data, whether it is whole genome sequencing (WGS), whole exome sequencing (WES) or targeted gene panels, is common practice for the study of genetic bases of biological processes. In biomedical research, NGS data are widely used to investigate the genetic causes of disease, allowing for the study of genomic variants including single nucleotide variations, small insertions or deletions of a few bases, as well as structural variants.
There are several practical challenges when processing NGS data. For example, 40x WGS data for one sample produced on the Illumina Hiseq 2000, one of the most popular sequencers, is about 400 gigabytes in its raw format (fastq format). Such big files are not easy to handle for the average non-specialised scientist or lab, since they require sophisticated tools, bioinformatics skills and high-performance computing for their analysis. Furthermore, with the increasing availability of next-generation sequencing data, non-specialists, including health care professionals and patients, are obtaining their genomic information without a corresponding ability to analyse and interpret it as the relevance of novel or existing variants in genes of interest is not always apparent. In this context, research and medical workers, with different scientific backgrounds and levels of bioinformatics skills, deal with NGS data constantly. Here we describe SNP-ReprteR, an extremely fast, accurate and computationally light bioinformatics pipeline for the analysis, annotation and visualisation of DNA NGS data. 

# Installation 

```
devtools::install_github("collaborativebioinformatics/snpReportR")
```


# Methods

# Implementation
The SNP-ReportR tool was developed to be implemented in a clinical setting to inform the interpretation of genetic variants. SNP-ReportR software is an open access, gene centric data browser for genetic analysis. SNP-ReportR is a package developed using R. After entering the gene name (HGNC, Ensembl gene (ENSG), or transcript SNP-ReportR will produce a comprehensive genetic report.  SNP-ReportR output includes a summary of expressed variants found in the gene, RIN score, allowing both clinicians and researchers to assess the accuracy of the report. Furthermore, SNP-ReportR provides visualisation of associated haplotype, gene location, and the reported chromosomal abnormalities.  
  
SNP-ReportR is available on GitHub(https://github.com/collaborativebioinformatics/expressed-variant-reporting). The repository provides detailed instructions for tool usage and installation. 

## Inputs
https://docs.google.com/spreadsheets/d/1pcB_bI_83B__sJ_Qw3tYDUhAYTz7Bh9SBvxjMzd8L4U/edit?usp=sharing

## Outputs
![](https://github.com/collaborativebioinformatics/expressed-variant-reporting/blob/main/Report_example.png)



# Operation
SNP-ReportR requires a VCF file that includes all the variants that should be genotyped (Flow chart).

SNP-ReportR is a web page application and an overview of the method development is demonstrated in Figure 1. SNPReprter creates two reports. The first report is aimed at patients, non-specialist clinicians. The second report is aimed for genetic researchers.

After uploading the VCF file or entering the gene name (HGNC, Ensembl gene (ENSG), or transcript (ENST) identifier) in the search box on the homepage, you will be directed to the gene-specific page containing: 1) Gene-level summary which includes , 2) Links to the gene's page on OMIM, GTEx, gnomAD, 3) A dynamic table with the annotated variants overlapping the gene, 4) A graph with the distribution of the allele frequency for variants matched with OMIM, GTEx, gnomAD. The profile of the genetic variants to consider, such as type and size range, can be specified on the side bar. Each column in the dynamic table can be "searched" into or reordered dynamically. All data used by the app will be available for download in tab-delimited files. By default, allele frequency is reported based on dbVar and gnomAD genomes and exomes. Furthermore, SNP-ReportR utilises multiple database and links variants to genes and annotate gene impact, allele frequency, and the overlap with clinically-relevant SNVs, SVs and indels. All data, including are available for download in a tab-delimited file. Each variant has been extensively annotated and aggregated in a customizable table using OMIM/OMIA Variant, dbSNP - Variant Allele origin and allele frequency known or predicted RNA-editing site, Repeat family from UCSC Genome Browser Repeatmasker Annotations Homopolymer adjacency Entropy around the variant Splice adjacency FATHMM pathogenicity prediction COSM.

## Flow Chart
![](https://github.com/collaborativebioinformatics/expressed-variant-reporting/blob/main/flowchart.v3.png)

Working flow chart link
https://docs.google.com/presentation/d/1GwtkY6gtXj9l-k23cuZHhgFGPFtRZ0RkTLCqfZa2p-o/edit?usp=sharing

SNPReprter creates two reports. The first report is aimed at patients, non-specialist clinicians. The second report is aimed for genetic researchers.

# Results

SNP-ReportR is a web page application. After uploading the VCF file or entering the gene name (HGNC, Ensembl gene (ENSG) in the search box on the homepage, you will be directed to the gene-specific page containing:
1.	Gene-level summary in addition to information about associated disease.  
2.	Links to the gene's page on OMIM, GTEx, gnomAD.
3.	A dynamic table with the annotated variants overlapping the gene.
4.	A graph with showing summary of all the variants within the gene.
5.	Chromosome location. 
6.	Table showing summary of gene expression analysis.
7.	Tissue specific gene expression summary table.
8.	Table with haplotype information, such as start, end and length. 
9.	Haplotype visualisation. 
10.	Link to 5 recent publications about the associated genes. 
11.	SNPReprter creates two reports. The first report is aimed at patients, non-specialist clinicians. The second report is aimed for genetic researchers




First report Patients, non-specialist clinicians
The report was designed to use a patient friendly language 

Example clinical/patient report
https://github.com/collaborativebioinformatics/expressed-variant-reporting/blob/main/Reports/genetic-report-3.pdf

### References:

* Smith RN, et al. InterMine: a flexible data warehouse system for the integration and analysis of heterogeneous biological data. Bioinformatics. 2012 Dec 1;28(23):3163-5.

* Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med. 2015;17(5):405-424. doi:10.1038/gim.2015.30

* [ACMG 59 Genes](https://www.coriell.org/1/NIGMS/Collections/ACMG-59-Genes?gclid=CjwKCAiA_9r_BRBZEiwAHZ_v1zXfYPL5eKXE2acAJmrzdjMurIU7y1GleMAJeoIkAjNCSzXHw20sDRoCaNsQAvD_BwE)



