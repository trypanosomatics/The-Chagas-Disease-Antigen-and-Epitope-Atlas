# The-Chagas-Disease-Antigen-and-Epitope-Atlas
A Trypanosoma cruzi Antigen and Epitope Atlas: deep characterization of individual antibody specificities in human Chagas Disease patients across the Americas.

This repository contains code that demonstrates the processing and analysis of peptide microarray data for the Atlas. You will need the following:

* The files included in this repository
* Raw data files from CHAGASTOPE-v1 arrays (such as *AR_PO_raw.tsv*) downloaded from Array Express [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11651/)
* Raw data files from CHAGASTOPE-v2 arrays (such as *AR_E1_PO_raw.tsv*) downloaded from Array Express [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11655/)
* The mapping file for CHAGASTOPE-v1 arrays (*Supplementary Table S09 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv*), available in our paper
* The mapping file for CHAGASTOPE-v2 arrays (*Supplementary Table S10 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv*), available in our paper

Once you have these, you should follow the instructions in each of the R code files in ascending order. In R you will need the following packages: **data.table**, 
**preprocessCore** and **zoo**. 

Keep in mind that if you are just interested in the outputs and not in the process, you can find the corresponding Normalized Data in the same Array Express links, and the Smoothed Data and Region Data as Supplementary Tables and Files in our paper.

Prefixes in the files and folders in this repository are there simply to indicate grouping and order.
