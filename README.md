# The-Chagas-Disease-Antigen-and-Epitope-Atlas
![LOGO](https://chagastope.org/images/home/chagastope-logo-letters-only-v3.png)
A Trypanosoma cruzi Antigen and Epitope Atlas: deep characterization of individual antibody specificities in human Chagas Disease patients across the Americas.

This repository contains code that demonstrates the processing and analysis of peptide microarray data for the Atlas. To analyze the entire dataset you will need the following:

* The files included in this repository
* The raw data files from CHAGASTOPE-v1 arrays (such as *AR_PO_raw.tsv*) downloaded from Array Express [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11651/). Place them in **chagastope_data/inputs/02_pools_raw_data**. 
* Raw data files from CHAGASTOPE-v2 arrays (such as *AR_E1_PO_raw.tsv*) downloaded from Array Express [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11655/). Place them in **chagastope_data/inputs/12_individual_serums_raw_data**.
* The mapping file for CHAGASTOPE-v1 arrays (*Supplementary Table S09 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv*), available in our paper. Place it in **chagastope_data/inputs/01_pools_array_design**.
* The mapping file for CHAGASTOPE-v2 arrays (*Supplementary Table S10 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv*), available in our paper. Place it in **chagastope_data/inputs/11_individual_serums_array_design**.

To run the R code you will need the following packages: **data.table**, **preprocessCore** and **zoo** for the antigenicity analysis, and **dplyr** and **reshape2** for the alanine scan analysis (and optional **pheatmap** and **colorblindr** for ploting heatmaps). All the outputs from running this code is already available; the Normalized Data can be found in the same Array Express links, and the Smoothed Data and Region Data can be found as Supplementary Tables and Files in our paper.

Prefixes in the files and folders in this repository are there simply to indicate grouping and order.

## Testing the code
By default this code will use a test dataset present in **test_data** and included in this repository.

The antigenicity analysis example for a subset of 20 proteins can be performed by running in UNIX:
```
$ Rscript 01_pools_normalize_data.R
$ Rscript 02_pools_smooth_data.R
$ Rscript 03_calculate_peaks.R
$ Rscript 04_calculate_regions.R
$ Rscript 11_individual_serums_normalize_data.R
$ Rscript 12_individual_serums_smooth_data.R
```

Alanine Scan example for **Ag2-antigen | TcCLB.511671.60** can be reproduced as in the example provided here by running in UNIX:
```
$ Rscript 21_alanine_scan_analisis.R
```
You can also run the code directly from Rstudio, in which case you want to set either the working directory (*setwd* function) or change the *main_folder* variable found in the CONFIG section in each of the codes to point towards the main folder of this repository (the folder containing all the codes).

## Running the code for the entire dataset
To run the entire dataset you will need to download all the data mentioned aboved and place it in the proper location. Afterwards you can either run each of the codes adding the following argument at the end:
```
-test F
```
for example
```
$ Rscript 01_pools_normalize_data.R -test F
```
or you can edit each of the codes and set:
```
testing <- FALSE
```
Also remember to change either the working directory or the *main_folder* variable if you are running this code directly from Rstudio.
