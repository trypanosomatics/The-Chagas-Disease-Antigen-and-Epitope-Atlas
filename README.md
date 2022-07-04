# The Chagas Disease Antigen and Epitope Atlas
<img src="https://chagastope.org/images/home/chagastope-logo-letters-only-v3.png" width="550px" alt="Chagastope Logo">

This is the web and electronic complement to the paper *"A Trypanosoma cruzi Antigen and Epitope Atlas: deep characterization of individual antibody specificities in human Chagas Disease patients across the Americas."*

The code demonstrates the processing and analysis of peptide microarray data as done for the Atlas. We make the code available so you can reproduce the results and play with the data. 

**BUT Remember** all the outputs from running this pipeline are already available: Normalized Data is available at the EMBL-EBI Array Express database, and the Smoothed Data and Region Data are available as Supplementary Tables and Files in the paper.


## Dependencies

The code is written in the R programming language, and depends on the following packages: **data.table**, **preprocessCore** and **zoo** for the antigenicity analysis, and **dplyr** and **reshape2** for the alanine scan analysis (and optional **pheatmap** and **colorblindr** for ploting heatmaps). 

## Testing the code

By default this code will use the test datasets present in the **test_data** directory included in this repository.

Prefixes in the filenames and folder names in the repository are there simply to indicate grouping and order of execution.

The antigenicity analysis example for a subset of 20 proteins can be performed by running in UNIX into the main folder of this repository (the folder containing all the ```.R``` scripts):
```
$ Rscript 01_pools_normalize_data.R
$ Rscript 02_pools_smooth_data.R
$ Rscript 03_calculate_peaks.R
$ Rscript 04_calculate_regions.R
$ Rscript 11_individual_serums_normalize_data.R
$ Rscript 12_individual_serums_smooth_data.R
```

And the analysis of single-residue mutagenesis data (aka AlanineScan) in the test data for **Ag2-antigen | TcCLB.511671.60** can be reproduced as in the example provided here by running in UNIX:
```
$ Rscript 21_alanine_scan_analisis.R
```

If you'd like to run the code from within Rstudio, you may want to set either the working directory (*setwd* function) or change the *main_folder* variable found in the CONFIG section in each of the script files to point towards the main folder of this repository.

## Analyzing Peptide Microarray Data from the Atlas

To analyze the entire dataset you will need the following:

* The files included in this repository
* The raw data files from CHAGASTOPE-v1 arrays (such as *AR_PO_raw.tsv*) downloaded from Array Express [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11651/). Place them in **chagastope_data/inputs/02_pools_raw_data**. 
* Raw data files from CHAGASTOPE-v2 arrays (such as *AR_E1_PO_raw.tsv*) downloaded from Array Express [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11655/). Place them in **chagastope_data/inputs/12_individual_serums_raw_data**.
* The mapping file for CHAGASTOPE-v1 arrays (*Supplementary Table S09 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv*), available in our paper. Place it in **chagastope_data/inputs/01_pools_array_design**.
* The mapping file for CHAGASTOPE-v2 arrays (*Supplementary Table S10 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv*), available in our paper. Place it in **chagastope_data/inputs/11_individual_serums_array_design**.

Download all data, place it in the proper location and run each script, adding ``` -test F``` at the end, for example:
```
$ Rscript 01_pools_normalize_data.R -test F
```
Alternatively, edit each script and set:

```
#### CONFIG ####
testing <- FALSE
```

Also remember to change either the working directory or the *main_folder* variable if you are running this code directly from Rstudio.
