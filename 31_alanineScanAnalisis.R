setwd("/home/leonel/Documents/paperChagastope/SUP_CODE/The-Chagas-Disease-Antigen-and-Epitope-Atlas/example/")
design_alanine <- read.csv("inputs/31_Design_AlanineScan_example.tsv",sep = "\t")

input_data_path <- "./inputs/32_AlanineScan_raw_data/"
files <- list.files(input_data_path)
files <- files[grepl("raw.tsv",files)]

for (i in 1:length(files)) {
  # i<-1
  Chips_roche_norm_dat <- read.csv(paste0(input_data_path,files[i],sep = ""),sep = "\t",stringsAsFactors = F)
  if (i == 1) {
    Chips_roche_norm <- Chips_roche_norm_dat
    next()
  }
  Chips_roche_norm <- rbind(Chips_roche_norm,Chips_roche_norm_dat)
}

AlanineScanData <- merge(Chips_roche_norm,design_alanine,by.x = "Reporter.Name",by.y = "array_express_id")
AlanineScanData$signal <- (AlanineScanData$allData.replica1 + AlanineScanData$allData.replica2) / 2
AlanineScanData$allData.replica1 <- NULL
AlanineScanData$allData.replica2 <- NULL


