library(dplyr)
library(reshape2)
setwd("/home/leonel/Documents/paperChagastope/SUP_CODE/The-Chagas-Disease-Antigen-and-Epitope-Atlas/example/")
design_alanine <- read.csv("inputs/31_Design_AlanineScan_example.tsv",sep = "\t")

input_data_path <- "./inputs/32_AlanineScan_raw_data/"
files <- list.files(input_data_path)
files <- files[grepl("raw.tsv",files)]

for (i in 1:length(files)) {
  # i<-1
  Chips_roche_norm_dat <- read.csv(paste0(input_data_path,files[i],sep = ""),sep = "\t",stringsAsFactors = F)
  Chips_roche_norm_dat$sample <- paste(strsplit(files[i],"_")[[1]][c(1:2)],collapse = "_")
  if (i == 1) {
    Chips_roche_norm <- Chips_roche_norm_dat
    next()
  }
  Chips_roche_norm <- rbind(Chips_roche_norm,Chips_roche_norm_dat)
}




estimateSignalChanges = function(selected_protein) {
  # selected_protein <- "TcCLB.511671.60"
  design_alanine_filtered <- filter(design_alanine, protein == selected_protein)
  AlanineScanData <- merge(Chips_roche_norm,design_alanine_filtered,by.x = "Reporter.Name",by.y = "array_express_id")
  AlanineScanData$signal <- (AlanineScanData$allData.replica1 + AlanineScanData$allData.replica2) / 2
  AlanineScanData$allData.replica1 <- NULL
  AlanineScanData$allData.replica2 <- NULL
  # head(AlanineScanData)
  orders <- group_by(filter(AlanineScanData,alanine_position == 0),sample) %>% summarise(signal = mean(signal))
  samples <- orders[order(orders$signal,decreasing = T),]$sample
  
  
  AlanineScanData$alanine_position_start <- AlanineScanData$peptide_to_scan_start + AlanineScanData$alanine_position - 1
  AlanineScanData$alanine_position_start[AlanineScanData$alanine_position == 0] <- 0
  
  temp_protein <- filter(AlanineScanData,!is.na(signal))
  if (length(temp_protein$Reporter.Name) < 2) {next()}
  temp_protein$change <- 0
  temp_originals <- select(filter(temp_protein,alanine_position_start == 0),peptide_to_scan,signal,sample) 
  temp_originals$original_signal <- temp_originals$signal
  temp_originals$signal <- NULL
  changeDT_total <- data.frame(change = 0,sample = "",stringsAsFactors = F,mutated_aa = "",protein = "")
  changeDT_total <- changeDT_total[0,]
  for (s in samples) {
    changeDT = filter(temp_protein,sample == s)
    temp_originals_filtered <- filter(temp_originals, sample == s)
    changeDT <- merge(changeDT,temp_originals_filtered,by = c("peptide_to_scan","sample"))
    changeDT$change <- changeDT$signal - changeDT$original_signal
    
    mutated_aa <- changeDT$alanine_position_start
    mutated_aa[changeDT$alanine_position == 0] <- 0
    mutated_aa <- unique(mutated_aa)

    changeDT$mutated_aa <- substr(changeDT$peptide_to_scan,changeDT$alanine_position,changeDT$alanine_position)
    changeDT_temp = dplyr::group_by(filter(changeDT,alanine_position != 0),alanine_position_start,sample) %>% dplyr::summarise(change = mean(change),
                                                                                                                               mutated_aa = unique(mutated_aa),
                                                                                                                               protein = unique(protein))
    changeDT_total <- rbind(changeDT_total,as.data.frame(changeDT_temp))
  }
  return(changeDT_total)
}

generateHeatmapMatrix = function(dt){
  temp_matrix <- acast(dt, paste(alanine_position_start,mutated_aa)~sample, value.var = "change",fun.aggregate = mean,na.rm = T)
  temp_matrix <- temp_matrix[order(as.numeric(gsub("[^0-9.]", "",rownames(temp_matrix)))),,drop = F]
  temp_matrix <- t(temp_matrix)
  temp_matrix
}

dt <- estimateSignalChanges("TcCLB.511671.60")
matrix <- generateHeatmapMatrix(dt)

write.table(dt,"./outputs/31_Alanine_scan_raw_data/Raw_data_signal_change_alanine_scan_TcCLB.511671.60.tsv",sep = "\t")
write.table(matrix,"./outputs/31_Alanine_scan_raw_data/Raw_data_signal_change_matrix_alanine_scan_TcCLB.511671.60.tsv",sep = "\t")
