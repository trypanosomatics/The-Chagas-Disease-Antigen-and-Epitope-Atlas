library(data.table)
library(plyr)

# Input files -----------------------------------------------------------------------------

TSSAs_singles_replicas = "Files/00_TSSA_variants_len16_allSingles_noReplicas.tsv"
TSSAs_singles_replicas <- fread(TSSAs_singles_replicas, header = T, sep = "\t", na.strings = NULL)

# Data treatment --------------------------------------------------------------------

TSSAs_replicas_full_signal = TSSAs_singles_replicas[,c("type","truncated", "raw_signal", "raw_signal_sd", "in_original_pool"):=NULL]
TSSAs_replicas_full_signal <- ddply(TSSAs_singles_replicas, c("protein","source", "parent_source"), summarise, sum(normalized_signal))
TSSAs_replicas_full_signal =  as.data.table(TSSAs_replicas_full_signal)

TSSAs_singles_replicas$total_signal = numeric()
for (i in 1:TSSAs_singles_replicas[,.N]){
    #i=2
    temp_protein = TSSAs_singles_replicas[i,"protein"]
    temp_source = TSSAs_singles_replicas[i,"source"]
    TSSAs_singles_replicas$total_signal[i] = TSSAs_replicas_full_signal[protein==temp_protein&source==temp_source,4]
}

TSSAs_singles_replicas$total_signal = as.numeric(TSSAs_singles_replicas$total_signal)
TSSAs_singles_replicas$graph_name = character()
for (i in 1:TSSAs_singles_replicas[,.N]){
    #i=2
    TSSAs_singles_replicas$graph_name[i] = paste(round(TSSAs_singles_replicas[i,total_signal]),"_",TSSAs_singles_replicas[i,source], sep = "")
}

# Outputs ---------------------------------------------------------------------------

TSSAs_mapped_output_file = "Files/09_singles_TSSAs_replicas.tsv"
write.table(TSSAs_singles_replicas, file = TSSAs_mapped_output_file, col.names = T, row.names = F, sep = "\t", quote = T)
