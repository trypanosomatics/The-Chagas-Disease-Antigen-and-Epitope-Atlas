library(data.table)

# Input files -----------------------------------------------------------------------------

complete_seq_file = "Files/chimeras.csv"
full_dataset_file = "Files//01_complete_dataset.tsv"

dt_complete_seq = fread(complete_seq_file, header = T, sep = ",", na.strings = NULL)
variants_full <- fread(full_dataset_file, header = T, sep = "\t", na.strings = NULL) #secuencia y variante

# Data treatment --------------------------------------------------------------------

sec_and_pep <- data.table(sec = as.character(),
                          pep = as.character(),
                          var = as.character(),
                          start = as.numeric(),
                          signal = as.numeric(),
                          source = as.character()
)

for (i in 1:dt_complete_seq[, .N]) {
    #all seqs
    #i=2
    temp_sec = dt_complete_seq[i]$secuencias
    
    for (j in 1:(nchar(dt_complete_seq[i]$secuencias)-14)) {
        #nchar - 14 = 15mers
        #j=12
        temp_pep = substr(temp_sec, start=j, stop= j+14)
        temp_var = dt_complete_seq[i]$variantes
        
        temp_dat = data.table(sec = temp_sec,
                              pep = temp_pep,
                              var = temp_var,
                              start = j,
                              signal = 0,
                              source = "")
        # data_map = rbind(data_map,temp_dat)
        sec_and_pep = rbindlist(list(sec_and_pep, temp_dat))
    }
}

#add SOURCE
all_sources = sec_and_pep[rep(seq_len(nrow(sec_and_pep)), times = 6), ]
sources = unique(variants_full$source)
all_sources$source[1:length(sec_and_pep$sec)]="AR"
all_sources$source[(length(sec_and_pep$sec)+1):(2*length(sec_and_pep$sec))]="BO"
all_sources$source[(2*length(sec_and_pep$sec)+1):(3*length(sec_and_pep$sec))]="BR"
all_sources$source[(3*length(sec_and_pep$sec)+1):(4*length(sec_and_pep$sec))]="CO"
all_sources$source[(4*length(sec_and_pep$sec)+1):(5*length(sec_and_pep$sec))]="MX"
all_sources$source[(5*length(sec_and_pep$sec)+1):(6*length(sec_and_pep$sec))]="US"

#add SIGNAL
signal_peptide <- data.table(sec = as.character(),
                             pep = as.character(),
                             var = as.character(),
                             start = as.numeric(),
                             signal = as.numeric(),
                             source = as.character()
)

for (i in 1:length(sources)){
    #i=3
    temp_source = sources[i]
    temp_df = all_sources[all_sources$source==temp_source,]
    temp_var_full = variants_full[variants_full$source==temp_source&variants_full$type=="PO",]
    temp_df$signal = temp_var_full$normalized_signal[match(temp_df$pep, temp_var_full$sequence)]
    
    signal_peptide = rbindlist(list(signal_peptide, temp_df))
}

# Outputs ---------------------------------------------------------------------------

data_map_output_file = "Files/02_data_map.tsv"
write.table(signal_peptide, file = data_map_output_file, col.names = T, row.names = F, sep = "\t", quote = T)

