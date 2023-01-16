library(data.table)

# Input files -----------------------------------------------------------------------------

tssas_file_complete = "Files/TSSA_variants_len15_allPools_noReplicas.tsv"
TSSAs_data_map_complete= fread(tssas_file_complete, header = T, sep = "\t", na.strings = NULL)

# Data treatment --------------------------------------------------------------------

#PO = positive
TSSAs_data_map_complete = TSSAs_data_map_complete[type=="PO",]
TSSAs_data_map_complete$normalized_signal = as.numeric(TSSAs_data_map_complete$normalized_signal)
TSSAs_data_map_complete = na.omit(TSSAs_data_map_complete, cols=seq_along(TSSAs_data_map_complete), invert=FALSE)

#Add sequence
TSSAs_data_map_complete$sec = TSSAs_data_map_complete$protein
TSSAs_data_map_complete[sec == "TSSA_I", sec :=  "MTTCRLLCALLALALCCCLSACTTANGGSTSSTPPGTDKKTAAGGTPSPSGASSGEAEASSNKNDGSLSSSAWVFAPLALAASVLAYTALG"]
TSSAs_data_map_complete[sec == "TSSA_II", sec := "MTTCRLLCALLALALCCCLTACTTANGGSTSSTPPSGTENKPATGEAPSQPGASSGEAEASSNKNDGSLSSSAWVSAPLALAASALAYTALG"]
TSSAs_data_map_complete[sec == "TSSA_III", sec := "MATCRLLCALLALALCCCLSACTTANGGSTISTPPSGTEKKAAAGEAPSPSGASSGEAEASSNKNDGSLSSSAWVFAPLALTASALAYTALG"]
TSSAs_data_map_complete[sec == "TSSA_IV", sec :=  "MTTCRLLCALLALALCCCLSACTTANGGSTISTPPSGTDKKTAAGEAPSPSGASSGEAEASSNKNDGSLSSSAWVFAPLALAASALAYTALG"]


#AA Signal
TSSAS_data_aa <- data.table(AA = as.character(),
                            Seq = as.character(),
                            Signal_mean = as.numeric(),
                            Pos = as.numeric(),
                            source = as.character()
)

sequences = unique(TSSAs_data_map_complete$sec)
sources = unique(TSSAs_data_map_complete$source)

for(k in 1:length(sources)){
    #k=1
    temp_source=sources[k]
    for (i in 1:length(sequences)) {
        #i=1
        temp_sec=sequences[i]
        for (j in 1:nchar(temp_sec)) {
            #j=3
            temp_aa = substr(temp_sec,j,j)
            
            #mean(all peptides from secuence I and source K that contain residue J)
            temp_signal = mean(TSSAs_data_map_complete$normalized_signal[TSSAs_data_map_complete$sec==temp_sec&TSSAs_data_map_complete$start<=j&TSSAs_data_map_complete$start+14>=j&TSSAs_data_map_complete$source==temp_source],na.rm = T)
            
            temp_tssas_aa = data.table(AA = temp_aa,
                                       Seq = temp_sec,
                                       Signal_mean = temp_signal,
                                       Pos = j,
                                       source = temp_source
            )
            
            TSSAS_data_aa  = rbindlist(list(TSSAS_data_aa, temp_tssas_aa))
        }
    }
}

TSSAS_data_aa = na.omit(TSSAS_data_aa, cols=seq_along(TSSAS_data_aa), invert=FALSE)

#Add isoform
TSSAS_data_aa$protein = TSSAS_data_aa$Seq
TSSAS_data_aa[protein == "MTTCRLLCALLALALCCCLSACTTANGGSTSSTPPGTDKKTAAGGTPSPSGASSGEAEASSNKNDGSLSSSAWVFAPLALAASVLAYTALG", protein := "TSSAI"]
TSSAS_data_aa[protein == "MTTCRLLCALLALALCCCLTACTTANGGSTSSTPPSGTENKPATGEAPSQPGASSGEAEASSNKNDGSLSSSAWVSAPLALAASALAYTALG", protein := "TSSAII"]
TSSAS_data_aa[protein == "MATCRLLCALLALALCCCLSACTTANGGSTISTPPSGTEKKAAAGEAPSPSGASSGEAEASSNKNDGSLSSSAWVFAPLALTASALAYTALG", protein := "TSSAIII"]
TSSAS_data_aa[protein == "MTTCRLLCALLALALCCCLSACTTANGGSTISTPPSGTDKKTAAGEAPSPSGASSGEAEASSNKNDGSLSSSAWVFAPLALAASALAYTALG", protein := "TSSAIV"]


#Fix TSSAI and TSSAIV positions
for (i  in 1:TSSAS_data_aa[, .N]){
    #i=2
    if (TSSAS_data_aa$Seq[i] == "LSCTTANGGSTSSTPPGTDKKTAAGGTPSPS"){
        TSSAS_data_aa$Pos[i] = TSSAS_data_aa$Pos[i]+1
    }
    if (TSSAS_data_aa$Seq[i] == "ANGGSTISTPPSGTDKKTAAGEAPSPSGASSG"){
        TSSAS_data_aa$Pos[i] = TSSAS_data_aa$Pos[i]+6
    }
}

# Outputs ---------------------------------------------------------------------------

tssas_aa_signal_output_file = "Files/04_TSSA_aa_signal"
write.table(TSSAS_data_aa, file = tssas_aa_signal_output_file, col.names = T, row.names = F, sep = "\t", quote = T)
