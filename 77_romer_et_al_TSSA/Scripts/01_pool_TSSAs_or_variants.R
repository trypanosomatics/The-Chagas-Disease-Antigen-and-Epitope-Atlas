library(data.table)

# Input files -----------------------------------------------------------------------------

full_dataset_file <- "Files/TSSA_variants_all_pools_averagedReplicas.tsv"
variants_full <- fread(full_dataset_file, header = T, sep = "\t", na.strings = NULL) #secuencia y variante

# Data treatment --------------------------------------------------------------------

#tssas_prot = "507511|508235|AAQ73324|AFF60282"

#Chimeras
Variants = variants_full[(variants_full$protein %like% "ACY54509.1"),]

#Variants$change = variable aa (nchar = 10)
for (i in 1:Variants[, .N]){
    #i=1
    temp_prot = Variants[i,1]
    temp_change = substr(temp_prot, 2, 11)
    Variants$change[i] = temp_change
}

#TSSA* DT (1 per isoform)
TSSAI = variants_full[variants_full$protein %like% "AFF60282",]
TSSAI$change = "SKDKTAGTPS"

TSSAII = variants_full[variants_full$protein %like% "507511",]
TSSAII$change = "STENPTEAQP"

TSSAIII = variants_full[variants_full$protein %like% "508235",]
TSSAIII$change = "ITEKAAEAPS"

TSSAIV = variants_full[variants_full$protein %like% "AAQ73324",]
TSSAIV$change = "SKDKTAGTPS"

TSSAs = rbindlist(list(TSSAI, TSSAII, TSSAIII, TSSAIV))

TSSAII_length = variants_full[variants_full$protein %like% "mapped",]
TSSAII_length$change = "STENPTEAQP"

full_dataset = rbindlist(list(TSSAs, TSSAII_length, Variants))

#Add TSSAII 15mer to TSSAII_length DT
TSSA_8_aa = TSSAII[substr(TSSAII$sequence, 1, 8) %in% TSSAII_length$sequence, ]
TSSAII_length = rbindlist(list(TSSAII_length, TSSA_8_aa))


# Outputs ---------------------------------------------------------------------------

tssas_output_file = "Files/01_tssas_complete.tsv"
tssaii_output_file = "Files/01_tssaii_length.tsv"
variants_output_file = "Files/01_variants_complete.tsv"
full_dataset_variants_output_file = "Files/01_complete_dataset.tsv"

write.table(TSSAII_length, file = tssaii_output_file, col.names = T, row.names = F, sep = "\t", quote = T)
write.table(full_dataset, file = full_dataset_variants_output_file, col.names = T, row.names = F, sep = "\t", quote = T)
write.table(TSSAs, file = tssas_output_file, col.names = T, row.names = F, sep = "\t", quote = T)
write.table(Variants, file = variants_output_file, col.names = T, row.names = F, sep = "\t", quote = T)

