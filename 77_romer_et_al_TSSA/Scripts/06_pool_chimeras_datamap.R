library(data.table)
library(stringdist)

# Input files -----------------------------------------------------------------------------
variants_file = "Files/02_data_map.tsv"
variants_data_map= fread(variants_file, header = T, sep = "\t", na.strings = NULL)

sec_original = ("TSSTPPSGTENKPATGEAPSQPGASS")

# Data treatment --------------------------------------------------------------------

#Number of changes vs reference TSSAII 
for (i in 1:variants_data_map[,.N]) {
    #i=2
    variants_data_map$Changes[i]=stringdist(sec_original, variants_data_map$sec[i], 
                                                method="hamming", p=0)  
    
}

#Treatment for different # changes
Data_plot1 = variants_data_map[variants_data_map$Changes<=1,] #Set here # changes

Values_plot1 = Values_plot1 = data.table(AA=as.character(),
                                           Seq=as.character(),
                                           Var=as.character(),
                                           Signal_mean=as.numeric(),
                                           Pos=as.numeric(),
                                           Source=as.character()
)


mapped_seq = as.character(unique(Data_plot1$sec))
sources = unique(Data_plot1$source)

for(k in 1:length(sources)){
    #k=1
    temp_source=sources[k]
    for (i in 1:length(mapped_seq)) {
        #i=1
        temp_sec=mapped_seq[i]
        for (j in 1:nchar(temp_sec)) {
            #j=1
            temp_aa = substr(temp_sec,j,j)
            temp_signal = mean(Data_plot1$signal[Data_plot1$sec==temp_sec&Data_plot1$start<=j&Data_plot1$start+14>=j&Data_plot1$source==temp_source],na.rm = T)
            temp_var = unique(Data_plot1[Data_plot1$sec==temp_sec,var])
            
            temp_Values_plot1 = data.table(AA=temp_aa,
                                            Seq=temp_sec,
                                            Var=temp_var,
                                            Signal_mean=temp_signal,
                                            Pos=j,
                                            Source=temp_source
            )
            Values_plot1 = rbind(Values_plot1,temp_Values_plot1)
            
        }
        
    }
    
}

# Outputs ---------------------------------------------------------------------------

chimera1_output_file = "Files/06_chimera_1change_aa_signal"
write.table(Values_plot1, file = chimera1_output_file, col.names = T, row.names = F, sep = "\t", quote = T)

variants_data_map_output_file = "Files/06_data_map_n_change"
write.table(variants_data_map, file = variants_data_map_output_file, col.names = T, row.names = F, sep = "\t", quote = T)

