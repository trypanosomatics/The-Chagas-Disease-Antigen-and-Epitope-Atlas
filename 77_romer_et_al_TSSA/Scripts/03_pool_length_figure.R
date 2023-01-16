library(ggplot2)
library(gridExtra)
library(data.table)

# Input files -----------------------------------------------------------------------------

tssa_length_file = "Files/01_tssaii_length.tsv"
TSSAII_mapped_length = fread(tssa_length_file, header = T, sep = "\t", na.strings = NULL)

TSSAII_len_positive = TSSAII_mapped_length[TSSAII_mapped_length$type=="PO",] #descarto sueros negativos 

#Final DT. Skip the data treatment section and load this.
temp_protein <- read.delim("Files/03_TSSAII_length_dt.tsv")

# Data treatment --------------------------------------------------------------------

temp_protein = TSSAII_len_positive #change start and end

for (i in 1:temp_protein[,.N]) {
    #i=1063
    if (temp_protein$protein[i]==">TcCLB.507511.91 | Trypanosoma cruzi CL Brener Esmeraldo-like | mucin TcMUCIII, putative | protein  |"){
        temp_protein$start[i] = temp_protein$start[i] - 5
        temp_protein$end[i] = temp_protein$end[i] - 5
    }
    
}

min_start = min(temp_protein$start)

temp_protein$ref_mat = temp_protein$start - min_start + 24 #start 1 = position 24

temp_protein$middle = as.numeric()

for (i in 1:temp_protein[,.N]){
    #i=2
    #even number criteria
    if(temp_protein$sequence_length[i]%%2==0) temp_protein$middle[i] = temp_protein$ref_mat[i] + temp_protein$sequence_length[i]/2 - 1
    #odd number criteria 
    else temp_protein$middle[i] = temp_protein$ref_mat[i] + temp_protein$sequence_length[i]/2 - 0.5
}


# Figures --------------------------------------------------------------------------

sources = unique(temp_protein$source)
plot_list_middle <- list() 
for (i in 1:length(sources)) {
    #i=1
    temp_source = sources[i]
    x <- temp_protein[temp_protein$source==temp_source,]
    
    # Plot
    plot_list_middle[[i]] = ggplot(x, aes(x = middle, y = (normalized_signal/10000), color = as.character(sequence_length))) +
        geom_point() +
        geom_line() +
        theme_minimal() +
        theme(panel.grid.minor = element_line(colour="gray87", size=0.5)) + 
        ylim(NA, max(temp_protein$normalized_signal/10000))+
        theme(legend.position = "bottom") +
        labs(x = "Position", y = "Signal", title = paste("Source:", temp_source))
}

grid.arrange(grobs=plot_list_middle,ncol=2)


#All kmers for one source, change SOURCE
#
kmers = unique(temp_protein$sequence_length)
kmers = sort(kmers)


color_plot = c("#374E55", "#80796B", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "white", "#6A6599")
plot_list_kmer <- list() 
#for (i in 1:length(largos)) { # Loop over loop.vector
for (i in 1:8) {
    #i=1
    temp_len = i+7
    x <- temp_protein[temp_protein$sequence_length==temp_len&temp_protein$source=="US",]
    
    # Plot
    plot_list_kmer[[i]] = ggplot(x, aes(x = middle, y = (normalized_signal/10000))) +
        geom_area(fill="lightgrey", alpha=0.3)+
        geom_point(color = color_plot[i]) +
        geom_line(color = color_plot[i]) +
        theme_minimal() +
        theme(panel.grid.minor = element_line(colour="gray87", size=0.5)) + 
        scale_x_continuous(minor_breaks = seq(min(temp_protein$middle), max(temp_protein$middle), 1), limits = c(min(temp_protein$middle) - 0.5, 0.5 + max(temp_protein$middle))) +
        ylim(NA, max(temp_protein$normalized_signal/10000))+
        labs(x = "Position", y = "Signal", title = paste("Length", temp_len))
}

grid.arrange(grobs=plot_list_kmer,ncol=2)


# Outputs ---------------------------------------------------------------------------

tssaII_length_outputfile = "Archivos/03_TSSAII_length_dt.tsv"
write.table(temp_protein, file = tssaII_length_outputfile, col.names = T, row.names = F, sep = "\t", quote = T)
