library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr) 
library(patchwork)
library(scales)

# Input files --------------------------------------------------------------------------

TSSAs_singles_replicas = "Files/09_TSSAs_singles_replicas.tsv"
TSSAs_singles_replicas <- fread(TSSAs_singles_replicas, header = T, sep = "\t", na.strings = NULL) 

# Data treatment----------------------------------------------------------------

TSSAs_mapped=TSSAs_singles_replicas

TSSAs_mapped$normalized_signal = TSSAs_mapped$normalized_signal/10000
TSSAs_mapped$normalized_signal_sd = TSSAs_mapped$normalized_signal_sd/10000
TSSAs_mapped$total_signal = TSSAs_mapped$total_signal/10000

TSSAs_mapped$middle = as.numeric()
for (i in 1:TSSAs_mapped[,.N]){
    #i=2
    TSSAs_mapped$middle[i] = TSSAs_mapped$start[i] + nchar(TSSAs_mapped$sequence[i])/2 - 1
    
}

country = as.vector(unique(TSSAs_mapped$parent_source))
sources = as.vector(unique(TSSAs_mapped$source))
prots = as.vector(unique(TSSAs_mapped$protein))


# Plots --------------------------------------------------------------------------

#Position vs signal, change sources and protein
sources_ar = unique(TSSAs_mapped[parent_source=="AR",]$source)
plot_list <- list() 
for (i in 1:length(sources_ar)) { # Loop over loop.vector
    #i=2
    temp_source = sources_ar[i]
    x <- TSSAs_mapped[source==temp_source&protein=="TSSA_II"]
    
    # Plot
    plot_list[[i]] = ggplot(x, aes(x = middle, y = normalized_signal)) +
        geom_point(alpha=0.6) +
        geom_line() +
        theme_minimal() +
        theme(panel.grid.minor = element_line(colour="gray87", size=0.2)) + 
        scale_x_continuous(minor_breaks=seq(1,85,1), breaks = seq(10, 85, 5), limits = c(20,65)) +
        ylim(NA, max(TSSAs_mapped$normalized_signal))+
        theme(legend.position = "bottom")+
        labs(x = "Posición", y = "Signal", title = sources_ar[i])
}
grid.arrange(grobs=plot_list,ncol=3)

#X = position, Y = source. Fill with mean signal. Change manually COUNTRY
#
low_color <- "#D55E00"
high_color <- "#56B4E9"
        
plot_list <- list() 
for (i in 1:length(prots)) {
    #i=1
    temp_prot = prots[i]
    x <- TSSAs_mapped[TSSAs_mapped$protein==temp_prot&parent_source=="AR",] #change source (AR, BO, BR, MX, CO or US)
    
    # Plot
    plot_list[[i]] = ggplot(x, aes(x = middle, y = graph_name))+
        geom_point(color="black", shape=22,size=4,aes(fill=normalized_signal),stroke = 0.1)+
        scale_fill_steps2(low = low_color,mid = "white", high = high_color, midpoint = 2, n.breaks = 10, limits=c(0,6))+
        labs(x = "Position", y = "Source", title = temp_prot)
    
}
        
plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + 
    plot_layout(ncol = 2, guides = "collect") + plot_annotation(title = paste("Source:",unique(x$parent_source)))

        

paletteLength <- 20
myColor <- function(paletteLength) {c(colorRampPalette(c(high_color, "white"))(paletteLength/2),colorRampPalette(c("white",low_color))(paletteLength/2))}

#X = protein, y = source. Fill = total signal
plot_list <- list() 
    for (i in 1:length(country)) { # Loop over loop.vector
        #i=1
        temp_country = country[i]
        x <- TSSAs_mapped[parent_source==temp_country,]
        
        # Plot
        plot_list[[i]] = ggplot(x, aes(x = protein, y = source))+
            geom_point(aes(fill=total_signal), color="black", shape=22, size=4, stroke = 0.2)+
            theme(legend.position="right")+
            scale_fill_stepsn(colors = myColor(paletteLength), n.breaks=12, limits=c(-1,84), 
                              values=rescale(c(-1,7,84)))+
            labs(x = "Protein", y = "Source", title = temp_country)
                    
}
                
grid.arrange(grobs=plot_list,ncol=3)
                
# Outputs ---------------------------------------------------------------------------
                
TSSAs_mapped_output_file = "Files/10_singles_TSSAs_mapped.tsv"
write.table(TSSAs_mapped, file = TSSAs_mapped_output_file, col.names = T, row.names = F, sep = "\t", quote = T)

                
                