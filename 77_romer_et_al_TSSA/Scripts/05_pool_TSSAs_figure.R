library(data.table)
library(ggplot2)
library(gridExtra)
library(plyr)
library(reshape2)

# Input files -----------------------------------------------------------------------------

tssas_file_complete = "Files/TSSA_variants_len15_allPools_noReplicas.tsv"
TSSAs_data_map_complete= fread(tssas_file_complete, header = T, sep = "\t", na.strings = NULL)
TSSAs_data_map_complete = TSSAs_data_map_complete[TSSAs_data_map_complete$type=="PO",]
TSSAs_data_map_complete$normalized_signal = as.numeric(TSSAs_data_map_complete$normalized_signal
                                                       )


tssas_file_complete = "Files/04_TSSA_aa_signal"
TSSAS_data_aa =  fread(tssas_file_complete, header = T, sep = "\t", na.strings = NULL)


# Plots --------------------------------------------------------------------------

sources = unique(TSSAS_data_aa$source)
sequences = unique(TSSAs_data_map_complete$protein)


#Heatmap TSSA complete sequence
plot_list <- list() 
for (i in 1:length(sources)) { # Loop over loop.vector
    #i=1
    temp_source = sources[i]
    x <- TSSAS_data_aa[TSSAS_data_aa$source==temp_source&TSSAS_data_aa$Pos>=1&TSSAS_data_aa$Pos<=80,]
    
    # Plot
    plot_list[[i]] = ggplot(x, aes(x = Pos, y = protein,label=AA)) + 
        geom_point(shape = 22, size = 5, aes(fill=Signal_mean/10000), color="grey", stroke = 0.2)+
        scale_fill_stepsn(
            colours=c("#f6f0f8","#fff8a6","#ffc268","#ff873a", "#e6211c", "#9e0023"),
            breaks=c(0.6,1,2,3,4,5), #0.567 = max. signal TSSAI in Arg.
            limits=c(0,6), #5.0591 es la señal maxima global
        )+
        geom_text(size=2.2,fontface = "bold")+
        theme_minimal()+
        theme(legend.position = "left")+
        theme(panel.grid.minor = element_line(colour="gray87", size=0.5)) + 
        #scale_fill_gradient2(midpoint = 2, low = "white", mid = "khaki2",
        #high = "firebrick2", limits = (c(0,max(TSSAS_data_aa$Signal_mean)/10000)))+
        labs(x = "Position", y = "Sequence", title = paste("Source:", temp_source))
}

grid.arrange(grobs=plot_list,ncol=1)



#Isoform signal, x = source
plot_list <- list() 
for (i in 1:length(sequences)) { 
    #i=1
    temp_protein = sequences[i]
    x <- TSSAs_data_map_complete[TSSAs_data_map_complete$protein==temp_protein,]
    x = na.omit(x, cols=seq_along(x), invert=FALSE)
    tssas_mean <- ddply(x, c("protein", "source", "type"), summarise, sum(normalized_signal))
    
    #Plot
    plot_list[[i]]=ggplot(tssas_mean, aes(x=source, y=sum/10000))+
        geom_jitter(aes(col=source),width = 0.25, height = 0)+ #jitter x point
        theme_minimal()+
        scale_y_continuous(name="Signal", limits=c(0, 100))+
        theme(axis.text.x = element_text(angle = 60, hjust = 1))+
        labs(x = "Source", y = "Signal", title = paste("TSSA", i))
    
}

grid.arrange(grobs=plot_list,ncol=2)


#X = source, y = signal, all isoforms
x = na.omit(TSSAS_data_aa, cols=seq_along(TSSAS_data_aa), invert=FALSE)
tssas_mean_full <- ddply(x, c("source","Seq", "protein") , summarise, sum(Signal_mean)) #tabla que tiene la media para cada variante
setnames(tssas_mean_full, "..1", "sum")

# Plot
ggplot(tssas_mean_full, aes(x=source, y=sum/10000))+
    geom_jitter(aes(col=protein),width = 0.05, height = 0)+ #jitter x point
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
