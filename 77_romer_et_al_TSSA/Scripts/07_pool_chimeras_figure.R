library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
library(data.table)

# Input files -----------------------------------------------------------------------------

Values_plot_1_file = "Files/06_chimera_1change_aa_signal"
Values_plot_1 <- fread(Values_plot_1_file, header = T, sep = "\t", na.strings = NULL)

Data_n_change_file = "Files/06_data_map_n_change"
Data_n_change <- fread(Data_n_change_file, header = T, sep = "\t", na.strings = NULL)
Data_plot_1 = Data_n_change[Data_n_change$Changes<=1,]

# Plots --------------------------------------------------------------------------

Values_plot_1[["Var"]] <- factor(as.character(Values_plot_1$Var), levels = c("STENPTEAQP", "ITENPTEAQP", "SKENPTEAQP", "STDNPTEAQP", "STEKPTEAQP", "STENATEAQP", 
                                                                             "STENTTEAQP", "STENPAEAQP", "STENPTGAQP", "STENPTETQP", "STENPTEAPP", "STENPTEAQS"))
sources = unique(Values_plot_1$Source)

low_color <- "#D55E00"
high_color <- "#56B4E9"

#AA Signal        
plot_list <- list() 
for (i in 1:length(sources)) { # Loop over loop.vector
#i=6
    temp_source = sources[i]
    x <- Values_plot_1[Values_plot_1$Source==temp_source,]

    # Plot
    plot_list[[i]] = ggplot(x, aes(x = Pos, y = Var,label=AA)) + 
    geom_point(shape = 22, size = 5, aes(fill=Signal_mean/10000), stroke = 0)+
    scale_fill_steps2(low = low_color,mid = "white", high = high_color, midpoint = 2, n.breaks = 10, limits=c(0,6))+
    geom_text(size=2,fontface = "bold")+
    theme_minimal()+
    theme(legend.position = "left")+
    theme(panel.grid.minor = element_line(colour="gray87", size=0.5)) + 
    labs(title = paste("Sueros", temp_source))
}
        
grid.arrange(grobs=plot_list,ncol=2)
        
#Mean of all peptides spanning the same chimeric sequence. X = Chimera, Y = signal
chimeras_1_mean <- ddply(Data_plot_1, c("sec","var", "source"), colwise(mean)) #tabla que tiene la media para cada variante
chimeras_1_mean[["var"]] <- factor(as.character(chimeras_1_mean$var), levels = c("STENPTEAQP", "ITENPTEAQP", "SKENPTEAQP", "STDNPTEAQP", "STEKPTEAQP", "STENATEAQP", 
                                                                                                           "STENTTEAQP", "STENPAEAQP", "STENPTGAQP", "STENPTETQP", "STENPTEAPP", "STENPTEAQS"))
plot_list <- list() 
for (i in 1:length(sources)) { # Loop over loop.vector
    #i=2
    temp_source = sources[i]
    x <- chimeras_1_mean[chimeras_1_mean$source==temp_source,]
        
    # Plot
    plot_list[[i]] = ggplot(x, aes(x=var, y=signal/10000))+
    geom_jitter(aes(col=var),width = 0, height = 0, show.legend = FALSE)+ #jitter x point
    theme_minimal()+
    scale_y_continuous(name="Signal", limits=c(0, 6), breaks = c(1:5))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    labs(title = paste("Source:", temp_source))
}
        
grid.arrange(grobs=plot_list,ncol=3)
        
#All sources
ggplot(chimeras_1_mean, aes(x=var, y=signal/10000))+
geom_point(aes(col=source))+
theme_minimal()+
scale_y_continuous(name="Señal", limits=c(0, 5.5), breaks = c(1:5))+
theme(axis.text.x = element_text(angle = 60, hjust = 1))
        
#Sum of all peptides spanning the same chimeric sequence. X = Chimera, Y = signal
chimeras_1_sum <- ddply(Data_plot_1, c("sec","var", "source"), summarise, sum(signal)) #tabla que tiene la señal total para cada variante
setnames(chimeras_1_sum, "..1", "sum")
chimeras_1_sum[["var"]] <- factor(as.character(chimeras_1_sum$var), levels = c("STENPTEAQP", "ITENPTEAQP", "SKENPTEAQP", "STDNPTEAQP", "STEKPTEAQP", "STENATEAQP", 
                                                                                                         "STENTTEAQP", "STENPAEAQP", "STENPTGAQP", "STENPTETQP", "STENPTEAPP", "STENPTEAQS"))
plot_list <- list() 
for (i in 1:length(sources)) { # Loop over loop.vector
    #i=2
    temp_source = sources[i]
    x <- chimeras_1_sum[chimeras_1_sum$source==temp_source,]
            
    # Plot
    plot_list[[i]] = ggplot(x, aes(x=var, y=sum/10000))+
    geom_jitter(aes(col=var),width = 0, height = 0, show.legend = FALSE)+ #jitter x point
    theme_minimal()+
    scale_y_continuous(name="Signal", limits=c(0, 61))+ #El maxino es 60.algo
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    labs(title = paste("Source:", temp_source))
}
        
grid.arrange(grobs=plot_list,ncol=3)
        
#All sources
ggplot(chimeras_1_sum, aes(x=var, y=sum/10000, label=source))+
geom_point(shape = 21, size = 2, aes(fill=source), color="grey", stroke = 0.2)+
scale_fill_manual(values=c("#7F3C8D","#E73F74","#11A579","#3969AC","#f97b72","#4b4b8f"))+
theme_minimal()+
scale_y_continuous(name="Signal", limits=c(0, 61))+
theme(axis.text.x = element_text(angle = 60, hjust = 1))

    
# Supplementary plots -------------------------------------------------------------------------

#S2
secs = c("STENPTEAQP", "SKENPTEAQP", "STEKPTEAQP", "SKEKPTEAQP", "STENTTEAQP", "STENPAEAQP", "STENPTGAQP",
                 "STENPTETQP", "STENTAGTQP", "SKENTTEAQP", "SKEKTAGTQP")
        
table_s2 = Data_n_change[Data_n_change$var%in%secs&Data_n_change$source=="BO",]
        
Values_s2 = data.table(AA=as.character(),
                        Seq=as.character(),
                        Var=as.character(),
                        Signal_mean=as.numeric(),
                        Pos=as.numeric(),
                        Source=as.character()
)
        
Mapped_seq = as.character(unique(table_s2$sec))
sources = unique(table_s2$source)
    
for(k in 1:length(sources)){
    #k=1
    temp_source=sources[k]
    for (i in 1:length(Mapped_seq)) {
        #i=1
        temp_sec=Mapped_seq[i]
        for (j in 1:nchar(temp_sec)) {
            #j=1
            temp_aa = substr(temp_sec,j,j)
            temp_signal = mean(table_s2$signal[table_s2$sec==temp_sec&table_s2$start<=j&table_s2$start+14>=j&table_s2$source==temp_source],na.rm = T)
            temp_var = as.character(unique(table_s2[table_s2$sec==temp_sec,"var"]))
        
            temp_Values_s2 = data.table(AA=temp_aa,
                                        Seq=temp_sec,
                                        Var=temp_var,
                                        Signal_mean=temp_signal,
                                        Pos=j,
                                        Source=temp_source
                    )
            Values_s2 = rbind(Values_s2,temp_Values_s2)
                    
}}}
        
ggplot(Values_s2, aes(x = Pos, y = Var,label=AA)) + 
    geom_point(shape = 22, size = 5, aes(fill=Signal_mean/10000), stroke = 0)+
    scale_fill_steps2(low = low_color,mid = "white", high = high_color, midpoint = 2, n.breaks = 10, limits=c(0,6))+
    geom_text(size=2,fontface = "bold")+
    theme_minimal()+
    theme(legend.position = "left")+
    theme(panel.grid.minor = element_line(colour="gray87", size=0.5)) + 
    labs(title = paste("Source:", temp_source))
