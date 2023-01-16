library(data.table)
library(stringdist)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

# Input files-------------------------------------------------------------

variants_file = "Files/06_data_map_n_change"
variants_data_map= fread(variants_file, header = T, sep = "\t", na.strings = NULL)

# Data treatment-------------------------------------------------------------

#Relevant and non-relevant changes (US)
variants_data_map = variants_data_map[variants_data_map$source=="US",]
variants_data_map$relevant_changes=0
variants_data_map$non_relevant_changes=0

#NR positions based on graphs obtained in 07_pool_chimeras_figure.R - NR positions may be different between regions.
original_var = ("STENPTEAQP")
var_NR_1 = ("ITENPTEAQP") 
var_NR_2 = ("STDNPTEAQP")
var_NR_3 = ("STENPTEAPP")
var_NR_4 = ("STENPTEAQS")

for (i in 1:variants_data_map[,.N]) {
    temp_change = 0
    #i=133
    
    #These ifs compare the sequence with each non-relevant sequence. 
    #If the aa is in the variant, the number of changes with respect to the original sequence
    #should be one less, so increase the non-relevant change by 1.
    if (stringdist(var_NR_1, variants_data_map$var[i], 
                   method="hamming", p=0)==stringdist(original_var, variants_data_map$var[i], method="hamming", p=0)-1){
        temp_change = temp_change + 1
    }
    if (stringdist(var_NR_2, variants_data_map$var[i], 
                   method="hamming", p=0)==stringdist(original_var, variants_data_map$var[i], method="hamming", p=0)-1){
        temp_change = temp_change + 1
    }
    if (stringdist(var_NR_3, variants_data_map$var[i], 
                   method="hamming", p=0)==stringdist(original_var, variants_data_map$var[i], method="hamming", p=0)-1){
        temp_change = temp_change + 1
    }
    if (stringdist(var_NR_4, variants_data_map$var[i], 
                   method="hamming", p=0)==stringdist(original_var, variants_data_map$var[i], method="hamming", p=0)-1){
        temp_change = temp_change + 1
    }
    
    variants_data_map$non_relevant_changes[i]=temp_change
    variants_data_map$relevant_changes[i]=stringdist(original_var, variants_data_map$var[i], method="hamming", p=0)-temp_change
}


# Mean plots-------------------------------------------------------------
mean_plots_US <- ddply(variants_data_map, c("sec","var", "source"), colwise(mean)) #mean for each variant
mean_plots_US <- mean_plots_US[,c(-4,-5)]
mean_plots_US$Changes = as.factor(mean_plots_US$Changes)

data_stack_mean <- melt(mean_plots_US, id = c("var", "signal", "Changes", "source", "sec"))

p1 = ggplot(data_stack_mean[data_stack_mean$source=="US",], aes(x=value, y=signal/10000, col = variable))+
    geom_point()+
    xlab("Changes")+
    ylab("Signal")+
    facet_wrap(vars(variable)) 
p1

p_NR_mean = ggplot(mean_plots_US[mean_plots_US$relevant_changes==0,], aes(x=Changes, y=signal/10000))+
    geom_violin(colour="violetred1", scale = "width", adjust = 0.8)+
    scale_y_continuous(name="Signal", limits=c(0, 6))+
    geom_jitter(aes(col=Changes),width = 0.1, height = 0.1)+
    labs(title = "Non-relevant changes")

p_R_mean = ggplot(mean_plots_US[mean_plots_US$non_relevant_changes==0,], aes(x=Changes, y=signal/10000))+
    geom_violin(colour="lightskyblue", scale = "width", adjust = 0.8)+ 
    scale_y_continuous(name="Signal", limits=c(0, 6))+
    geom_jitter(aes(col=Changes),width = 0.1, height = 0.1)+
    labs(title = "Relevant changes")

full_US_mean = ggarrange(p_NR_mean, p_R_mean, align = "h")
full_US_mean

mean_plots_US$total_changes = 0
for (i in 1:nrow(mean_plots_US)){
    #i=6
    mean_plots_US$total_changes[i] = paste(mean_plots_US[i,"relevant_changes"],"R","+",mean_plots_US[i,"non_relevant_changes"],"NR")
}

mean_plots_US$total_changes = as.factor(mean_plots_US$total_changes)

ggplot(mean_plots_US, aes(x=total_changes, y=signal/10000))+
    geom_jitter(colour="gray38", shape=21, aes(fill = Changes), size = 2,
                width = 0.25, height = 0) +
    theme_minimal()+
    scale_y_continuous(name="Señal", limits=c(0, 5), breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    scale_fill_brewer(palette = "Paired")

# Sum plots-------------------------------------------------------------

sum_plots_US <- ddply(variants_data_map, c("sec","var", "source", "Changes", "non_relevant_changes", "relevant_changes"), summarise, sum(signal))
setnames(sum_plots_US, "..1", "sum")
sum_plots_US$Changes = as.factor(sum_plots_US$Changes)

data_stack_sum <- melt(sum_plots_US, id = c("var", "sum", "Changes", "source", "sec"))
p2 = ggplot(data_stack_sum[data_stack_sum$source=="US",], aes(x=value, y=sum/10000, col = variable))+
    geom_point()+
    xlab("Changes")+
    ylab("Signal")+
    facet_wrap(vars(variable)) 
p2

p_NR_sum = ggplot(sum_plots_US[sum_plots_US$relevant_changes==0,], aes(x=Changes, y=sum/10000))+
    geom_violin(colour="violetred1", scale = "width", adjust = 0.8)+
    scale_y_continuous(name="Signal", limits=c(10, 60))+
    geom_jitter(aes(col=Changes),width = 0.1, height = 0.1)+
    labs(title = "Non-relevant changes")

p_R_sum = ggplot(sum_plots_US[sum_plots_US$non_relevant_changes==0,], aes(x=Changes, y=sum/10000))+
    geom_violin(colour="violetred1", scale = "width", adjust = 0.8)+
    scale_y_continuous(name="Signal", limits=c(10, 60))+
    geom_jitter(aes(col=Changes),width = 0.1, height = 0.1)+
    labs(title = "Relevant changes")

full_US_sum = ggarrange(p_NR_sum, p_R_sum, align = "h")
full_US_sum

sum_plots_US$total_changes = 0
for (i in 1:nrow(sum_plots_US)){
    #i=6
    sum_plots_US$total_changes[i] = paste(sum_plots_US[i,"relevant_changes"],"R","+",sum_plots_US[i,"non_relevant_changes"],"NR")
}

sum_plots_US$total_changes = as.factor(sum_plots_US$total_changes)

ggplot(sum_plots_US, aes(x=total_changes, y=sum/10000))+
    geom_jitter(colour="gray38", shape=21, aes(fill = Changes), size = 2,
                width = 0.25, height = 0) +
    theme_minimal()+
    scale_y_continuous(name="Señal", limits = c(0,65))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    scale_fill_brewer(palette = "Paired")

# Outputs ---------------------------------------------------------------------------

variants_data_map_output_file = "Files/08_chimeras_US.tsv"
write.table(variants_data_map, file = variants_data_map_output_file, col.names = T, row.names = F, sep = "\t", quote = T)
