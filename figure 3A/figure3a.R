library(tidyverse)    
library(fst)  
library(ggfittext) 
library(data.table)              
library(RColorBrewer)              
library(grid)              
library(gridExtra)              
library(cowplot)              
library(ggpubr)              
library(ggsci)              
library(gtable) 
data_all <- read.csv("RA_baseline_clin.csv")
colnames(data_all)[25]<-"early.disease"
colnames(data_all)[26]<-"established.disease"

data_all$NEUTROPHIL <- as.numeric(data_all$NEUTROPHIL)
data_all$LYMPHOCYTE <- as.numeric(data_all$LYMPHOCYTE)
data_all$MONOCYTE <- as.numeric(data_all$MONOCYTE)
data_all$ESOPHIL <- as.numeric(data_all$ESOPHIL)
data_all$BASOPHIL <- as.numeric(data_all$BASOPHIL)
data_all$WBC <- as.numeric(data_all$WBC)

cell_total <- data_all$NEUTROPHIL + data_all$LYMPHOCYTE + data_all$MONOCYTE + data_all$ESOPHIL + data_all$BASOPHIL

data_all$NEUTROPHIL_ratio <- data_all$NEUTROPHIL / cell_total
data_all$LYMPHOCYTE_ratio <- data_all$LYMPHOCYTE / cell_total
data_all$MONOCYTE_ratio <- data_all$MONOCYTE / cell_total
data_all$ESOPHIL_ratio <- data_all$ESOPHIL / cell_total
data_all$BASOPHIL_ratio <- data_all$BASOPHIL / cell_total

#Extract the legend from a ggplot object with this function:              
g_legend <- function(a.gplot){              
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))              
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")              
  if(length(leg) == 0){              
    return(NULL)              
  }              
  legend <- tmp$grobs[[leg]]              
  return(legend)              
}              



#Then use the following function on a list of legends to return a single grob that will contain all the legends aligned:              
align.legends.fun  <- function(legend.list){              
  aligned.leg <- legend.list[[1]]              
  for(i in 2:length(legend.list)){              
    leg1        <- legend.list[[i]]$grobs[[1]]              
    leg_a       <- gtable_add_rows(aligned.leg, pos = nrow(aligned.leg) - 1, heights = sum(leg1$heights))              
    leg_final   <- gtable_add_grob(leg_a, leg1, t = nrow(leg_a) - 1, l = 3)              
    aligned.leg <- leg_final              
  }              
  return(aligned.leg)              
}  
## set gab_size              
gab_size    <- 0.2              
text_size   <- 18              
legend_text_size  <- 13              
legend_title_size <- 15 

WBC_list <- c("NEUTROPHIL_ratio", "LYMPHOCYTE_ratio", "MONOCYTE_ratio", "ESOPHIL_ratio", "BASOPHIL_ratio")
DAS28_list <- c("share_VAS", "share_TJC", "share_SJC", "share_CRP")
Sub_list <- c("Clinical Remission", "low disease activity", "moderate disease activity", "high disease activity")

order_samples <- c("Clinical Remission", "low disease activity", "moderate disease activity", "high disease activity")

data_all$disease.activity <- factor(data_all$disease.activity, levels = order_samples)

#####plot!
p_combined_list <- lapply(Sub_list, function(x){              
  if(x == Sub_list[1]){              
    data <- data_all[data_all$disease.activity == Sub_list[1],]              
  } else if (x == Sub_list[2]){              
    data <- data_all[data_all$disease.activity == Sub_list[2],]              
  } else if (x == Sub_list[3]){              
    data <- data_all[data_all$disease.activity == Sub_list[3],]              
  } else if (x == Sub_list[4]){              
    data <- data_all[data_all$disease.activity == Sub_list[4],]              
  }         

data <- data %>%
  arrange(disease.activity,share_TJC, share_SJC, share_VAS, share_CRP)
ordered_samples <- data$sample
data$sample <- factor(data$sample, levels = ordered_samples)

p_subtypes <- ggplot(data, aes(x = sample, y = 1, fill = disease.activity)) +               
  geom_tile(colour = 'white', size = gab_size) +              
  scale_fill_manual(name = 'disease activity',               
                    values = c("Clinical Remission" = "#FDBF6F","low disease activity" ="#A6CEE3", "moderate disease activity"="#B2DF8A","high disease activity"="#FB9A99" ), drop = F) +               
  scale_y_continuous(expand = c(0,0)) +                
  scale_x_discrete(expand = c(0,0), drop = F) +              
  ylab('disease activity') +              
  theme_bw() +              
  theme(axis.text = element_blank(), axis.ticks = element_blank(),               
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
        text = element_text(size = 20),              
        legend.text = element_text(size = 15), legend.title = element_text(size = 18))        

plot_data <- data %>%
  gather(key = share_DAS28, value = proportion, share_VAS, share_TJC, share_SJC, share_CRP) %>%
  mutate(share_DAS28 = factor(share_DAS28, levels = c("share_VAS", "share_TJC", "share_SJC", "share_CRP")))

p_prop.subtype <- ggplot(plot_data, aes(x = sample, y = proportion, fill = share_DAS28)) +
  geom_bar(stat = 'identity', position = position_stack(reverse = F)) +
  scale_fill_manual(name = 'DAS28-CRP Component', values = c("share_VAS" = "#EFC000FF", 
                                                         "share_TJC" = "#1F78B4", 
                                                         "share_SJC" = "#B2DF8A", 
                                                         "share_CRP" = "#A73030FF"), drop = F) +
  scale_x_discrete(drop = F) +
  scale_y_continuous(expand = c(0,0)) + ylab('% Proportion') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18), panel.grid = element_blank()) 

plot_data2 <- data %>%
  gather(key = cell_type, value = ratio, NEUTROPHIL_ratio, LYMPHOCYTE_ratio, MONOCYTE_ratio, ESOPHIL_ratio, BASOPHIL_ratio) %>%
  mutate(cell_type = factor(cell_type))


p_wbc_ratios <- ggplot(plot_data2, aes(x = sample, y = ratio, fill = cell_type)) +
  geom_bar(stat = 'identity', position = position_stack(reverse = F)) +
  scale_color_lancet(name = 'WBC Cell Type Ratio') +
  scale_x_discrete(drop = F) +
  scale_y_continuous(expand = c(0,0)) + ylab('cell ratio') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18), panel.grid = element_blank())


p_age <- ggplot(data, aes(x = sample, y = 0.6 , fill = age)) +              
  geom_tile(colour = 'white', size = gab_size) +              
  scale_fill_gradient2(name = 'Age',  low = "#4575B4" , mid = "#FFFFBF", high = "#D73027" , midpoint = 50, limits = c(16,77), na.value = 'white') +              
  scale_y_continuous(expand = c(0,0)) +              
  scale_x_discrete(expand = c(0,0), drop = F) +              
  ylab('Age') +              
  theme_bw() +              
  theme(axis.text = element_blank(), axis.ticks = element_blank(),              
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
        text = element_text(size = 20),              
        legend.text = element_text(size = 15), legend.title = element_text(size = 18))              

# plot Gender 
data <- data %>%
  mutate(GENDER = ifelse(GENDER == "F", "female", 
                         ifelse(GENDER == "M", "male", GENDER)))
             
p_sex <- ggplot(data, aes(x = sample, y = 0.6, fill = GENDER)) +              
  geom_tile(colour = 'white', size = gab_size) +              
  scale_fill_manual(name = 'Gender',              
                    values = c("#FB9A99","#A6CEE3"), na.value = 'white', drop = F) +              
  scale_y_continuous(expand = c(0,0)) +              
  scale_x_discrete(expand = c(0,0), drop = F) +              
  ylab('Gender') +              
  theme_bw() +              
  theme(axis.text = element_blank(), axis.ticks = element_blank(),              
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
        text = element_text(size = 20),              
        legend.text = element_text(size = 15), legend.title = element_text(size = 18))
# plot ACPA              
data$ACPA.subtype <- ifelse(data$ACPA.subtype == "", NA, data$ACPA.subtype)
p_ACPA <- ggplot(data, aes(x = sample, y = 0.6, fill = ACPA.subtype)) +              
  geom_tile(colour = 'white', size = gab_size) +              
  scale_fill_manual(name = 'ACPA',              
                    values = c("#4575B4","#EFC000FF"), na.value = 'white', drop = F,na.translate = FALSE) +              
  scale_y_continuous(expand = c(0,0)) +              
  scale_x_discrete(expand = c(0,0), drop = F) +              
  ylab('ACPA') +              
  theme_bw() +              
  theme(axis.text = element_blank(), axis.ticks = element_blank(),              
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
        text = element_text(size = 20),              
        legend.text = element_text(size = 15), legend.title = element_text(size = 18))
#####RF
data$RF <- ifelse(data$RF == "", NA, data$RF)
p_RF <- ggplot(data, aes(x = sample, y = 0.6, fill = RF)) +              
  geom_tile(colour = 'white', size = gab_size) +              
  scale_fill_manual(name = 'RF',              
                    values = c("#B2DF8A","#FB9A99"), na.value = 'white', drop = F,na.translate = FALSE) +              
  scale_y_continuous(expand = c(0,0)) +              
  scale_x_discrete(expand = c(0,0), drop = F) +              
  ylab('RF') +              
  theme_bw() +              
  theme(axis.text = element_blank(), axis.ticks = element_blank(),              
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
        text = element_text(size = 20),              
        legend.text = element_text(size = 15), legend.title = element_text(size = 18))
###HAQ
p_HAQ <- ggplot(data, aes(x = sample, y = 0.6 , fill = log2(HAQ+1))) +              
  geom_tile(colour = 'white', size = gab_size) +              
  scale_fill_gradient2(name = 'log2(HAQ)',  low = "#B2DF8A" , mid = "#FFFFBF", high = "orange" , midpoint=log2(6),limits = c(0,log2(81)), na.value = 'white') +              
  scale_y_continuous(expand = c(0,0)) +              
  scale_x_discrete(expand = c(0,0), drop = F) +              
  ylab('log2(HAQ)') +              
  theme_bw() +              
  theme(axis.text = element_blank(), axis.ticks = element_blank(),              
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
        text = element_text(size = 20),              
        legend.text = element_text(size = 15), legend.title = element_text(size = 18)) 
###duration
data$early.disease <- ifelse(data$early.disease == "0", "early disease",
                                   ifelse(data$early.disease == "Y", "established disease", 
                                    data$early.disease))
###duration
p_duration <- ggplot(data, aes(x = sample, y = 0.6, fill = early.disease)) +              
  geom_tile(colour = 'white', size = gab_size) +              
  scale_fill_manual(name = 'disease duration',              
                    values = c("early disease"="#BDD7E7", "established disease"="#08519C"), na.value = 'white', drop = F) +              
  scale_y_continuous(expand = c(0,0)) +              
  scale_x_discrete(expand = c(0,0), drop = F) +              
  ylab('disease duration') +              
  theme_bw() +              
  theme(axis.text = element_blank(), axis.ticks = element_blank(),              
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
        text = element_text(size = 20),              
        legend.text = element_text(size = 15), legend.title = element_text(size = 18))              
#####BMI
p_BMI <- ggplot(data, aes(x = sample, y = 0.6 , fill = BMI)) +              
  geom_tile(colour = 'white', size = gab_size) +              
  scale_fill_gradient2(low = "darkgreen", high = "red", mid = "yellow",midpoint=22,limits = c(14,47), na.value = 'white') +              
  scale_y_continuous(expand = c(0,0)) +              
  scale_x_discrete(expand = c(0,0), drop = F) +              
  ylab('BMI') +              
  theme_bw() +              
  theme(axis.text = element_blank(), axis.ticks = element_blank(),              
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
        text = element_text(size = 20),              
        legend.text = element_text(size = 15), legend.title = element_text(size = 18))  
#########################################              
###          combine plots             ###              
##########################################              
#align legends              
legends <- list(g_legend(p_subtypes),       
                g_legend(p_prop.subtype),   
                g_legend(p_wbc_ratios),     
                g_legend(p_age),          
                g_legend(p_sex),            
                g_legend(p_ACPA),           
                g_legend(p_RF),            
                g_legend(p_HAQ),         
                g_legend(p_duration),      
                g_legend(p_BMI))                     



aligend.legends1 <- align.legends.fun(legends[c(1:10)]) #               



#exclude legends from original plots              
p_subtypes <- p_subtypes + theme(legend.position = 'none')
p_prop.subtype <- p_prop.subtype + theme(legend.position = 'none')
p_wbc_ratios <- p_wbc_ratios + theme(legend.position = 'none')
p_age <- p_age + theme(legend.position = 'none')
p_sex <- p_sex + theme(legend.position = 'none')
p_ACPA <- p_ACPA + theme(legend.position = 'none')
p_RF <- p_RF + theme(legend.position = 'none')
p_HAQ <- p_HAQ + theme(legend.position = 'none')
p_duration <- p_duration + theme(legend.position = 'none')
p_BMI <- p_BMI + theme(legend.position = 'none')         



#delete y-labels except the one on the very left              
if(x != Sub_list[1]){ 
  p_subtypes <- p_subtypes + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_prop.subtype <- p_prop.subtype + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_wbc_ratios <- p_wbc_ratios + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_age <- p_age + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_sex <- p_sex + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_ACPA <- p_ACPA + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_RF <- p_RF + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_HAQ <- p_HAQ + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_duration <- p_duration + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_BMI <- p_BMI + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}              


#####combine both
p_subtypes <- ggplotGrob(p_subtypes + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
p_prop.subtype <- ggplotGrob(p_prop.subtype + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
p_wbc_ratios <- ggplotGrob(p_wbc_ratios + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
p_age <- ggplotGrob(p_age + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
p_sex <- ggplotGrob(p_sex + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
p_ACPA <- ggplotGrob(p_ACPA + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
p_RF <- ggplotGrob(p_RF + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
p_HAQ <- ggplotGrob(p_HAQ + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
p_duration <- ggplotGrob(p_duration + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
p_BMI <- ggplotGrob(p_BMI + theme(plot.margin = unit(c(0.1, 0.1, 1, 0.1), "cm")))          



grid_unit_pmax <- grid::unit.pmax(p_subtypes$widths, p_prop.subtype$widths, p_wbc_ratios$widths,
                                  p_age$widths, p_sex$widths, p_ACPA$widths, p_RF$widths, 
                                  p_HAQ$widths, p_duration$widths, p_BMI$widths)

p_subtypes$widths <- grid_unit_pmax
p_prop.subtype$widths <- grid_unit_pmax
p_wbc_ratios$widths <- grid_unit_pmax
p_age$widths <- grid_unit_pmax
p_sex$widths <- grid_unit_pmax
p_ACPA$widths <- grid_unit_pmax
p_RF$widths <- grid_unit_pmax
p_HAQ$widths <- grid_unit_pmax
p_duration$widths <- grid_unit_pmax
p_BMI$widths <- grid_unit_pmax


g <- gtable_rbind(p_subtypes,
                  p_prop.subtype,
                  p_HAQ,
                  p_age,
                  p_sex,
                  p_wbc_ratios,
                  p_ACPA,
                  p_RF,
                  p_duration,
                  p_BMI)


id_panels_h <- unique(g$layout[g$layout$name == "panel", "t"])

g$heights[id_panels_h] <- grid::unit(c(0.6, 2,0.5, 0.5, 0.5,2, 0.5, 0.5, 0.5, 0.5), "null")

g_combined <- arrangeGrob(grobs = list(g, aligend.legends1),               
                          nrow = 1,              
                          ncol = 2, #3,              
                          layout_matrix = matrix(c(rep(1, 4), rep(2,4)), nrow = 1))

g_combined<-g
return(g_combined)              

})

    
pdf(paste0('Fig.3a.sample_overview.pdf'), width = 20, height = 6)              
grid.arrange(p_combined_list[[1]], p_combined_list[[2]], p_combined_list[[3]],               
             p_combined_list[[4]],             
             ncol = 4, widths = c(65,45,95,28)   # add extra margin to first column               
)              
dev.off()

pdf( 'Fig.3a.sample_overview_legend.pdf', width = 10, height = 30)              
grid.arrange(p_combined_list[[1]],                   
             ncol = 1, widths = c(5) #,  # add extra margin to first column               
)  
dev.off()
