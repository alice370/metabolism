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
data_all <- read_excel("RA information.xlsx")
names(data_all)
data_all$`Smoke status` <- ifelse(data_all$`Smoke status` == 1, "yes",
                                  ifelse(data_all$`Smoke status` == 0, "no", NA))

data_all$Hypertension <- ifelse(data_all$Hypertension == 1, "yes",
                                ifelse(data_all$Hypertension == 0, "no", NA))

data_all$Hyperglycemia <- ifelse(data_all$Hyperglycemia == 1, "yes",
                                 ifelse(data_all$Hyperglycemia == 0, "no", NA))

data_all$Hyperlipidemia <- ifelse(data_all$Hyperlipidemia == 1, "yes",
                                  ifelse(data_all$Hyperlipidemia == 0, "no", NA))




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
data_all$DAS28_CRP <- 0.56 * sqrt(data_all$TJC) +
  0.28 * sqrt(data_all$SJC) +
  0.36 * log(data_all$CRP + 1) +
  0.014 * data_all$VAS +
  0.96

# 确保以下列名在你的数据框中存在：TJC28, SJC28, CRP, VAS, DAS28_CRP
data_all$share_TJC <- (0.56 * sqrt(data_all$TJC)) / (data_all$`DAS28_CRP` - 0.96)
data_all$share_SJC <- (0.28 * sqrt(data_all$SJC)) / (data_all$`DAS28_CRP` - 0.96)
data_all$share_CRP <- (0.36 * log(data_all$CRP + 1)) / (data_all$`DAS28_CRP` - 0.96)
data_all$share_VAS <- (0.014 * data_all$VAS) / (data_all$`DAS28_CRP` - 0.96)
write.csv(data_all,'data_all.csv')
table(data_all$`Disease Activity Class`)
DAS28_list <- c("share_VAS", "share_TJC", "share_SJC", "share_CRP")
Sub_list <- c("Clinical Remission", "Low disease activity", "moderate disease activity", "high disease activity")

order_samples <- c("Clinical Remission", "Low disease activity", "moderate disease activity", "high disease activity")

data_all$`Disease Activity Class` <- factor(data_all$`Disease Activity Class`, levels = order_samples)


p_combined_list<-list()
####plot####
p_combined_list <- lapply(Sub_list, function(x){              
  if(x == Sub_list[1]){              
    data <- data_all[data_all$`Disease Activity Class` == Sub_list[1],]              
  } else if (x == Sub_list[2]){              
    data <- data_all[data_all$`Disease Activity Class` == Sub_list[2],]              
  } else if (x == Sub_list[3]){              
    data <- data_all[data_all$`Disease Activity Class` == Sub_list[3],]              
  } else if (x == Sub_list[4]){              
    data <- data_all[data_all$`Disease Activity Class` == Sub_list[4],]              
  }         
  
  data <- data %>%
    arrange(`Disease Activity Class`,share_TJC, share_SJC, share_VAS, share_CRP)
  ordered_samples <- data$Sample
  data$Sample <- factor(data$Sample, levels = ordered_samples)
  ####disease activity####
  p_subtypes <- ggplot(data, aes(x = Sample, y = 1, fill = `Disease Activity Class`)) +               
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Disease activity',               
                      values = c("Clinical Remission" = "#FDBF6F","Low disease activity" ="#A6CEE3", "moderate disease activity"="#B2DF8A","high disease activity"="#FB9A99" ), drop = F) +               
    scale_y_continuous(expand = c(0,0)) +                
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Disease activity') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),               
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))        
  ####shared_DAS28####
  plot_data <- data %>%
    gather(key = share_DAS28, value = proportion, share_VAS, share_TJC, share_SJC, share_CRP) %>%
    mutate(share_DAS28 = factor(share_DAS28, levels = c("share_VAS", "share_TJC", "share_SJC", "share_CRP")))
  
  p_prop.subtype <- ggplot(plot_data, aes(x = Sample, y = proportion, fill = share_DAS28)) +
    geom_bar(stat = 'identity', position = position_stack(reverse = F)) +
    scale_fill_manual(name = '%DAS28-CRP Parameter', values = c("share_VAS" = "#EFC000FF", 
                                                                "share_TJC" = "#1F78B4", 
                                                                "share_SJC" = "#B2DF8A", 
                                                                "share_CRP" = "#A73030FF"), drop = F) +
    scale_x_discrete(drop = F) +
    scale_y_continuous(expand = c(0,0)) + ylab('% Proportion') +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18), panel.grid = element_blank()) 
  
  ####age####
  p_age <- ggplot(data, aes(x = Sample, y = 0.6 , fill = Age)) +              
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
  
  
  ####sex####
  p_sex <- ggplot(data, aes(x = Sample, y = 0.6, fill = Gender)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Sex',              
                      values = c("#FB9A99","#A6CEE3"), na.value = 'white', drop = F) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Sex') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  ####ACPA####             
  data$`ACPA status` <- ifelse(data$`ACPA status` == "NA", NA, data$`ACPA status`)
  p_ACPA <- ggplot(data, aes(x = Sample, y = 0.6, fill = `ACPA status`)) +              
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
  ####smoke####
  p_smoke <- ggplot(data, aes(x = Sample, y = 0.6, fill = data$`Smoke status`)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Smoke status',              
                      values = c("#B2DF8A","#FB9A99"), na.value = 'white', drop = F,na.translate = FALSE) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Smoke Status') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  
  
  ####hypertension####
  p_hypertension <- ggplot(data, aes(x = Sample, y = 0.6, fill = data$Hypertension)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Hypertension',              
                      values = c("yes"="#D8A7CA", "no"="#9CBACF"), na.value = 'white', drop = F) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Hypertension') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))   
  ####hyperglucose####
  p_Hyperglycemia <- ggplot(data, aes(x = Sample, y = 0.6, fill = data$Hyperglycemia)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Hyperglycemia',              
                      values = c("yes"="#D6BFAF", "no"="#A8CBB7"), na.value = 'white', drop = F) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Hyperglycemia') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18)) 
  ####hyperlipidemia####
  p_Hyperlipidemia <- ggplot(data, aes(x = Sample, y = 0.6, fill = data$Hyperlipidemia)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Hyperlipidemia',              
                      values = c("yes"="#F9E79F", "no"="#A8C686"), na.value = 'white', drop = F) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Hyperlipidemia') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18)) 
  ####BMI####
  p_BMI <- ggplot(data, aes(x = Sample, y = 0.6 , fill = data$`Body Mass Index`)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_gradient2(name='BMI',low = "darkgreen", high = "red", mid = "yellow",midpoint=22,limits = c(14,47), na.value = 'white') +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('BMI') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))  
  
  #exclude legends from original plots              
  p_subtypes <- p_subtypes + theme(legend.position = 'none')
  p_prop.subtype <- p_prop.subtype + theme(legend.position = 'none')
  p_age <- p_age + theme(legend.position = 'none')
  p_sex <- p_sex + theme(legend.position = 'none')
  p_ACPA <- p_ACPA + theme(legend.position = 'none')
  p_smoke <- p_smoke + theme(legend.position = 'none')
  p_hypertension <- p_hypertension + theme(legend.position = 'none')
  p_Hyperglycemia <- p_Hyperglycemia + theme(legend.position = 'none')
  p_Hyperlipidemia <- p_Hyperlipidemia + theme(legend.position = 'none')
  p_BMI <- p_BMI + theme(legend.position = 'none')         
  
  
  
  #delete y-labels except the one on the very left              
  if(x != Sub_list[1]){ 
    p_subtypes <- p_subtypes + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_prop.subtype <- p_prop.subtype + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_smoke <- p_smoke + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_age <- p_age + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_sex <- p_sex + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_ACPA <- p_ACPA + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_hypertension <- p_hypertension + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_Hyperglycemia <- p_Hyperglycemia + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_Hyperlipidemia <- p_Hyperlipidemia + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_BMI <- p_BMI + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }              
  
  
  #####combine both
  p_subtypes <- ggplotGrob(p_subtypes + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_prop.subtype <- ggplotGrob(p_prop.subtype + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_smoke <- ggplotGrob(p_smoke + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_age <- ggplotGrob(p_age + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_sex <- ggplotGrob(p_sex + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_ACPA <- ggplotGrob(p_ACPA + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_hypertension <- ggplotGrob(p_hypertension + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_Hyperglycemia <- ggplotGrob(p_Hyperglycemia + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_Hyperlipidemia <- ggplotGrob(p_Hyperlipidemia + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_BMI <- ggplotGrob(p_BMI + theme(plot.margin = unit(c(0.1, 0.1, 1, 0.1), "cm")))          
  
  
  
  grid_unit_pmax <- grid::unit.pmax(p_subtypes$widths, p_prop.subtype$widths, p_smoke$widths,
                                    p_age$widths, p_sex$widths, p_ACPA$widths, p_hypertension$widths, 
                                    p_Hyperglycemia$widths, p_Hyperlipidemia$widths, p_BMI$widths)
  
  p_subtypes$widths <- grid_unit_pmax
  p_prop.subtype$widths <- grid_unit_pmax
  p_smoke$widths <- grid_unit_pmax
  p_age$widths <- grid_unit_pmax
  p_sex$widths <- grid_unit_pmax
  p_ACPA$widths <- grid_unit_pmax
  p_hypertension$widths <- grid_unit_pmax
  p_Hyperglycemia$widths <- grid_unit_pmax
  p_Hyperlipidemia$widths <- grid_unit_pmax
  p_BMI$widths <- grid_unit_pmax
  
  
  g <- gtable_rbind(p_subtypes,
                    p_prop.subtype,
                    p_ACPA,
                    p_age,
                    p_sex,
                    p_BMI,
                    p_smoke,
                    p_Hyperglycemia,
                    p_hypertension,
                    p_Hyperlipidemia)
  
  
  id_panels_h <- unique(g$layout[g$layout$name == "panel", "t"])
  
  g$heights[id_panels_h] <- grid::unit(c(0.6, 2,0.5, 0.5, 0.5,0.5, 0.5, 0.5, 0.5, 0.5), "null")
  
  g_combined <- arrangeGrob(grobs = list(g, aligend.legends1),               
                            nrow = 1,              
                            ncol = 2, #3,              
                            layout_matrix = matrix(c(rep(1, 4), rep(2,4)), nrow = 1))
  
  g_combined<-g
  return(g_combined)              
  
})


pdf(paste0('Fig.3A.pdf'), width = 20, height = 6)              
grid.arrange(p_combined_list[[1]], p_combined_list[[2]], p_combined_list[[3]],               
             p_combined_list[[4]],             
             ncol = 4, widths = c(65,45,95,28)   # add extra margin to first column               
)              
dev.off()


####for legend####
  data <- data_all %>%
    arrange(`Disease Activity Class`,share_TJC, share_SJC, share_VAS, share_CRP)
  ordered_samples <- data$Sample
  data$Sample <- factor(data$Sample, levels = ordered_samples)
####disease activity####
  p_subtypes <- ggplot(data, aes(x = Sample, y = 1, fill = `Disease Activity Class`)) +               
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Disease activity',               
                      values = c("Clinical Remission" = "#FDBF6F","Low disease activity" ="#A6CEE3", "moderate disease activity"="#B2DF8A","high disease activity"="#FB9A99" ), drop = F) +               
    scale_y_continuous(expand = c(0,0)) +                
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Disease activity') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),               
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))        
  ####shared_DAS28####
  plot_data <- data %>%
    gather(key = share_DAS28, value = proportion, share_VAS, share_TJC, share_SJC, share_CRP) %>%
    mutate(share_DAS28 = factor(share_DAS28, levels = c("share_VAS", "share_TJC", "share_SJC", "share_CRP")))
  
  p_prop.subtype <- ggplot(plot_data, aes(x = Sample, y = proportion, fill = share_DAS28)) +
    geom_bar(stat = 'identity', position = position_stack(reverse = F)) +
    scale_fill_manual(name = '%DAS28-CRP Parameter', values = c("share_VAS" = "#EFC000FF", 
                                                               "share_TJC" = "#1F78B4", 
                                                               "share_SJC" = "#B2DF8A", 
                                                               "share_CRP" = "#A73030FF"), drop = F) +
    scale_x_discrete(drop = F) +
    scale_y_continuous(expand = c(0,0)) + ylab('% Proportion') +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18), panel.grid = element_blank()) 
  
  ####age####
  p_age <- ggplot(data, aes(x = Sample, y = 0.6 , fill = Age)) +              
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
  
  
  ####sex####
  p_sex <- ggplot(data, aes(x = Sample, y = 0.6, fill = Gender)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Sex',              
                      values = c("#FB9A99","#A6CEE3"), na.value = 'white', drop = F) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Sex') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  ####ACPA####             
  data$`ACPA status` <- ifelse(data$`ACPA status` == "NA", NA, data$`ACPA status`)
  p_ACPA <- ggplot(data, aes(x = Sample, y = 0.6, fill = `ACPA status`)) +              
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
  ####smoke####
  p_smoke <- ggplot(data, aes(x = Sample, y = 0.6, fill = data$`Smoke status`)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Smoke status',              
                      values = c("#B2DF8A","#FB9A99"), na.value = 'white', drop = F,na.translate = FALSE) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Smoke Status') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  

  ####hypertension####
  p_hypertension <- ggplot(data, aes(x = Sample, y = 0.6, fill = data$Hypertension)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Hypertension',              
                      values = c("yes"="#D8A7CA", "no"="#9CBACF"), na.value = 'white', drop = F) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Hypertension') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))   
  ####hyperglucose####
  p_Hyperglycemia <- ggplot(data, aes(x = Sample, y = 0.6, fill = data$Hyperglycemia)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Hypertension',              
                      values = c("yes"="#D6BFAF", "no"="#A8CBB7"), na.value = 'white', drop = F) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Hypertension') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18)) 
  ####hyperlipidemia####
  p_Hyperlipidemia <- ggplot(data, aes(x = Sample, y = 0.6, fill = data$Hyperlipidemia)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_manual(name = 'Hypertension',              
                      values = c("yes"="#F9E79F", "no"="#A8C686"), na.value = 'white', drop = F) +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('Hypertension') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18)) 
  ####BMI####
  p_BMI <- ggplot(data, aes(x = Sample, y = 0.6 , fill = data$`Body Mass Index`)) +              
    geom_tile(colour = 'white', size = gab_size) +              
    scale_fill_gradient2(name='BMI',low = "darkgreen", high = "red", mid = "yellow",midpoint=22,limits = c(14,47), na.value = 'white') +              
    scale_y_continuous(expand = c(0,0)) +              
    scale_x_discrete(expand = c(0,0), drop = F) +              
    ylab('BMI') +              
    theme_bw() +              
    theme(axis.text = element_blank(), axis.ticks = element_blank(),              
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),              
          text = element_text(size = 20),              
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))  

  legends <- list(g_legend(p_subtypes),       
                  g_legend(p_prop.subtype),   
                  g_legend(p_age),          
                  g_legend(p_sex),            
                  g_legend(p_ACPA),           
                  g_legend(p_smoke),            
                  g_legend(p_hypertension),         
                  g_legend(p_Hyperglycemia),
                  g_legend(p_Hyperlipidemia),
                  g_legend(p_BMI))                     
  
  
  
  aligend.legends1 <- align.legends.fun(legends[c(1:10)]) #               
  
  legend_row <- arrangeGrob(grobs = legends,
                            nrow = 1)  # 一行多个图例


pdf( 'Fig.3A legend.pdf', width = 30, height = 10)              
grid.draw(legend_row)
dev.off()
