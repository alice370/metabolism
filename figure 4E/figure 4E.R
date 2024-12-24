library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(eoffice)
rm(list=ls())
data<-read.csv('364sample-188metabolites-log10-auto.csv',row.names = 1)
clin<-read.csv('sample_information.csv',check.names = F)
activity_test<-read.csv('disease activity class test.csv',row.names = 1)
clinical_test<-read.csv('metabolite_comparison_results_scaled.csv')
activity_sig<-activity_test$metabolites[which(activity_test$P_DL<0.05|activity_test$P_LM<0.05|activity_test$P_MH<0.05|activity_test$P_LH<0.05|activity_test$P_DM<0.05|activity_test$P_DH<0.05)]
clinical_sig<-clinical_test$Metabolite[clinical_test$PValue_RA_HC<0.05]
con<-intersect(activity_sig,clinical_sig)
colnames(clinical_test)[1]<-colnames(activity_test)[1]
con_test<-merge(clinical_test[clinical_test$metabolites%in%con,c(1:3)],
                activity_test[activity_test$metabolites%in%con,c(1,3:8,9:14)],by="metabolites")

con_up <- con_test[con_test[,2] > 0 & 
                     ((con_test[,4] < 0.05 & con_test[,10] < 0) | 
                        (con_test[,5] < 0.05 & con_test[,11] < 0) | 
                        (con_test[,6] < 0.05 & con_test[,12] < 0) | 
                        (con_test[,7] < 0.05 & con_test[,13] < 0) | 
                        (con_test[,8] < 0.05 & con_test[,14] < 0) | 
                        (con_test[,9] < 0.05 & con_test[,15] < 0)), ]
con_down <- con_test[con_test[,2] < 0 & 
                       ((con_test[,4] < 0.05 & con_test[,10] > 0) | 
                          (con_test[,5] < 0.05 & con_test[,11] >0) | 
                          (con_test[,6] < 0.05 & con_test[,12] > 0) | 
                          (con_test[,7] < 0.05 & con_test[,13] > 0) | 
                          (con_test[,8] < 0.05 & con_test[,14] >0) | 
                          (con_test[,9] < 0.05 & con_test[,15] >0)), ]

write.csv(con_test,'con_test_clin_activity.csv')

con<-rbind(con_down[1:4,],con_up[c(1,2,3,5),])

R<-data[con$metabolites,clin$Sample[clin$`disease activity class` == "Clinical Remission"]]
L<-data[con$metabolites,clin$Sample[clin$`disease activity class` == " Low disease activity"]]
M<-data[con$metabolites,clin$Sample[clin$`disease activity class`== "moderate disease activity"]]
H<-data[con$metabolites,clin$Sample[clin$`disease activity class` == "high disease activity"]]
RA<-data[con$metabolites,clin$Sample[clin$Group == "RA"]]
HC<-data[con$metabolites,clin$Sample[clin$Group == "Health"]]


group_labels <- rep(c("RA","HC","D","L", "M", "H"), 
                    times = c(length(RA), length(HC),length(R), length(L), length(M),length(H)))

filtered_data<-cbind(RA,HC,R,L,M,H)
filtered_data$name<-row.names(filtered_data)
long_data <- filtered_data %>%
  pivot_longer(-name, names_to = "sample", values_to = "Value")


long_data$group<- rep(group_labels,8)
long_data$group <- factor(long_data$group, levels = c('HC', 'RA', 'D','L', 'M', 'H'))



get_signif_data <- function(con, meta, plot, step) {
  
  signif_data <- con %>%
    filter(metabolites == meta) %>%
    
    select(starts_with("P_")) %>%
  
    pivot_longer(cols = everything(), names_to = "comparison", values_to = "p.value") %>%

    filter(p.value < 0.05) %>%

    mutate(
      # 假设列名格式为 "P_Group1Group2"，如 "P_MH"
      group1 = str_sub(comparison, 3, 3),      
      group2 = str_sub(comparison, 4, 4),    
    
      label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    select(-comparison)
  
  

  max_values <- max(plot$Value[plot$group%in%c('L','M','D','H')])

  signif_data <- signif_data %>%
    rowwise() %>%
    mutate(
      y.start = max_values + step,
      y.end = y.start + step
    ) %>%
    ungroup()
  
  signif_data <- signif_data %>%
    arrange(y.start) %>%
    group_by(y.start) %>%
    mutate(y.position = y.end + (row_number() - 1) * step) %>%
    ungroup()
  
  return(signif_data)
}
max_values <- max(plot$Value[plot$group%in%c('L','M','D','H')])
####循环绘图####
for (n in 1:nrow(con)) {
  meta<-con$metabolites[n]
  plot<-long_data[long_data$name==meta,]
  sig_data <- get_signif_data(con, meta, plot, step = 0.5)
  p <- ggplot(plot, aes(x = group, y = Value, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.8, size = 0.2, alpha=0.8) +
    geom_jitter(aes(color=group), position = position_jitterdodge(1), size=0.8, alpha=0.4) +
    scale_color_manual(values = c(
      "HC" = "#4E79A7", "RA" = "#E15759", "M" = "#59A14F",
      "H" = "#F28E2B", "L" = "#76B7B2", "D" = "#9C755F")) +
    scale_fill_manual(values = c("HC" = "#A0CBE8", "M" = "#8CD17D", "RA" = "#FF9D9A",
                                 "H" = "#FFBE7D", "L" = "#86BCB6", "D" = "#D4A76A")) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill = 'transparent', color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.2), 
      axis.ticks = element_line(colour = "black", linewidth = 0.2), 
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10, color = "black")
    )+
    ylab('') +
    labs(title = meta) +
    ylim(NA, max_values + (nrow(sig_data) + 2) * 0.5)  
  
  line_to_star_gap <- 0.05  
  
  for (i in 1:nrow(sig_data)) {
    group1 <- sig_data$group1[i]
    group2 <- sig_data$group2[i]
    
    y1 <- sig_data$y.position[i] - 0.5
    y2 <- y1 + line_to_star_gap  
    
    label <- sig_data$label[i]
    
    x1 <- which(levels(plot$group) == group1)
    x2 <- which(levels(plot$group) == group2)

    dt <- data.frame(x = (x1 + x2) / 2, y = y2, label = label)

    p <- p +
      geom_segment(x = x1, xend = x2, y = y1, yend = y1, color = "black",linewidth=0.2) +
      geom_text(data = dt, 
                aes(x = x, y = y, label = label),
                size =4 , 
                inherit.aes = FALSE, 
                show.legend = FALSE)
  }
  max_RH <- max(plot$Value[plot$group%in%c('RA','HC')])
  p_RH<-con$PValue_RA_HC[con$metabolites==meta]
  label_RH<-case_when(
    p_RH < 0.001 ~ "***",
    p_RH < 0.01 ~ "**",
    p_RH < 0.05 ~ "*",
    TRUE ~ ""
  )
  dt2 <- data.frame(x = 1.5, y = max_RH+0.8, label = label_RH)
  p<- p +
    geom_segment(x = 1, xend = 2, y = max_RH+0.5, yend = max_RH+0.5, color = "black",linewidth=0.2) +
    geom_text(data = dt2, 
              aes(x = x, y = y, label = label),
              size = 4,  
              inherit.aes = FALSE, 
              show.legend = FALSE)
  assign(paste0("p",n),p)
}
a<-ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2,ncol =4,legend = NULL)
a
cairo_pdf("plot.pdf", width = 210/25.4, height = 100/25.4)
print(a)
dev.off()
