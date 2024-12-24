rm(list=ls())
library("dplyr")
library(tidyr)
library(ggplot2)
library(eoffice)
library(ggpubr)
data<-read.csv("209sample_188metabolites_auto_log10.csv",row.names = 1)
clin<-read.csv('RA_baseline_clin.csv',check.names = F)
sig<-read.csv('significance_count.csv')
test<-read.csv('M-U-ACPA.csv',row.names = 1)
meta<-intersect(sig$Metabolite[sig$Significant_Count>25],test$metabolic[test$p.value<0.05])
merge1 <- as.data.frame(t(data[meta, clin$sample[clin$ACPA=='positive'|clin$ACPA=='negative']]))

merge1$sample<-row.names(merge1)
merge1<-merge(clin[,c(1,6)],merge1,by="sample")

row.names(merge1)<-merge1$sample
merge1<-merge1[,-1]
long_data <- merge1 %>%
  pivot_longer(-ACPA, names_to = "metabolites", values_to = "Value")

p <- ggplot(long_data, aes(x = ACPA, y = Value, fill = ACPA)) + 
  geom_violin(trim = TRUE, alpha = 0.8) +  
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("#a1b8e1", "#df7f7f")) +  
  theme_classic2() + 
  theme(legend.position = "none",
        strip.background = element_blank()) +  
  labs(x = "", y = "") +  
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("positive", "negative")),  
                     label.y = max(long_data$Value, na.rm = TRUE) + 0.5) +  
  facet_wrap(~metabolites,nrow = 1,scales = "free_y")  
topptx(p,'ACPA-4meta-violin.pptx')
