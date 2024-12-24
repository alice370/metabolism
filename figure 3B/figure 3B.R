data<-read.csv("209sample_188metabolites_auto_log10.csv",row.names = 1)
clin<-read.csv("RA_baseline_clin.csv",check.names = F)
library(ggplot2)
library(ggpubr)
library(eoffice)

p <- ggplot(clin, aes(x = GENDER, y = `DAS28-CRP`, fill = GENDER)) + 
  geom_violin(trim = FALSE, alpha = 0.8) +  
  geom_boxplot(width = 0.1, fill = "white") +  
  scale_fill_manual(values = c("#a1b8e1", "#df7f7f")) + 
  theme_classic2() + 
  theme(legend.position = "right") +  
  labs(x = "", y = "DAS28-CRP") +  
  stat_compare_means(method = "t.test", label = "p.format", 
                     comparisons = list(c("F", "M")), 
                     label.y = max(clin$`DAS28-CRP`, na.rm = TRUE) + 0.5)

print(p)
topptx(p,filename = "gender_DAS28_violin.pptx")