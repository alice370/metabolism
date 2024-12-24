rm(list=ls())
data<-read.csv("209sample_188metabolites_auto_log10.csv",row.names = 1)
clin<-read.csv('RA_baseline_clin.csv',check.names = F)
merge1<-as.data.frame(t(data))
merge1$sample<-row.names(merge1)
merge1<-merge(clin[,c(1,3)],merge1,by="sample")
result_age = data.frame(ID = colnames(merge1)[3:ncol(merge1)],
                        p_age_lm = sapply(merge1[,3:ncol(merge1)], function(x){summary(lm(merge1$age~x))$coefficients[2,4]}),
                        beta_age_lm = sapply(merge1[,3:ncol(merge1)], function(x){t(coefficients(lm(merge1$age~x)))[[2]]}))

result_age$age_lm=ifelse(result_age$p_age_lm<0.05 & result_age$beta_age_lm>0,"positive",
                         ifelse(result_age$p_age_lm<0.05 & result_age$beta_age_lm<0,"negative","n.s."))
write.csv(result_age,'age-meta-ln.csv')

data1<-result_age

top_10_indices <- order(data1$p_age_lm)[1:10]
significant_points <- data.frame(
  beta.lm = data1$beta_age_lm[top_10_indices],
  p.value = -log10(data1$p_age_lm[top_10_indices]),
  metabolic = data1$ID[top_10_indices]
)
library(ggplot2)
library(ggrepel)
library(eoffice)
gg <- ggplot(data1, aes(x = beta_age_lm, y = -log10(data1$p_age_lm), color = data1$age_lm)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('positive' = '#DE3024', 'down' = 'negative', 'n.s.' = 'grey')) + 
  geom_vline(xintercept = 0, linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(size = 0.5), 
    axis.line.y = element_line(size = 0.5),
    axis.text = element_text(size = 6), 
    axis.ticks = element_line(size = 1),

    legend.position = "none" 
  ) +
  labs(
    x = "beta lm",
    y = "-log10(p.value)",
    color = "",
    size=8
  )
gg

gg <- gg + geom_text_repel(
  data = significant_points,
  aes(x = beta.lm, y = p.value, label = metabolic),
  size = 2,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5 
)
gg
topptx(gg,filename = "age_lm_plot.pptx")
