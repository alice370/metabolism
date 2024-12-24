rm(list = ls())
data<-read.csv("364sample_188metabolites_mean_log10.csv",row.names = 1)
clin<-read.csv('sample_information.csv')
merge1<-as.data.frame(t(data[,clin$sample[clin$Group=='At-risk of RA']]))
merge1$sample<-row.names(merge1)
colnames(clin)[1]<-'sample'
merge1<-merge(clin[,c(1,4)],merge1,by="sample")
result_age = data.frame(ID = colnames(merge1)[3:ncol(merge1)],
                        p_age_lm = sapply(merge1[,3:ncol(merge1)], function(x){summary(lm(merge1$Age~x))$coefficients[2,4]}),
                        beta_age_lm = sapply(merge1[,3:ncol(merge1)], function(x){t(coefficients(lm(merge1$Age~x)))[[2]]}))

result_age$age_lm=ifelse(result_age$p_age_lm<0.05 & result_age$beta_age_lm>0,"positive",
                         ifelse(result_age$p_age_lm<0.05 & result_age$beta_age_lm<0,"negative","n.s."))
write.csv(result_age,'IAR-age-meta-ln.csv')

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
    axis.ticks = element_line(size = 0.3),
    axis.ticks.length = unit(5, "pt"),
    legend.position = "none"  
  ) +
  labs(
    x = "Coefficient",
    y = "-log10(p value)",
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
topptx(gg,filename = "IAR_age_lm_plot.pptx")
######HC####
merge1<-as.data.frame(t(data[,clin$sample[clin$Group=='Health']]))
merge1$sample<-row.names(merge1)
colnames(clin)[1]<-'sample'
merge1<-merge(clin[,c(1,4)],merge1,by="sample")
result_age = data.frame(ID = colnames(merge1)[3:ncol(merge1)],
                        p_age_lm = sapply(merge1[,3:ncol(merge1)], function(x){summary(lm(merge1$Age~x))$coefficients[2,4]}),
                        beta_age_lm = sapply(merge1[,3:ncol(merge1)], function(x){t(coefficients(lm(merge1$Age~x)))[[2]]}))

result_age$age_lm=ifelse(result_age$p_age_lm<0.05 & result_age$beta_age_lm>0,"positive",
                         ifelse(result_age$p_age_lm<0.05 & result_age$beta_age_lm<0,"negative","n.s."))
write.csv(result_age,'HC-age-meta-ln.csv')

data1<-result_age

top_10_indices <- data1$p_age_lm<0.05
significant_points <- data.frame(
  beta.lm = data1$beta_age_lm[top_10_indices],
  p.value = -log10(data1$p_age_lm[top_10_indices]),
  metabolic = data1$ID[top_10_indices]
)
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
    axis.ticks = element_line(size = 0.3),
    axis.ticks.length = unit(5, "pt"),
   
    legend.position = "none"  
  ) +
  labs(
    x = "Coefficient",
    y = "-log10(p value)",
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
topptx(gg,filename = "HC_age_lm_plot.pptx")

A<-IAR$ID[IAR$age_lm=='positive']
B<-HC$ID[HC$age_lm=='positive']
C<-RA$ID[RA$age_lm=='positive']
intersect(A,C)
venn_list <- list(
  "RA" = A,
  "at risk individual" =  C
)
p1<-ggvenn(venn_list,show_percentage =FALSE,show_elements = FALSE,label_sep = ",",
           digits = 1,stroke_color = "white",
           fill_color = c("#DC00007F", "#F39B7F7F"),
           set_name_color = "black")
p1
topptx(p1,filename = "age-intersect.pptx")

