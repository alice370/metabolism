library(ggplot2)
library(eoffice)
library(ggrepel)
rm(list=ls())
data<- read.csv('209sample_188metabolites_auto_log10.csv',row.names = 1,check.names = F)
clin<-read.csv("RA_baseline_clin.csv",check.names = F)
merge<-as.data.frame(t(data))
merge$sample<-row.names(merge)
merge1<-merge(clin[,c(1:2,3,6:8)],merge,by="sample")
result = data.frame(
  ID = colnames(merge1)[7:ncol(merge1)],
  p_DAS28CRP_lm = sapply(merge1[,7:ncol(merge1)], function(x) {
    summary(lm(merge1$`DAS28-CRP` ~ x))$coefficients[2,4]  # p-value for the predictor variable x
  }),
  beta_DAS28CRP_lm = sapply(merge1[,7:ncol(merge1)], function(x) {
    coefficients(lm(merge1$`DAS28-CRP` ~ x))[2]  # beta coefficient for the predictor variable x
  })
)

result$DAS28CRP_lm=ifelse(result$p_DAS28CRP_lm<0.05 & result$beta_DAS28CRP_lm>0,"positive",
                          ifelse(result$p_DAS28CRP_lm<0.05 & result$beta_DAS28CRP_lm<0,"negative","n.s."))
write.csv(result,'auto_log10_lm_DAS28.csv')                          
table(result$DAS28CRP_lm)


# Assuming 'result' is your data frame
# Select rows where DAS28CRP_lm is either 'positive' or 'negative'
positive <- result[order(result$beta_DAS28CRP_lm,decreasing = T)[1:2],]
# For the 'negative' subset, select the top five rows with the smallest values in the second column
negative <- result[result$DAS28CRP_lm == 'negative', ]
top_five_negative <- negative[order(negative$p_DAS28CRP_lm)[1:5], ]
significant_points<-rbind(positive,top_five_negative)


gg<-ggplot(result, aes(x = result$beta_DAS28CRP_lm, y = -log10(p_DAS28CRP_lm), color = result$beta_DAS28CRP_lm, size = -log10(p_DAS28CRP_lm))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  geom_point() +
  scale_color_gradient2(low = "darkblue", mid='grey', high = "darkred") + 
  scale_size_continuous(range = c(0, 4)) + 
  theme(panel.background = element_rect(fill = 'transparent', color = 'black')) +
  labs(
    x = "beta-lm-DAS28-CRP",
    y = "-log10(FDR)"
  )+ theme_classic()+
  theme(axis.title = element_text(size = 8),
           axis.text = element_text(size = 8),
           legend.position = 'none')

gg
gg <- gg + geom_text_repel(
  data = significant_points,
  aes(x = beta_DAS28CRP_lm, y = -log10(p_DAS28CRP_lm), label = ID),
  size = 5,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5,
  max.overlaps = 5
)

gg
topptx(gg,filename = "lm-DAS28CRP-auto.pptx")