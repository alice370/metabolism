library(readr)
library(readxl)
rm(list=ls())
origin <- read.csv("364sample_185metabolites.csv",check.names = F,row.names = 1)
data <- read_csv("364sample_185metabolites_mean_log10.csv")
clin <- read_excel("RA information.xlsx")
table(clin$`Response Classification`)
MHY <- clin$Sample[clin$`Response Classification` %in% c('Response') &
                     clin$Drug == 'MTX+HCQ']
MHN<- clin$Sample[clin$`Response Classification` %in% c('No Response') &
                    clin$Drug == 'MTX+HCQ']
MLY<- clin$Sample[clin$`Response Classification` %in% c('Response') &
                    clin$Drug == 'MTX+LEF']
MLN<-clin$Sample[clin$`Response Classification` %in% c('No Response') &
                   clin$Drug == 'MTX+LEF']


result_response <- data.frame()


for (metabolite in row.names(origin)) {
  group_MHY <- as.numeric(origin[metabolite, MHY])
  group_MHN <- as.numeric(origin[metabolite, MHN])
  group_MLY <- as.numeric(origin[metabolite, MLY])
  group_MLN <- as.numeric(origin[metabolite, MLN])
  
  # 统计非NA数量
  nonNA_MHY <- sum(!is.na(group_MHY))
  nonNA_MHN <- sum(!is.na(group_MHN))
  nonNA_MLY <- sum(!is.na(group_MLY))
  nonNA_MLN <- sum(!is.na(group_MLN))
  
  # NA数量
  NA_MHY <- sum(is.na(group_MHY))
  NA_MHN <- sum(is.na(group_MHN))
  NA_MLY <- sum(is.na(group_MLY))
  NA_MLN <- sum(is.na(group_MLN))
  
  # MHY vs MHN 比较
  contingency_table_MH <- matrix(c(nonNA_MHY, NA_MHY,
                                   nonNA_MHN, NA_MHN),
                                 nrow = 2,
                                 dimnames = list(Status = c("NonNA", "NA"), 
                                                 Group = c("MHY", "MHN")))
  
  if (any(contingency_table_MH < 5)) {
    test_result_MH <- fisher.test(contingency_table_MH)
    test_used_MH <- "Fisher"
  } else {
    test_result_MH <- chisq.test(contingency_table_MH)
    test_used_MH <- "Chi-square"
  }
  
  # MLY vs MLN 比较
  contingency_table_ML <- matrix(c(nonNA_MLY, NA_MLY,
                                   nonNA_MLN, NA_MLN),
                                 nrow = 2,
                                 dimnames = list(Status = c("NonNA", "NA"), 
                                                 Group = c("MLY", "MLN")))
  
  if (any(contingency_table_ML < 5)) {
    test_result_ML <- fisher.test(contingency_table_ML)
    test_used_ML <- "Fisher"
  } else {
    test_result_ML <- chisq.test(contingency_table_ML)
    test_used_ML <- "Chi-square"
  }
  
  # 存储结果
  result_row <- data.frame(
    metabolic = metabolite,
    nonNA_MHY = nonNA_MHY,
    detected_MHY = nonNA_MHY / length(group_MHY),
    nonNA_MHN = nonNA_MHN,
    detected_MHN = nonNA_MHN / length(group_MHN),
    nonNA_MLY = nonNA_MLY,
    detected_MLY = nonNA_MLY / length(group_MLY),
    nonNA_MLN = nonNA_MLN,
    detected_MLN = nonNA_MLN / length(group_MLN),
    MH_pvalue = test_result_MH$p.value,
    MH_test_used = test_used_MH,
    ML_pvalue = test_result_ML$p.value,
    ML_test_used = test_used_ML
  )
  
  result_response <- rbind(result_response, result_row)
}

write.csv(result_response, 'response_new卡方频率.csv', row.names = FALSE)

library(ggplot2)
library(eoffice)
library(ggrepel)
library(tidyr)
rm(list=ls())
data <- read_csv("364sample_185metabolites_mean_log10.csv")
data<-as.data.frame(data)
row.names(data)<-data$metabolites
data<-data[,-1]
clin <- read_excel("RA information.xlsx")
merge<-as.data.frame(t(data))
merge$Sample<-row.names(merge)
merge1<-merge(clin[,c(1,20,21)],merge,by='Sample')

MHY <- clin$Sample[clin$`Response Classification` %in% c('Response') &
                     clin$Drug == 'MTX+HCQ']
MHN<- clin$Sample[clin$`Response Classification` %in% c('No Response') &
                    clin$Drug == 'MTX+HCQ']
MLY<- clin$Sample[clin$`Response Classification` %in% c('Response') &
                    clin$Drug == 'MTX+LEF']
MLN<-clin$Sample[clin$`Response Classification` %in% c('No Response') &
                   clin$Drug == 'MTX+LEF']
MHY<-merge1[merge1$Sample%in%MHY,-c(1:3)]
MHN<-merge1[merge1$Sample%in%MHN,-c(1:3)]
MLY<-merge1[merge1$Sample%in%MLY,-c(1:3)]
MLN<-merge1[merge1$Sample%in%MLN,-c(1:3)]

MHY_median<- apply(MHY,2,median,na.rm=T)
MHN_median<- apply(MHN,2,median,na.rm=T)
MLY_median<- apply(MLY,2,median,na.rm=T)
MLN_median<- apply(MLN,2,median,na.rm=T)
log10_FC_MH<-MHY_median-MHN_median
log10_FC_ML<-MLY_median-MLN_median
P_MH<-rep(NA,185)
P_ML<-rep(NA,185)

for(i in 1:185) try({
  P_MH[i] = wilcox.test(as.numeric(MHY[,i]), as.numeric(MHN[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_ML[i] = wilcox.test(as.numeric(MLY[,i]), as.numeric(MLN[,i]), alternative = "two.sided", paired = FALSE)$p.value
})
test<-data.frame(colnames(merge1)[-c(1:3)],P_MH,P_ML,
                 log10_FC_MH,log10_FC_ML,MHY_median,MHN_median,MLY_median,MLN_median)
colnames(test)[1]<-'metabolites'
write.csv(test,"response test.csv")
test$MH_sig<-'no'
test$MH_sig[test$P_MH<0.05&test$log10_FC_MH>0]<-'up'
test$MH_sig[test$P_MH<0.05&test$log10_FC_MH<0]<-'down'
table(test$MH_sig)
test$ML_sig<-'no'
test$ML_sig[test$P_ML<0.05&test$log10_FC_ML>0]<-'up'
test$ML_sig[test$P_ML<0.05&test$log10_FC_ML<0]<-'down'
table(test$ML_sig)

####MTX+HCQ####
MH.gg <- ggplot(test, aes(x =`log10_FC_MH` , y = -log10(P_MH), color = MH_sig)) +
  geom_point(size = 4) +
  scale_color_manual(values = c('up' = '#DE3024', 'down' = '#25108f', 'no' = 'grey'))+
  geom_vline(xintercept = log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_vline(xintercept = -log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  theme(
    plot.margin = margin(),
    panel.background   = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(size = 0.7, color = "black"),
    axis.line.y        = element_line(size = 0.7, color = "black"),
    axis.text          = element_text(size = 14, color = "black"),   # 刻度文字
    axis.title         = element_text(size = 15, color = "black"),   # 坐标轴标题
    axis.ticks         = element_line(size = 0.7, color = "black"),
    legend.position    = "none")+
  labs(
    x = "Log10(responder/non-responder)",
    y = "-Log10(p value)",
    color = "black",
    size=15
  )
MH.gg <- MH.gg + geom_text_repel(data=test[test$MH_sig!='no',],
                                 aes(x =`log10_FC_MH` , y = -log10(P_MH), label = metabolites),
                                 size = 5,
                                 color = "black",
                                 box.padding = unit(0.4, "lines"),
                                 segment.color = "black",
                                 segment.size = 0.6,
                                 min.segment.length = 0.5  
)
MH.gg 
ggsave("MH_YN_DEM.pdf", width = 5, height = 5) 
####MTX+LEF####
ML.gg <- ggplot(test, aes(x =`log10_FC_ML` , y = -log10(P_ML), color = ML_sig)) +
  geom_point(size = 4) +
  scale_color_manual(values = c('up' = '#DE3024', 'down' = '#25108f', 'no' = 'grey'))+
  geom_vline(xintercept = log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_vline(xintercept = -log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  theme(
    plot.margin = margin(),
    panel.background   = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(size = 0.7, color = "black"),
    axis.line.y        = element_line(size = 0.7, color = "black"),
    axis.text          = element_text(size = 14, color = "black"),   # 刻度文字
    axis.title         = element_text(size = 15, color = "black"),   # 坐标轴标题
    axis.ticks         = element_line(size = 0.7, color = "black"),
    legend.position    = "none")+
  labs(
    x = "Log10(responder/non-responder)",
    y = "-Log10(p value)",
    color = "black",
    size=15
  )
ML.gg <- ML.gg + geom_text_repel(data=test[test$ML_sig!='no',],
                                 aes(x =`log10_FC_ML` , y = -log10(P_ML), label = metabolites),
                                 size = 5,
                                 color = "black",
                                 box.padding = unit(0.4, "lines"),
                                 segment.color = "black",
                                 segment.size = 0.6,
                                 min.segment.length = 0.5  
)
ML.gg 
ggsave("ML_YN_DEM.pdf", width = 5, height = 5) 


