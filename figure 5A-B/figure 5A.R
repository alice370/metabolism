library(ggplot2)
library(eoffice)
library(ggrepel)
library(tidyr)
rm(list=ls())
data<-read.csv('209sample_188metabolites_auto_log10.csv',row.names = 1)
clin<-read.csv('sample_information.csv',check.names = F)
merge<-as.data.frame(t(data))
merge$sample<-row.names(merge)
merge1<-merge(clin[,c(1,3,4,27,28)],merge,by='sample')

MHY <- merge1[merge1$`response classification` %in% c('Good Response', 'Moderate Response') &
    merge1$`DRUG` == 'MTQ+HCQ',  -c(1:5)]
MHN<-merge1[merge1$`response classification` %in% c('No Response') &
                    merge1$`DRUG` == 'MTQ+HCQ',  -c(1:5)]
MLY<-merge1[merge1$`response classification` %in% c('Good Response', 'Moderate Response') &
              merge1$`DRUG` == 'MTQ+LEF',  -c(1:5)]
MLN<-merge1[merge1$`response classification` %in% c('No Response') &
              merge1$`DRUG` == 'MTQ+LEF',  -c(1:5)]
MHY_median<- apply(MHY,2,median,na.rm=T)
MHN_median<- apply(MHN,2,median,na.rm=T)
MLY_median<- apply(MLY,2,median,na.rm=T)
MLN_median<- apply(MLN,2,median,na.rm=T)
log10_FC_MH<-MHY_median-MHN_median
log10_FC_ML<-MLY_median-MLN_median
P_MH<-rep(NA,188)
P_ML<-rep(NA,188)

for(i in 1:188) try({
  P_MH[i] = wilcox.test(as.numeric(MHY[,i]), as.numeric(MHN[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_ML[i] = wilcox.test(as.numeric(MLY[,i]), as.numeric(MLN[,i]), alternative = "two.sided", paired = FALSE)$p.value
  })
test<-data.frame(colnames(merge1)[-c(1:5)],P_MH,P_ML,
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
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(size = 0.5),
    axis.line.y = element_line(size = 0.5), 
    axis.text = element_text(size = 12), 
    axis.ticks = element_line(linewidth = 0.5), 
    legend.position = "none" 
  ) +
  labs(
    x = "log10(Y/N)",
    y = "-log10(p value)",
    color = "",
    size=12
  )
MH.gg <- MH.gg + geom_text_repel(data=test[test$MH_sig!='no',],
  aes(x =`log10_FC_MH` , y = -log10(P_MH), label = metabolites),
  size = 4,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5  
)
MH.gg 
topptx(MH.gg,filename = "MH_YN_DEM.pptx")
####MTX+LEF####
ML.gg <- ggplot(test, aes(x =`log10_FC_ML` , y = -log10(P_ML), color = ML_sig)) +
  geom_point(size = 4) +
  scale_color_manual(values = c('up' = '#DE3024', 'down' = '#25108f', 'no' = 'grey'))+
  geom_vline(xintercept = log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_vline(xintercept = -log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(size = 0.5),
    axis.line.y = element_line(size = 0.5), 
    axis.text = element_text(size = 12), 
    axis.ticks = element_line(linewidth = 0.5), 
    legend.position = "none" 
  ) +
  labs(
    x = "log10(Y/N)",
    y = "-log10(p value)",
    color = "",
    size=12
  )
ML.gg <- ML.gg + geom_text_repel(data=test[test$ML_sig!='no',],
                                 aes(x =`log10_FC_ML` , y = -log10(P_ML), label = metabolites),
                                 size = 4,
                                 color = "black",
                                 box.padding = unit(0.4, "lines"),
                                 segment.color = "black",
                                 segment.size = 0.6,
                                 min.segment.length = 0.5  
)
ML.gg 
topptx(ML.gg,filename = "ML_YN_DEM.pptx")


