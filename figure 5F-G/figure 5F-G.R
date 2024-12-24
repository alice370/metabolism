library(ggrepel)
library(ggplot2)
library(ggpubr)
library(eoffice)
rm(list = ls())
data<-read.csv('561sample_188meta_log10_auto.csv',row.names = 1)
clin<-read.csv('sample_information.csv',check.names = F)

data <- as.data.frame(t(data))
MLYA<-clin$sample[which(clin$DRUG =="MTQ+LEF"&
                          clin$`response classification` %in% c('Good Response', 'Moderate Response'))]
MLYA<-data[MLYA,]
MLYB<-clin$`follow up-sample`[which(clin$DRUG =="MTQ+LEF"&
                                      clin$`response classification` %in% c('Good Response', 'Moderate Response'))]
MLYB<-data[MLYB,]
MLNA<-clin$sample[which(clin$DRUG =="MTQ+LEF"&
                          clin$`response classification` %in% c('No Response'))]
MLNA<-data[MLNA,]
MLNB<-clin$`follow up-sample`[which(clin$DRUG =="MTQ+LEF"&
                                      clin$`response classification` %in% c('No Response'))]
MLNB<-data[MLNB,]
MHYA<-clin$sample[which(clin$DRUG =="MTQ+HCQ"&
                          clin$`response classification` %in% c('Good Response', 'Moderate Response'))]
MHYA<-data[MHYA,]
MHYB<-clin$`follow up-sample`[which(clin$DRUG =="MTQ+HCQ"&
                                      clin$`response classification` %in% c('Good Response', 'Moderate Response'))]
MHYB<-data[MHYB,]
MHNA<-clin$sample[which(clin$DRUG =="MTQ+HCQ"&
                          clin$`response classification` %in% c('No Response'))]
MHNA<-data[MHNA,]
MHNB<-clin$`follow up-sample`[which(clin$DRUG =="MTQ+HCQ"&
                                      clin$`response classification` %in% c('No Response'))]
MHNB<-data[MHNB,]

HC<-clin$sample[clin$Group=='Health']
HC<-data[HC,]
MHYA_median<- apply(MHYA,2,median,na.rm=T)
MHYB_median<- apply(MHYB,2,median,na.rm=T)
MHNA_median<- apply(MHNA,2,median,na.rm=T)
MHNB_median<- apply(MHNB,2,median,na.rm=T)
MLYA_median<- apply(MLYA,2,median,na.rm=T)
MLYB_median<- apply(MLYB,2,median,na.rm=T)
MLNA_median<- apply(MLNA,2,median,na.rm=T)
MLNB_median<- apply(MLNB,2,median,na.rm=T)
HC_median<-apply(HC,2,median,na.rm=T)
log10_FC_MHY<-MHYB_median-MHYA_median
log10_FC_MHN<-MHNB_median-MHNA_median
log10_FC_MLY<-MLYB_median-MLYA_median
log10_FC_MLN<-MLNB_median-MLNA_median
log10_FC_MHY_HC<-MHYA_median-HC_median
log10_FC_MLY_HC<-MLYA_median-HC_median

P_MHY<-rep(NA,188)
P_MHN<-rep(NA,188)
P_MLY<-rep(NA,188)
P_MLN<-rep(NA,188)
P_MHY_H<-rep(NA,188)
P_MLY_H<-rep(NA,188)

for(i in 1:188) try({
  P_MHY[i] = wilcox.test(as.numeric(MHYA[,i]), as.numeric(MHYB[,i]), alternative = "two.sided", paired = TRUE)$p.value
  P_MHN[i] = wilcox.test(as.numeric(MHNA[,i]), as.numeric(MHNB[,i]), alternative = "two.sided", paired = TRUE)$p.value
  P_MLY[i] = wilcox.test(as.numeric(MLYA[,i]), as.numeric(MLYB[,i]), alternative = "two.sided", paired = TRUE)$p.value
  P_MLN[i] = wilcox.test(as.numeric(MLNA[,i]), as.numeric(MLNB[,i]), alternative = "two.sided", paired = TRUE)$p.value
  P_MLY_H[i] = wilcox.test(as.numeric(MLYA[,i]), as.numeric(HC[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_MHY_H[i] = wilcox.test(as.numeric(MHYA[,i]), as.numeric(HC[,i]), alternative = "two.sided", paired = FALSE)$p.value
  
})
test<-data.frame(colnames(data),P_MHY,P_MHN,P_MLY,P_MLN,P_MHY_H,P_MLY_H,
                 log10_FC_MHY,log10_FC_MHN,log10_FC_MLY,log10_FC_MLN,log10_FC_MHY_HC,log10_FC_MLY_HC)
colnames(test)[1]<-'metabolites'
write.csv(test,"drug-effect-test.csv")
test<-read.csv("drug-effect-test.csv",row.names = 1)
data1<-test
data1$ML<- "no"
data1$ML[(data1$P_MLY<0.05)&(data1$log10_FC_MLY>0)] <- "up"
data1$ML[(data1$P_MLY<0.05)&(data1$log10_FC_MLY<0)]<-"down"
table(data1$ML)

data1$MH<- "no"
data1$MH[(data1$P_MHY<0.05)&(data1$log10_FC_MHY>0)] <- "up"
data1$MH[(data1$P_MHY<0.05)&(data1$log10_FC_MHY<0)]<-"down"
table(data1$MH)

library(ggplot2)
library(ggrepel)
library(eoffice)
ML_gg <- ggplot(data1, aes(x = log10_FC_MLY, y = -log10(P_MLY), color = ML)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('up' = '#DE3024', 'down' = '#009ACD', 'no' = 'grey')) + 
  geom_vline(xintercept = log10(1), linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(size = 0.5), 
    axis.line.y = element_line(size = 0.5), 
    axis.text = element_text(size = 10),
    axis.ticks = element_line(linewidth =0.5),
    axis.ticks.length = unit(0.15, "cm"),
    legend.position = "none"  
  ) +
  labs(
    x = "log10(After/Before)",
    y = "-log10(p value)",
    ML = "",
    size=10
  )
ML_gg

ML_gg <- ML_gg + geom_text_repel(
  data = data1[data1$ML!='no',],
  aes(x = log10_FC_MLY, y = -log10(P_MLY), label = metabolites),
  size = 3,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5  
)
ML_gg
topptx(ML_gg,filename = "ML_effect_wixon_plot.pptx")
####MH####
MH_gg <- ggplot(data1, aes(x = log10_FC_MHY, y = -log10(P_MHY), color = MH)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('up' = '#DE3024', 'down' = '#009ACD', 'no' = 'grey')) + 
  geom_vline(xintercept = log10(1), linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(size = 0.5), 
    axis.line.y = element_line(size = 0.5),
    axis.text = element_text(size = 10), 
    axis.ticks = element_line(linewidth =0.5),
    axis.ticks.length = unit(0.15, "cm"),
    legend.position = "none"  
  ) +
  labs(
    x = "log10(After/Before)",
    y = "-log10(p value)",
    MH = "",
    size=10
  )
MH_gg

MH_gg <- MH_gg + geom_text_repel(
  data = data1[data1$MH!='no',],
  aes(x = log10_FC_MHY, y = -log10(P_MHY), label = metabolites),
  size = 3,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5 
)
MH_gg
topptx(MH_gg,filename = "MH_effect_wixon_plot.pptx")
