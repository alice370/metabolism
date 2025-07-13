rm(list=ls())
library(tidyr)
library(ggplot2)
library(eoffice)
library(GSVA)
library(ggpubr)
rm(list=ls())
file_content <- readLines("gsea_geneset.txt")

gsea_geneset <- list()

for (line in file_content) {

  elements <- strsplit(line, "\t")[[1]]

  elements <- gsub('"', '', elements)
  elements <- elements[elements != ""]

  list_name <- gsub("\\(.*", "", elements[1])

  if (!list_name %in% names(gsea_geneset)) {
    list_elements <- elements[-1]
   
    gsea_geneset[[list_name]] <- list_elements
  }
}
data<-read.csv('561sample_185meta_log10_auto.csv',row.names = 1)
clin<-read.csv('sample_information.csv',check.names = F)

param <- gsvaParam(as.matrix(data), gsea_geneset,kcdf = c("Gaussian"),minSize = 3)

gsva_matrix <- gsva(param, verbose = TRUE)
gsva_matrix<-as.data.frame(gsva_matrix)
write.csv(gsva_matrix,"gsva_pathway.csv")

pathway_names<-row.names(gsva_matrix)


data <- as.data.frame(t(gsva_matrix))
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

P_MHY<-rep(NA,447)
P_MHN<-rep(NA,447)
P_MLY<-rep(NA,447)
P_MLN<-rep(NA,447)
P_MHY_H<-rep(NA,447)
P_MLY_H<-rep(NA,447)

for(i in 1:447) try({
  P_MHY[i] = wilcox.test(as.numeric(MHYA[,i]), as.numeric(MHYB[,i]), alternative = "two.sided", paired = TRUE)$p.value
  P_MHN[i] = wilcox.test(as.numeric(MHNA[,i]), as.numeric(MHNB[,i]), alternative = "two.sided", paired = TRUE)$p.value
  P_MLY[i] = wilcox.test(as.numeric(MLYA[,i]), as.numeric(MLYB[,i]), alternative = "two.sided", paired = TRUE)$p.value
  P_MLN[i] = wilcox.test(as.numeric(MLNA[,i]), as.numeric(MLNB[,i]), alternative = "two.sided", paired = TRUE)$p.value
  P_MLY_H[i] = wilcox.test(as.numeric(MLYA[,i]), as.numeric(HC[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_MHY_H[i] = wilcox.test(as.numeric(MHYA[,i]), as.numeric(HC[,i]), alternative = "two.sided", paired = FALSE)$p.value
  
})
test<-data.frame(colnames(data),P_MHY,P_MHN,P_MLY,P_MLN,P_MHY_H,P_MLY_H,
                 log10_FC_MHY,log10_FC_MHN,log10_FC_MLY,log10_FC_MLN,log10_FC_MHY_HC,log10_FC_MLY_HC)
colnames(test)[1]<-'pathway'
write.csv(test,"drug-effect-gsva-test.csv")
####plot on selected pathway£¨HC-BEFORE-AFTER<0.05£©####
test<-read.csv("drug-effect-gsva-test.csv",row.names = 1)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(eoffice)
library(dplyr)
library(tidyr)
rm(list = ls())
data<-read.csv('gsva_pathway.csv',row.names = 1)
sig<-read.table('sig_meta_gsva.txt',sep=',',header = T)
sig_data<-data[sig$pathway,]
clin<-read.csv('sample_information.csv',check.names = F)
test<-read.csv('drug-effect-test.csv',row.names = 1)

MLYA<-clin$sample[which(clin$DRUG =="MTQ+LEF"&
                          clin$`response classification` %in% c('Good Response', 'Moderate Response'))]
MLYB<-clin$`follow up-sample`[which(clin$DRUG =="MTQ+LEF"&
                             clin$`response classification` %in% c('Good Response', 'Moderate Response'))]
HC<-clin$sample[clin$Group=='Health']
####ML####

ML_dt<-sig_data[,c(HC,MLYA,MLYB)]
group_labels <- rep(c("HC","MLYA", "MLYB"), 
                    times = c(length(HC),length(MLYA), length(MLYB)))
ML_dt$name<-row.names(ML_dt)
long_data <- ML_dt %>%
  pivot_longer(-name,names_to = "sample", values_to = "Value")
long_data$group <-rep(group_labels,6)
long_data$group <- factor(long_data$group, levels = c("HC","MLYA", "MLYB"))
ML<-row.names(sig_data)
####ML####
for (i in 1:length(ML)) {
  meta<-ML[i]
  plot<-long_data[long_data$name==meta,]
  max_values<-max(plot$Value)
  p <- ggplot(plot, aes(x = group, y = Value, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.5,width=1, size = 0.1) + 
    geom_boxplot(width = 0.1, fill = "white",size = 0.1, outlier.size = 0.1) +  
    scale_fill_manual(values = c(
      "HC" = "#4E79A7", "MLYA" = "#E15759", "MLYB" = "#F28E2B")) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill = 'transparent', color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.1),
      axis.ticks = element_line(colour = "black", linewidth = 0.1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 3),
      axis.text.y = element_text(size = 3, color = "black")
    ) +
    ylab('') +
    labs(title = meta) +

    stat_compare_means(comparisons = list(c("MLYA", "MLYB")),
                       method = "wilcox.test", 
                       paired = TRUE,  # Åä¶Ô±È½Ï
                       label = "p.signif", size = 1,label.y = max_values+0.5,tip.length=0, bracket.size = 0.1) +
    

    stat_compare_means(comparisons = list(c("HC", "MLYA")),
                       method = "wilcox.test", 
                       paired = FALSE,  
                       label = "p.signif", size = 1,label.y = max_values+0.3,tip.length=0,bracket.size = 0.1) +

    stat_compare_means(comparisons = list(c("HC", "MLYB")),
                       method = "wilcox.test", 
                       paired = FALSE,  
                       label = "p.signif", size = 1,label.y = max_values+0.05,tip.length=0,bracket.size = 0.1)+
    ylim(NA,max_values+0.7)
  
  assign(paste0("p",i),p)
}

a<-ggarrange(p1,p2,p3,p4,p5,p6,nrow=2,ncol =3,legend = NULL)
a
cairo_pdf("plot.pdf", width = 130/25.4, height = 42/25.4)
print(a)
dev.off()
