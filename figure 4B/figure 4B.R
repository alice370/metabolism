library(ggplot2)
library(tidyr)
library(eoffice)
library(dplyr)
library(reshape)
rm(list=ls())
origin <- read.csv("364sample_185metabolites.csv",check.names = F,row.names = 1)
data <- read.csv("364sample_185metabolites_mean_log10.csv", header = TRUE,row.names = 1)
clin<-read.csv("RA_baseline_clin.csv",check.names = F)
H<-clin$sample[clin$`disease activity class`=='high disease activity']
M<-clin$sample[clin$`disease activity class`=='moderate disease activity']
L<-clin$sample[clin$`disease activity class`==" Low disease activity"]
D<-clin$sample[clin$`disease activity class`=="Clinical Remission"]
result_activity <- data.frame()

for (metabolite in row.names(origin)) {
  group_H <- as.numeric(origin[metabolite, H])
  group_M <- as.numeric(origin[metabolite, M])
  group_L <- as.numeric(origin[metabolite, L])
  group_D <- as.numeric(origin[metabolite, D])
  
  nonNA_H <- sum(!is.na(group_H))
  nonNA_M <- sum(!is.na(group_M))
  nonNA_L <- sum(!is.na(group_L))
  nonNA_D <- sum(!is.na(group_D))
 
  NA_H <- sum(is.na(group_H))
  NA_M <- sum(is.na(group_M))
  NA_L <- sum(is.na(group_L))
  NA_D <- sum(is.na(group_D))

  contingency_table <- matrix(c(nonNA_H, NA_H,
                                nonNA_M, NA_M,
                                nonNA_L, NA_L,
                                nonNA_D, NA_D),
                              nrow = 2,
                              dimnames = list(Status = c("NonNA", "NA"), 
                                              Activity = c("H", "M", "L", "D")))

  if (any(contingency_table < 5)) {
    test_result <- fisher.test(contingency_table)
    test_used <- "Fisher"
  } else {
    test_result <- chisq.test(contingency_table)
    test_used <- "Chi-square"
  }

  result_row <- data.frame(
    metabolic = metabolite,
    nonNA_H = nonNA_H,
    detected_H = nonNA_H / length(group_H),
    nonNA_M = nonNA_M,
    detected_M = nonNA_M / length(group_M),
    nonNA_L = nonNA_L,
    detected_L = nonNA_L / length(group_L),
    nonNA_D = nonNA_D,
    detected_D = nonNA_D / length(group_D),
    NA_test_pvalue = test_result$p.value,
    NA_test_used = test_used
  )
  
  result_activity <- rbind(result_activity, result_row)
}

write.csv(result_activity, 'disease activity-frequency.csv', row.names = FALSE)

clin<-read.csv("RA_baseline_clin.csv",check.names = F)
data<-as.data.frame(t(data))
data$sample<-row.names(data)
data<-merge(clin[,c(1,7,8)],data,by='sample')
table(data$`disease activity class`)
H<-data[which(data$`disease activity class` =="high disease activity"),-c(1:3)]
M<-data[which(data$`disease activity class`=="moderate disease activity"),-c(1:3)]
L<-data[which(data$`disease activity class`==" Low disease activity"),-c(1:3)]
D<-data[which(data$`disease activity class`=="Clinical Remission" ),-c(1:3)]
H_median<- apply(H,2,median,na.rm=T)
M_median<- apply(M,2,median,na.rm=T)
L_median<- apply(L,2,median,na.rm=T)
D_median<- apply(D,2,median,na.rm=T)
log10_FC_DH<-D_median-H_median
log10_FC_DM<-D_median-M_median
log10_FC_DL<-D_median-L_median
log10_FC_LH<-L_median-H_median
log10_FC_LM<-L_median-M_median
log10_FC_MH<-M_median-H_median
P_ANOVE<-rep(NA,185)
P_DH<-rep(NA,185)
P_DM<-rep(NA,185)
P_DL<-rep(NA,185)
P_LH<-rep(NA,185)
P_LM<-rep(NA,185)
P_MH<-rep(NA,185)

for(i in 1:185) try({
  P_DH[i] = wilcox.test(as.numeric(D[,i]), as.numeric(H[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_DM[i] = wilcox.test(as.numeric(D[,i]), as.numeric(M[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_DL[i] = wilcox.test(as.numeric(D[,i]), as.numeric(L[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_LH[i] = wilcox.test(as.numeric(L[,i]), as.numeric(H[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_LM[i] = wilcox.test(as.numeric(L[,i]), as.numeric(M[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_MH[i] = wilcox.test(as.numeric(M[,i]), as.numeric(H[,i]), alternative = "two.sided", paired = FALSE)$p.value
  P_ANOVE[i] = summary(aov(data[,(i+3)] ~ data$`disease activity class`, data = data))[[1]][,5][1]
})
test<-data.frame(colnames(data)[-c(1:3)],P_ANOVE,P_DH,P_DM,P_DL,P_LH,P_LM,P_MH,
                 log10_FC_DH,log10_FC_DM,log10_FC_DL,log10_FC_LH,log10_FC_LM,log10_FC_MH,D_median,L_median,M_median,H_median)
colnames(test)[1]<-'metabolites'
write.csv(test,"disease activity class test.csv")

########################Fig4A:plot of disease activity levels in different groups################
library(pheatmap)
data<-read.csv("disease activity class test.csv",row.names = 1,check.names = F)
data1<-data[which(data$P_DL<0.05|data$P_LM<0.05|data$P_MH<0.05|data$P_LH<0.05|data$P_DM<0.05|data$P_DH<0.05),]
data2<-data1[,c(15:18)]
data2<-data.frame(round(t(apply(data2, 1, scale)),2))
colnames(data2)<-c("Disease remission","Low disease activity","Moderate disease activity","High disease activity")
row.names(data2)<-row.names(data1)
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
p <-pheatmap(data2,scale='none',show_rownames = T, cluster_cols = F,clustering_method = "ward",border_color = NA,cutree_rows=4, 
             color = c(colorRampPalette(colors = c("#3C82B9","white"))(length(bk)/2),
                       colorRampPalette(colors = c("white","#EE3434"))(length(bk)/2)),
             cellwidth = 25,cellheight = 8,
             breaks=bk)
ggsave("disease-activity-cluster-phetmap.pdf", plot=p,width = 10, height = 10)

row_cluster = cutree(p$tree_row,k=4)
{
  newOrder = data2[p$tree_row$order,]
  newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder),names(row_cluster))]
  colnames(newOrder)[ncol(newOrder)]="Cluster"
  head(newOrder)
  unique(newOrder$Cluster)
  newOrder$Cluster = paste0("cluster",newOrder$Cluster)
  newOrder$gene = rownames(newOrder)
  head(newOrder)
  write.csv(newOrder,"disease-activity-4-cluster.csv")
  library(reshape2)
  data_new = melt(newOrder)
  head(data_new)
  library(stringr)
  data_new$group1<-data_new$variable
  data_new$group2<-str_c(data_new$gene)
  data_new$group3<-str_c(data_new$Cluster,data_new$variable)
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  medians<-aggregate(data_new$value,list(data_new$group3),median)
  medians$variable<-str_sub(medians$Group.1,9,-1)
  medians$new<-str_sub(medians$Group.1,1,8)
  medians$new1<-str_sub(medians$Group.1,1,8)
  
  medians$new
}
{p1<-ggplot(data_new[data_new$Cluster == "cluster1", ], aes(x = variable, y = value, group = variable)) +
    # Violin plot with grouping and coloring by variable
    geom_violin(aes(fill = variable), alpha = 0.7, size = 0.8) +
    geom_line(data=medians[medians$new1=="cluster1",],aes(variable,x,group=new),size=0.75,color="black")+
    # Boxplot overlaid on violin plots
    geom_boxplot(aes(group = variable), outlier.shape = NA, width = 0.42, size = 0.5, fill = "white") +
    # Manual color scale for variable groups
    scale_fill_manual(values = c("Disease remission" = "#638759", "Low disease activity" = "#A1B4D8", "Moderate disease activity" = "#eca75e","High disease activity"="#D8ACAC"), guide = FALSE) +
    # Labels for axes
    labs(x = "", y = "") +
    # Classic theme with custom modifications
    theme_classic() +
    theme(
      strip.background.x = element_rect(color = "white", fill = "white"),
      strip.text.x = element_text(size = 15, color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 15, color = "black", vjust = 0.5, hjust = 1),
      axis.text.x = element_blank(),
      title = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5)
    ) 
  
  p2<-ggplot(data_new[data_new$Cluster == "cluster2", ], aes(x = variable, y = value, group = variable)) +
    # Violin plot with grouping and coloring by variable
    geom_violin(aes(fill = variable), alpha = 0.7, size = 0.8) +
    geom_line(data=medians[medians$new1=="cluster2",],aes(variable,x,group=new),size=0.75,color="black")+
    # Boxplot overlaid on violin plots
    geom_boxplot(aes(group = variable), outlier.shape = NA, width = 0.42, size = 0.5, fill = "white") +
    # Manual color scale for variable groups
    scale_fill_manual(values = c("Disease remission" = "#638759", "Low disease activity" = "#A1B4D8", "Moderate disease activity" = "#eca75e","High disease activity"="#D8ACAC"), guide = FALSE) +
    # Labels for axes
    labs(x = "", y = "") +
    # Classic theme with custom modifications
    theme_classic() +
    theme(
      strip.background.x = element_rect(color = "white", fill = "white"),
      strip.text.x = element_text(size = 15, color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 15, color = "black", vjust = 0.5, hjust = 1),
      axis.text.x = element_blank(),
      title = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5)
    ) 
  
  p3<-ggplot(data_new[data_new$Cluster == "cluster3", ], aes(x = variable, y = value, group = variable)) +
    # Violin plot with grouping and coloring by variable
    geom_violin(aes(fill = variable), alpha = 0.7, size = 0.8) +
    geom_line(data=medians[medians$new1=="cluster3",],aes(variable,x,group=new),size=0.75,color="black")+
    # Boxplot overlaid on violin plots
    geom_boxplot(aes(group = variable), outlier.shape = NA, width = 0.42, size = 0.5, fill = "white") +
    # Manual color scale for variable groups
    scale_fill_manual(values = c("Disease remission" = "#638759", "Low disease activity" = "#A1B4D8", "Moderate disease activity" = "#eca75e","High disease activity"="#D8ACAC"), guide = FALSE) +
    # Labels for axes
    labs(x = "", y = "") +
    # Classic theme with custom modifications
    theme_classic() +
    theme(
      strip.background.x = element_rect(color = "white", fill = "white"),
      strip.text.x = element_text(size = 15, color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 15, color = "black", vjust = 0.5, hjust = 1),
      axis.text.x = element_blank(),
      title = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5)
    ) 
  p4<-ggplot(data_new[data_new$Cluster == "cluster4", ], aes(x = variable, y = value, group = variable)) +
    # Violin plot with grouping and coloring by variable
    geom_violin(aes(fill = variable), alpha = 0.7, size = 0.8) +
    geom_line(data=medians[medians$new1=="cluster4",],aes(variable,x,group=new),size=0.75,color="black")+
    # Boxplot overlaid on violin plots
    geom_boxplot(aes(group = variable), outlier.shape = NA, width = 0.42, size = 0.5, fill = "white") +
    # Manual color scale for variable groups
    scale_fill_manual(values = c("Disease remission" = "#638759", "Low disease activity" = "#A1B4D8", "Moderate disease activity" = "#eca75e","High disease activity"="#D8ACAC"), guide = FALSE) +
    # Labels for axes
    labs(x = "", y = "") +
    # Classic theme with custom modifications
    theme_classic() +
    theme(
      strip.background.x = element_rect(color = "white", fill = "white"),
      strip.text.x = element_text(size = 15, color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 15, color = "black", vjust = 0.5, hjust = 1),
      axis.text.x = element_blank(),
      title = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5)
    ) 
}
a<-ggarrange(p1,p2,p3,p4,nrow=4,ncol =1,legend = NULL)
a
ggsave("disease-activity-class-violin.pdf", width = 3.3, height = 8.35)
