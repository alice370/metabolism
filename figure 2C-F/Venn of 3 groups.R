setwd("C:/Users/alice/Desktop/实验数据/代谢组/数据处理/Code/figure 2C-E")
library(readxl)
rm(list = ls())
library('ggvenn')
RH<-read.csv('RH_GEEglm_meta.csv',row.names = 1)
RH$color<- "no"
RH$color[(RH$p_value<0.05)&(RH$Estimate >0)] <- "positive"
RH$color[(RH$p_value<0.05)&(RH$Estimate <0)]  <- "negative"
table(RH$color)

RI<-read.csv('RI_GEEglm_meta.csv',row.names = 1)
RI$color<- "no"
RI$color[(RI$p_value<0.05)&(RI$Estimate >0)] <- "positive"
RI$color[(RI$p_value<0.05)&(RI$Estimate <0)]  <- "negative"
table(RI$color)

IH<-read.csv('IH_GEEglm_meta.csv',row.names = 1)
IH$color<- "no"
IH$color[(IH$p_value<0.05)&(IH$Estimate >0)] <- "positive"
IH$color[(IH$p_value<0.05)&(IH$Estimate <0)]  <- "negative"
table(IH$color)

set1 <- RH$Metabolite[RH$color == 'positive']
set2 <- RI$Metabolite[RI$color == 'positive']
set3 <- IH$Metabolite[IH$color == 'positive']

venn_list <- list(
  "RA vs HC" = set1,
  "RA vs PRA" = set2,
  "PRA vs HC" = set3
)
p1<-ggvenn(venn_list,show_percentage = F,show_elements = F,label_sep = ",",
           digits = 1,stroke_color = "white",
           fill_color = c("#DC00007F", "#F39B7F7F","#D7A246"),
           set_name_color = "black")
p1
topptx(p1,filename = "up-intersect.pptx")
positive <- Reduce(intersect, list(set1, set2, set3))
write.csv(positive,'12 positive-common-metabolites.csv')


set1 <- RH$Metabolite[RH$color == 'negative']
set2 <- RI$Metabolite[RI$color == 'negative']
set3 <- IH$Metabolite[IH$color == 'negative']

venn_list <- list(
  "RA vs HC" = set1,
  "RA vs PRA" = set2,
  "PRA vs HC" = set3
)
p2<-ggvenn(venn_list,show_percentage = F,show_elements = F,label_sep = ",",
           digits = 1,stroke_color = "white",
           fill_color = c("#4DBBD57F", "#00A0877F","#91D1C27F"),
           set_name_color = "black")
p2
topptx(p2,filename = "down-intersect.pptx")
negative <- Reduce(intersect, list(set1, set2, set3))
write.csv(negative,'4 positive-common-metabolites.csv')
