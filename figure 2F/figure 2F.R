library(ggplot2)
library(eoffice)
rm(list=ls())
results_df<-read.csv("metabolite_comparison_results_scaled.csv")
str(results_df)
results_df$RA_HC.color_ra<- "no"
results_df$RA_HC.color_ra[(results_df$PValue_RA_HC<0.05)&(results_df$Coef_RA_HC>0)] <- "positive"
results_df$RA_HC.color_ra[(results_df$PValue_RA_HC<0.05)&(results_df$Coef_RA_HC<0)]  <- "negative"
table(results_df$RA_HC.color_ra)

results_df$PRA_HC.color_ra<- "no"
results_df$PRA_HC.color_ra[(results_df$PValue_HC_atRisk<0.05)&(results_df$Coef_HC_atRisk<0)] <- "positive"
results_df$PRA_HC.color_ra[(results_df$PValue_HC_atRisk<0.05)&(results_df$Coef_HC_atRisk>0)]  <- "negative"
table(results_df$PRA_HC.color_ra)

results_df$RA_PRA.color_ra<- "no"
results_df$RA_PRA.color_ra[(results_df$PValue_RA_atRisk<0.05)&(results_df$Coef_RA_atRisk>0)] <- "positive"
results_df$RA_PRA.color_ra[(results_df$PValue_RA_atRisk<0.05)&(results_df$Coef_RA_atRisk<0)]  <- "negative"
table(results_df$RA_PRA.color_ra)

library('ggvenn')

set1 <- results_df$Metabolite[results_df$RA_HC.color_ra == 'positive']
set2 <- results_df$Metabolite[results_df$PRA_HC.color_ra == 'positive']
set3 <- results_df$Metabolite[results_df$RA_PRA.color_ra == 'positive']

venn_list <- list(
  "RA vs HC" = set1,
  "PRA vs HC" = set2,
  "RA vs PRA" = set3
)
p1<-ggvenn(venn_list,show_percentage = F,show_elements = F,label_sep = ",",
           digits = 1,stroke_color = "white",
           fill_color = c("#DC00007F", "#F39B7F7F","#D7A246"),
           set_name_color = "black")
p1
topptx(p1,filename = "up-intersect.pptx")
positive <- Reduce(intersect, list(set1, set2, set3))


set1 <- results_df$Metabolite[results_df$RA_HC.color_ra == 'negative']
set2 <- results_df$Metabolite[results_df$PRA_HC.color_ra == 'negative']
set3 <- results_df$Metabolite[results_df$RA_PRA.color_ra == 'negative']

venn_list <- list(
  "RA vs HC" = set1,
  "PRA vs HC" = set2,
  "RA vs PRA" = set3
)
p2<-ggvenn(venn_list,show_percentage = F,show_elements = F,label_sep = ",",
           digits = 1,stroke_color = "white",
           fill_color = c("#4DBBD57F", "#00A0877F","#91D1C27F"),
           set_name_color = "black")
p2
topptx(p2,filename = "down-intersect.pptx")
negative <- Reduce(intersect, list(set1, set2, set3))
