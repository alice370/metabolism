library(GSVA)
file_content <- readLines("gsea_geneset.txt")
rm(list=ls())

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
data<-read.csv('364sample-188metabolites-log10-auto.csv',row.names = 1)
clin<-read.csv('sample_information.csv')

param <- gsvaParam(as.matrix(data), gsea_geneset,kcdf = c("Gaussian"),minSize = 3)

gsva_matrix <- gsva(param, verbose = TRUE)
gsva_matrix<-as.data.frame(gsva_matrix)
write.csv(gsva_matrix,"gsva_pathway.csv")

pathway_names<-row.names(gsva_matrix)

normality_results <- data.frame(
  Metabolite = row.names(gsva_matrix),
  P_Value = apply(gsva_matrix, 1, function(x) shapiro.test(x)$p.value)
)
write.csv(normality_results,"shapiro.test_gsva.csv")


data <- as.data.frame(t(gsva_matrix))
data$Sample <- row.names(data)


metabolite_data <- merge(clin[, c(1:4)], data, by = "Sample")

results_df <- data.frame()

for (i in 5:ncol(metabolite_data)) {

  metabolite_name <- colnames(metabolite_data)[i]

  metabolite_abundance <- metabolite_data[, i]

  temp_data <- data.frame(
    metabolite_abundance = metabolite_abundance,
    group = metabolite_data$Group,
    sex = as.factor(metabolite_data$Gender),
    age = as.numeric(metabolite_data$Age)
  )

  model_RA_HC <- lm(metabolite_abundance ~ group + sex + age, 
                    data = subset(temp_data, group %in% c("RA", "Health")))
  coef_summary_RA_HC <- summary(model_RA_HC)$coefficients
  confint_RA_HC <- confint(model_RA_HC, "groupRA", level = 0.95)
  coef_RA_HC <- coef_summary_RA_HC["groupRA", c("Estimate", "Pr(>|t|)")]
  

  model_RA_atRisk <- lm(metabolite_abundance ~ group + sex + age, 
                        data = subset(temp_data, group %in% c("RA", "At-risk of RA")))
  coef_summary_RA_atRisk <- summary(model_RA_atRisk)$coefficients
  confint_RA_atRisk <- confint(model_RA_atRisk, "groupRA", level = 0.95)
  coef_RA_atRisk <- coef_summary_RA_atRisk["groupRA", c("Estimate", "Pr(>|t|)")]

  model_HC_atRisk <- lm(metabolite_abundance ~ group + sex + age, 
                        data = subset(temp_data, group %in% c("Health", "At-risk of RA")))
  coef_summary_HC_atRisk <- summary(model_HC_atRisk)$coefficients
  confint_HC_atRisk <- confint(model_HC_atRisk, "groupHealth", level = 0.95)
  coef_HC_atRisk <- coef_summary_HC_atRisk["groupHealth", c("Estimate", "Pr(>|t|)")]

  coef_df <- data.frame(
    Metabolite = metabolite_name,
    
    # RA vs HC
    Coef_RA_HC = coef_RA_HC["Estimate"],
    PValue_RA_HC = coef_RA_HC["Pr(>|t|)"],
    CI_Lower_RA_HC = confint_RA_HC[1],
    CI_Upper_RA_HC = confint_RA_HC[2],
    
    # RA vs At-risk
    Coef_RA_atRisk = coef_RA_atRisk["Estimate"],
    PValue_RA_atRisk = coef_RA_atRisk["Pr(>|t|)"],
    CI_Lower_RA_atRisk = confint_RA_atRisk[1],
    CI_Upper_RA_atRisk = confint_RA_atRisk[2],
    
    # HC vs At-risk
    Coef_HC_atRisk = coef_HC_atRisk["Estimate"],
    PValue_HC_atRisk = coef_HC_atRisk["Pr(>|t|)"],
    CI_Lower_HC_atRisk = confint_HC_atRisk[1],
    CI_Upper_HC_atRisk = confint_HC_atRisk[2]
  )

  results_df <- rbind(results_df, coef_df)
}


data1<-results_df
str(data1)
data1$RA_HC.color_ra<- "no"
data1$RA_HC.color_ra[(data1$PValue_RA_HC<0.01)&(data1$Coef_RA_HC>0)] <- "positive"
data1$RA_HC.color_ra[(data1$PValue_RA_HC<0.01)&(data1$Coef_RA_HC<0)]  <- "negative"
table(data1$RA_HC.color_ra)

data1$PRA_HC.color_ra<- "no"
data1$PRA_HC.color_ra[(data1$PValue_HC_atRisk<0.01)&(data1$Coef_HC_atRisk<0)] <- "positive"
data1$PRA_HC.color_ra[(data1$PValue_HC_atRisk<0.01)&(data1$Coef_HC_atRisk>0)]  <- "negative"
table(data1$PRA_HC.color_ra)

data1$RA_PRA.color_ra<- "no"
data1$RA_PRA.color_ra[(data1$PValue_RA_atRisk<0.01)&(data1$Coef_RA_atRisk>0)] <- "positive"
data1$RA_PRA.color_ra[(data1$PValue_RA_atRisk<0.01)&(data1$Coef_RA_atRisk<0)]  <- "negative"
table(data1$RA_PRA.color_ra)

write.csv(data1, "gsva_comparison_results_scaled.csv", row.names = FALSE)
sig_gsva<-read.csv('selected_sig_gsva_pathway.csv')
plot_data<-results_df[results_df$Metabolite%in%sig_gsva$Metabolite,]
write.csv(plot_data,"top_GSVA_lm.csv")
test<-plot_data[,c(1,2,6,10)]
test$Coef_HC_atRisk<--test$Coef_HC_atRisk
colnames(test)[4]<-'Coef_atRisk_HC'
row.names(test)<-test$Metabolite
test<-test[,-1]
library(pheatmap)
library(eoffice)

breaks <- c(seq(min(test), 0, length.out = 128), seq(0, max(test), length.out = 129)[-1])

p <- pheatmap(
  test,
  scale = "none",  
  color = colorRampPalette(c("#3980b8", "white", "#ef3b3c"))(256),
  breaks = breaks, 
  cluster_cols = FALSE,
  cellwidth = 20,
  cellheight = 12,
  border_color = "white",
  fontsize = 8,
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 0,
  treeheight_col = 0,
  cutree_rows = 1,
  clustering_method = "average",
  legend = TRUE
)
topptx(p,filename = "gsva-pheatmap.pptx")
