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
data <- read_csv("364sample_185metabolites_mean_log10.csv")

param <- gsvaParam(as.matrix(data), gsea_geneset,kcdf = c("Gaussian"),minSize = 3)

gsva_matrix <- gsva(param, verbose = TRUE)
gsva_matrix<-as.data.frame(gsva_matrix)
write.csv(gsva_matrix,"gsva_pathway.csv")
####geeglm####
library(geepack)  
library(readxl)
rm(list=ls())
gsva_matrix<-read.csv("gsva_pathway.csv",row.names = 1)
IAR_information <- read_excel("IAR information.xlsx")
RA_information <- read_excel("RA information.xlsx")
HC_information <- read_excel("HC information.xlsx")
RA_information$group <- "RA"
IAR_information$group <- "IAR"
HC_information$group <- "HC"

gsva_data<-as.data.frame(t(gsva_matrix))
gsva_data$Sample<-row.names(gsva_data)
# 只保留需要的列（确保列名一致）
RA_sub  <- RA_information[, c("Sample", 'Gender','Age','smoke',"BMI", "BP", "BG", "BL", "group")]
colnames(RA_sub)[2:3]<-c('gender','age')
IAR_sub <- IAR_information[, c("Sample", 'gender','age','smoke',"BMI", "BP", "BG", "BL", "group")]
HC_sub  <- HC_information[, c("Sample", 'gender','age','smoke',"BMI", "BP", "BG", "BL", "group")]
# 合并
combined_data <- rbind(RA_sub, IAR_sub, HC_sub)
####RA vs HC####
RH_df<-combined_data[combined_data$group%in%c('RA','HC'),]
RH_df<-merge(RH_df,gsva_data,by='Sample',all.x = TRUE)
RH_df$group <- as.numeric(ifelse(RH_df$group == "RA", 1,
                                 ifelse(RH_df$group == "HC", 0, RH_df$group)))
RH_df$gender <- as.numeric( ifelse(RH_df$gender == "Female", 1,
                                   ifelse(RH_df$group == "Male", 0, RH_df$group)))
table(RH_df$group)
table(RH_df$gender)
RH_df$ID<-c(1:nrow(RH_df))
gsva_list<-colnames(RH_df)[10:ncol(RH_df)]

RH_results_list <- list()

for (meta in gsva_list) {
  formula <- as.formula(paste0("`", meta, "`~ group  + gender + age + BMI + smoke + BP + BG + BL"))
  
  model <- geeglm(formula, data = RH_df, id = ID, family = gaussian, corstr = "exchangeable")
  coef_summary <- summary(model)$coefficients
  
  est <- coef_summary[2, "Estimate"]
  se  <- coef_summary[2, "Std.err"]
  ci_lower <- if (!is.na(se)) est - 1.96 * se else NA
  ci_upper <- if (!is.na(se)) est + 1.96 * se else NA
  p_value <- coef_summary[2, "Pr(>|W|)"]
  
  # 无论有无值，都保存
  RH_results_list[[meta]] <- data.frame(
    Metabolite = meta,
    Estimate.RH = est,
    Std.err.RH = se,
    CI_lower.RH = ci_lower,
    CI_upper.RH = ci_upper,
    p_value.RH = p_value,  # 用格式化后的值
    stringsAsFactors = FALSE
  )
}
# 合并所有结果为一个数据框
RH_final_results <- do.call(rbind, RH_results_list)
RH_final_results$sig <-ifelse(RH_final_results$p_value<0.05, 'sig','not sig')
table(RH_final_results$sig)

write.csv(RH_final_results,'RH_GEEglm_gsva.csv')

####RA vs IAR####
RI_df<-combined_data[combined_data$group%in%c('RA','IAR'),]
RI_df<-merge(RI_df,gsva_data,by='Sample',all.x = TRUE)
RI_df$group <- as.numeric(ifelse(RI_df$group == "RA", 1,
                                 ifelse(RI_df$group == "IAR", 0, RI_df$group)))
RI_df$gender <- as.numeric( ifelse(RI_df$gender == "Female", 1,
                                   ifelse(RI_df$group == "Male", 0, RI_df$group)))
table(RI_df$group)
table(RI_df$gender)
RI_df$ID<-c(1:nrow(RI_df))
gsva_list<-colnames(RI_df)[10:setwd("C:/Users/alice/Desktop/实验数据/代谢组/数据处理/Code/figure 2H")
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
data <- read_csv("364sample_185metabolites_mean_log10.csv")
clin<-read.csv('sample_information.csv')

param <- gsvaParam(as.matrix(data), gsea_geneset,kcdf = c("Gaussian"),minSize = 3)

gsva_matrix <- gsva(param, verbose = TRUE)
gsva_matrix<-as.data.frame(gsva_matrix)
write.csv(gsva_matrix,"gsva_pathway.csv")
####geeglm####
library(geepack)  
library(readxl)
rm(list=ls())
gsva_matrix<-read.csv("gsva_pathway.csv",row.names = 1)
IAR_information <- read_excel("IAR information.xlsx")
RA_information <- read_excel("RA information.xlsx")
HC_information <- read_excel("HC information.xlsx")
RA_information$group <- "RA"
IAR_information$group <- "IAR"
HC_information$group <- "HC"

gsva_data<-as.data.frame(t(gsva_matrix))
gsva_data$Sample<-row.names(gsva_data)
# 只保留需要的列（确保列名一致）
RA_sub  <- RA_information[, c("Sample", 'Gender','Age','smoke',"BMI", "BP", "BG", "BL", "group")]
colnames(RA_sub)[2:3]<-c('gender','age')
IAR_sub <- IAR_information[, c("Sample", 'gender','age','smoke',"BMI", "BP", "BG", "BL", "group")]
HC_sub  <- HC_information[, c("Sample", 'gender','age','smoke',"BMI", "BP", "BG", "BL", "group")]
# 合并
combined_data <- rbind(RA_sub, IAR_sub, HC_sub)
####RA vs HC####
RH_df<-combined_data[combined_data$group%in%c('RA','HC'),]
RH_df<-merge(RH_df,gsva_data,by='Sample',all.x = TRUE)
RH_df$group <- as.numeric(ifelse(RH_df$group == "RA", 1,
                                 ifelse(RH_df$group == "HC", 0, RH_df$group)))
RH_df$gender <- as.numeric( ifelse(RH_df$gender == "Female", 1,
                                   ifelse(RH_df$group == "Male", 0, RH_df$group)))
table(RH_df$group)
table(RH_df$gender)
RH_df$ID<-c(1:nrow(RH_df))
gsva_list<-colnames(RH_df)[10:ncol(RH_df)]

RH_results_list <- list()

for (meta in gsva_list) {
  formula <- as.formula(paste0("`", meta, "`~ group  + gender + age + BMI + smoke + BP + BG + BL"))
  
  model <- geeglm(formula, data = RH_df, id = ID, family = gaussian, corstr = "exchangeable")
  coef_summary <- summary(model)$coefficients
  
  est <- coef_summary[2, "Estimate"]
  se  <- coef_summary[2, "Std.err"]
  ci_lower <- if (!is.na(se)) est - 1.96 * se else NA
  ci_upper <- if (!is.na(se)) est + 1.96 * se else NA
  p_value <- coef_summary[2, "Pr(>|W|)"]
  
  # 无论有无值，都保存
  RH_results_list[[meta]] <- data.frame(
    Metabolite = meta,
    Estimate.RH = est,
    Std.err.RH = se,
    CI_lower.RH = ci_lower,
    CI_upper.RH = ci_upper,
    p_value.RH = p_value,  # 用格式化后的值
    stringsAsFactors = FALSE
  )
}
# 合并所有结果为一个数据框
RH_final_results <- do.call(rbind, RH_results_list)
RH_final_results$sig <-ifelse(RH_final_results$p_value<0.05, 'sig','not sig')
table(RH_final_results$sig)

write.csv(RH_final_results,'RH_GEEglm_gsva.csv')

####RA vs IAR####
RI_df<-combined_data[combined_data$group%in%c('RA','IAR'),]
RI_df<-merge(RI_df,gsva_data,by='Sample',all.x = TRUE)
RI_df$group <- as.numeric(ifelse(RI_df$group == "RA", 1,
                                 ifelse(RI_df$group == "IAR", 0, RI_df$group)))
RI_df$gender <- as.numeric( ifelse(RI_df$gender == "Female", 1,
                                   ifelse(RI_df$group == "Male", 0, RI_df$group)))
table(RI_df$group)
table(RI_df$gender)
RI_df$ID<-c(1:nrow(RI_df))
gsva_list<-colnames(RI_df)[10:ncol(RI_df)]

RI_results_list <- list()

for (meta in gsva_list) {
  formula <- as.formula(paste0("`", meta, "`~ group  + gender + age + BMI + smoke + BP + BG + BL"))
  
  model <- geeglm(formula, data = RI_df, id = ID, family = gaussian, corstr = "exchangeable")
  coef_summary <- summary(model)$coefficients
  
  est <- coef_summary[2, "Estimate"]
  se  <- coef_summary[2, "Std.err"]
  ci_lower <- if (!is.na(se)) est - 1.96 * se else NA
  ci_upper <- if (!is.na(se)) est + 1.96 * se else NA
  p_value <- coef_summary[2, "Pr(>|W|)"]
  
  # 无论有无值，都保存
  RI_results_list[[meta]] <- data.frame(
    Metabolite = meta,
    Estimate.RI = est,
    Std.err.RI = se,
    CI_lower.RI = ci_lower,
    CI_upper.RI = ci_upper,
    p_value.RI = p_value,
    stringsAsFactors = FALSE
  )
}
# 合并所有结果为一个数据框
RI_final_results <- do.call(rbind, RI_results_list)
RI_final_results$sig <-ifelse(RI_final_results$p_value<0.05, 'sig','not sig')
table(RI_final_results$sig)
write.csv(RI_final_results,'RI_GEEglm_gsva.csv')

####IAR VS HC####
IH_df<-combined_data[combined_data$group%in%c('IAR','HC'),]
IH_df<-merge(IH_df,gsva_data,by='Sample',all.x = TRUE)
IH_df$group <- as.numeric(ifelse(IH_df$group == "IAR", 1,
                                 ifelse(IH_df$group == "HC", 0, IH_df$group)))
IH_df$gender <- as.numeric( ifelse(IH_df$gender == "Female", 1,
                                   ifelse(IH_df$group == "Male", 0, IH_df$group)))
table(IH_df$group)
table(IH_df$gender)
IH_df$ID<-c(1:nrow(IH_df))
####move the metabolites####
gsva_list<-colnames(RI_df)[10:ncol(IH_df)]

IH_results_list <- list()

for (meta in gsva_list) {
  formula <- as.formula(paste0("`", meta, "`~ group  + gender + age + BMI + smoke + BP + BG + BL"))
  
  model <- geeglm(formula, data = IH_df, id = ID, family = gaussian, corstr = "exchangeable")
  coef_summary <- summary(model)$coefficients
  
  est <- coef_summary[2, "Estimate"]
  se  <- coef_summary[2, "Std.err"]
  ci_lower <- if (!is.na(se)) est - 1.96 * se else NA
  ci_upper <- if (!is.na(se)) est + 1.96 * se else NA
  p_value <- coef_summary[2, "Pr(>|W|)"]
  
  # 无论有无值，都保存
  IH_results_list[[meta]] <- data.frame(
    Metabolite = meta,
    Estimate.IH = est,
    Std.err.IH = se,
    CI_lower.IH = ci_lower,
    CI_upper.IH = ci_upper,
    p_value.IH = p_value,
    stringsAsFactors = FALSE
  )
}
# 合并所有结果为一个数据框
IH_final_results <- do.call(rbind, IH_results_list)
IH_final_results$sig <-ifelse(IH_final_results$p_value<0.05, 'sig','not sig')
table(IH_final_results$sig)
write.csv(IH_final_results,'IH_GEEglm_gsva.csv')


Combine1<-merge(RH_final_results,RI_final_results,by='Metabolite',all.x=TRUE)
Combine_final<-merge(Combine1,IH_final_results,by='Metabolite',all.x=TRUE)

write.csv(Combine_final, "Combined_gsva_geeglm.csv", row.names = FALSE)
####select all sig pathway as top pathway####
metabolites<-data$metabolites
top_pathway <- read_excel("top_pathway.xlsx")
top_pathway$overlap <- sapply(top_pathway$Metabolite, function(p) {
  if (! p %in% names(gsea_geneset)) return(NA_character_)
  ov <- intersect(gsea_geneset[[p]], metabolites)
  if (length(ov)==0) return(NA_character_)
  paste(ov, collapse = ";")  # 用分号连接，或改成 "," 
})
write.csv(top_pathway,'top_pathway.csv')
plot_data <- read_excel("top_pathway.xlsx")
test<-as.data.frame(plot_data[,c(1,2,8,14)])
row.names(test)<-test$pathway
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

RI_results_list <- list()

for (meta in gsva_list) {
  formula <- as.formula(paste0("`", meta, "`~ group  + gender + age + BMI + smoke + BP + BG + BL"))
  
  model <- geeglm(formula, data = RI_df, id = ID, family = gaussian, corstr = "exchangeable")
  coef_summary <- summary(model)$coefficients
  
  est <- coef_summary[2, "Estimate"]
  se  <- coef_summary[2, "Std.err"]
  ci_lower <- if (!is.na(se)) est - 1.96 * se else NA
  ci_upper <- if (!is.na(se)) est + 1.96 * se else NA
  p_value <- coef_summary[2, "Pr(>|W|)"]
  
  # 无论有无值，都保存
  RI_results_list[[meta]] <- data.frame(
    Metabolite = meta,
    Estimate.RI = est,
    Std.err.RI = se,
    CI_lower.RI = ci_lower,
    CI_upper.RI = ci_upper,
    p_value.RI = p_value,
    stringsAsFactors = FALSE
  )
}
# 合并所有结果为一个数据框
RI_final_results <- do.call(rbind, RI_results_list)
RI_final_results$sig <-ifelse(RI_final_results$p_value<0.05, 'sig','not sig')
table(RI_final_results$sig)
write.csv(RI_final_results,'RI_GEEglm_gsva.csv')

####IAR VS HC####
IH_df<-combined_data[combined_data$group%in%c('IAR','HC'),]
IH_df<-merge(IH_df,gsva_data,by='Sample',all.x = TRUE)
IH_df$group <- as.numeric(ifelse(IH_df$group == "IAR", 1,
                                 ifelse(IH_df$group == "HC", 0, IH_df$group)))
IH_df$gender <- as.numeric( ifelse(IH_df$gender == "Female", 1,
                                   ifelse(IH_df$group == "Male", 0, IH_df$group)))
table(IH_df$group)
table(IH_df$gender)
IH_df$ID<-c(1:nrow(IH_df))
####move the metabolites####
gsva_list<-colnames(RI_df)[10:115]

IH_results_list <- list()

for (meta in gsva_list) {
  formula <- as.formula(paste0("`", meta, "`~ group  + gender + age + BMI + smoke + BP + BG + BL"))
  
  model <- geeglm(formula, data = IH_df, id = ID, family = gaussian, corstr = "exchangeable")
  coef_summary <- summary(model)$coefficients
  
  est <- coef_summary[2, "Estimate"]
  se  <- coef_summary[2, "Std.err"]
  ci_lower <- if (!is.na(se)) est - 1.96 * se else NA
  ci_upper <- if (!is.na(se)) est + 1.96 * se else NA
  p_value <- coef_summary[2, "Pr(>|W|)"]
  
  # 无论有无值，都保存
  IH_results_list[[meta]] <- data.frame(
    Metabolite = meta,
    Estimate.IH = est,
    Std.err.IH = se,
    CI_lower.IH = ci_lower,
    CI_upper.IH = ci_upper,
    p_value.IH = p_value,
    stringsAsFactors = FALSE
  )
}
# 合并所有结果为一个数据框
IH_final_results <- do.call(rbind, IH_results_list)
IH_final_results$sig <-ifelse(IH_final_results$p_value<0.05, 'sig','not sig')
table(IH_final_results$sig)
write.csv(IH_final_results,'IH_GEEglm_gsva.csv')


Combine1<-merge(RH_final_results,RI_final_results,by='Metabolite',all.x=TRUE)
Combine_final<-merge(Combine1,IH_final_results,by='Metabolite',all.x=TRUE)

write.csv(Combine_final, "Combined_gsva_geeglm.csv", row.names = FALSE)
####select all sig pathway as top pathway####
top_pathway <- read_excel("top_pathway.xlsx")
top_pathway$RA_median <- apply(gsva_matrix[top_pathway$Metabolite, RA_sub$Sample], 1, median)
top_pathway$IAR_median <- apply(gsva_matrix[top_pathway$Metabolite, IAR_sub$Sample], 1, median)
top_pathway$HC_median <- apply(gsva_matrix[top_pathway$Metabolite, HC_sub$Sample], 1, median)

test<-as.data.frame(top_pathway[,c(1,20:22)])
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


