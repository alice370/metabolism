setwd("C:/Users/alice/Desktop/实验数据/代谢组/数据处理/Code/figure 3I-J")
rm(list=ls())
library("dplyr")
library(tidyr)
library(ggplot2)
library(eoffice)
library(ggpubr)
library(readxl)
clin <- read_excel("RA information.xlsx")
origin <- read.csv("364sample_185metabolites.csv",check.names = F,row.names = 1)

F <- clin$Sample[clin$`ACPA status` == 'positive']
M <- clin$Sample[clin$`ACPA status` == 'negative']

result_ACPA <- data.frame()
for (metabolite in row.names(origin)) {
  group_F <- as.numeric(origin[metabolite, F])
  group_M <- as.numeric(origin[metabolite, M])
  # NA值统计
  NA_F <- sum(is.na(group_F))
  NA_M <- sum(is.na(group_M))
  nonNA_F <- length(group_F) - NA_F
  nonNA_M <- length(group_M) - NA_M
  
  contingency_table <- matrix(c(nonNA_F, NA_F, nonNA_M, NA_M), nrow = 2, 
                              dimnames = list(Status = c("NonNA", "NA"), Gender = c("F", "M")))
  
  # 检验方法选择
  if (any(contingency_table < 5)) {
    test_result <- fisher.test(contingency_table)
    test_used <- "Fisher"
  } else {
    test_result <- chisq.test(contingency_table)
    test_used <- "Chi-square"
  }
  
  result_row <- data.frame(
    metabolic = metabolite,
    N_P_nonNA = nonNA_F,
    N_N_nonNA = nonNA_M,
    detected_P = nonNA_F/length(group_F),
    detected_N = nonNA_M/length(group_M),
    NA_test_pvalue = test_result$p.value,
    NA_test_used = test_used
  )
  
  result_ACPA <- rbind(result_ACPA, result_row)
}
write.csv(result_ACPA, 'ACPA-new-frequency.csv', row.names = FALSE)
