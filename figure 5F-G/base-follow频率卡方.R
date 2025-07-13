rm(list=ls())
origin <- read.csv("POS-NEG_all.csv",check.names = F,row.names = 1)
clin<-read.csv('sample_information.csv',check.names = F)
data<-read.csv('561sample_185meta_log10_auto.csv',row.names = 1)

MLYA<-clin$sample[which(clin$DRUG =="MTQ+LEF"&
                          clin$`response classification` %in% c('Good Response', 'Moderate Response'))]

MLYB<-clin$`follow up-sample`[which(clin$DRUG =="MTQ+LEF"&
                                      clin$`response classification` %in% c('Good Response', 'Moderate Response'))]

MLNA<-clin$sample[which(clin$DRUG =="MTQ+LEF"&
                          clin$`response classification` %in% c('No Response'))]

MLNB<-clin$`follow up-sample`[which(clin$DRUG =="MTQ+LEF"&
                                      clin$`response classification` %in% c('No Response'))]

MHYA<-clin$sample[which(clin$DRUG =="MTQ+HCQ"&
                          clin$`response classification` %in% c('Good Response', 'Moderate Response'))]

MHYB<-clin$`follow up-sample`[which(clin$DRUG =="MTQ+HCQ"&
                                      clin$`response classification` %in% c('Good Response', 'Moderate Response'))]

MHNA<-clin$sample[which(clin$DRUG =="MTQ+HCQ"&
                          clin$`response classification` %in% c('No Response'))]

MHNB<-clin$`follow up-sample`[which(clin$DRUG =="MTQ+HCQ"&
                                      clin$`response classification` %in% c('No Response'))]

result_group <- data.frame()

for (metabolite in row.names(data)) {
  comparisons <- list(
    list("MHN", MHNA, MHNB),
    list("MHY", MHYA, MHYB),
    list("MLN", MLNA, MLNB),
    list("MLY", MLYA, MLYB)
  )
  
  result_row <- data.frame(metabolic = metabolite)
  
  for (comp in comparisons) {
    name <- comp[[1]]
    group_A <- as.numeric(origin[origin$name==metabolite, comp[[2]]])
    group_B <- as.numeric(origin[origin$name==metabolite, comp[[3]]])
    
    nonNA_A <- sum(!is.na(group_A))
    nonNA_B <- sum(!is.na(group_B))
    NA_A <- sum(is.na(group_A))
    NA_B <- sum(is.na(group_B))
    
    contingency_table <- matrix(c(nonNA_A, NA_A, nonNA_B, NA_B),
                                nrow = 2,
                                dimnames = list(Status = c("NonNA", "NA"), 
                                                Group = c(paste0(name, "A"), paste0(name, "B"))))
    
    if (any(contingency_table < 5)) {
      test_result <- fisher.test(contingency_table)
      test_used <- "Fisher"
    } else {
      test_result <- chisq.test(contingency_table)
      test_used <- "Chi-square"
    }
    
    result_row[[paste0(name, "A_nonNA")]] <- nonNA_A
    result_row[[paste0(name, "A_detected")]] <- nonNA_A / length(group_A)
    result_row[[paste0(name, "B_nonNA")]] <- nonNA_B
    result_row[[paste0(name, "B_detected")]] <- nonNA_B / length(group_B)
    result_row[[paste0(name, "_pvalue")]] <- test_result$p.value
    result_row[[paste0(name, "_test_used")]] <- test_used
  }
  
  result_group <- rbind(result_group, result_row)
}

write.csv(result_group, 'base-followÆµÂÊ¿¨·½.csv', row.names = FALSE)
