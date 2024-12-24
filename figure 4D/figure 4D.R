rm(list=ls())
library("dplyr")
library(tidyr)
library(ggplot2)
library(eoffice)
library(ggrepel)
library(scales)
squash_axis <- function(from, to, factor) {
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  
  trans <- function(x) {    
    # Initialize a logical vector to identify NAs
    is_na <- is.na(x)
    
    # get indices for the relevant regions
    isq <- !is_na & x > from & x < to
    ito <- !is_na & x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    # Initialize a logical vector to identify NAs
    is_na <- is.na(x)
    
    # get indices for the relevant regions
    isq <- !is_na & x > from & x < from + (to - from)/factor
    ito <- !is_na & x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation and inverse transformation
  return(scales::trans_new("squash_axis", trans, inv))
}


data<-read.csv("209sample_188metabolites_auto_log10.csv",row.names = 1)
clin<-read.csv('RA_baseline_clin.csv',check.names = F)
test<-read.csv('disease activity class test.csv')


merge1<-as.data.frame(t(data))
merge1$sample<-row.names(merge1)
merge1<-merge(clin[,c(1,8:11,13)],merge1,by="sample")

# 创建一个空的B表格
sig <- data.frame(metabolites = character(), compare=character(),p = numeric(), FC = numeric(), stringsAsFactors = FALSE)

# 循环遍历A表格的每一行
for (i in 1:nrow(test)) {
  for (j in 3:8) {
    # 检查值是否小于0.05
    if (test[i, j] < 0.05) {
      # 将A表格中的对应值添加到B表格
      sig <- rbind(sig, data.frame(
        metabolites = test[i, 1],
        compare=colnames(test)[j],
        p = test[i, j],
        FC = test[i, j + 6]
      ))
    }
  }
}
table(sig$compare)

H<-merge1[which(merge1$`disease activity class` =="high disease activity"),]
M<-merge1[which(merge1$`disease activity class`=="moderate disease activity"),]
L<-merge1[which(merge1$`disease activity class`==" Low disease activity"),]
D<-merge1[which(merge1$`disease activity class`=="Clinical Remission" ),]

for (i in 1:nrow(sig)) {
  # 获取meta值（当前行第一列）
  meta <- sig[i, 1]
  
  # 分解compare列单元格字符串的最后两个字符
  compare_str <- sig[i, "compare"]
  A<- get(substr(compare_str, nchar(compare_str)-1,nchar(compare_str)-1))
  B<-get(substr(compare_str, nchar(compare_str),nchar(compare_str)))
  
  combined_df <- rbind(
      A[, c(3:6, which(colnames(A) == meta))],
      B[, c(3:6, which(colnames(B) == meta))]
    )
  

  for (j in 1:4) {

    cor_test <- cor.test(combined_df[, j], combined_df[, 5])
    

    var_name <- colnames(combined_df)[j]
    

    sig[i, paste0(var_name, "_correlation")] <- cor_test$estimate  
    sig[i, paste0(var_name, "_pvalue")] <- cor_test$p.value       
  }
} 
write.csv(sig,'correlation analysis.csv')


####prepare for plot####
cor_range<-range(sig[,c(5,7,9,11)])
sig_split <- split(sig, sig$compare)
for (name in names(sig_split)) {
  small_table <- sig_split[[name]]
  long_table <- small_table %>%
    pivot_longer(
      cols = 5:12,
      names_to = "var_type",
      values_to = "value"
    ) %>%
    mutate(
      var = sub("_.*", "", var_type),          # Extracting the first part of column name 'VAS_correlation' as 'VAS'
      metric_type = sub(".*_", "", var_type)   # Extracting 'pvalue' or 'correlation' from the second part
    ) %>%
    select(-var_type) %>% 
    pivot_wider(
      names_from = metric_type,  
      values_from = value      
    )
  long_table <- long_table %>%
    mutate(
      color = scales::col_numeric(
        palette = c("blue", "white", "red"),  
        domain = range(cor_range) 
      )(correlation)
    )
  
  
  top5_metabolites <- long_table %>%
    filter(pvalue < 0.05) %>% 
    arrange(pvalue) %>%        
    slice(1:5)     
  
  
  plot <-long_table %>%
    ggplot(aes(x = -log10(pvalue), y = var, color = correlation)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey")+
    geom_point(size=4) +
    scale_color_gradient2(
      low = "#24108e", mid = "white", high = "#de3024", midpoint = 0,
      limits =cor_range)+
    labs(
      title = name,
      x = "-log10(pvalue)",
      y = "",
      color = "Correlation"
    ) +
    theme_minimal()+
    geom_text_repel(
      data = top5_metabolites, 
      aes(label = metabolites), 
      size = 3.5, 
      box.padding = 0.2,
      point.padding = 0.3,color='black'
    ) 
assign(paste0(name,'_p'),plot)
}
library(patchwork)
combined_plot <- (P_DH_p | P_DM_p | P_DL_p) / (P_LH_p | P_LM_p | P_MH_p)
combined_plot
ggsave("combined_plot.pdf", plot = combined_plot, width = 12, height = 4)
topptx(plot,'plot_legend.pptx')
