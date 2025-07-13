rm(list=ls())
library("dplyr")
library(tidyr)
library(ggplot2)
library(eoffice)
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


data<-read.csv("IAR-185metabolites-log10-auto.csv",row.names = 1)
clin<-read.csv('at risk of RA_follow up.csv',check.names = F)
data<-data[,clin$Sample[clin$`Follow-up visits`=='yes']]
merge1<-as.data.frame(t(data))
merge1$Sample<-row.names(merge1)
merge1<-merge(clin[,c(1,5)],merge1,by="Sample")
merge1$`Developed RA`

converter <- merge1 %>% filter(`Developed RA` == "yes")
other <- merge1 %>% filter(`Developed RA` == "no")


significance_matrix <- matrix(0, nrow = nrow(data), ncol = 100)
rownames(significance_matrix) <- row.names(data)
####
FC_matrix <- matrix(0, nrow = nrow(data), ncol = 100)
rownames(FC_matrix) <- row.names(data)


n_iterations <- 100


n_sample <- 4

set.seed(123)
metabolite_cols <- colnames(merge1)[3:187]
for (i in 1:n_iterations) {
  
  
  sampled_other <-other %>% sample_n(n_sample)
  
  
  sampled_data <- bind_rows(sampled_other, converter)
  
  
  for (met in metabolite_cols) {
    
    
    converter_met <- sampled_data %>% filter(`Developed RA` == "yes") %>% pull(met)
    other_met <- sampled_data %>% filter(`Developed RA` == "no") %>% pull(met)
    
    
    test_result <- try(wilcox.test(converter_met, other_met, exact = FALSE), silent = TRUE)
    
    
    if(class(test_result) == "try-error") {
      next  
    }
    
    p_value <- test_result$p.value
    FC<-mean(converter_met)-mean(other_met)
    
    significance_matrix[met, i] <- p_value
    FC_matrix[met, i]<-FC
    
  }
  
  if(i %% 10 == 0) {
    cat("complete", i, "sampling")
  }
}

significant_counts <- apply(significance_matrix, 1, function(x) sum(x < 0.05))


result_df <- data.frame(
  Metabolite = rownames(significance_matrix),
  Significant_Count = significant_counts
)

significance_df <- as.data.frame(significance_matrix) %>%
  mutate(Metabolite = rownames(significance_matrix))


write.csv(significance_matrix, "significance_100_Iteration.csv", row.names = TRUE)
write.csv(FC_matrix, "FC_100_Iteration.csv", row.names = TRUE)
write.csv(result_df, "significance_count.csv", row.names = FALSE)
sig_meta<-result_df$Metabolite[result_df$Significant_Count>=10]
# Sort data by count in descending order
df <- result_df %>% arrange(desc(Significant_Count))
sig_meta<-df$Metabolite[df$Significant_Count>=10]
# Create color gradient: from orange to gold for counts >= 10, light gray for < 10
df$color <- ifelse(df$Significant_Count >= 10, 
                   scales::gradient_n_pal(c("orange", "gold"))(scales::rescale(df$Significant_Count[df$Significant_Count >= 10])), 
                   "lightgray")

nrow(df[df$Significant_Count>0,])
# Plot using ggplot2
p <- ggplot(df[df$Significant_Count>0,], aes(y = reorder(Metabolite, Significant_Count), x = Significant_Count, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +  # Use the colors directly
  geom_vline(xintercept = 10, linetype = "dashed", color = "gray") +  # Add horizontal dashed line at y = 10
  geom_text(aes(x=0,label = ifelse(Significant_Count >= 10, Metabolite, "")), 
            hjust = 0) +  
  labs(y = "", x = "Count", title = "") +
  theme_classic()+theme(axis.text.y = element_blank())+
  coord_trans(y = squash_axis(0, 35, 10))
print(p)
topptx(p,'significant_count.pptx')

result <- data.frame(Element = character(), pvalue=numeric(),FC = numeric(), stringsAsFactors = FALSE)

for (row_name in sig_meta) {
  row_index <- which(rownames(significance_matrix) == row_name)
  col_indices <- which(significance_matrix[row_index, ] < 0.05)
  
  for (col_index in col_indices) {
    result <- rbind(result, data.frame(Element = row_name, 
                                       pvalue=significance_matrix[row_index, col_index],
                                       Value = FC_matrix[row_index, col_index]))
  }
}
write.csv(result, "sig_meta_FC_in sig_result.csv", row.names = FALSE)
####plot####
result$Element <- factor(result$Element, levels = rev(unique(result$Element)))

plot <- ggplot(result, aes(x = Value, y = Element)) +
  geom_boxplot(aes(group = Element), width=0.5,outlier.shape = NA,lwd = 0.2,fatten= 0.2) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6,color='#154599') +
  theme_classic() +
  xlab("") +
  ylab("") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8b0000")+
  theme(axis.text = element_blank(),
        axis.line = element_line(linewidth = 0.35))
print(plot)
ggsave('sig_meta_FC.pdf',plot=plot,width=3,height=4)
