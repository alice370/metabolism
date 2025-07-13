setwd("C:/Users/alice/Desktop/实验数据/代谢组/数据处理/Code/figure 3I-J")
rm(list=ls())
library("dplyr")
library(tidyr)
library(ggplot2)
library(eoffice)
library(geepack)  
library(readxl)
data <- read.csv("364sample-185metabolites-log10-auto.csv",row.names = 1)
clin <- read_excel("RA information.xlsx")
merge1<-as.data.frame(t(data))
merge1$Sample<-row.names(merge1)
merge1<-merge(clin[,c(1,3,4,5:9,14,22)],merge1,by="Sample",all.x = TRUE)
colnames(merge1)[9]<-'ACPA'
merge1<-merge1[!merge1$ACPA=='NA',]
table(merge1$ACPA)
merge1$ACPA <- as.numeric(ifelse(merge1$ACPA == "positive", 1,
                                 ifelse(merge1$ACPA == "negative", 0, merge1$ACPA)))
merge1$Gender <- as.numeric( ifelse(merge1$Gender == "Female", 1,
                                    ifelse(merge1$Gender == "Male", 0, merge1$Gender)))
table(merge1$Gender)
ACPA_positive <- merge1 %>% filter(`ACPA` == "1")
ACPA_negative <- merge1 %>% filter(`ACPA`== "0")


significance_matrix <- matrix(0, nrow = nrow(data), ncol = 100)
rownames(significance_matrix) <- row.names(data)
coef_matrix <- matrix(0, nrow = nrow(data), ncol = 100)
rownames(coef_matrix) <- row.names(data)

n_iterations <- 100


n_sample <- 14

set.seed(123)  
metabolite_cols <- colnames(merge1)[11:195]

for (i in 1:n_iterations) {
  
  set.seed(i)
  sampled_ACPA_pos <- ACPA_positive %>% sample_n(n_sample)
  
  sampled_data <- bind_rows(sampled_ACPA_pos, ACPA_negative)
  
  
  sampled_data$ACPA<- factor(sampled_data$ACPA, levels = c("0", "1"))
  sampled_data$ID<-c(1:nrow(sampled_data))
  meta<-metabolite_cols[1]
  for (meta in metabolite_cols) {
    
    
    formula <- as.formula(paste0("`", meta, "`~ ACPA + Gender + Age + CRP"))
    
    model <- geeglm(formula, data = sampled_data, id = ID, family = gaussian, corstr = "exchangeable")
    coef_summary <- summary(model)$coefficients
    
    p_value <- coef_summary[2, "Pr(>|W|)"]
    coef<-coef_summary[2, "Estimate"]
    significance_matrix[meta, i] <- p_value
    coef_matrix[meta, i] <- coef
  }
  # print progress#
  if(i %% 10 == 0) {
    cat("complete", i, "iteration\n")
  }
}

significant_counts <- apply(significance_matrix, 1, function(x) sum(x < 0.05))

result_df <- data.frame(
  Metabolite = rownames(significance_matrix),
  Significant_Count = significant_counts
)

significance_df <- as.data.frame(significance_matrix) %>%
  mutate(Metabolite = rownames(significance_matrix))

significance_long <- significance_df %>%
  pivot_longer(
    cols = -Metabolite,      
    names_to = "Iteration",  
    values_to = "P_Value"     
  )
target_metabolites <- c("Fumarate", "Phosphoenolpyruvic acid",'Argininosuccinic acid','Xanthosine')

significance_long <- significance_long %>%
  mutate(Color = ifelse(Metabolite %in% target_metabolites,
                        "target",
                        "Other"))

median_p_values <- significance_long %>%
  filter(Metabolite %in% target_metabolites) %>%
  group_by(Metabolite) %>%
  summarize(Median_P = median(P_Value, na.rm = TRUE))

color_mapping<-c("target"="#8b0000","Other"="lightgrey")

p<-ggplot(significance_long, aes(x = -log10(P_Value), group = Metabolite, fill = Color)) +
  
  geom_density(data = subset(significance_long, !Metabolite %in% target_metabolites), alpha = 0.3,color='transparent') +
  
  geom_density(data = subset(significance_long, Metabolite %in% target_metabolites), alpha = 0.3,color='transparent') +
  
  geom_vline(xintercept = -log10(0.05), color = "grey", linetype = "dashed") +
  
  geom_vline(data = median_p_values,
             aes(xintercept = -log10(Median_P)),
             linetype = "dashed",
             color = "#5e5ea1") +
  
  scale_fill_manual(values = color_mapping) +
  labs(title = " ",
       x = "-log10(p value)",
       y = "Density") +
  theme_classic() +
  theme(
    legend.position = 'none',
    axis.line = element_line(size = 0.75),
    axis.ticks = element_line(size = 0.75),        
    axis.ticks.length = unit(0.2, "cm") 
  )+scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
p
write.csv(significance_matrix, "significance_100_Iteration.csv", row.names = FALSE)
write.csv(coef_matrix, "coef_100_Iteration.csv", row.names = FALSE)
result_df$p_median<- apply(significance_matrix[result_df$Metabolite, ], 1, median)
write.csv(result_df, "significance_count.csv", row.names = FALSE)
ggsave("significance-density.pdf", width = 6, height = 4)

####violin####
rm(list=ls())
library("dplyr")
library(tidyr)
library(ggplot2)
library(eoffice)
library(ggpubr)
data<-read.csv("364sample-185metabolites-log10-auto.csv",row.names = 1)
clin <- read_excel("RA information.xlsx")
sig<-read.csv('significance_count.csv')
meta<-sig$Metabolite[sig$Significant_Count>50]
merge1 <- as.data.frame(t(data[meta, clin$Sample[clin$`ACPA status`=='positive'|clin$`ACPA status`=='negative']]))

merge1$Sample<-row.names(merge1)
merge1<-merge(clin[,c(1,14)],merge1,by="Sample")
custom_order <- c("Fumarate", "Phosphoenolpyruvic acid", "Argininosuccinic acid", "Xanthosine")

row.names(merge1)<-merge1$Sample
merge1<-merge1[,-1]
long_data <- merge1 %>%
  pivot_longer(-`ACPA status`, names_to = "metabolites", values_to = "Value")

long_data$metabolites <- factor(long_data$metabolites, levels = custom_order)
p <- ggplot(long_data, aes(x = `ACPA status`, y = Value, fill = `ACPA status`)) + 
  geom_violin(trim = TRUE, alpha = 0.8) +  
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("#a1b8e1", "#df7f7f")) +  
  theme_classic2() + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.line.x = element_line(linewidth = 0.25), 
        axis.line.y = element_line(linewidth = 0.25),
        axis.ticks = element_line(size = 0.25),
        axis.text  = element_text(color = 'black')) +  
  labs(x = "", y = "") +  
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("positive", "negative")),  
                     label.y = max(long_data$Value, na.rm = TRUE) + 0.5) +  
  facet_wrap(~metabolites,nrow = 1,scales = "free_y")  
p
topptx(p,'ACPA-4meta-violin.pptx')


