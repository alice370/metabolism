rm(list=ls())
library("dplyr")
library(tidyr)
library(ggplot2)
library(eoffice)
data<-read.csv("209sample_188metabolites_auto_log10.csv",row.names = 1)
clin<-read.csv('RA_baseline_clin.csv',check.names = F)
merge1<-as.data.frame(t(data))
merge1$sample<-row.names(merge1)
merge1<-merge(clin[,c(1,2,3,6,7)],merge1,by="sample")


ACPA_positive <- merge1 %>% filter(ACPA == "positive")
ACPA_negative <- merge1 %>% filter(ACPA == "negative")


significance_matrix <- matrix(0, nrow = nrow(data), ncol = 100)
rownames(significance_matrix) <- row.names(data)
coef_matrix <- matrix(0, nrow = nrow(data), ncol = 100)
rownames(coef_matrix) <- row.names(data)

n_iterations <- 100


n_sample <- 14

set.seed(123)  
metabolite_cols <- colnames(merge1)[6:193]
for (i in 1:n_iterations) {
  
  
  sampled_ACPA_pos <- ACPA_positive %>% sample_n(n_sample)

  sampled_data <- bind_rows(sampled_ACPA_pos, ACPA_negative)
  

  sampled_data$ACPA<- factor(sampled_data$ACPA, levels = c("negative", "positive"))
  

  for (met in metabolite_cols) {
    
   
    formula <- as.formula(paste0("`", met, "` ~ `ACPA` + `GENDER` + `age` + `CRP`"))
    
 
    model <- try(lm(formula, data = sampled_data), silent = TRUE)
    

    if(class(model) == "try-error") {
      next  
    }
    
 
    summary_model <- summary(model)
    
   
    coef_summary <- summary_model$coefficients
    if("ACPApositive" %in% rownames(coef_summary)) {
      p_value <- coef_summary["ACPApositive", "Pr(>|t|)"]
      coef<-coef_summary["ACPApositive", "Estimate"]
        significance_matrix[met, i] <- p_value
        coef_matrix[met, i] <- coef
      }
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
target_metabolites <- c("Dimethyl fumarate", "Phosphoenolpyruvic acid")

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
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1),        
    axis.ticks.length = unit(0.3, "cm") 
  )+scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
p
write.csv(significance_matrix, "significance_100_Iteration.csv", row.names = FALSE)
write.csv(coef_matrix, "coef_100_Iteration.csv", row.names = FALSE)
write.csv(result_df, "significance_count.csv", row.names = FALSE)
ggsave("significance-density.pdf", width = 6, height = 4)
