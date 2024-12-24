library(ggplot2)
library(ggrepel)
library(eoffice)
rm(list=ls())

data<-read.csv('364sample-188metabolites-log10-auto.csv',row.names = 1)
clin<-read.csv('sample_information.csv')
normality_results <- data.frame(
  Metabolite = row.names(data),
  P_Value = apply(data, 1, function(x) shapiro.test(x)$p.value)
)
write.csv(normality_results,"shapiro.test_metabolites.csv")




data <- as.data.frame(t(data))
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

write.csv(results_df, "metabolite_comparison_results_scaled.csv", row.names = FALSE)

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

write.csv(results_df, "metabolite_comparison_results_scaled.csv", row.names = FALSE)



RA_HC.top_10_indices <- order(data1$PValue_RA_HC)[1:10]
RA_HC.significant_points <- data.frame(
  estimate = data1$Coef_RA_HC[RA_HC.top_10_indices],
  p.value = data1$PValue_RA_HC[RA_HC.top_10_indices],
  metabolic = data1$Metabolite[RA_HC.top_10_indices]
)


PRA_HC.top_10_indices <- order(data1$PValue_HC_atRisk)[1:10]
PRA_HC.significant_points <- data.frame(
  estimate = data1$Coef_HC_atRisk[PRA_HC.top_10_indices],
  p.value = data1$PValue_HC_atRisk[PRA_HC.top_10_indices],
  metabolic = data1$Metabolite[PRA_HC.top_10_indices]
)


RA_PRA.top_10_indices <- order(data1$PValue_RA_atRisk)[1:10]
RA_PRA.significant_points <- data.frame(
  estimate = data1$Coef_RA_atRisk[RA_PRA.top_10_indices],
  p.value = data1$PValue_RA_atRisk[RA_PRA.top_10_indices],
  metabolic = data1$Metabolite[RA_PRA.top_10_indices]
)


RA_HC.gg<-ggplot(data1, aes(y = -log10(data1$PValue_RA_HC) , x = data1$Coef_RA_HC, color = data1$RA_HC.color_ra)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('positive' = '#DE3024', 'negative' = '#25108f', 'no' = 'grey')) + 
  geom_point(data = RA_HC.significant_points, aes(y = -log10(p.value), x = estimate), colour = "yellow", size = 4) +
  geom_hline(yintercept = -log10(0.01), linetype = 2, color = "grey", size = 1) +
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(size = 0.5),
    axis.line.y = element_line(size = 0.5), 
    axis.text = element_text(size = 12), 
    axis.ticks = element_line(size = 0.5), 
    legend.position = "none"  
  ) +
  labs(
    x = "coefficient",
    y = "-log10(p value)",
    color = "",
    size=8
  )

RA_HC.gg <- RA_HC.gg + geom_text_repel(
  data = RA_HC.significant_points,
  aes(y = -log10(p.value), x = estimate, label = metabolic),
  size = 4,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5 
)
RA_HC.gg
topptx(RA_HC.gg,filename = "RA_HC_coeff.pptx")


PRA_HC.gg<-ggplot(data1, aes(y = -log10(data1$PValue_HC_atRisk) , x = -data1$Coef_HC_atRisk, color = data1$PRA_HC.color_ra)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('positive' = '#DE3024', 'negative' = '#25108f', 'no' = 'grey')) + 
  geom_point(data = PRA_HC.significant_points, aes(y = -log10(p.value), x = -estimate), colour = "yellow", size = 4) +
  geom_hline(yintercept = -log10(0.01), linetype = 2, color = "grey", size = 1) +
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(size = 0.5), 
    axis.line.y = element_line(size = 0.5), 
    axis.text = element_text(size = 12), 
    axis.ticks = element_line(size = 0.5), 
    legend.position = "none"  
  ) +
  labs(
    x = "coefficient",
    y = "-log10(p value)",
    color = "",
    size=8
  )


PRA_HC.gg <- PRA_HC.gg + geom_text_repel(
  data = PRA_HC.significant_points,
  aes(y = -log10(p.value), x = -estimate, label = metabolic),
  size = 4,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5  
)
PRA_HC.gg
topptx(PRA_HC.gg,filename = "PRA_HC_coeff.pptx")


RA_PRA.gg<-ggplot(data1, aes(y = -log10(data1$PValue_RA_atRisk) , x = data1$Coef_RA_atRisk, color = data1$RA_PRA.color_ra)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('positive' = '#DE3024', 'negative' = '#25108f', 'no' = 'grey')) + 
  geom_point(data = RA_PRA.significant_points, aes(y = -log10(p.value), x = estimate), colour = "yellow", size = 4) +
  geom_hline(yintercept = -log10(0.01), linetype = 2, color = "grey", size = 1) +
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(size = 0.5), 
    axis.line.y = element_line(size = 0.5), 
    axis.text = element_text(size = 12),
    axis.ticks = element_line(size = 0.5), 
    legend.position = "none" 
  ) +
  labs(
    x = "coefficient",
    y = "-log10(p value)",
    color = "",
    size=8
  )

RA_PRA.gg <- RA_PRA.gg + geom_text_repel(
  data = RA_PRA.significant_points,
  aes(y = -log10(p.value), x = estimate, label = metabolic),
  size = 4,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5 
)
RA_PRA.gg
topptx(RA_PRA.gg,filename = "RA_PRA_coeff.pptx")

