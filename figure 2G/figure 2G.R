library(ggplot2)
library(tidyr)
library(dplyr)
library(eoffice)
rm(list=ls())
results_df<-read.csv("metabolite_comparison_results_scaled.csv")
results_df<-results_df
str(results_df)
results_df$RA_HC.color_ra<- "no"
results_df$RA_HC.color_ra[(results_df$PValue_RA_HC<0.01)&(results_df$Coef_RA_HC>0)] <- "positive"
results_df$RA_HC.color_ra[(results_df$PValue_RA_HC<0.01)&(results_df$Coef_RA_HC<0)]  <- "negative"
table(results_df$RA_HC.color_ra)

results_df$PRA_HC.color_ra<- "no"
results_df$PRA_HC.color_ra[(results_df$PValue_HC_atRisk<0.01)&(results_df$Coef_HC_atRisk<0)] <- "positive"
results_df$PRA_HC.color_ra[(results_df$PValue_HC_atRisk<0.01)&(results_df$Coef_HC_atRisk>0)]  <- "negative"
table(results_df$PRA_HC.color_ra)

results_df$RA_PRA.color_ra<- "no"
results_df$RA_PRA.color_ra[(results_df$PValue_RA_atRisk<0.01)&(results_df$Coef_RA_atRisk>0)] <- "positive"
results_df$RA_PRA.color_ra[(results_df$PValue_RA_atRisk<0.01)&(results_df$Coef_RA_atRisk<0)]  <- "negative"
table(results_df$RA_PRA.color_ra)

up<-results_df$Metabolite[results_df$RA_HC.color_ra=='positive'&results_df$PRA_HC.color_ra=='positive'&results_df$RA_PRA.color_ra=='positive']
down<-results_df$Metabolite[results_df$RA_HC.color_ra=='negative'&results_df$PRA_HC.color_ra=='negative'&results_df$RA_PRA.color_ra=='negative']
selected_results_df <- results_df %>% filter(Metabolite %in% c(up,down))

long_df <- selected_results_df[,-c(14:16)] %>%
  pivot_longer(
    cols = -Metabolite,
    names_to = c("Variable", "Comparison"),
    names_pattern = "([^_]+)_(RA_HC|RA_atRisk|HC_atRisk)"
  ) %>%
  pivot_wider(
    names_from = "Variable",
    values_from = "value"
  ) %>%
  mutate(
    Comparison = case_when(
      Comparison == "RA_HC" ~ "RA vs HC",
      Comparison == "RA_atRisk" ~ "RA vs At-risk",
      Comparison == "HC_atRisk" ~ "HC vs At-risk",
      TRUE ~ Comparison
    )
  ) 

long_df <- long_df %>%
  mutate(
   
    ToInvert = Comparison == "HC vs At-risk",
    

    Coef = if_else(ToInvert, -Coef, Coef),
    Upper = if_else(ToInvert, -Upper, Upper),
    Lower = if_else(ToInvert, -Lower, Lower),
 
    Comparison = if_else(ToInvert, "At-risk vs HC", Comparison)
  ) %>%

  select(-ToInvert)

long_df <- long_df %>%
  mutate(
    Significance = case_when(
      PValue < 0.001 ~ "***",
      PValue < 0.01 ~ "**",
      PValue < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

long_df$Comparison <- factor(long_df$Comparison, levels = rev(c("RA vs HC", "At-risk vs HC", "RA vs At-risk")))

p <- ggplot(long_df[long_df$Metabolite!='NAD',], aes(y = Comparison)) +
  geom_vline(xintercept = 0,linetype = "dashed", color = "grey") +
  geom_hline(yintercept = long_df$Comparison,linetype = "dashed", color = "grey")+
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, color = "black")+
  geom_point(aes(x = Coef), size = 1.2, color = "#2c5155") +
  theme_classic() +
  labs(
    x = "Coefficients",
    y = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, family  = "sans"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.ticks.y = element_blank(), 
    axis.line.y = element_blank(),
    legend.position = "none",
    strip.background = element_blank()
  ) +

  xlim(-0.1, 2.5)+
  facet_wrap(~ Metabolite, scales = "free_x")

print(p)
topptx(p,'forest_plot_9_up.pptx')
