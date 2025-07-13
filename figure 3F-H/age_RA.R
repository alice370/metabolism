library(mgcv)
library(readxl)
library(readr)
library(gratia)
library(tidyr)
library(dplyr)
library(eoffice)
library(ggplot2)
rm(list=ls())
clin <- read_excel("209RA-age.xlsx")
meta_dt <- read_csv("364sample_185metabolites_mean_log10.csv")
meta_dt<-as.data.frame(meta_dt)
row.names(meta_dt)<-meta_dt[,1]
meta_dt<-meta_dt[,-1]
meta_dt<-as.data.frame(t(meta_dt))
meta_dt$Sample<-row.names(meta_dt)
data<-merge(clin[,c(1,4,12)],meta_dt,by='Sample',all.x = TRUE)
df<-na.omit(data)
colnames(df)[3]<-'duration'

original_names <- colnames(df)[4:188]  


colnames(df)[4:188] <- paste0("meta", 1:185)

metabolite_list<-colnames(df)[4:188]
results_list <- list()

for (meta in metabolite_list) {

  formula_smooth <- as.formula(paste0("`", meta, "` ~ s(Age) + s(duration)"))
  formula_linear <- as.formula(paste0("`", meta, "` ~ Age + s(duration)"))
  

  df_sub <- na.omit(df[, c(meta, "Age", "duration")])

  fit <- tryCatch({
    gam(formula_smooth, data = df_sub, method = "REML")
  }, error = function(e) NULL)

  fit_linear <- tryCatch({
    gam(formula_linear, data = df_sub, method = "REML")
  }, error = function(e) NULL)
  
  if (!is.null(fit) && !is.null(fit_linear)) {
    s_tab <- summary(fit)$s.table
    edf <- s_tab["s(Age)", "edf"]
    p_smooth <- s_tab["s(Age)", "p-value"]
    
    aic_smooth <- AIC(fit)
    aic_linear <- AIC(fit_linear)
    
    p_tab <- summary(fit_linear)$p.table
    if ("Age" %in% rownames(p_tab)) {
      est_age <- p_tab["Age", "Estimate"]
      se_age <- p_tab["Age", "Std. Error"]
      p_age <- p_tab["Age", "Pr(>|t|)"]
    } else {
      est_age <- se_age <- p_age <- NA
    }
  } else {
    edf <- p_smooth <- aic_smooth <- aic_linear <- est_age <- se_age <- p_age <- NA
  }

  results_list[[meta]] <- data.frame(
    Metabolite = meta,
    edf = edf,
    p_smooth = p_smooth,
    AIC_smooth = aic_smooth,
    AIC_linear = aic_linear,
    Estimate_age = est_age,
    StdErr_age = se_age,
    p_age = p_age,
    stringsAsFactors = FALSE
  )
}


results_df <- do.call(rbind, results_list)

results_df$original_name <- original_names[match(results_df$Metabolite, metabolite_list)]

results_df_notlm <- read.csv("results_df_notlm.csv",row.names = 1)

meta_non_linear <- results_df_notlm$Metabolite
age_grid <- seq(16, 77, by = 1)
deriv_matrix <- matrix(NA, nrow = length(meta_non_linear), ncol = length(age_grid))
rownames(deriv_matrix) <- results_df_notlm$original_name
colnames(deriv_matrix) <- age_grid

trajectory_matrix <- matrix(NA, nrow = length(meta_non_linear), ncol = length(age_grid))
rownames(trajectory_matrix) <- results_df_notlm$original_name  
colnames(trajectory_matrix) <- age_grid  
for (i in seq_along(meta_non_linear)) {
  meta <- meta_non_linear[i]
  formula <- as.formula(paste0("`", meta, "` ~ s(Age) + s(duration)"))
  fit <- gam(formula, data = df, method = "REML")
  
  newdata <- data.frame(Age = age_grid, duration = median(df$duration, na.rm = TRUE))
  

  deriv_df <- derivatives(fit, select = "s(Age)", data = newdata)
  

  deriv_matrix[i, ] <- deriv_df$.derivative
  
  pred_vals <- predict(fit, newdata = newdata)
  trajectory_matrix[i, ] <- pred_vals
}
deriv_matrix<-as.data.frame(deriv_matrix)
trajectory_matrix<-as.data.frame(trajectory_matrix)

dist_deriv <- dist(scale(deriv_matrix))
clust_deriv <- hclust(dist_deriv, method = "ward.D2")
plot(clust_deriv)

cluster_assignments <- cutree(clust_deriv, k = 4) 
trajectory_matrix$Cluster <- factor(cluster_assignments[rownames(trajectory_matrix)])
####output trajectory_matrix and dist_deriv####
####plot clusters of trajectory####
trajectory_fitted_clustered <- read_csv("trajectory_fitted_clustered.csv")
deriv_matrix <- read.csv("deriv_matrix.csv",row.names = 1,check.names = F)
results_df_notlm <- read_csv("results_df_notlm.csv")
plot_df <- trajectory_fitted_clustered
colnames(plot_df)[1]<-'Metabolite'
# 转为 tidy 格式
plot_df_long <- plot_df %>%
  pivot_longer(
    cols = -c(Metabolite, Cluster),  
    names_to = "Age",
    values_to = "Fitted"
  ) %>%
  mutate(Age = as.numeric(Age))

cluster_colors <- c("1" = "#E64B35FF", 
                    "4" = "#4DBBD5FF", 
                    "3" = "#00A087FF", 
                    "2" = "#3C5488FF")
age_grid <- seq(16, 77, by = 1)
####plot####

cluster_metas <- unique(plot_df_long$Metabolite[plot_df_long$Cluster == k])

min_age_points <- apply(deriv_matrix[cluster_metas, , drop=FALSE], 1, function(x) age_grid[which.min(abs(x))])
age_range <- quantile(min_age_points, probs = c(0.25, 0.75))
age_median <- median(min_age_points)

n_meta <- length(cluster_metas)

y_range <- range(plot_df_long$Fitted[plot_df_long$Cluster == k], na.rm = TRUE)

top5_names <- results_df_notlm %>%
  filter(original_name %in% cluster_metas) %>%
  arrange(p_smooth) %>%
  slice_head(n = 5) %>%
  pull(original_name)

label_data <- data.frame(
  Metabolite = top5_names,
  Age = min(plot_df_long$Age, na.rm = TRUE),  
  Fitted = seq(min(y_range) + 0.05 * diff(y_range),
               min(y_range) + 0.30 * diff(y_range),
               length.out = length(top5_names)) 
)
color_k <- cluster_colors[as.character(4)]
  
p <- ggplot(plot_df_long %>% filter(Cluster == k),
              aes(x = Age, y = Fitted, group = Metabolite)) +
    geom_line(alpha = 0.3, color = "grey50") +
    stat_summary(fun = mean, geom = "line", color = color_k, size = 1.2) +
    geom_vline(xintercept = age_range, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = age_median, linetype = "dashed", color = "grey20") +
    geom_text(data = label_data,
            aes(x = Age, y = Fitted, label = Metabolite),
            hjust = 0,    
            vjust = 0,   
            size = 5,
            color = "black") +
    
    labs(
      title = paste0("Cluster ", k,"(n=",n_meta,")"),
      x = "Age",
      y = "Fitted GAM value"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(color = "black", size = 16, hjust = 0.5),
      axis.title = element_text(color = "black", size = 16),
      axis.ticks  = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = 14),
      strip.text = element_text(color = "black", size = 16),
      legend.text = element_text(color = "black", size = 14),
      legend.title = element_text(color = "black", size = 16)
    ) +
annotate("text",
         x = age_range[1]-2,
         y = max(y_range, na.rm = TRUE) - 0.05 * diff(y_range),
         label = paste0( round(age_range[1])),
         angle = 0, vjust = -0.5, size = 5, color = "black") +
  
  annotate("text",
           x = age_range[2]-2,
           y = max(y_range, na.rm = TRUE) - 0.05 * diff(y_range),
           label = paste0( round(age_range[2])),
           angle = 0, vjust = -0.5, size = 5, color = "black") +
  
  annotate("text",
           x = age_median-2,
           y = max(y_range, na.rm = TRUE) - 0.05 * diff(y_range),
           label = paste0(round(age_median)),
           angle = 0, vjust = -0.5, size = 5, color = "black")+ 
  scale_x_continuous(
    limits = c(16, 77),
    breaks = seq(20, 70, by = 10)
  )
p
ggsave(paste0("cluster", k, ".pdf"), plot = p, height = 4, width = 5)

####p2 for special####
p2 <- ggplot(plot_df_long %>% filter(Cluster == k),
             aes(x = Age, y = Fitted, group = Metabolite)) +
  geom_line(alpha = 0.3, color = "grey50") +
  stat_summary(fun = mean, geom = "line", color = color_k, size = 1.2) +

geom_text(data = label_data,
          aes(x = Age, y = Fitted, label = Metabolite),
          hjust = 0,    
          vjust = 0,    
          size = 5,
          color = "black") +
  
  labs(
    title = paste0("Cluster ", k,"(n=",n_meta,")"),
    x = "Age",
    y = "Fitted GAM value"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(color = "black", size = 16, hjust = 0.5),
    axis.title = element_text(color = "black", size = 16),
    axis.ticks  = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 14),
    strip.text = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14),
    legend.title = element_text(color = "black", size = 16)
  ) +
  scale_x_continuous(
    limits = c(16, 77),
    breaks = seq(20, 70, by = 10)
  )
p2
ggsave(paste0("cluster", k, ".pdf"), plot = p2, height = 4, width = 5)

####plot of lm####
rm(list=ls())

result_age<-read.csv('results_df_all.csv',row.names = 1)
data1<-result_age[result_age$suit.for.lm==1|result_age$edf<1.5,]
data1$age_lm <- ifelse(data1$p_age < 0.05 & data1$Estimate_age > 0, "positive",
                       ifelse(data1$p_age < 0.05 & data1$Estimate_age < 0, "negative", "n.s."))

top_10_indices <- order(data1$p_age)[1:9]
top_10_indices<-c(top_10_indices,which(data1$Metabolite=='meta22'))
significant_points <- data.frame(
  beta.lm = data1$Estimate_age[top_10_indices],
  p.value = -log10(data1$p_age[top_10_indices]),
  metabolic = data1$original_name[top_10_indices]
)
library(ggplot2)
library(ggrepel)
library(eoffice)
gg <- ggplot(data1, aes(x = Estimate_age, y = -log10(data1$p_age), color = data1$age_lm)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c('positive' = '#DE3024', 'negative' = '#009acd', 'n.s.' = 'grey')) + 
  labs(
    x = "Coefficients",
    y = "-Log10(p value)",
    color = "",
    size=8
  ) + geom_text_repel(
  data = significant_points,
  aes(x = beta.lm, y = p.value, label = metabolic),
  size = 5,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5 
)+
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(color = "black", size = 16),
    axis.ticks  = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 15),
    legend.position = 'none'
  )
gg
ggsave("age_lm_plot.pdf", plot = gg, height = 4.5, width = 4.5)
topptx(gg,filename = "age_lm_plot.pptx")

####intersect RH and age####
result_age<-read.csv('results_df_all.csv',row.names = 1)
data1<-result_age[result_age$suitble.for.lm==1,]
data1$age_lm <- ifelse(data1$p_age < 0.05 & data1$Estimate_age > 0, "positive",
                       ifelse(data1$p_age < 0.05 & data1$Estimate_age < 0, "negative", "n.s."))
table(data1$age_lm)
RH_GEEglm_meta <- read_csv("RH_GEEglm_meta.csv")
RH_pos<-RH_GEEglm_meta$Metabolite[RH_GEEglm_meta$Estimate>0&RH_GEEglm_meta$p_value<0.05]
common<-intersect(age_positive,RH_pos)

write.csv(common,'increased in age and RA.csv')



