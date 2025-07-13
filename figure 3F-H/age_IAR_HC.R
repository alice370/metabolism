library(mgcv)
library(readxl)
library(readr)
library(gratia)
library(tidyr)
library(dplyr)
library(eoffice)
library(ggplot2)
rm(list=ls())
meta_dt <- read_csv("364sample_185metabolites_mean_log10.csv")
meta_dt<-as.data.frame(meta_dt)
row.names(meta_dt)<-meta_dt[,1]
meta_dt<-meta_dt[,-1]
meta_dt<-as.data.frame(t(meta_dt))
meta_dt$Sample<-row.names(meta_dt)
sample_information <- read_csv("sample_information.csv")
table(sample_information$Group)
IAR<-sample_information$Sample[sample_information$Group=='At-risk of RA']
HC<-sample_information$Sample[sample_information$Group=='Health']
####IAR####
data<-meta_dt[row.names(meta_dt)%in%IAR,]
data<-merge(sample_information[,c(1,4)],data,by='Sample',all.y  = TRUE)

original_names <- colnames(data)[3:187] 

colnames(data)[3:187] <- paste0("meta", 1:185)

metabolite_list<-colnames(data)[3:187]
results_list <- list()

for (meta in metabolite_list) {
  formula_smooth <- as.formula(paste0("`", meta, "` ~ s(Age)"))
  formula_linear <- as.formula(paste0("`", meta, "` ~ Age "))
  
 
  data_sub <- na.omit(data[, c(meta, "Age")])
  

  fit <- tryCatch({
    gam(formula_smooth, data = data_sub, method = "REML")
  }, error = function(e) NULL)
  
 
  fit_linear <- tryCatch({
    gam(formula_linear, data = data_sub, method = "REML")
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


results_data <- do.call(rbind, results_list)

results_data$original_name <- original_names[match(results_data$Metabolite, metabolite_list)]
write.csv(results_data,'IAR_results_data_all2.csv')
####Health####
data<-meta_dt[row.names(meta_dt)%in%HC,]
data<-merge(sample_information[,c(1,4)],data,by='Sample',all.y  = TRUE)


original_names <- colnames(data)[3:187]  

colnames(data)[3:187] <- paste0("meta", 1:185)

metabolite_list<-colnames(data)[3:187]
results_list <- list()

for (meta in metabolite_list) {
 
  formula_smooth <- as.formula(paste0("`", meta, "` ~ s(Age)"))
  formula_linear <- as.formula(paste0("`", meta, "` ~ Age "))
  

  data_sub <- na.omit(data[, c(meta, "Age")])
  

  fit <- tryCatch({
    gam(formula_smooth, data = data_sub, method = "REML")
  }, error = function(e) NULL)

  fit_linear <- tryCatch({
    gam(formula_linear, data = data_sub, method = "REML")
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

results_data <- do.call(rbind, results_list)

results_data$original_name <- original_names[match(results_data$Metabolite, metabolite_list)]
write.csv(results_data,'HC_results_data_all.csv')
####plot of lm####
rm(list=ls())

result_age<-read.csv('HC_results_data_all.csv',row.names = 1)
data1<-result_age[result_age$suitble.for.lm==1,]
data1$age_lm <- ifelse(data1$p_age < 0.05 & data1$Estimate_age > 0, "positive",
                       ifelse(data1$p_age < 0.05 & data1$Estimate_age < 0, "negative", "n.s."))

top_10_indices <- order(data1$p_age)[1:10]
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
ggsave("HC_age_lm_plot.pdf", plot = gg, height = 4.5, width = 4.5)
