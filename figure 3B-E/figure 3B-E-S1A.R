setwd("C:/Users/alice/Desktop/实验数据/代谢组/数据处理/Code/figure 3B-E-S1A")
data <- read.csv("364sample-185metabolites-log10-auto.csv", row.names = 1)
clin <- read.csv('RA_baseline_clin.csv')
rm(list=ls())
origin <- read.csv("364sample_185metabolites.csv",check.names = F,row.names = 1)
F <- clin$sample[clin$GENDER == 'F']
M <- clin$sample[clin$GENDER == 'M']

F <- intersect(F, colnames(data))
M <- intersect(M, colnames(data))


result_gender <- data.frame()

for (metabolite in row.names(data)) {
  group_F <- as.numeric(data[metabolite, F])
  group_M <- as.numeric(data[metabolite, M])
  result <- wilcox.test(group_F, group_M)
  p <- result$p.value
  log10.fc <- mean(group_F)-mean(group_M)
  
  group_F <- as.numeric(origin[metabolite, F])
  group_M <- as.numeric(origin[metabolite, M])

  NA_F <- sum(is.na(group_F))
  NA_M <- sum(is.na(group_M))
  nonNA_F <- length(group_F) - NA_F
  nonNA_M <- length(group_M) - NA_M
  
  contingency_table <- matrix(c(nonNA_F, NA_F, nonNA_M, NA_M), nrow = 2, 
                              dimnames = list(Status = c("NonNA", "NA"), Gender = c("F", "M")))
  

  if (any(contingency_table < 5)) {
    test_result <- fisher.test(contingency_table)
    test_used <- "Fisher"
  } else {
    test_result <- chisq.test(contingency_table)
    test_used <- "Chi-square"
  }
  
  result_row <- data.frame(
    metabolic = metabolite,
    p.value = p,
    log10.fc = log10.fc,
    N_F_nonNA = nonNA_F,
    N_M_nonNA = nonNA_M,
    detected_F = nonNA_F/length(group_F),
    detected_M = nonNA_M/length(group_M),
    NA_test_pvalue = test_result$p.value,
    NA_test_used = test_used
  )
  
  result_gender <- rbind(result_gender, result_row)
}

write.csv(result_gender, 'gender-wilcox-new.csv', row.names = FALSE)

####plot RA####
data1<-read.csv('gender-wilcox-new.csv')
data1$color<- "no"
data1$color[(data1$p.value<0.05)&(data1$log10.fc>log10(1.2))] <- "up"
data1$color[(data1$p.value<0.05)&(data1$log10.fc<(-log10(1.2)))] <-"down"
table(data1$color)
top_10_indices <- order(data1$p.value)[1:10]
significant_points <- data.frame(
  log10.fc = data1$log10.fc[top_10_indices],
  p.value = -log10(data1$p.value[top_10_indices]),
  metabolic = data1$metabolites[top_10_indices]
)
gg <- ggplot(data1, aes(x = -log10.fc, y = -log10(data1$p.value), color = data1$color)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('down' = '#DE3024', 'up' = '#009ACD', 'no' = 'grey')) + 
  geom_vline(xintercept = log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_vline(xintercept = -log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), #
    panel.grid.minor = element_blank(), # 
    axis.line.x = element_line(size = 0.5), #
    axis.line.y = element_line(size = 0.5), # 
    axis.text = element_text(size = 6), # 
    axis.ticks = element_line(size = 1),
    
    legend.position = "none"  
  ) +
  labs(
    x = "Log10(male/female)",
    y = "-Log10(p value)",
    color = "",
    size=8
  )+ geom_text_repel(
    data = significant_points,
    aes(x = -log10.fc, y = p.value, label = metabolic),
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
ggsave("RA_gender_wilcox.pdf", plot = gg, height = 4.5, width = 4.5)
####plot IAR####
data1<-read.csv('IAR-gender-wilcox.csv')
data1$color<- "no"
data1$color[(data1$p.value<0.05)&(data1$log10.fc>log10(1.2))] <- "up"
data1$color[(data1$p.value<0.05)&(data1$log10.fc<(-log10(1.2)))] <-"down"
table(data1$color)
top_10_indices <- order(data1$p.value)[1:10]
significant_points <- data.frame(
  log10.fc = data1$log10.fc[top_10_indices],
  p.value = -log10(data1$p.value[top_10_indices]),
  metabolic = data1$metabolites[top_10_indices]
)
gg <- ggplot(data1, aes(x = -log10.fc, y = -log10(data1$p.value), color = data1$color)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('down' = '#DE3024', 'up' = '#009ACD', 'no' = 'grey')) + 
  geom_vline(xintercept = log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_vline(xintercept = -log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), #
    panel.grid.minor = element_blank(), # 
    axis.line.x = element_line(size = 0.5), #
    axis.line.y = element_line(size = 0.5), # 
    axis.text = element_text(size = 6), # 
    axis.ticks = element_line(size = 1),
    
    legend.position = "none"  
  ) +
  labs(
    x = "Log10(male/female)",
    y = "-Log10(p value)",
    color = "",
    size=8
  )+ geom_text_repel(
    data = significant_points,
    aes(x = -log10.fc, y = p.value, label = metabolic),
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
ggsave("IAR_gender_wilcox.pdf", plot = gg, height = 4.5, width = 4.5)
####plot IAR####
data1<-read.csv('HC-gender-wilcox.csv')
data1$color<- "no"
data1$color[(data1$p.value<0.05)&(data1$log10.fc>log10(1.2))] <- "up"
data1$color[(data1$p.value<0.05)&(data1$log10.fc<(-log10(1.2)))] <-"down"
table(data1$color)
top_10_indices <- order(data1$p.value)[1:10]
significant_points <- data.frame(
  log10.fc = data1$log10.fc[top_10_indices],
  p.value = -log10(data1$p.value[top_10_indices]),
  metabolic = data1$metabolites[top_10_indices]
)
gg <- ggplot(data1, aes(x = -log10.fc, y = -log10(data1$p.value), color = data1$color)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('down' = '#DE3024', 'up' = '#009ACD', 'no' = 'grey')) + 
  geom_vline(xintercept = log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_vline(xintercept = -log10(1.2), linetype = 2, color = "grey", size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey", size = 1) +
  
  theme(
    panel.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.grid.major = element_blank(), #
    panel.grid.minor = element_blank(), # 
    axis.line.x = element_line(size = 0.5), #
    axis.line.y = element_line(size = 0.5), # 
    axis.text = element_text(size = 6), # 
    axis.ticks = element_line(size = 1),
    
    legend.position = "none"  
  ) +
  labs(
    x = "Log10(male/female)",
    y = "-Log10(p value)",
    color = "",
    size=8
  )+ geom_text_repel(
    data = significant_points,
    aes(x = -log10.fc, y = p.value, label = metabolic),
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
ggsave("HC_gender_wilcox.pdf", plot = gg, height = 4.5, width = 4.5)


library('ggvenn')
RA<-read.csv('gender-wilcox-new.csv')
RA$color<- "no"
RA$color[(RA$p.value<0.05)&(RA$log10.fc>log10(1.2))] <- "up"
RA$color[(RA$p.value<0.05)&(RA$log10.fc<(-log10(1.2)))] <-"down"
table(RA$color)
RA_negative <- RA$metabolites[RA$color == 'down']

IAR<-read.csv('IAR-gender-wilcox.csv')
IAR$color<- "no"
IAR$color[(IAR$p.value<0.05)&(IAR$log10.fc>log10(1.2))] <- "up"
IAR$color[(IAR$p.value<0.05)&(IAR$log10.fc<(-log10(1.2)))] <-"down"
table(IAR$color)
IAR_negative <- IAR$metabolites[IAR$color == 'down']

HC<-read.csv('HC-gender-wilcox.csv')
HC$color<- "no"
HC$color[(HC$p.value<0.05)&(HC$log10.fc>log10(1.2))] <- "up"
HC$color[(HC$p.value<0.05)&(HC$log10.fc<(-log10(1.2)))] <-"down"
table(HC$color)
HC_negative <- HC$metabolites[HC$color == 'down']
venn_list <- list(
  "RA" = RA_negative,
  "IAR" = IAR_negative,
  'HC'=HC_negative
)
library(purrr)
intersect(RA_negative,IAR_negative)
common_metabolites <- reduce(venn_list, intersect)
topptx(p1,filename = "intersect_lm_test.pptx")
####intersect RH and gender####
RH_GEEglm_meta <- read_csv("RH_GEEglm_meta.csv")
RH_pos<-RH_GEEglm_meta$Metabolite[RH_GEEglm_meta$Estimate>0&RH_GEEglm_meta$p_value<0.05]
common<-intersect(RA_negative,RH_pos)

write.csv(common,'increased in male and RA.csv')