data<-read.csv("209sample_188metabolites_auto_log10.csv",row.names = 1)
clin<-read.csv('RA_baseline_clin.csv')
F<-clin$sample[clin$GENDER=='F']
M<-clin$sample[clin$GENDER=='M']
result_gender<-data.frame()
for (i in row.names(data)) {
  group_F <- as.numeric(as.character(data[i, F]))
  group_M <- as.numeric(as.character(data[i, M]))
  result <- wilcox.test(group_F, group_M)
  p <- result$p.value
  log10.fc <- mean(group_F)-mean(group_M)
  result_row <- data.frame(
    metabolic = i,
    p.value = p,
    log10.fc = log10.fc
  )
  result_gender <- rbind(result_gender, result_row)
}
write.csv(result_gender,'gender-wilcox.csv')

data1<-result_gender
data1$color<- "no"
data1$color[(data1$p.value<0.05)&(data1$log10.fc>log10(1.2))] <- "up"
data1$color[(data1$p.value<0.05)&(data1$log10.fc<(-log10(1.2)))] <-"down"
table(data1$color)
top_10_indices <- order(data1$p.value)[1:10]
significant_points <- data.frame(
  log10.fc = data1$log10.fc[top_10_indices],
  p.value = -log10(data1$p.value[top_10_indices]),
  metabolic = data1$metabolic[top_10_indices]
)
library(ggplot2)
library(ggrepel)
library(eoffice)
gg <- ggplot(data1, aes(x = log10.fc, y = -log10(data1$p.value), color = data1$color)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('up' = '#DE3024', 'down' = '#009ACD', 'no' = 'grey')) + 
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
    x = "log10(fold.change)",
    y = "-log10(p.value)",
    color = "",
    size=8
  )
gg

gg <- gg + geom_text_repel(
  data = significant_points,
  aes(x = log10.fc, y = p.value, label = metabolic),
  size = 2,
  color = "black",
  box.padding = unit(0.4, "lines"),
  segment.color = "black",
  segment.size = 0.6,
  min.segment.length = 0.5 
)
gg
topptx(gg,filename = "gender_wixon_plot.pptx")
