library(GGally)
library(rstatix)
library(dplyr)
df<-read.csv('QC13_ALL.csv',row.names = 1)
p1<-ggpairs(
  df[,-1],
  aes(color = df$group),
  upper = list(continuous = wrap("cor", size = 3)), 
  lower = list(continuous = wrap("points", alpha = 0.5, size = 1)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5)) 
)
a<-cor(df[,-1], method = "pearson",use="pairwise" )
mean(a)
d<-range(a)
p2<-pheatmap::pheatmap(a,cluster_cols = FALSE,cluster_rows = F,scale = "row",
                   color = c(colorRampPalette(colors = c("#ffff00","#ff8050"))(100)),
                   legend_breaks=seq(d[1],d[2]), border=FALSE)
####p1 and p2 were combined using AI####