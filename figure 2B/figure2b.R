library(mixOmics)
library(ggplot2)
library(ggrepel)
library(eoffice)
metabolite_data <- read.csv("364sample_185metabolites_mean_log10.csv", row.names = 1)
clin<-read.csv('clin_new.csv')
metabolite_data<-metabolite_data[,clin$Sample]
group_vector<- clin$Group[match(colnames(metabolite_data),clin$Sample)]
dt<-as.data.frame(t(metabolite_data))

plsda_result <- plsda(dt, group_vector,ncomp = 10)

df <- unclass(plsda_result)

df1 = as.data.frame(df$variates$X)
df1$samples = rownames(df1)
df1$group=group_vector

col=c("Health"="royalblue","At-risk of RA"="#FFC24B","RA"='firebrick2')

explain = df$prop_expl_var$X
x_lable <- round(explain[1],digits=3)
y_lable <- round(explain[2],digits=3)

p1 <- ggplot(df1, aes(x = comp1, y = comp2, color = group)) +
  theme_bw() +
  scale_color_manual(values = col) +
  geom_point(size = 1.8) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 0, lty = "dashed") +
  geom_hline(yintercept = 0, lty = "dashed") +
  guides(color = guide_legend(title = NULL)) +
  labs(x = paste0("t(1) (", x_lable * 100, "%)"),
       y = paste0("t(2) (", y_lable * 100, "%)")) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "none") 

p1<-p1+stat_ellipse(data=df1,
                    geom = "polygon",level = 0.95,
                    linetype = 2,size=0.5,
                    aes(fill=group),
                    alpha=0.2,
                    show.legend = T)+
  scale_fill_manual(values = col)
p1
write.csv(df1,'3group-PLSDA-score.csv')
topptx(p1,filename = "RA-PRA-HC-PLSDA.pptx")