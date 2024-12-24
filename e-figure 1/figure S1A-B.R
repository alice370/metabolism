library(dplyr)
library(ggplot2)
library(ggpubr)
library(eoffice)
rm(list = ls())
data=read.csv("sample_information.csv",check.names = F)
colnames(data)=data[1,];data=data[-1,];data[,7:ncol(data)] = lapply(data[,7:ncol(data)], as.numeric)
####gender-ML####
dat=data%>%filter(DRUG=="MTQ+LEF")
dat=dat[!is.na(dat$`response classification`),]
dat <- dat %>%
  mutate(`response classification` = case_when(
    `response classification` %in% c("Good Response", "Moderate Response") ~ "response",
    TRUE ~ `response classification`
  ))

gender=as.data.frame(table(dat$Gender,dat$`response classification`))
gender_table <- xtabs(Freq ~ Var1 + Var2, data = gender)
chisq_result <- chisq.test(gender_table)
chisq_result$p.value


p<-ggplot(gender, aes(Var1, Freq, fill=Var2))+
  geom_bar(stat="identity", position = position_dodge(.8), width = .6)+
  geom_text(aes(label=Freq), position = position_dodge(.8))+
  scale_fill_manual(values = c("#21908c","#fde725"))+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,40), breaks=c(0,20,40))
topptx(p,'gender-response-ML.pptx')
####gender-MH####
dat=data%>%filter(DRUG=="MTQ+HCQ")
dat=dat[!is.na(dat$`response classification`),]
dat <- dat %>%
  mutate(`response classification` = case_when(
    `response classification` %in% c("Good Response", "Moderate Response") ~ "response",
    TRUE ~ `response classification`
  ))
gender=as.data.frame(table(dat$Gender,dat$`response classification`))
gender_table <- xtabs(Freq ~ Var1 + Var2, data = gender)
chisq_result <- chisq.test(gender_table)
chisq_result$p.value


p<-ggplot(gender, aes(Var1, Freq, fill=Var2))+
  geom_bar(stat="identity", position = position_dodge(.8), width = .6)+
  geom_text(aes(label=Freq), position = position_dodge(.8))+
  scale_fill_manual(values = c("#21908c","#fde725"))+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,40), breaks=c(0,20,40))
topptx(p,'gender-response-MH.pptx')
####age-MH####
p<-ggplot(dat[dat$DRUG=="MTQ+HCQ",], aes(`response classification`, Age, fill=`response classification`))+
  geom_violin(width=.6, alpha=.3)+
  geom_boxplot(width = .25, outlier.colour = NA, alpha=1, fill="white")+
  geom_jitter(shape=21, width=.1 , color="transparent")+
  scale_fill_manual(values = c("#21908c","#fde725"), name="")+
  theme_classic()+
  stat_compare_means(method="wilcox.test", label.y = 90, size=5)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks=c(0,25,50,75,100))+
  theme(axis.text = element_text(color="black",size=15),
        axis.line = element_line(color="black"),
        axis.ticks =  element_line(color="black"),
        axis.title = element_text(color="black",size=15))
topptx(p,'age-response-MH.pptx')
####age-ML####
p<-ggplot(dat[dat$DRUG=="MTQ+LEF",], aes(`response classification`, Age, fill=`response classification`))+
  geom_violin(width=.6, alpha=.3)+
  geom_boxplot(width = .25, outlier.colour = NA, alpha=1, fill="white")+
  geom_jitter(shape=21, width=.1 , color="transparent")+
  scale_fill_manual(values = c("#21908c","#fde725"), name="")+
  theme_classic()+
  stat_compare_means(method="wilcox.test", label.y = 90, size=5)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks=c(0,25,50,75,100))+
  theme(axis.text = element_text(color="black",size=15),
        axis.line = element_line(color="black"),
        axis.ticks =  element_line(color="black"),
        axis.title = element_text(color="black",size=15))
topptx(p,'age-response-ML.pptx')
####DAS28-ML####
dat1<-dat[dat$DRUG=="MTQ+LEF",c(1,27,9:11,13)]
dat1[,c(3:6)]= lapply(dat1[,c(3:6)], scale)
dat2=reshape2::melt(dat1, id.vars=c(1,2), measure.vars=c(3:6))

p<-ggplot(dat2, aes(variable, value, color=`response classification`, fill=`response classification`))+
  geom_boxplot(outlier.colour = NA, width=.5, alpha=.01, position = position_dodge(width = .7),color="black")+
  geom_jitter(position = position_jitterdodge(jitter.width = .15, dodge.width = .7), alpha=1)+
  stat_compare_means(label="p.signif",method="wilcox.test", label.y = 5.2,size=5)+
  scale_color_manual(values = c("#21908c","#fde725"), name="")+ 
  scale_fill_manual(values = c("#21908c","#fde725"), name="")+
  theme_classic()+
  labs(x="",y="Scaled value")+
  theme(axis.text = element_text(color="black",size=15),
        axis.line = element_line(color="black"),
        axis.ticks =  element_line(color="black"),
        axis.title = element_text(color="black",size=15),
        legend.position="none")
p
topptx(p,'DAS28-ML.pptx')
####DAS28-MH####
dat1<-dat[dat$DRUG=="MTQ+HCQ",c(1,27,9:11,13)]
dat1[,c(3:6)]= lapply(dat1[,c(3:6)], scale)
dat2=reshape2::melt(dat1, id.vars=c(1,2), measure.vars=c(3:6))

p<-ggplot(dat2, aes(variable, value, color=`response classification`, fill=`response classification`))+
  geom_boxplot(outlier.colour = NA, width=.5, alpha=.01, position = position_dodge(width = .7),color="black")+
  geom_jitter(position = position_jitterdodge(jitter.width = .15, dodge.width = .7), alpha=1)+
  stat_compare_means(label="p.signif",method="wilcox.test", label.y = 5.2,size=5)+
  scale_color_manual(values = c("#21908c","#fde725"), name="")+ 
  scale_fill_manual(values = c("#21908c","#fde725"), name="")+
  coord_cartesian(ylim=c(-2.5,5.5))+
  theme_classic()+
  labs(x="",y="Scaled value")+
  theme(axis.text = element_text(color="black",size=15),
        axis.line = element_line(color="black"),
        axis.ticks =  element_line(color="black"),
        axis.title = element_text(color="black",size=15),
        legend.position="none")
topptx(p,'DAS28-MH.pptx')
