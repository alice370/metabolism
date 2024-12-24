library(ggplot2)
library(eoffice)
data=read.csv("pos-neg_all.csv",row.names = 1)
clin<-read.csv('sample_information.csv',check.names = F)
follow<-clin$`follow up-sample`[!clin$`follow up-sample`=='']
follow<-follow[!is.na(follow)]
A<-c(clin$Sample[clin$Group=='RA'],follow)
P<-clin$Sample[clin$Group=='At-risk of RA' ]         
H<-clin$Sample[clin$Group=='Health' ]
sample<-c(A,P,H)
dt<-data[,sample]
data<-as.data.frame(t(dt))
data$group<-c(rep("RA",406),rep("PRE",56),rep("HP",99))
count<- apply(data,2,function(c)sum(c!=0,na.rm = T))
a<-cbind(names(data),count)
write.csv(a,"sample_count.csv")

HP<-data[which(data$group=="HP"),-238]

count=c()
data1<-HP[1,]
sum<-colSums(data1,na.rm = T)
count<-sum(sum!=0,na.rm = T)
for (i in 2:nrow(HP)) {
  data1<-HP[c(1:i),]
  sum<-colSums(data1,na.rm = T)
  c<-sum(sum!=0,na.rm = T)
  count<-c(count,c)
}
count<-as.data.frame(count)
count$group<-rep('HP')
count$number<-as.numeric(row.names(count)) 
count$number<-count$number /99


RA<-data[which(data$group=="RA"),-238]
data1<-RA[1,]
sum<-colSums(data1,na.rm = T)
count_A<-sum(sum!=0,na.rm = T)
for (i in 2:nrow(RA)) {
  data1<-RA[c(1:i),]
  sum<-colSums(data1,na.rm = T)
  c<-sum(sum!=0,na.rm = T)
  count_A<-c(count_A,c)
}
count_A<-as.data.frame(count_A)
count_A$group<-rep('RA')
colnames(count_A)[1]<-'count'
count_A$number<-as.numeric(row.names(count_A)) 
count_A$number<-count_A$number / 406

PRE<-data[which(data$group=="PRE"),-238]
count_P=c()
data1<-PRE[1,]
sum<-colSums(data1,na.rm = T)
count_P<-sum(sum!=0,na.rm = T)
for (i in 2:nrow(PRE)) {
  data1<-PRE[c(1:i),]
  sum<-colSums(data1,na.rm = T)
  c<-sum(sum!=0,na.rm = T)
  count_P<-c(count_P,c)
}
count_P<-as.data.frame(count_P)
count_P$group<-rep('PRA')
count_P$number<-as.numeric(row.names(count_P)) 
count_P$number<-count_P$number / 56
colnames(count_P)[1]<-'count'
count_dt<-rbind(count,count_P,count_A)
write.csv(count_dt,'Ord_count.csv')

data<-read.csv('Ord_count.csv')

p<-ggplot(data,aes(x=data$number,y=count,col=group,shape=group))+
  geom_point(aes(x=data$number,y=count,col=group,shape=group),size=4,alpha=0.9)+
  geom_line(size=1.5,alpha=0.7)+
  scale_shape_manual(values=c(17,16,18))+
  scale_color_manual(values = c("#3578AC","#EAC450","#D5231F"))+
  theme_classic()+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())

topptx(p, 'sample_metaboltes_count.pptx')