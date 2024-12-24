
library(dplyr)
library(glmnet)
library(multipleROC)
library(ggridges)

rm(list = ls())
####MTX+LEF####
clin=read.csv("sample_information.csv",check.names = F)
data<-read.csv("209sample_188metabolites_auto_log10.csv",check.names = F,row.names = 1)
test<-read.csv("response test.csv",check.names = F)
meta<-test$metabolites[test$P_ML<0.05]
sample<-clin$sample[clin$DRUG=="MTQ+LEF"&!is.na(clin$`response classification`)]
dat1<-data[meta,sample]
dat1<-as.data.frame(t(dat1))
dat1$sample<-row.names(dat1)
dat1<-merge(clin[,c(1,9:11,13,27)],dat1,by='sample')
dat1 <- dat1 %>%
  mutate(`response classification` = case_when(
    `response classification` %in% c("Good Response", "Moderate Response") ~ "1",
    `response classification` %in% c("No Response") ~ "0",
    TRUE ~ `response classification`
  ))
Y=dat1[dat1$`response classification`==1, ]
N=dat1[dat1$`response classification`==0, ]
result=NULL 
list_train=list()
list_test=list()
coef_all=NULL
####ridge####
for(i in 1:100){
  set.seed(i)
  train=rbind(Y[sample(1:nrow(Y),round(nrow(Y)*0.5),replace = F),],N[sample(1:nrow(N),round(nrow(N)*0.5),replace = F),])
  test=dat1[!c(dat1$sample%in%train$sample),]
  y=as.numeric(train$`response classification`)
  y_test <- as.numeric(test$`response classification`)
  
  x <- cbind(train[,c(2:4,7:13)])
  x_test <- cbind(test[,c(2:4,7:13)])
  
  cvfit=cv.glmnet(as.matrix(x), y, nfolds = 10,family="binomial",alpha=0)
  ridge <- glmnet(as.matrix(x),y, family="binomial", lambda=cvfit$lambda.min, alpha=0)
  coef=data.frame(ID=colnames(x),i=rep(i,ncol(x)),coef=coef(ridge)@x[2:(ncol(x)+1)])
  
  train_y <- predict(ridge,as.matrix(x), type="response")
  df=data.frame(y,train_y=as.numeric(train_y))
  train_auc=multipleROC(y~train_y,data=df,plot =F)[["auc"]]
  
  
  test_y <- predict(ridge,as.matrix(x_test), type="response") 
  df1=data.frame(y_test,test_y=as.numeric(test_y))
  test_auc=multipleROC(y_test~test_y,data=df1,plot =F)[["auc"]]
  
  result=rbind(result, data.frame(i,train_auc, test_auc))
  
  coef_all=rbind(coef_all,coef)
  
  list_train[[i]]=multipleROC(y~train_y,data=df,plot =F)
  list_test[[i]]=multipleROC(y_test~test_y,data=df1,plot =F)
}

re=reshape2::melt(result,id.vars=1)
ggplot(re,aes(value, fill=variable))+
  geom_density(alpha=.5)+
  geom_vline(xintercept = c(median(result$train_auc), median(result$test_auc)), linetype=2, color=c("#AF322F","#5D669F"))+
  scale_fill_manual(values = c("#AF322F","#5D669F"),name="", label=c("Train","Test"))+
  theme_minimal()+
  labs(x="AUROC",y="Density")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        legend.position = "bottom")
ggsave("ML-7meta+3clin-AUC.pdf", width = 3, height = 3.3)
range(result$train_auc);range(result$test_auc)
median(result$train_auc);median(result$test_auc)
range1=paste0("AUROC range:",round(range(result$train_auc)[1],2),'-',
              round(range(result$train_auc)[2],2))
median1=paste0("Median AUROC:",round(median(result$train_auc),2))
range2=paste0("AUROC range:",round(range(result$test_auc)[1],2),'-',
              round(range(result$test_auc)[2],2))
median2=paste0("Median AUROC:",round(median(result$test_auc),2))

plot_ROC(list_train,show.AUC = F,show.points = F,
         show.eta=FALSE,
         show.sens=FALSE 
)+scale_color_manual(values = rep(rgb(93,102,159,20, maxColorValue = 255),100))+
  theme(axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        panel.grid = element_blank())+
  annotate("text",x=0.65,y=0.2, label=range1, size=4)+
  annotate("text",x=0.65,y=0.1, label=median1, size=4)

ggsave("ML-7meta+3clin-train.pdf.pdf", width = 3.2, height = 3)


plot_ROC(list_test,show.AUC = F,show.points = F,
         show.eta=FALSE,
         show.sens=FALSE 
)+scale_color_manual(values = rep(rgb(175,50,47,20, maxColorValue = 255),100))+
  theme(axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        panel.grid = element_blank())+
  annotate("text",x=0.65,y=0.2, label=range2, size=4)+
  annotate("text",x=0.65,y=0.1, label=median2, size=4)
ggsave("ML-7meta+3clin-test.pdf", width = 3.2, height = 3)


ggplot(coef_all, aes(x = coef, y = ID)) +
  geom_density_ridges_gradient(aes(fill = stat(x))) +  
  geom_vline(xintercept = 0, linetype=2)+
  scale_fill_gradient2(high="#AF322F", low="#5D669F")+
  theme_classic()+
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        legend.position = "bottom")
ggsave("ML-7meta+3clin-coef.pdf", width = 6, height = 5)
