library(eoffice)
data<-read.csv("364sample_185metabolites_mean_log10.csv",row.names = 1)
clin<-read.csv("clin_new.csv")
range(data)
selected_columns <- data[, clin$Sample]
set.seed(1234)
boxplot(selected_columns[,sample(1:209)],axes=F,col="white",
        border = "#B22F36",lwd=1.5, cex=.3, 
        xlim=c(1,364),ylim=c(-8,5),xlab = "Samples",ylab = "log10 (metabolites)",
        cex.lab=1, font.lab=1) 
boxplot(selected_columns[,sample(210:265)],add = T,axes=F,at=c(210:265),col="white",
        border="#EFC000FF",lwd=1.5, cex=.3,ylim=c(-8,5))
boxplot(selected_columns[,sample(266:364)],add = T,axes=F,at=c(266:364),col="white",
        border="#006BBE",lwd=1.5, cex=.3,ylim=c(-8,5))
axis(2,at=c(-8,0,5), 
     label=c("-8","0","5"),lwd=0.8,
     lwd.ticks = 0.7,
     font.axis=1,
     cex.axis=0.8)
####save it as png using R plots device it self####
