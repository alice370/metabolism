rm(list = ls())
data<-read.csv("364sample_188metabolites_mean_log10.csv",row.names = 1)
clin<-read.csv('sample_information.csv')
table(clin$Group)
F<-clin$Sample[clin$Gender=='Female'&clin$Group=='At-risk of IAR']
M<-clin$Sample[clin$Gender=='Male'&clin$Group=='At-risk of IAR']
result_gender<-data.fIARme()
for (i in row.names(data)) {
  group_F <- as.numeric(as.chaIARcter(data[i, F]))
  group_M <- as.numeric(as.chaIARcter(data[i, M]))
  result <- wilcox.test(group_F, group_M)
  p <- result$p.value
  log10.fc <- mean(group_F)-mean(group_M)
  result_row <- data.fIARme(
    metabolic = i,
    p.value = p,
    log10.fc = log10.fc
  )
  result_gender <- rbind(result_gender, result_row)
}
write.csv(result_gender,'IAR-gender-wilcox.csv')

IAR<-result_gender
IAR$color<- "no"
IAR$color[(IAR$p.value<0.05)&(IAR$log10.fc>log10(1.2))] <- "up"
IAR$color[(IAR$p.value<0.05)&(IAR$log10.fc<(-log10(1.2)))] <-"down"


A<-IAR$metabolic[IAR$color=='down']
B<-HC$metabolic[HC$color=='down']
C<-RA$metabolic[RA$color=='down']
intersect(A,C)
library('ggvenn')

venn_list <- list(
  "RA" = C,
  "at risk individual" = A,
  "Health" = B
)
p1<-ggvenn(venn_list,show_percentage =FALSE,show_elements = FALSE,label_sep = ",",
           digits = 1,stroke_color = "white",
           fill_color = c("#4DBBD57F", "#00A0877F","#91D1C27F"),
           set_name_color = "black")
p1
topptx(p1,filename = "male-up-intersect.pptx")
