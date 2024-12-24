result_gender<-read.csv('gender-wilcox.csv',row.names = 1)
result_RH<-read.csv('metabolite_comparison_results_scaled.csv')
set1<-result_RH$Metabolite[(result_RH$Coef_RA_HC<0)&result_RH$PValue_RA_HC<0.05]
set2<-result_gender$metabolic[(result_gender$p.value<0.05)&(result_gender$log10.fc<(-log10(1.2)))]
set3<-result_RH$Metabolite[(result_RH$Coef_RA_HC>0)&result_RH$PValue_RA_HC<0.05]
con<-intersect(set3,set2)
write.csv(con,"male-RA-intersect.csv")
####the plot draw in ppt####