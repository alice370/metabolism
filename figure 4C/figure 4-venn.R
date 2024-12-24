library('ggvenn')
result<-read.csv('auto_log10_lm_DAS28.csv',row.names = 1)
positive <- result$ID[order(result$beta_DAS28CRP_lm,decreasing = T)[1:2]]
negative <- result$ID[result$DAS28CRP_lm == 'negative']
lm_sig<-c(positive,negative)
M_U_result<-read.csv('disease activity class test.csv',row.names = 1)
M_U_sig<-M_U_result$metabolites[which(M_U_result$P_DL<0.05|M_U_result$P_LM<0.05|M_U_result$P_MH<0.05|M_U_result$P_LH<0.05|M_U_result$P_DM<0.05|M_U_result$P_DH<0.05)]
con<-intersect(M_U_sig,lm_sig)
write.csv(con,'intersect_lm_test.csv')
library(eoffice)

venn_list <- list(
  "significant in lm" = lm_sig,
  "significant in comparison between disease activity groups" = M_U_sig
)
p1<-ggvenn(venn_list,show_percentage = F,show_elements = F,label_sep = ",",
           digits = 1,stroke_color = "white",
           fill_color = c("#DC00007F", "#D7A246"),
           set_name_color = "black")
p1
topptx(p1,filename = "intersect_lm_test.pptx")
