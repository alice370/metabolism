library(dplyr)
library(glmnet)
library(multipleROC)
library(ggridges)
library(ggplot2)
library(ggpubr)
library(pROC)
library(scales)
library(PredictABEL)
library(tidyr)
setwd("C:/Users/alice/Desktop/实验数据/代谢组/返修/Revised 0525/药物反应code修改/NRI and IDI-MH")
rm(list = ls())
####MTX+LEF####
#################################################################### 5-fold feature selection
##################################### input data
d1 = read.csv("MTX+HCQ_dataset.csv",row.names = 1,check.names = F)
features_dt=read.csv("features by obs coef.csv",check.names = F)
####initial####
{ result_clin    <- NULL 
result_full    <- NULL 
coef_all_clin  <- NULL 
coef_all_full  <- NULL   
list_train_full=list()
list_test_full=list()
list_train_clinical=list()
list_test_clinical=list()  

}
features <- features_dt$feature[1:6]
dat1 <- d1[, c(1, 2, 5, which(colnames(d1) %in% features))]
dat1$sample<-row.names(dat1)
Y=dat1[dat1$`response`==1, ]
N=dat1[dat1$`response`==0, ]

####construct####
for (i in 1:100) {
  set.seed(i)
  # —— 划分 train / test —— 
  train <- rbind(
    Y[sample(1:nrow(Y), round(nrow(Y)*0.7), replace = FALSE), ],
    N[sample(1:nrow(N), round(nrow(N)*0.7), replace = FALSE), ]
  )
  test  <- dat1[!dat1$sample %in% train$sample, ]
  y      <- as.numeric(train$`response`)
  y_test <- as.numeric(test$`response`)
  
  # —— 特征提取 —— 

  x_clinical      <- train[, 1:2]
  x_test_clinical <- test[,  1:2]

  x_full      <- train[, c(1:2,4:9)]
  x_test_full <- test[,  c(1:2,4:9)]
  
  # —— 基线模型（临床）—— 
  cv_clinical       <- cv.glmnet(as.matrix(x_clinical), y,
                                 nfolds = 10, family = "binomial", alpha = 0)
  ridge_clinical    <- glmnet(as.matrix(x_clinical), y,
                              family = "binomial", lambda = cv_clinical$lambda.min, alpha = 0)
  coef_clin=data.frame(ID=colnames(x_clinical),i=rep(i,ncol(x_clinical)),coef=coef(ridge_clinical)@x[2:(ncol(x_clinical)+1)])
  train_y_clinical  <- predict(ridge_clinical, as.matrix(x_clinical), type = "response")
  df_train_clin     <- data.frame(y, train_y_clinical = as.numeric(train_y_clinical))
  train_auc_clin    <- multipleROC(y ~ train_y_clinical, data = df_train_clin, plot = FALSE)[["auc"]]
  train_cutoff_clin <- unique(multipleROC(y ~ train_y_clinical, data = df_train_clin, plot = FALSE)[["cutoff"]]$train_y_clinical)
  pred_train_clin   <- ifelse(train_y_clinical >= train_cutoff_clin, 1, 0)
  acc_train_clin    <- mean(pred_train_clin == y)
  
  test_y_clinical   <- predict(ridge_clinical, as.matrix(x_test_clinical), type = "response")
  df_test_clin      <- data.frame(y_test, test_y_clinical = as.numeric(test_y_clinical))
  test_auc_clin     <- multipleROC(y_test ~ test_y_clinical, data = df_test_clin, plot = FALSE)[["auc"]]
  test_cutoff_clin  <- unique(multipleROC(y_test ~ test_y_clinical, data = df_test_clin, plot = FALSE)[["cutoff"]]$test_y_clinical)
  pred_test_clin    <- ifelse(test_y_clinical >= test_cutoff_clin, 1, 0)
  acc_test_clin     <- mean(pred_test_clin == y_test)
  
  # —— 尝试跑训练集，如果触发 warning，直接跳过这一轮 —— 
  skip_iteration <- FALSE
  # —— 增强模型（临床+代谢物）—— 
  tryCatch({
  cv_full        <- cv.glmnet(as.matrix(x_full), y,
                              nfolds = 10, family = "binomial", alpha = 0)
  ridge_full     <- glmnet(as.matrix(x_full), y,
                           family = "binomial", lambda = cv_full$lambda.min, alpha = 0)
  coef_full=data.frame(ID=colnames(x_full),i=rep(i,ncol(x_full)),coef=coef(ridge_full)@x[2:(ncol(x_full)+1)])
  train_y_full   <- predict(ridge_full, as.matrix(x_full), type = "response")
  df_train_full  <- data.frame(y, train_y_full = as.numeric(train_y_full))
  train_auc_full <- multipleROC(y ~ train_y_full, data = df_train_full, plot = FALSE)[["auc"]]
  train_cutoff_full <- unique(multipleROC(y ~ train_y_full, data = df_train_full, plot = FALSE)[["cutoff"]]$train_y_full)
  pred_train_full   <- ifelse(train_y_full >= train_cutoff_full, 1, 0)
  acc_train_full    <- mean(pred_train_full == y)
  
  test_y_full    <- predict(ridge_full, as.matrix(x_test_full), type = "response")
  df_test_full   <- data.frame(y_test, test_y_full = as.numeric(test_y_full))
  test_auc_full  <- multipleROC(y_test ~ test_y_full, data = df_test_full, plot = FALSE)[["auc"]]
  test_cutoff_full <-unique(multipleROC(y_test ~ test_y_full, data = df_test_full, plot = FALSE)[["cutoff"]]$test_y_full)
  pred_test_full   <- ifelse(test_y_full >= test_cutoff_full, 1, 0)
  acc_test_full    <- mean(pred_test_full == y_test)
  
  }, warning = function(w) {
    # 如果捕获到 warning，标记跳过
    skip_iteration <<- TRUE
  }, error = function(e) {
    # 如果捕获到 error，也标记跳过
    skip_iteration <<- TRUE
  })
  
  if (skip_iteration) {
    next  # 跳过后续代码，进入下一轮
  }
  list_train_full[[i]] <- multipleROC(y ~ train_y_full, data = df_train_full, plot = FALSE)
  list_test_full[[i]]  <- multipleROC(y_test ~ test_y_full, data = df_test_full, plot = FALSE)
  
  list_train_clinical[[i]] <- multipleROC(y ~ train_y_clinical, data = df_train_clin, plot = FALSE)
  list_test_clinical[[i]]  <- multipleROC(y_test ~ test_y_clinical, data = df_test_clin, plot = FALSE)
  
  result_clin=rbind(result_clin, data.frame(i,train_auc_clin, acc_train_clin,test_auc_clin,acc_test_clin))
  result_full=rbind(result_full, data.frame(i,train_auc_full, acc_train_full,test_auc_full,acc_test_full))
  coef_all_clin=rbind(coef_all_clin,coef_clin)
  coef_all_full=rbind(coef_all_full,coef_full)
  
  # —— 训练集比较 —— 
  df_cmp_train <- data.frame(
    event    = y,
    baseline = as.numeric(train_y_clinical),
    newmodel = as.numeric(train_y_full)
  )
  # —— 训练集比较 —— 
  df_cmp_test <- data.frame(
    event    = y_test,
    baseline = as.numeric(test_y_clinical),
    newmodel = as.numeric(test_y_full)
  )
  # 分类 NRI/IDI：使用训练集的动态 cutoff
  ### 保存四个 reclassification 输出到文件
  # 1?? 训练集，临床模型
  sink(paste0("iter", i, "_train_clin.txt"))
  rc_train_clin <- reclassification(
    data      = df_cmp_train,
    cOutcome  = 1,
    predrisk1 = df_cmp_train$baseline,
    predrisk2 = df_cmp_train$newmodel,
    cutoff    = c(0, train_cutoff_clin, 1)
  )
  sink()
  
  # 2?? 训练集，全特征模型
  sink(paste0("iter", i, "_train_full.txt"))
  rc_train_full <- reclassification(
    data      = df_cmp_train,
    cOutcome  = 1,
    predrisk1 = df_cmp_train$baseline,
    predrisk2 = df_cmp_train$newmodel,
    cutoff    = c(0, train_cutoff_full, 1)
  )
  sink()
  
  # 3?? 测试集，临床模型
  sink(paste0("iter", i, "_test_clin.txt"))
  rc_test_clin <- reclassification(
    data      = df_cmp_test,
    cOutcome  = 1,
    predrisk1 = df_cmp_test$baseline,
    predrisk2 = df_cmp_test$newmodel,
    cutoff    = c(0, test_cutoff_clin, 1)
  )
  sink()
  
  # 4?? 测试集，全特征模型
  sink(paste0("iter", i, "_test_full.txt"))
  rc_test_full <- reclassification(
    data      = df_cmp_test,
    cOutcome  = 1,
    predrisk1 = df_cmp_test$baseline,
    predrisk2 = df_cmp_test$newmodel,
    cutoff    = c(0, test_cutoff_full, 1)
  )
  sink()
}

# 将所有模型提取出来并合并成数据框
extract_metrics <- function(model, iteration) {
  sens_text <- model['sens']
  
  # 提取数字
  values <- as.numeric(unlist(regmatches(sens_text, gregexpr("[0-9.]+", sens_text))))
  
  # 组合成一行数据
  c(i = iteration, Sens = values[1], Spec = values[2], PPV = values[3], NPV = values[4])
}
####merge the result of clin####
result_df1 <- do.call(rbind, mapply(extract_metrics, list_train_clinical, seq_along(list_train_clinical), SIMPLIFY = FALSE))
result_df2<-do.call(rbind, mapply(extract_metrics, list_test_clinical, seq_along(list_test_clinical), SIMPLIFY = FALSE))
result_df1 <- as.data.frame(result_df1)
result_df2 <- as.data.frame(result_df2)
colnames(result_df1) <- c("i", "Sens_train", "Spec_train", "PPV_train", "NPV_train")
colnames(result_df2) <- c("i", "Sens_test", "Spec_test", "PPV_test", "NPV_test")
result_merged_clin <- merge(result_df1, result_df2, by = "i", suffixes = c("_train", "_test"))
result_final_clin<-merge(result_clin, result_merged_clin, by = "i")
####merge the result of full####
result_df1 <- do.call(rbind, mapply(extract_metrics, list_train_full, seq_along(list_train_full), SIMPLIFY = FALSE))
result_df2<-do.call(rbind, mapply(extract_metrics, list_test_full, seq_along(list_test_full), SIMPLIFY = FALSE))
result_df1 <- as.data.frame(result_df1)
result_df2 <- as.data.frame(result_df2)
colnames(result_df1) <- c("i", "Sens_train", "Spec_train", "PPV_train", "NPV_train")
colnames(result_df2) <- c("i", "Sens_test", "Spec_test", "PPV_test", "NPV_test")
result_merged_full <- merge(result_df1, result_df2, by = "i", suffixes = c("_train", "_test"))
result_final_full<-merge(result_full, result_merged_full, by = "i")

####output NRI####
files <- list.files(pattern = "iter\\d+_.*\\.txt$")

# 初始化结果表
summary_list <- list()

# 正则表达式提取器
extract_metric <- function(lines, pattern) {
  line <- grep(pattern, lines, value = TRUE)
  if (length(line) == 0) return(c(NA, NA, NA))
  
  # 主值
  number <- sub(".*: ([0-9\\.\\-]+) \\[.*", "\\1", line)
  
  # 非贪婪提取 CI
  ci     <- sub(".*\\[ (.*?) \\].*", "\\1", line)
  
  # 支持 e-记法的 p 值
  pval   <- sub(".*p-value: ([0-9\\.eE\\-]+).*", "\\1", line)
  
  return(c(as.numeric(number), ci, as.numeric(pval)))
}


for (f in files) {
  lines <- readLines(f)
  
  # 提取迭代号和模型类型
  iter  <- sub("iter(\\d+)_.*\\.txt", "\\1", f)
  model <- sub("iter\\d+_(.*)\\.txt", "\\1", f)
  
  # 提取 NRI (Categorical)
  nri_cat <- extract_metric(lines, "NRI\\(Categorical\\)")
  
  if (any(is.na(nri_cat))) {
    message("?? Warning: NA detected in NRI (Categorical) for file: ", f)
  }
  
  # 提取 IDI
  idi <- extract_metric(lines, "IDI ")
  
  if (any(is.na(idi))) {
    message("?? Warning: NA detected in IDI for file: ", f)
  }
  
  # 保存一行结果
  summary_list[[f]] <- data.frame(
    iter        = as.numeric(iter),
    model       = model,
    nri_value   = nri_cat[1],
    nri_ci      = nri_cat[2],
    nri_p       = nri_cat[3],
    idi_value   = idi[1],
    idi_ci      = idi[2],
    idi_p       = idi[3],
    stringsAsFactors = FALSE
  )
}
# 合并所有行
summary_table <- do.call(rbind, summary_list)

# 分四个表格保存（按模型）
for (m in unique(summary_table$model)) {
  subset <- summary_table[summary_table$model == m, ]
  write.csv(subset, paste0("summary_", m, ".csv"), row.names = FALSE)
}

# 总表保存
write.csv(summary_table, "summary_all.csv", row.names = FALSE)

####plot and output of meta+clin####
write.csv(result_final_full,'MH-6feature+clin.csv')
write.csv(result_final_clin,'MH-only-clin.csv')
re=reshape2::melt(result_final_full[,c(1,2,4)],id.vars=1)
ggplot(re,aes(value, fill=variable))+
  geom_density(alpha=.5)+
  geom_vline(xintercept = c(median(result_final_full$train_auc), median(result_final_full$test_auc)), linetype=2, color=c("#AF322F","#5D669F"))+
  scale_fill_manual(values = c("#AF322F","#5D669F"),name="", label=c("Train","Test"))+
  theme_minimal()+
  labs(x="AUROC",y="Density")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        legend.position = "bottom")
ggsave("MH-6feature+clin-AUC.pdf", width = 3, height = 3.3)

range1=paste0("AUROC range:",round(range(result_final_full$train_auc)[1],2),'-',
              round(range(result_final_full$train_auc)[2],2))
median1=paste0("Median AUROC:",round(median(result_final_full$train_auc),2))
range2=paste0("AUROC range:",round(range(result_final_full$test_auc)[1],2),'-',
              round(range(result_final_full$test_auc)[2],2))
median2=paste0("Median AUROC:",round(median(result_final_full$test_auc),2))


list_train_full <- list_train_full[!sapply(list_train_full, is.null)]

plot_ROC(list_train_full,show.AUC = F,show.points = F,
         show.eta=FALSE,
         show.sens=FALSE 
)+scale_color_manual(values = rep(rgb(93,102,159,20, maxColorValue = 255),100))+
  theme(axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        panel.grid = element_blank())+
  annotate("text",x=0.65,y=0.2, label=range1, size=4)+
  annotate("text",x=0.65,y=0.1, label=median1, size=4)

ggsave("ML-6feature+clin-train.pdf.pdf", width = 3.2, height = 3)

list_test_full <- list_test_full[!sapply(list_test_full, is.null)]
plot_ROC(list_test_full,show.AUC = F,show.points = F,
         show.eta=FALSE,
         show.sens=FALSE 
)+scale_color_manual(values = rep(rgb(175,50,47,20, maxColorValue = 255),100))+
  theme(axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        panel.grid = element_blank())+
  annotate("text",x=0.65,y=0.2, label=range2, size=4)+
  annotate("text",x=0.65,y=0.1, label=median2, size=4)
ggsave("MH-6feature+clin-test.pdf", width = 3.2, height = 3)


ggplot(coef_all_full, aes(x = coef, y = ID)) +
  geom_density_ridges_gradient(aes(fill = stat(x))) +  
  geom_vline(xintercept = 0, linetype=2)+
  scale_fill_gradient2(high="#AF322F", low="#5D669F")+
  theme_classic()+
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        legend.position = "bottom")
ggsave("MH-6feature+clin-coef.pdf", width = 6, height = 5)  

####plot and output clin####
re=reshape2::melt(result_final_clin[,c(1,2,4)],id.vars=1)
ggplot(re,aes(value, fill=variable))+
  geom_density(alpha=.5)+
  geom_vline(xintercept = c(median(result_final_clin$train_auc), median(result_final_clin$test_auc)), linetype=2, color=c("#AF322F","#5D669F"))+
  scale_fill_manual(values = c("#AF322F","#5D669F"),name="", label=c("Train","Test"))+
  theme_minimal()+
  labs(x="AUROC",y="Density")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        legend.position = "bottom")
ggsave("MH-2clin-AUC.pdf", width = 3, height = 3.3)

range1=paste0("AUROC range:",round(range(result_final_clin$train_auc)[1],2),'-',
              round(range(result_final_clin$train_auc)[2],2))
median1=paste0("Median AUROC:",round(median(result_final_clin$train_auc),2))
range2=paste0("AUROC range:",round(range(result_final_clin$test_auc)[1],2),'-',
              round(range(result_final_clin$test_auc)[2],2))
median2=paste0("Median AUROC:",round(median(result_final_clin$test_auc),2))

list_train_clinical <- list_train_clinical[!sapply(list_train_clinical, is.null)]
plot_ROC(list_train_clinical,show.AUC = F,show.points = F,
         show.eta=FALSE,
         show.sens=FALSE 
)+scale_color_manual(values = rep(rgb(93,102,159,20, maxColorValue = 255),100))+
  theme(axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        panel.grid = element_blank())+
  annotate("text",x=0.65,y=0.2, label=range1, size=4)+
  annotate("text",x=0.65,y=0.1, label=median1, size=4)

ggsave("MH-2clin-train.pdf.pdf", width = 3.2, height = 3)

list_test_clinical <- list_test_clinical[!sapply(list_test_clinical, is.null)]
plot_ROC(list_test_clinical,show.AUC = F,show.points = F,
         show.eta=FALSE,
         show.sens=FALSE 
)+scale_color_manual(values = rep(rgb(175,50,47,20, maxColorValue = 255),100))+
  theme(axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        panel.grid = element_blank())+
  annotate("text",x=0.65,y=0.2, label=range2, size=4)+
  annotate("text",x=0.65,y=0.1, label=median2, size=4)
ggsave("MH-2clin-test.pdf", width = 3.2, height = 3)

coef_all_clin<-unique(coef_all_clin)
ggplot(coef_all_clin, aes(x = coef, y = ID)) +
  geom_density_ridges_gradient(aes(fill = stat(x))) +  
  geom_vline(xintercept = 0, linetype=2)+
  scale_fill_gradient2(high="#AF322F", low="#5D669F")+
  theme_classic()+
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text = element_text(color="black", size=13),
        axis.title = element_text(color="black", size=13),
        legend.position = "bottom")
ggsave("MH-2clin-coef.pdf", width = 6, height = 5)  

####plot IDI of test#####
df <- read.csv("summary_test_full.csv")
df <- df %>%
  separate(idi_ci, into = c("ci_lower", "ci_upper"), sep = " - ") %>%
  mutate(
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    significance = ifelse(idi_p < 0.05, "Significant", "Not Significant")
  )
# 计算各自 median
median_test <- median(df$idi_value, na.rm = TRUE)
ggplot(df, aes(x = as.numeric(iter), y = idi_value, color = significance)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, size = 0.8) +
  scale_color_manual(values = c("Significant" = "#AF322F", "Not Significant" = "gray")) +
  theme_classic(base_size = 14) +
  labs(
    title = " ",
    x = "Iteration",
    y = "IDI Value",
    color = "Significance"
  ) +
  theme(
    legend.position = "top",
    text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 18, face = "bold", color = "black")
  )+
  annotate(
    "text",
    x = Inf,  # 右边界
    y = -Inf,  # 下方（你可调整高度）
    label = paste0("Median IDI of test set: ", round(median_test, 3)),
    color = "black",
    size = 5,
    hjust = 1 , # 右对齐
    vjust = -1 
  )
ggsave("MH-IDI-test.pdf", width = 6, height = 4)
####plot nri of train#####
df1 <- read.csv("summary_train_full.csv")
df1 <- df1 %>%
  separate(idi_ci, into = c("ci_lower", "ci_upper"), sep = " - ") %>%
  mutate(
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    significance = ifelse(idi_p < 0.05, "Significant", "Not Significant")
  )
# 计算各自 median
median_train <- median(df1$idi_value, na.rm = TRUE)
ggplot(df1, aes(x = as.numeric(iter), y = idi_value, color = significance)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, size = 0.8) +
  scale_color_manual(values = c("Significant" = "#5D669F", "Not Significant" = "gray")) +
  theme_classic(base_size = 14) +
  labs(
    title = " ",
    x = "Iteration",
    y = "IDI Value",
    color = "Significance"
  ) +
  theme(
    legend.position = "top",
    text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 18, face = "bold", color = "black")
  )+
  annotate(
    "text",
    x = Inf,  # 右边界
    y = -Inf,  # 下方（你可调整高度）
    label = paste0("Median IDI of train set: ", round(median_train, 3)),
    color = "black",
    size = 5,
    hjust = 1 , # 右对齐
    vjust = -1 
  )
ggsave("MH-IDI-train.pdf", width = 6, height = 4)


