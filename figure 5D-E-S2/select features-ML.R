library(dplyr)
library(glmnet)
library(multipleROC)
library(ggridges)
library(ggplot2)
library(ggpubr)
library(pROC)
library(scales)
library(eoffice)
library(PredictABEL)
setwd("C:/Users/alice/Desktop/实验数据/代谢组/数据处理/Code/figure 5D-E-S2C-S3")
rm(list = ls())
####MTX+LEF####
#################################################################### 5-fold feature selection
##################################### input data
d1 = read.csv("MTX+LEF_dataset.csv",row.names = 1,check.names = F)
table(d1$response)
df = d1[,-c(1:4,6)]; colnames(df)[1]="target"
names(d1)
yes = df %>% filter(target==1)
no = df %>% filter(target==0)
################################################################
Features = NULL
sf <- list()
data <- df
x = as.matrix(data[, 2:ncol(data)])
y = as.numeric(data$target)
# feature selection
set.seed(123); fitcv<-cv.glmnet(x,y,family = "binomial",alpha = 1,type.measure = "auc", nfolds=4)
set.seed(123); lasso_model <- glmnet(x, y, alpha=1, family="binomial", lambda=fitcv$lambda.min)
sf_lasso <- colnames(x)[which(coef(lasso_model) != 0)-1]  
coef_lasso <- coef(lasso_model)

coef_lasso <- data.frame(
  feature = colnames(x),
  coef    = as.numeric(coef_lasso[-1]) 
)
coef_lasso <- coef_lasso[coef_lasso$feature %in% sf_lasso, ]

coef_lasso <- coef_lasso[order(-abs(coef_lasso$coef)), ]
write.csv(coef_lasso,'features by obs coef.csv')
# ½«ËùÓÐÄ£ÐÍÌáÈ¡³öÀ´²¢ºÏ²¢³ÉÊý¾Ý¿ò
extract_metrics <- function(model, iteration) {
  sens_text <- model['sens']
  

  values <- as.numeric(unlist(regmatches(sens_text, gregexpr("[0-9.]+", sens_text))))

  c(i = iteration, Sens = values[1], Spec = values[2], PPV = values[3], NPV = values[4])
}
####select the optimal numbers of features####
setwd("C:/Users/alice/Desktop/实验数据/代谢组/数据处理/Code/figure 5D-E-S2C-S3/select features-ML")

for (n in 1:9) {
  top_features <- head(coef_lasso$feature, n)
  
  dat1 <- d1[, c(1, 2, 3,5, which(colnames(d1) %in% top_features))]
  dat1$sample <- row.names(dat1)
  

  Y <- dat1[dat1$`response` == 1, ]
  N <- dat1[dat1$`response` == 0, ]
  
  ####initial output####
  result_full <- NULL
  list_train_full <- list()
  list_test_full <- list()

  ####construct####
  for (i in 1:100) {
    set.seed(i)
    
    # ³õÊ¼»¯Ä¬ÈÏ NA
    this_result <- data.frame(
      i = i,
      train_auc_full = NA,
      acc_train_full = NA,
      test_auc_full = NA,
      acc_test_full = NA,
      train_cutoff_full = NA,
      test_cutoff_full = NA
    )
    
    list_train_full[[i]] <- NA
    list_test_full[[i]] <- NA
    
    tryCatch({
      #### divide according to 7:3 ####
      train <- rbind(
        Y[sample(1:nrow(Y), round(nrow(Y) * 0.7), replace = FALSE), ],
        N[sample(1:nrow(N), round(nrow(N) * 0.7), replace = FALSE), ]
      )
      test <- dat1[!dat1$sample %in% train$sample, ]
      y <- as.numeric(train$`response`)
      y_test <- as.numeric(test$`response`)
      
      x_full <- train[, c(1, 2, 3,5:(ncol(train) - 1))]
      x_test_full <- test[, c(1, 2, 3,5:(ncol(train) - 1))]
      
      # ½¨Ä££¨Áë»Ø¹é£©
      cv_full <- cv.glmnet(as.matrix(x_full), y, nfolds = 10, family = "binomial", alpha = 0)
      ridge_full <- glmnet(as.matrix(x_full), y, family = "binomial", lambda = cv_full$lambda.min, alpha = 0)
      
      # ¡ª¡ª ÑµÁ·¼¯Ô¤²â ¡ª¡ª 
      train_y_full <- predict(ridge_full, as.matrix(x_full), type = "response")
      df_train_full <- data.frame(y, train_y_full = as.numeric(train_y_full))
      roc_train <- multipleROC(y ~ train_y_full, data = df_train_full, plot = FALSE)
      train_auc_full <- roc_train[["auc"]]
      train_cutoff_full <- roc_train[["cutoff"]]$train_y_full
      pred_train_full <- ifelse(train_y_full >= train_cutoff_full, 1, 0)
      acc_train_full <- mean(pred_train_full == y)
      
      # ¡ª¡ª ²âÊÔ¼¯Ô¤²â ¡ª¡ª 
      test_y_full <- predict(ridge_full, as.matrix(x_test_full), type = "response")
      df_test_full <- data.frame(y_test, test_y_full = as.numeric(test_y_full))
      roc_test <- multipleROC(y_test ~ test_y_full, data = df_test_full, plot = FALSE)
      test_auc_full <- roc_test[["auc"]]
      test_cutoff_full <- roc_test[["cutoff"]]$test_y_full
      pred_test_full <- ifelse(test_y_full >= test_cutoff_full, 1, 0)
      acc_test_full <- mean(pred_test_full == y_test)
      
      # Èç¹û¶¼Ã»±¨´í£¬ÌîÈëÊµ¼Ê½á¹û
      this_result <- data.frame(
        i = i,
        train_auc_full = train_auc_full,
        acc_train_full = acc_train_full,
        test_auc_full = test_auc_full,
        acc_test_full = acc_test_full,
        train_cutoff_full = train_cutoff_full,
        test_cutoff_full = test_cutoff_full
      )
      list_train_full[[i]] <- roc_train
      list_test_full[[i]] <- roc_test
    }, warning = function(w) {
      # ¾²Ä¬Ìø¹ý£¬±£Áô NA
    }, error = function(e) {
      # ¾²Ä¬Ìø¹ý£¬±£Áô NA
    })
    
    # Ã¿ÂÖ½á¹û»ã×Ü
    result_full <- rbind(result_full, this_result)
  }
  
  # ºÏ²¢Ãô¸Ð¶È/ÌØÒì¶ÈµÈÖ¸±ê
  result_df1 <- do.call(rbind, mapply(extract_metrics, list_train_full, seq_along(list_train_full), SIMPLIFY = FALSE))
  result_df2 <- do.call(rbind, mapply(extract_metrics, list_test_full, seq_along(list_test_full), SIMPLIFY = FALSE))
  result_df1 <- as.data.frame(result_df1)
  result_df2 <- as.data.frame(result_df2)
  colnames(result_df1) <- c("i", "Sens_train", "Spec_train", "PPV_train", "NPV_train")
  colnames(result_df2) <- c("i", "Sens_test", "Spec_test", "PPV_test", "NPV_test")
  
  result_merged_full <- merge(result_df1, result_df2, by = "i", suffixes = c("_train", "_test"))
  result_final_full <- merge(result_full, result_merged_full, by = "i")
  
  # ±£´æÃ¿Ò»ÂÖµÄ½á¹û
  write.csv(result_final_full, paste0('top', n, '_feature_parameters.csv'), row.names = FALSE)
}


summary_list <- list()
n=1
for (n in 1:30) {
  filename <- paste0("top", n, "_feature_parameters.csv")
  
  if (file.exists(filename)) {
    df <- read.csv(filename)
    
    na_rows <- sum(apply(df[, -1], 1, function(row) all(is.na(row))))
    
    if (na_rows < 50) {
      df_numeric <- df[, -1]
      
      metric_stats <- lapply(df_numeric, function(col) {
        c(median = median(col, na.rm = TRUE),
          min = min(col, na.rm = TRUE),
          max = max(col, na.rm = TRUE))
      })
      
      metric_flat <- unlist(metric_stats)
      
      summary_row <- data.frame(
        FeatureCount = n,
        t(metric_flat),
        row.names = NULL,
        check.names = FALSE
      )
      
      summary_list[[length(summary_list) + 1]] <- summary_row
    }
  }
}

summary_table <- do.call(rbind, summary_list)

# ±£´æ½á¹û
write.csv(summary_table, "summary_over50_features_median_range.csv", row.names = FALSE)


# ¶ÁÈëÄãµÄ±í¸ñ
df <- read.csv("summary_over50_features_median_range.csv")

# ÕûÀíÑµÁ·¼¯Êý¾Ý
train_df <- data.frame(
  FeatureCount = df$FeatureCount,
  Median = df$train_auc_full.median,
  Lower = df$train_auc_full.min,
  Upper = df$train_auc_full.max,
  Set = "Training"
)

# ÕûÀí²âÊÔ¼¯Êý¾Ý
test_df <- data.frame(
  FeatureCount = df$FeatureCount,
  Median = df$test_auc_full.median,
  Lower = df$test_auc_full.min,
  Upper = df$test_auc_full.max,
  Set = "Testing"
)

# ºÏ²¢
plot_df <- rbind(train_df, test_df)

# ×Ô¶¨ÒåÑÕÉ«
color_map <- c("Training" = "royalblue", "Testing" = "#bd5553")

ggplot(plot_df, aes(x = FeatureCount, y = Median, color = Set)) +
  geom_line(size = 1.5) +  # ¼Ó´ÖÏß
  geom_point(size = 3) +   # ¼Ó´Öµã
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.4, size = 1) +  # ¼Ó´ÖÎó²î°ô
  scale_color_manual(values = color_map) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(limits = c(0.25, 1)) +
  theme_classic(base_size = 16) +
  labs(
    title = "",
    x = "Number of Features",
    y = "AUC",
    color = "Dataset"
  ) +
  theme(
    legend.position = "right",
    text = element_text(color = "black"),            
    axis.text.x = element_text(angle = 0, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black"),
  )
topptx(p,'Model Performance vs. Number of Features.pptx')

ggsave("Model Performance vs. Number of Features.pdf", width = 8, height = 5)  

