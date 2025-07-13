library(ggplot2)
library(eoffice)
data<-read.csv('185meta_QC13.csv')
cv <- apply(data[,4:16], 1, function(x) {
   mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  if (mean_x != 0) {
    return(sd_x / mean_x)
  } else {
    return(NA) 
  }
})
data$CV<-cv


p<-ggplot(data , aes(x = X, y = CV, fill = CV)) +
  geom_bar(stat = "identity", color = "transparent") +
  scale_fill_gradient(high = "yellow", low = "orange") + 
  theme_minimal() +
  labs(
    x = "Metabolites",
    y = "CV VALUE"
  ) +theme_classic()+
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
topptx(p,'QC_13cv.pptx')
