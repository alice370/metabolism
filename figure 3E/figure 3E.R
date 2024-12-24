clin<-read.csv('RA_baseline_clin.csv',check.names = F)
# Fit separate linear models for each gender
model_male <- lm(`DAS28-CRP` ~ age, data = clin[clin$GENDER == "M", ])
model_female <- lm(`DAS28-CRP` ~ age, data = clin[clin$GENDER == "F", ])

# Extract R-squared values
r2_male <- summary(model_male)$r.squared
r2_female <- summary(model_female)$r.squared

# Extract p-values
p_value_male <- summary(model_male)$coefficients[2,4] # p-value for age in male model
p_value_female <- summary(model_female)$coefficients[2,4] # p-value for age in female model

# Assuming your dataframe is named merge1 with 'GENDER', 'age', and the variable of interest
p1<-ggplot(clin, aes(x = age, y = `DAS28-CRP`, color = GENDER)) +
  geom_point(size=2) +
  stat_smooth(method = "lm", se = TRUE, aes(label = paste("R2 =", round(summary(lm(`DAS28-CRP` ~ age + GENDER, data = clin))$r.squared, 2))), parse = TRUE) +
  theme_minimal() +
  theme(legend.position = 'none',
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("M" = "#154599", "F" = "darkred"))
p1
topptx(p1,filename = "age-DAS28CRP.pptx") 