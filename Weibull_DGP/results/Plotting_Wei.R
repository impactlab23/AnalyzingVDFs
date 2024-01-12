setwd("/AnalyzingVDFs/Weibull_DGP")
source(paste(getwd(),"/results/1_ttest_results.R", sep = ""))
source(paste(getwd(),"/results/2_wilcox_results.R", sep = ""))
source(paste(getwd(),"/results/3_KJL_results.R", sep = ""))
source(paste(getwd(),"/results/4_hurdle_results.R", sep = ""))
source(paste(getwd(),"/results/5_ordinal1_results.R", sep = ""))
source(paste(getwd(),"/results/6_ordinal2_results.R", sep = ""))
source(paste(getwd(),"/results/7_CR_results.R", sep = ""))
source(paste(getwd(),"/results/8_binSurv_results.R", sep = ""))

ttest_1[, 1:2]
wilcox_1[, 1:2]
LR_1[, 1:2]
binSurv_1[, 1:2]
CR_1[, 1:2]
Ordinal1_1[, 1:2]
Ordinal3_1[, 1:2]
dat <- data.frame(scenario = ttest_1[, 1], 
                  ttest = ttest_1[, 2],
                  wilcox = wilcox_1[, 2],
                  KJL = LR_1[, 2],
                  binSurv = binSurv_1[, 2],
                  CR = CR_1[, 2],
                  Cat3 = Ordinal1_1[, 2],
                  Cat10 = Ordinal3_1[, 2]) %>% 
  pivot_longer( 2:8) %>%
  mutate(scenario = factor(scenario, levels = c("imp_lib_no_mort",
                                                "imp_mort_no_lib",
                                                "worsen_mort_imp_lib",
                                                "worsen_lib_imp_mort",
                                                "worsen_both",
                                                "imp_both",
                                                "no_eff"
  )))



library(ggplot2)

ggplot(dat, aes(x = value, y = scenario, color = name)) +
  geom_point(size = 2.5) + 
  theme_bw()  + 
  scale_y_discrete(labels = c('Shorter Ventilation Only',
                              'Benefit for Death Alone',
                              "Shorter Ventilation \n and Harm for Death",
                              "Benefit for Death and \n Longer Ventilation",
                              "Harm for Both",
                              'Benefit for Both',
                              "No Treatment Effect"
  )) + 
  labs(x = "Probability of Positive Treatment Effect",
       y = "Treatment Effect Scenario",
       color = "Analysis Model") + 
  scale_color_discrete(label = c(
    "Binary Death with \n Time to Liberation",
    "10 Category Ordinal",
                                  "3 Category Ordinal",
    "Competing Risk Model",
                                  "KJL Test",
                                  "T-Test",
                                  "Wilcoxon Test"))
