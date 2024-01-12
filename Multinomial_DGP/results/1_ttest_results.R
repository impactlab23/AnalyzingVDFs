library(tidyverse)
library(dplyr)

setwd("/GitHub/AnalyzingVDFs/Multinomial_DGP")
wdname <- paste(getwd(),"/1 - t test", sep = "")
setwd(wdname)

ttest_noeffect <- read.csv("ttest_noeffect.csv")
ttest_imp_both <- read.csv("ttest_imp_both.csv")
ttest_imp_lib_worse_mort <- read.csv("ttest_imp_lib_worse_mort.csv")
ttest_imp_mort_worse_lib <- read.csv("ttest_imp_mort_worse_lib.csv")
ttest_worsen_both <- read.csv("ttest_worsen_both.csv")
ttest_imp_lib_no_mort <- read.csv("ttest_imp_lib_no_mort.csv")
ttest_imp_mort_no_lib <- read.csv("ttest_imp_mort_no_lib.csv")

#-------calculating lambda that controls type one error----------

ttest_noeffect <- na.omit(ttest_noeffect) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2,
         trt_effect = V3)

#type 1 error

(ttest_type1_error <- nrow(ttest_noeffect[ttest_noeffect$positive_difference == "Yes",])/nrow(ttest_noeffect))


#--------based on the lambda, calculating "power" for each scenario----------
ttest_worsen_both <- na.omit(ttest_worsen_both) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2,
         trt_effect = V3)
(ttest_power_worsen_both <- nrow(ttest_worsen_both[ttest_worsen_both$positive_difference == "Yes",])/nrow(ttest_worsen_both))

ttest_imp_both <- na.omit(ttest_imp_both) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2,
         trt_effect = V3)
(ttest_power_imp_both <- nrow(ttest_imp_both[ttest_imp_both$positive_difference == "Yes",])/nrow(ttest_imp_both))

ttest_imp_lib_worse_mort <- na.omit(ttest_imp_lib_worse_mort) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2,
         trt_effect = V3)
(ttest_post_imp_lib_worse_mort <- nrow(ttest_imp_lib_worse_mort[ttest_imp_lib_worse_mort$positive_difference == "Yes",])/nrow(ttest_imp_lib_worse_mort))

ttest_imp_mort_worse_lib <- na.omit(ttest_imp_mort_worse_lib) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2,
         trt_effect = V3)
(ttest_post_imp_mort_worse_lib <- nrow(ttest_imp_mort_worse_lib[ttest_imp_mort_worse_lib$positive_difference == "Yes",])/nrow(ttest_imp_mort_worse_lib))

ttest_imp_lib_no_mort <- na.omit(ttest_imp_lib_no_mort) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2,
         trt_effect = V3)
(ttest_post_imp_lib_no_mort <- nrow(ttest_imp_lib_no_mort[ttest_imp_lib_no_mort$positive_difference == "Yes",])/nrow(ttest_imp_lib_no_mort))

ttest_imp_mort_no_lib <- na.omit(ttest_imp_mort_no_lib) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2,
         trt_effect = V3)
(ttest_post_imp_mort_no_lib <- nrow(ttest_imp_mort_no_lib[ttest_imp_mort_no_lib$positive_difference == "Yes",])/nrow(ttest_imp_mort_no_lib))

ttest_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "worsen_lib_imp_mort", "imp_lib_no_mort", "imp_mort_no_lib"),
                      power = c(ttest_type1_error, ttest_power_worsen_both, ttest_power_imp_both, ttest_post_imp_lib_worse_mort, ttest_post_imp_mort_worse_lib, ttest_post_imp_lib_no_mort, ttest_post_imp_mort_no_lib),
                      trt_effect_est_mean = c(mean(ttest_noeffect$trt_effect), mean(ttest_worsen_both$trt_effect), mean(ttest_imp_both$trt_effect), mean(ttest_imp_lib_worse_mort$trt_effect), mean(ttest_imp_mort_worse_lib$trt_effect), mean(ttest_imp_lib_no_mort$trt_effect), mean(ttest_imp_mort_no_lib$trt_effect)))

sqrt(ttest_1$power*(1 - ttest_1$power)/nrow(ttest_imp_mort_no_lib))[c(1,3,2,5,4,7,6)]
  