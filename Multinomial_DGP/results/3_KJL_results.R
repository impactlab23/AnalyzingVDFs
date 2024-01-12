library(tidyverse)
library(dplyr)


setwd("/GitHub/AnalyzingVDFs/Multinomial_DGP")
wdname <- paste(getwd(),"/3 - KJL", sep = "")
setwd(wdname)

LR_noeffect <- read.csv("LR_noeffect.csv")
LR_imp_both <- read.csv("LR_imp_both.csv")
LR_imp_lib_worse_mort <- read.csv("LR_imp_lib_worse_mort.csv")
LR_imp_mort_worse_lib <- read.csv("LR_imp_mort_worse_lib.csv")
LR_worsen_both <- read.csv("LR_worsen_both.csv")
LR_imp_lib_no_mort <- read.csv("LR_imp_lib_no_mort.csv")
LR_imp_mort_no_lib <- read.csv("LR_imp_mort_no_lib.csv")

#-------calculating lambda that controls type one error----------

LR_noeffect <- LR_noeffect %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(LR_noeffect) %>%
  filter(V1 != "Error") %>%
  mutate(trt_effect_YN = V1,
         p_val = as.numeric(V2),
         trt_est = as.numeric(V3),
         trt_sd = as.numeric(V4)) 

(LR_type1_error <- mean(LR_noeffect$p_val/2 <= 0.05 & LR_noeffect$trt_est >= 0))

#--------based on the lambda, calculating "power" for each scenario----------

################## worsening both ################################
LR_worsen_both <- LR_worsen_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(LR_worsen_both) %>%
  filter(V1 != "Error") %>%
  mutate(trt_effect_YN = V1,
         p_val = as.numeric(V2),
         trt_est = as.numeric(V3),
         trt_sd = as.numeric(V4)) 

(LR_worsen_both_pow <- mean(LR_worsen_both$p_val/2 <= 0.05 & LR_worsen_both$trt_est >= 0))

################## improving both ################################
LR_imp_both <- LR_imp_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(LR_imp_both) %>%
  filter(V1 != "Error") %>%
  mutate(trt_effect_YN = V1,
         p_val = as.numeric(V2),
         trt_est = as.numeric(V3),
         trt_sd = as.numeric(V4)) 

(LR_imp_both_pow <- mean(LR_imp_both$p_val/2 <= 0.05 & LR_imp_both$trt_est >= 0))


################## improving liberation worsening mortality ################################
LR_imp_lib_worse_mort <- LR_imp_lib_worse_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(LR_imp_lib_worse_mort) %>%
  filter(V1 != "Error") %>%
  mutate(trt_effect_YN = V1,
         p_val = as.numeric(V2),
         trt_est = as.numeric(V3),
         trt_sd = as.numeric(V4)) 

(LR_imp_lib_worse_mort_pow <- mean(LR_imp_lib_worse_mort$p_val/2 <= 0.05 & LR_imp_lib_worse_mort$trt_est >= 0))


################## improving mortality worsening liberation ################################
LR_imp_mort_worse_lib <- LR_imp_mort_worse_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(LR_imp_mort_worse_lib) %>%
  filter(V1 != "Error") %>%
  mutate(trt_effect_YN = V1,
         p_val = as.numeric(V2),
         trt_est = as.numeric(V3),
         trt_sd = as.numeric(V4)) 

(LR_imp_mort_worse_lib_pow <- mean(LR_imp_mort_worse_lib$p_val/2 <= 0.05 & LR_imp_mort_worse_lib$trt_est >= 0))

################## improving liberation no mortality ################################
LR_imp_lib_no_mort <- LR_imp_lib_no_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(LR_imp_lib_no_mort) %>%
  filter(V1 != "Error") %>%
  mutate(trt_effect_YN = V1,
         p_val = as.numeric(V2),
         trt_est = as.numeric(V3),
         trt_sd = as.numeric(V4)) 

(LR_imp_lib_no_mort_pow <- mean(LR_imp_lib_no_mort$p_val/2 <= 0.05 & LR_imp_lib_no_mort$trt_est >= 0))

################## improving mortality no liberation ################################
LR_imp_mort_no_lib <- LR_imp_mort_no_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(LR_imp_mort_no_lib) %>%
  filter(V1 != "Error") %>%
  mutate(trt_effect_YN = V1,
         p_val = as.numeric(V2),
         trt_est = as.numeric(V3),
         trt_sd = as.numeric(V4)) 

(LR_imp_mort_no_lib_pow <- mean(LR_imp_mort_no_lib$p_val/2 <= 0.05 & LR_imp_mort_no_lib$trt_est >= 0))

LR_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "worsen_lib_imp_mort", "imp_lib_no_mort", "imp_mort_no_lib"),
                   post_p = c(LR_type1_error, LR_worsen_both_pow, LR_imp_both_pow, LR_imp_lib_worse_mort_pow, LR_imp_mort_worse_lib_pow, LR_imp_lib_no_mort_pow, LR_imp_mort_no_lib_pow),
                   trt_est = c(mean(LR_noeffect$trt_est), mean(LR_worsen_both$trt_est), mean(LR_imp_both$trt_est), mean(LR_imp_lib_worse_mort$trt_est), mean(LR_imp_mort_worse_lib$trt_est), mean(LR_imp_lib_no_mort$trt_est), mean(LR_imp_mort_no_lib$trt_est)),
                   sd_est = c(mean(LR_noeffect$trt_sd), mean(LR_worsen_both$trt_sd), mean(LR_imp_both$trt_sd), mean(LR_imp_lib_worse_mort$trt_sd), mean(LR_imp_mort_worse_lib$trt_sd), mean(LR_imp_lib_no_mort$trt_sd), mean(LR_imp_mort_no_lib$trt_sd)))

sqrt(LR_1$post_p*(1 - LR_1$post_p)/nrow(LR_imp_lib_no_mort))[c(1,3,2,5,4,7,6)]
