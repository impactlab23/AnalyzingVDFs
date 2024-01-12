library(tidyverse)
library(dplyr)


setwd("/GitHub/AnalyzingVDFs/Weibull_DGP")
wdname <- paste(getwd(),"/8 - binSurv", sep = "")
setwd(wdname)

binSurv_noeffect <- read.csv("binSurv_noeffect.csv")
binSurv_imp_both <- read.csv("binSurv_imp_both.csv")
binSurv_imp_lib_worse_mort <- read.csv("binSurv_imp_lib_worse_mort.csv")
binSurv_imp_mort_worse_lib <- read.csv("binSurv_imp_mort_worse_lib.csv")
binSurv_worsen_both <- read.csv("binSurv_worsen_both.csv")
binSurv_imp_lib_no_mort <- read.csv("binSurv_imp_lib_no_mort.csv")
binSurv_imp_mort_no_lib <- read.csv("binSurv_imp_mort_no_lib.csv")

#-------calculating lambda that controls type one error----------

binSurv_noeffect <- binSurv_noeffect %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(binSurv_noeffect) %>%
  filter(post_prob_death_BM != "Error") %>%
  mutate(post_p_death = as.numeric(post_prob_death_BM),
         post_p_extubation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_extubation_surv = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

binSurv_noeffect <- binSurv_noeffect %>%
  mutate(max_post = pmax(post_p_death, post_p_extubation))

binSurv_threshold <- binSurv_noeffect %>%
  filter(post_p_death >= 0.2 & post_p_extubation >= 0.2) %>%
  arrange(desc(max_post))

first_100 <- head(binSurv_threshold, nrow(binSurv_noeffect)*0.05)

#0.9666774
lambda_binSurv <- as.numeric(first_100[nrow(first_100),]$max_post)
lambda_binSurv <- 0.975

(binSurv_post_p_null <- sum(binSurv_threshold$max_post > lambda_binSurv)/nrow(binSurv_noeffect))


#--------based on the lambda, calculating "power" for each scenario----------

################## worsening both ################################

binSurv_worsen_both <- binSurv_worsen_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(binSurv_noeffect) %>%
  filter(post_prob_death_BM != "Error") %>%
  mutate(post_p_death = as.numeric(post_prob_death_BM),
         post_p_extubation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_extubation_surv = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv),
         max_post = pmax(post_p_death, post_p_extubation)) 

binSurv_worsen_both_power <- binSurv_worsen_both %>%
  filter(post_p_death >= 0.2 & post_p_extubation >= 0.2) %>%
  arrange(desc(max_post))

(binSurv_post_p_worsen_both <- nrow(binSurv_worsen_both_power[binSurv_worsen_both_power$max_post > lambda_binSurv,])/nrow(binSurv_worsen_both))

################## improving both ################################

binSurv_imp_both <- binSurv_imp_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(binSurv_noeffect) %>%
  filter(post_prob_death_BM != "Error") %>%
  mutate(post_p_death = as.numeric(post_prob_death_BM),
         post_p_extubation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_extubation_surv = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv),
         max_post = pmax(post_p_death, post_p_extubation)) 

binSurv_imp_both_power <- binSurv_imp_both %>%
  filter(post_p_death >= 0.2 & post_p_extubation >= 0.2) %>%
  arrange(desc(max_post))

(binSurv_post_p_imp_both <- nrow(binSurv_imp_both_power[binSurv_imp_both_power$max_post > lambda_binSurv,])/nrow(binSurv_imp_both))


################## improving liberation worsening mortality ################################
binSurv_imp_lib_worse_mort <- binSurv_imp_lib_worse_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(binSurv_noeffect) %>%
  filter(post_prob_death_BM != "Error") %>%
  mutate(post_p_death = as.numeric(post_prob_death_BM),
         post_p_extubation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_extubation_surv = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv),
         max_post = pmax(post_p_death, post_p_extubation)) 

binSurv_imp_lib_worse_mort_power <- binSurv_imp_lib_worse_mort %>%
  filter(post_p_death >= 0.2 & post_p_extubation >= 0.2) %>%
  arrange(desc(max_post))

(binSurv_post_p_imp_lib_worse_mort <- nrow(binSurv_imp_lib_worse_mort_power[binSurv_imp_lib_worse_mort_power$max_post > lambda_binSurv,])/nrow(binSurv_imp_lib_worse_mort))


################## improving mortality worsening liberation ################################
binSurv_imp_mort_worse_lib <- binSurv_imp_mort_worse_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(binSurv_noeffect) %>%
  filter(post_prob_death_BM != "Error") %>%
  mutate(post_p_death = as.numeric(post_prob_death_BM),
         post_p_extubation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_extubation_surv = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv),
         max_post = pmax(post_p_death, post_p_extubation)) 

binSurv_imp_mort_worse_lib_power <- binSurv_imp_mort_worse_lib %>%
  filter(post_p_death >= 0.2 & post_p_extubation >= 0.2) %>%
  arrange(desc(max_post))

(binSurv_post_p_imp_mort_worse_lib <- nrow(binSurv_imp_mort_worse_lib_power[binSurv_imp_mort_worse_lib_power$max_post > lambda_binSurv,])/nrow(binSurv_imp_mort_worse_lib))

################## improving liberation noning mortality ################################
binSurv_imp_lib_no_mort <- binSurv_imp_lib_no_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(binSurv_noeffect) %>%
  filter(post_prob_death_BM != "Error") %>%
  mutate(post_p_death = as.numeric(post_prob_death_BM),
         post_p_extubation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_extubation_surv = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv),
         max_post = pmax(post_p_death, post_p_extubation)) 

binSurv_imp_lib_no_mort_power <- binSurv_imp_lib_no_mort %>%
  filter(post_p_death >= 0.2 & post_p_extubation >= 0.2) %>%
  arrange(desc(max_post))

(binSurv_post_p_imp_lib_no_mort <- nrow(binSurv_imp_lib_no_mort_power[binSurv_imp_lib_no_mort_power$max_post > lambda_binSurv,])/nrow(binSurv_imp_lib_no_mort))


################## improving mortality noning liberation ################################
binSurv_imp_mort_no_lib <- binSurv_imp_mort_no_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(binSurv_noeffect) %>%
  filter(post_prob_death_BM != "Error") %>%
  mutate(post_p_death = as.numeric(post_prob_death_BM),
         post_p_extubation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_extubation_surv = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv),
         max_post = pmax(post_p_death, post_p_extubation)) 

binSurv_imp_mort_no_lib_power <- binSurv_imp_mort_no_lib %>%
  filter(post_p_death >= 0.2 & post_p_extubation >= 0.2) %>%
  arrange(desc(max_post))

(binSurv_post_p_imp_mort_no_lib <- nrow(binSurv_imp_mort_no_lib_power[binSurv_imp_mort_no_lib_power$max_post > lambda_binSurv,])/nrow(binSurv_imp_mort_no_lib))



binSurv_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "worsen_lib_imp_mort", "imp_lib_no_mort", "imp_mort_no_lib"),
                        post_p = c(binSurv_post_p_null, binSurv_post_p_worsen_both, binSurv_post_p_imp_both, binSurv_post_p_imp_lib_worse_mort, binSurv_post_p_imp_mort_worse_lib, binSurv_post_p_imp_lib_no_mort, binSurv_post_p_imp_mort_no_lib),
                        trt_death_est_mean = c(mean(binSurv_noeffect$trt_death), mean(binSurv_worsen_both$trt_death), mean(binSurv_imp_both$trt_death), mean(binSurv_imp_lib_worse_mort$trt_death), mean(binSurv_imp_mort_worse_lib$trt_death), mean(binSurv_imp_lib_no_mort$trt_death), mean(binSurv_imp_mort_no_lib$trt_death)),
                        trt_extubation_est_mean = c(mean(binSurv_noeffect$trt_extubation_surv), mean(binSurv_worsen_both$trt_extubation_surv), mean(binSurv_imp_both$trt_extubation_surv), mean(binSurv_imp_lib_worse_mort$trt_extubation_surv), mean(binSurv_imp_mort_worse_lib$trt_extubation_surv), mean(binSurv_imp_lib_no_mort$trt_extubation_surv), mean(binSurv_imp_mort_no_lib$trt_extubation_surv)),
                        sd_death_est_mean = c(mean(binSurv_noeffect$sd_death), mean(binSurv_worsen_both$sd_death), mean(binSurv_imp_both$sd_death), mean(binSurv_imp_lib_worse_mort$sd_death), mean(binSurv_imp_mort_worse_lib$sd_death), mean(binSurv_imp_lib_no_mort$sd_death), mean(binSurv_imp_mort_no_lib$sd_death)),
                        sd_extubation_est_mean = c(mean(binSurv_noeffect$sd_extubation), mean(binSurv_worsen_both$sd_extubation), mean(binSurv_imp_both$sd_extubation), mean(binSurv_imp_lib_worse_mort$sd_extubation), mean(binSurv_imp_mort_worse_lib$sd_extubation), mean(binSurv_imp_lib_no_mort$sd_extubation), mean(binSurv_imp_mort_no_lib$sd_extubation)))
binSurv_1[, 1:2]
sqrt(binSurv_1$post_p*(1 - binSurv_1$post_p)/nrow(binSurv_imp_mort_no_lib))[c(1,3,2,5,4,7,6)]
