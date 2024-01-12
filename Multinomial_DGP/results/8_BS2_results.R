library(tidyverse)
library(dplyr)


setwd("/GitHub/AnalyzingVDFs/Multinomial_DGP")
wdname <- paste(getwd(),"/8 - BS2", sep = "")
setwd(wdname)

BS2_noeffect <- read.csv("BS2_noeffect.csv")
BS2_imp_both <- read.csv("BS2_imp_both.csv")
BS2_imp_lib_worse_mort <- read.csv("BS2_imp_lib_worse_mort.csv")
BS2_imp_mort_worse_lib <- read.csv("BS2_imp_mort_worse_lib.csv")
BS2_worsen_both <- read.csv("BS2_worsen_both.csv")
BS2_imp_mort_no_lib <- read.csv("BS2_imp_mort_no_lib.csv")
BS2_imp_lib_no_mort <- read.csv("BS2_imp_lib_no_mort.csv")

#-------calculating lambda that controls type one error----------

BS2_noeffect <- BS2_noeffect %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS2_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS2_noeffect <- BS2_noeffect %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS2_threshold <- BS2_noeffect %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

first_100 <- head(BS2_threshold, nrow(BS2_noeffect)*0.05)

#0.9666774
lambda_BS2 <- as.numeric(first_100[nrow(first_100),]$max_post)

lambda_BS2 <- 0.96

(BS2_post_p_null <- sum(BS2_threshold$max_post > lambda_BS2)/nrow(BS2_noeffect))

#--------based on the lambda, calculating "power" for each scenario----------

################## worsening both ################################

BS2_worsen_both <- BS2_worsen_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS2_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS2_worsen_both <- BS2_worsen_both %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS2_worsen_both_power <- BS2_worsen_both %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS2_post_p_worsen_both <- nrow(BS2_worsen_both_power[BS2_worsen_both_power$max_post > lambda_BS2,])/nrow(BS2_worsen_both))

################## improving both ################################

BS2_imp_both <- BS2_imp_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS2_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS2_imp_both <- BS2_imp_both %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS2_imp_both_power <- BS2_imp_both %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS2_post_p_imp_both <- nrow(BS2_imp_both_power[BS2_imp_both_power$max_post > lambda_BS2,])/nrow(BS2_imp_both))


################## improving liberation worsening mortality ################################
BS2_imp_lib_worse_mort <- BS2_imp_lib_worse_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS2_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS2_imp_lib_worse_mort <- BS2_imp_lib_worse_mort %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS2_imp_lib_worse_mort_power <- BS2_imp_lib_worse_mort %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS2_post_p_imp_lib_worse_mort <- nrow(BS2_imp_lib_worse_mort_power[BS2_imp_lib_worse_mort_power$max_post > lambda_BS2,])/nrow(BS2_imp_lib_worse_mort))

################## improving mortality worsening liberation ################################
BS2_imp_mort_worse_lib <- BS2_imp_mort_worse_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS2_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS2_imp_mort_worse_lib <- BS2_imp_mort_worse_lib %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS2_imp_mort_worse_lib_power <- BS2_imp_mort_worse_lib %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS2_post_p_imp_mort_worse_lib <- nrow(BS2_imp_mort_worse_lib_power[BS2_imp_mort_worse_lib_power$max_post > lambda_BS2,])/nrow(BS2_imp_mort_worse_lib))


################## improving mortality no liberation ################################
BS2_imp_mort_no_lib <- BS2_imp_mort_no_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS2_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS2_imp_mort_no_lib <- BS2_imp_mort_no_lib %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS2_imp_mort_no_lib_power <- BS2_imp_mort_no_lib %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS2_post_p_imp_mort_no_lib <- nrow(BS2_imp_mort_no_lib_power[BS2_imp_mort_no_lib_power$max_post > lambda_BS2,])/nrow(BS2_imp_mort_no_lib))

################## improving liberation no mortality ################################
BS2_imp_lib_no_mort <- BS2_imp_lib_no_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS2_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS2_imp_lib_no_mort <- BS2_imp_lib_no_mort %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS2_imp_lib_no_mort_power <- BS2_imp_lib_no_mort %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS2_post_p_imp_lib_no_mort <- nrow(BS2_imp_lib_no_mort_power[BS2_imp_lib_no_mort_power$max_post > lambda_BS2,])/nrow(BS2_imp_lib_no_mort))


BS2_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "imp_mort_worse_lib", "imp_mort_no_lib", "imp_lib_no_mort"),
                    post_p = c(BS2_post_p_null, BS2_post_p_worsen_both, BS2_post_p_imp_both, BS2_post_p_imp_lib_worse_mort, BS2_post_p_imp_mort_worse_lib, BS2_post_p_imp_mort_no_lib, BS2_post_p_imp_lib_no_mort),
                    trt_death_est_mean = c(mean(BS2_noeffect$trt_death), mean(BS2_worsen_both$trt_death), mean(BS2_imp_both$trt_death), mean(BS2_imp_lib_worse_mort$trt_death), mean(BS2_imp_mort_worse_lib$trt_death), mean(BS2_imp_mort_no_lib$trt_death), mean(BS2_imp_lib_no_mort$trt_death)),
                    trt_liberation_est_mean = c(mean(BS2_noeffect$trt_liberation), mean(BS2_worsen_both$trt_liberation), mean(BS2_imp_both$trt_liberation), mean(BS2_imp_lib_worse_mort$trt_liberation), mean(BS2_imp_mort_worse_lib$trt_liberation), mean(BS2_imp_mort_no_lib$trt_liberation), mean(BS2_imp_lib_no_mort$trt_liberation)),
                    sd_death_est_mean = c(mean(BS2_noeffect$sd_death), mean(BS2_worsen_both$sd_death), mean(BS2_imp_both$sd_death), mean(BS2_imp_lib_worse_mort$sd_death), mean(BS2_imp_mort_worse_lib$sd_death), mean(BS2_imp_mort_no_lib$sd_death), mean(BS2_imp_lib_no_mort$sd_death)),
                    sd_liberation_est_mean = c(mean(BS2_noeffect$sd_extubation), mean(BS2_worsen_both$sd_extubation), mean(BS2_imp_both$sd_extubation), mean(BS2_imp_lib_worse_mort$sd_extubation), mean(BS2_imp_mort_worse_lib$sd_extubation), mean(BS2_imp_mort_no_lib$sd_extubation), mean(BS2_imp_lib_no_mort$sd_extubation)))

cbind(BS2_1[, 1:2], se = sqrt(BS2_1$post_p*(1 - BS2_1$post_p)/nrow(BS2_imp_lib_no_mort)))
