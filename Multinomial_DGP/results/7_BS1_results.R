library(tidyverse)
library(dplyr)

setwd("/GitHub/AnalyzingVDFs/Multinomial_DGP")
wdname <- paste(getwd(),"/7 - BS1", sep = "")
setwd(wdname)

BS1_noeffect <- read.csv("BS1_noeffect.csv")
BS1_imp_both <- read.csv("BS1_imp_both.csv")
BS1_imp_lib_worse_mort <- read.csv("BS1_imp_lib_worse_mort.csv")
BS1_imp_mort_worse_lib <- read.csv("BS1_imp_mort_worse_lib.csv")
BS1_worsen_both <- read.csv("BS1_worsen_both.csv")
BS1_imp_mort_no_lib <- read.csv("BS1_imp_mort_no_lib.csv")
BS1_imp_lib_no_mort <- read.csv("BS1_imp_lib_no_mort.csv")

#-------calculating lambda that controls type one error----------

BS1_noeffect <- BS1_noeffect %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS1_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS1_noeffect <- BS1_noeffect %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS1_threshold <- BS1_noeffect %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

first_100 <- head(BS1_threshold, nrow(BS1_noeffect)*0.05)

#0.9666774
lambda_BS1 <- as.numeric(first_100[nrow(first_100),]$max_post)
lambda_BS1 <- 0.965

(BS1_post_p_null <- sum(BS1_threshold$max_post > lambda_BS1)/nrow(BS1_noeffect))


#--------based on the lambda, calculating "power" for each scenario----------

################## worsening both ################################

BS1_worsen_both <- BS1_worsen_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS1_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS1_worsen_both <- BS1_worsen_both %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS1_worsen_both_power <- BS1_worsen_both %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS1_post_p_worsen_both <- nrow(BS1_worsen_both_power[BS1_worsen_both_power$max_post > lambda_BS1,])/nrow(BS1_worsen_both))

################## improving both ################################

BS1_imp_both <- BS1_imp_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS1_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS1_imp_both <- BS1_imp_both %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS1_imp_both_power <- BS1_imp_both %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS1_post_p_imp_both <- nrow(BS1_imp_both_power[BS1_imp_both_power$max_post > lambda_BS1,])/nrow(BS1_imp_both))


################## improving liberation worsening mortality ################################
BS1_imp_lib_worse_mort <- BS1_imp_lib_worse_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS1_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS1_imp_lib_worse_mort <- BS1_imp_lib_worse_mort %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS1_imp_lib_worse_mort_power <- BS1_imp_lib_worse_mort %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS1_post_p_imp_lib_worse_mort <- nrow(BS1_imp_lib_worse_mort_power[BS1_imp_lib_worse_mort_power$max_post > lambda_BS1,])/nrow(BS1_imp_lib_worse_mort))


################## improving mortality no liberation ################################
BS1_imp_mort_worse_lib <- BS1_imp_mort_worse_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS1_worseeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS1_imp_mort_worse_lib <- BS1_imp_mort_worse_lib %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS1_imp_mort_worse_lib_power <- BS1_imp_mort_worse_lib %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS1_post_p_imp_mort_worse_lib <- nrow(BS1_imp_mort_worse_lib_power[BS1_imp_mort_worse_lib_power$max_post > lambda_BS1,])/nrow(BS1_imp_mort_worse_lib))

################## improving mortality no liberation ################################
BS1_imp_mort_no_lib <- BS1_imp_mort_no_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS1_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS1_imp_mort_no_lib <- BS1_imp_mort_no_lib %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS1_imp_mort_no_lib_power <- BS1_imp_mort_no_lib %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS1_post_p_imp_mort_no_lib <- nrow(BS1_imp_mort_no_lib_power[BS1_imp_mort_no_lib_power$max_post > lambda_BS1,])/nrow(BS1_imp_mort_no_lib))

################## improving mortality no liberation ################################
BS1_imp_lib_no_mort <- BS1_imp_lib_no_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(BS1_noeffect) %>%
  filter(trt_death_BM != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_BM),
         post_prob_liberation = as.numeric(post_prob_extubation_surv),
         trt_death = as.numeric(trt_death_BM),
         trt_liberation = as.numeric(trt_extubation_surv),
         sd_death = as.numeric(sd_death_BM),
         sd_extubation = as.numeric(sd_extubation_surv)) 

BS1_imp_lib_no_mort <- BS1_imp_lib_no_mort %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

BS1_imp_lib_no_mort_power <- BS1_imp_lib_no_mort %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(BS1_post_p_imp_lib_no_mort <- nrow(BS1_imp_lib_no_mort_power[BS1_imp_lib_no_mort_power$max_post > lambda_BS1,])/nrow(BS1_imp_lib_no_mort))

BS1_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "imp_mort_worse_lib", "imp_mort_no_lib", "imp_lib_no_mort"),
                    post_p = c(BS1_post_p_null, BS1_post_p_worsen_both, BS1_post_p_imp_both, BS1_post_p_imp_lib_worse_mort, BS1_post_p_imp_mort_worse_lib, BS1_post_p_imp_mort_no_lib, BS1_post_p_imp_lib_no_mort),
                    trt_death_est_mean = c(mean(BS1_noeffect$trt_death), mean(BS1_worsen_both$trt_death), mean(BS1_imp_both$trt_death), mean(BS1_imp_lib_worse_mort$trt_death), mean(BS1_imp_mort_worse_lib$trt_death), mean(BS1_imp_mort_no_lib$trt_death), mean(BS1_imp_lib_no_mort$trt_death)),
                    trt_liberation_est_mean = c(mean(BS1_noeffect$trt_liberation), mean(BS1_worsen_both$trt_liberation), mean(BS1_imp_both$trt_liberation), mean(BS1_imp_lib_worse_mort$trt_liberation), mean(BS1_imp_mort_worse_lib$trt_liberation), mean(BS1_imp_mort_no_lib$trt_liberation), mean(BS1_imp_lib_no_mort$trt_liberation)),
                    sd_death_est_mean = c(mean(BS1_noeffect$sd_death), mean(BS1_worsen_both$sd_death), mean(BS1_imp_both$sd_death), mean(BS1_imp_lib_worse_mort$sd_death), mean(BS1_imp_mort_worse_lib$sd_death), mean(BS1_imp_mort_no_lib$sd_death), mean(BS1_imp_lib_no_mort$sd_death)),
                    sd_liberation_est_mean = c(mean(BS1_noeffect$sd_extubation), mean(BS1_worsen_both$sd_extubation), mean(BS1_imp_both$sd_extubation), mean(BS1_imp_lib_worse_mort$sd_extubation), mean(BS1_imp_mort_worse_lib$sd_extubation), mean(BS1_imp_mort_no_lib$sd_extubation), mean(BS1_imp_lib_no_mort$sd_extubation)))
BS1_1[, 1:2]
sqrt(BS1_1$post_p*(1 - BS1_1$post_p)/nrow(BS1_imp_lib_no_mort))[c(1,3,2,5,4,6,7)]
