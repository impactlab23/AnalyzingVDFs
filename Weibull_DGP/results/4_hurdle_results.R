library(tidyverse)
library(dplyr)


setwd("/GitHub/AnalyzingVDFs/Weibull_DGP")
wdname <- paste(getwd(),"/4 - hurdle", sep = "")
setwd(wdname)

hurdle_noeffect <- read.csv("hurdle_noeffect.csv")
hurdle_imp_both <- read.csv("hurdle_imp_both.csv")
hurdle_imp_lib_worse_mort <- read.csv("hurdle_imp_lib_worse_mort.csv")
hurdle_imp_mort_worse_lib <- read.csv("hurdle_imp_mort_worse_lib.csv")
hurdle_worsen_both <- read.csv("hurdle_worsen_both.csv")

hurdle_imp_lib_no_mort <- read.csv("hurdle_imp_lib_no_mort.csv")
hurdle_imp_mort_no_lib <- read.csv("hurdle_imp_mort_no_lib.csv")

#-------calculating lambda that controls type one error----------

hurdle_noeffect <- hurdle_noeffect %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(hurdle_noeffect) %>%
  filter(V1 != "Error") %>%
  mutate(post_prob_death = as.numeric(V1),
         post_prob_liberation = as.numeric(V2),
         trt_death = as.numeric(V3),
         trt_liberation = as.numeric(V4),
         sd_death = as.numeric(V5),
         sd_extubation = as.numeric(V6)) 

hurdle_noeffect <- hurdle_noeffect %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

hurdle_threshold <- hurdle_noeffect %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

first_100 <- head(hurdle_threshold, nrow(hurdle_noeffect)*0.05)

#0.9666774
lambda_hurdle <- as.numeric(first_100[nrow(first_100),]$max_post)

lambda_hurdle <- 0.998

hurdle_noeffect_power <- hurdle_noeffect %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(hurdle_p_noeffect <- nrow(hurdle_noeffect_power[
  hurdle_noeffect_power$max_post >= lambda_hurdle,])/
    nrow(hurdle_noeffect))

#--------based on the lambda, calculating "power" for each scenario----------

################## worsening both ################################

hurdle_worsen_both <- hurdle_worsen_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(hurdle_worsen_both) %>%
  filter(V1 != "Error") %>%
  mutate(post_prob_death = as.numeric(V1),
         post_prob_liberation = as.numeric(V2),
         trt_death = as.numeric(V3),
         trt_liberation = as.numeric(V4),
         sd_death = as.numeric(V5),
         sd_extubation = as.numeric(V6)) 

hurdle_worsen_both <- hurdle_worsen_both %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

hurdle_worsen_both_power <- hurdle_worsen_both %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(hurdle_post_p_worsen_both <- nrow(hurdle_worsen_both_power[hurdle_worsen_both_power$max_post >= lambda_hurdle,])/nrow(hurdle_worsen_both))

################## improving both ################################

hurdle_imp_both <- hurdle_imp_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(hurdle_imp_both) %>%
  filter(V1 != "Error") %>%
  mutate(post_prob_death = as.numeric(V1),
         post_prob_liberation = as.numeric(V2),
         trt_death = as.numeric(V3),
         trt_liberation = as.numeric(V4),
         sd_death = as.numeric(V5),
         sd_extubation = as.numeric(V6)) 

hurdle_imp_both <- hurdle_imp_both %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

hurdle_imp_both_power <- hurdle_imp_both %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(hurdle_post_p_imp_both <- nrow(hurdle_imp_both_power[hurdle_imp_both_power$max_post >= lambda_hurdle,])/nrow(hurdle_imp_both))


################## improving liberation worsening mortality ################################
hurdle_imp_lib_worse_mort <- hurdle_imp_lib_worse_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(hurdle_imp_lib_worse_mort) %>%
  filter(V1 != "Error") %>%
  mutate(post_prob_death = as.numeric(V1),
         post_prob_liberation = as.numeric(V2),
         trt_death = as.numeric(V3),
         trt_liberation = as.numeric(V4),
         sd_death = as.numeric(V5),
         sd_extubation = as.numeric(V6))

hurdle_imp_lib_worse_mort <- hurdle_imp_lib_worse_mort %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

hurdle_imp_lib_worse_mort_power <- hurdle_imp_lib_worse_mort %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(hurdle_post_p_imp_lib_worse_mort <- nrow(hurdle_imp_lib_worse_mort_power[hurdle_imp_lib_worse_mort_power$max_post >= lambda_hurdle,])/nrow(hurdle_imp_lib_worse_mort))


################## improving mortality worsening liberation ################################
hurdle_imp_mort_worse_lib <- hurdle_imp_mort_worse_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(hurdle_imp_mort_worse_lib) %>%
  filter(V1 != "Error") %>%
  mutate(post_prob_death = as.numeric(V1),
         post_prob_liberation = as.numeric(V2),
         trt_death = as.numeric(V3),
         trt_liberation = as.numeric(V4),
         sd_death = as.numeric(V5),
         sd_extubation = as.numeric(V6)) 

hurdle_imp_mort_worse_lib <- hurdle_imp_mort_worse_lib %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

hurdle_imp_mort_worse_lib_power <- hurdle_imp_mort_worse_lib %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(hurdle_post_p_imp_mort_worse_lib <- nrow(hurdle_imp_mort_worse_lib_power[hurdle_imp_mort_worse_lib_power$max_post >= lambda_hurdle,])/nrow(hurdle_imp_mort_worse_lib))

################## improving liberation noning mortality ################################
hurdle_imp_lib_no_mort <- hurdle_imp_lib_no_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(hurdle_imp_lib_no_mort) %>%
  filter(V1 != "Error") %>%
  mutate(post_prob_death = as.numeric(V1),
         post_prob_liberation = as.numeric(V2),
         trt_death = as.numeric(V3),
         trt_liberation = as.numeric(V4),
         sd_death = as.numeric(V5),
         sd_extubation = as.numeric(V6))

hurdle_imp_lib_no_mort <- hurdle_imp_lib_no_mort %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

hurdle_imp_lib_no_mort_power <- hurdle_imp_lib_no_mort %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(hurdle_post_p_imp_lib_no_mort <- nrow(hurdle_imp_lib_no_mort_power[hurdle_imp_lib_no_mort_power$max_post >= lambda_hurdle,])/nrow(hurdle_imp_lib_no_mort))


################## improving mortality noning liberation ################################
hurdle_imp_mort_no_lib <- hurdle_imp_mort_no_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(hurdle_imp_mort_no_lib) %>%
  filter(V1 != "Error") %>%
  mutate(post_prob_death = as.numeric(V1),
         post_prob_liberation = as.numeric(V2),
         trt_death = as.numeric(V3),
         trt_liberation = as.numeric(V4),
         sd_death = as.numeric(V5),
         sd_extubation = as.numeric(V6)) 

hurdle_imp_mort_no_lib <- hurdle_imp_mort_no_lib %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))

hurdle_imp_mort_no_lib_power <- hurdle_imp_mort_no_lib %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(hurdle_post_p_imp_mort_no_lib <- nrow(hurdle_imp_mort_no_lib_power[hurdle_imp_mort_no_lib_power$max_post >= lambda_hurdle,])/nrow(hurdle_imp_mort_no_lib))


hurdle_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "worsen_lib_imp_mort", "imp_lib_no_mort", "imp_mort_no_lib"),
                       post_p = c(hurdle_p_noeffect, hurdle_post_p_worsen_both, hurdle_post_p_imp_both, hurdle_post_p_imp_lib_worse_mort, hurdle_post_p_imp_mort_worse_lib, hurdle_post_p_imp_lib_no_mort, hurdle_post_p_imp_mort_no_lib),
                       trt_death_est_mean = c(mean(hurdle_noeffect$trt_death), mean(hurdle_worsen_both$trt_death), mean(hurdle_imp_both$trt_death), mean(hurdle_imp_lib_worse_mort$trt_death), mean(hurdle_imp_mort_worse_lib$trt_death), mean(hurdle_imp_lib_no_mort$trt_death), mean(hurdle_imp_mort_no_lib$trt_death)),
                       trt_liberation_est_mean = c(mean(hurdle_noeffect$trt_liberation), mean(hurdle_worsen_both$trt_liberation), mean(hurdle_imp_both$trt_liberation), mean(hurdle_imp_lib_worse_mort$trt_liberation), mean(hurdle_imp_mort_worse_lib$trt_liberation), mean(hurdle_imp_lib_no_mort$trt_liberation), mean(hurdle_imp_mort_no_lib$trt_liberation)),
                       sd_death_est_mean = c(mean(hurdle_noeffect$sd_death), mean(hurdle_worsen_both$sd_death), mean(hurdle_imp_both$sd_death), mean(hurdle_imp_lib_worse_mort$sd_death), mean(hurdle_imp_mort_worse_lib$sd_death), mean(hurdle_imp_lib_no_mort$sd_death), mean(hurdle_imp_mort_no_lib$sd_death)),
                       sd_liberation_est_mean = c(mean(hurdle_noeffect$sd_extubation), mean(hurdle_worsen_both$sd_extubation), mean(hurdle_imp_both$sd_extubation), mean(hurdle_imp_lib_worse_mort$sd_extubation), mean(hurdle_imp_mort_worse_lib$sd_extubation), mean(hurdle_imp_lib_no_mort$sd_extubation), mean(hurdle_imp_mort_no_lib$sd_extubation)))
hurdle_1[, 1:2]
sqrt(hurdle_1$post_p*(1 - hurdle_1$post_p)/nrow(hurdle_imp_mort_no_lib))[c(1,3,2,5,4,7,6)]
