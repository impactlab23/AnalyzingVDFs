library(tidyverse)
library(dplyr)


setwd("/GitHub/AnalyzingVDFs/Weilbull_DGP")
wdname <- paste(getwd(),"/7 - CR", sep = "")
setwd(wdname)

CR_noeffect <- read.csv("CR_noeffect.csv")
CR_imp_both <- read.csv("CR_imp_both.csv")
CR_imp_lib_worse_mort <- read.csv("CR_imp_lib_worse_mort.csv")
CR_imp_mort_worse_lib <- read.csv("CR_imp_mort_worse_lib.csv")
CR_worsen_both <- read.csv("CR_worsen_both.csv")

CR_imp_lib_no_mort <- read.csv("CR_imp_lib_no_mort.csv")
CR_imp_mort_no_lib <- read.csv("CR_imp_mort_no_lib.csv")

#-------calculating lambda that controls type one error----------

CR_noeffect <- CR_noeffect %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(CR_noeffect) %>%
  filter(trt_death_CR != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_CR),
         post_prob_liberation = as.numeric(post_prob_liberation_CR),
         trt_death = as.numeric(trt_death_CR),
         trt_liberation = as.numeric(trt_liberation_CR),
         sd_death = as.numeric(sd_death_CR),
         sd_extubation = as.numeric(sd_liberation_CR)) 

CR_noeffect <- CR_noeffect %>%
  mutate(max_post = pmax(post_prob_death, post_prob_liberation))
n_sim <- nrow(CR_noeffect)

CR_threshold <- CR_noeffect %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

first_100 <- head(CR_threshold, nrow(CR_noeffect)*0.05)

#0.9666774
lambda_CR <- as.numeric(first_100[nrow(first_100),]$max_post)
lambda_CR <- 0.96

(CR_post_p_null <- sum(CR_threshold$max_post > lambda_CR)/n_sim)


#--------based on the lambda, calculating "power" for each scenario----------

################## worsening both ################################

CR_worsen_both <- CR_worsen_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(CR_noeffect) %>%
  filter(trt_death_CR != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_CR),
         post_prob_liberation = as.numeric(post_prob_liberation_CR),
         trt_death = as.numeric(trt_death_CR),
         trt_liberation = as.numeric(trt_liberation_CR),
         sd_death = as.numeric(sd_death_CR),
         sd_extubation = as.numeric(sd_liberation_CR),
         max_post = pmax(post_prob_death, post_prob_liberation)) 

CR_worsen_both_power <- CR_worsen_both %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(CR_post_p_worsen_both <- nrow(CR_worsen_both_power[CR_worsen_both_power$max_post > lambda_CR,])/nrow(CR_worsen_both))

################## improving both ################################

CR_imp_both <- CR_imp_both %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(CR_noeffect) %>%
  filter(trt_death_CR != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_CR),
         post_prob_liberation = as.numeric(post_prob_liberation_CR),
         trt_death = as.numeric(trt_death_CR),
         trt_liberation = as.numeric(trt_liberation_CR),
         sd_death = as.numeric(sd_death_CR),
         sd_extubation = as.numeric(sd_liberation_CR),
         max_post = pmax(post_prob_death, post_prob_liberation)) 

CR_imp_both_power <- CR_imp_both %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(CR_post_p_imp_both <- nrow(CR_imp_both_power[CR_imp_both_power$max_post > lambda_CR,])/nrow(CR_imp_both))


################## improving liberation worsening mortality ################################
CR_imp_lib_worse_mort <- CR_imp_lib_worse_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(CR_noeffect) %>%
  filter(trt_death_CR != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_CR),
         post_prob_liberation = as.numeric(post_prob_liberation_CR),
         trt_death = as.numeric(trt_death_CR),
         trt_liberation = as.numeric(trt_liberation_CR),
         sd_death = as.numeric(sd_death_CR),
         sd_extubation = as.numeric(sd_liberation_CR),
         max_post = pmax(post_prob_death, post_prob_liberation)) 

CR_imp_lib_worse_mort_power <- CR_imp_lib_worse_mort %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(CR_post_p_imp_lib_worse_mort <- nrow(CR_imp_lib_worse_mort_power[CR_imp_lib_worse_mort_power$max_post > lambda_CR,])/nrow(CR_imp_lib_worse_mort))


################## improving mortality worsening liberation ################################
CR_imp_mort_worse_lib <- CR_imp_mort_worse_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(CR_noeffect) %>%
  filter(trt_death_CR != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_CR),
         post_prob_liberation = as.numeric(post_prob_liberation_CR),
         trt_death = as.numeric(trt_death_CR),
         trt_liberation = as.numeric(trt_liberation_CR),
         sd_death = as.numeric(sd_death_CR),
         sd_extubation = as.numeric(sd_liberation_CR),
         max_post = pmax(post_prob_death, post_prob_liberation)) 

CR_imp_mort_worse_lib_power <- CR_imp_mort_worse_lib %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(CR_post_p_imp_mort_worse_lib <- nrow(CR_imp_mort_worse_lib_power[CR_imp_mort_worse_lib_power$max_post > lambda_CR,])/nrow(CR_imp_mort_worse_lib))


################## improving liberation noning mortality ################################
CR_imp_lib_no_mort <- CR_imp_lib_no_mort %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(CR_noeffect) %>%
  filter(trt_death_CR != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_CR),
         post_prob_liberation = as.numeric(post_prob_liberation_CR),
         trt_death = as.numeric(trt_death_CR),
         trt_liberation = as.numeric(trt_liberation_CR),
         sd_death = as.numeric(sd_death_CR),
         sd_extubation = as.numeric(sd_liberation_CR),
         max_post = pmax(post_prob_death, post_prob_liberation)) 

CR_imp_lib_no_mort_power <- CR_imp_lib_no_mort %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(CR_post_p_imp_lib_no_mort <- nrow(CR_imp_lib_no_mort_power[CR_imp_lib_no_mort_power$max_post > lambda_CR,])/nrow(CR_imp_lib_no_mort))


################## improving mortality noning liberation ################################
CR_imp_mort_no_lib <- CR_imp_mort_no_lib %>% 
  mutate(X = seq(1, 2000, 1)) %>% 
  na.omit(CR_noeffect) %>%
  filter(trt_death_CR != "Error") %>%
  mutate(post_prob_death = as.numeric(post_prob_death_CR),
         post_prob_liberation = as.numeric(post_prob_liberation_CR),
         trt_death = as.numeric(trt_death_CR),
         trt_liberation = as.numeric(trt_liberation_CR),
         sd_death = as.numeric(sd_death_CR),
         sd_extubation = as.numeric(sd_liberation_CR),
         max_post = pmax(post_prob_death, post_prob_liberation)) 

CR_imp_mort_no_lib_power <- CR_imp_mort_no_lib %>%
  filter(post_prob_death >= 0.2 & post_prob_liberation >= 0.2) %>%
  arrange(desc(max_post))

(CR_post_p_imp_mort_no_lib <- nrow(CR_imp_mort_no_lib_power[CR_imp_mort_no_lib_power$max_post > lambda_CR,])/nrow(CR_imp_mort_no_lib))



CR_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "worsen_lib_imp_mort", "imp_lib_no_mort", "imp_mort_no_lib"),
                   post_p = c(CR_post_p_null, CR_post_p_worsen_both, CR_post_p_imp_both, CR_post_p_imp_lib_worse_mort, CR_post_p_imp_mort_worse_lib, CR_post_p_imp_lib_no_mort, CR_post_p_imp_mort_no_lib),
                   trt_death_est_mean = c(mean(CR_noeffect$trt_death), mean(CR_worsen_both$trt_death), mean(CR_imp_both$trt_death), mean(CR_imp_lib_worse_mort$trt_death), mean(CR_imp_mort_worse_lib$trt_death), mean(CR_imp_lib_no_mort$trt_death), mean(CR_imp_mort_no_lib$trt_death)),
                   trt_liberation_est_mean = c(mean(CR_noeffect$trt_liberation), mean(CR_worsen_both$trt_liberation), mean(CR_imp_both$trt_liberation), mean(CR_imp_lib_worse_mort$trt_liberation), mean(CR_imp_mort_worse_lib$trt_liberation), mean(CR_imp_lib_no_mort$trt_liberation), mean(CR_imp_mort_no_lib$trt_liberation)),
                   sd_death_est_mean = c(mean(CR_noeffect$sd_death), mean(CR_worsen_both$sd_death), mean(CR_imp_both$sd_death), mean(CR_imp_lib_worse_mort$sd_death), mean(CR_imp_mort_worse_lib$sd_death), mean(CR_imp_lib_no_mort$sd_death), mean(CR_imp_mort_no_lib$sd_death)),
                   sd_liberation_est_mean = c(mean(CR_noeffect$sd_extubation), mean(CR_worsen_both$sd_extubation), mean(CR_imp_both$sd_extubation), mean(CR_imp_lib_worse_mort$sd_extubation), mean(CR_imp_mort_worse_lib$sd_extubation), mean(CR_imp_lib_no_mort$sd_extubation), mean(CR_imp_mort_no_lib$sd_extubation)))
CR_1[, 1:2]
sqrt(CR_1$post_p*(1 - CR_1$post_p)/nrow(CR_imp_mort_no_lib))[c(1,3,2,5,4,7,6)]
