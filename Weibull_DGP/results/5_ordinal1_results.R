library(tidyverse)
library(dplyr)


setwd("/GitHub/AnalyzingVDFs/Weibull_DGP")
wdname <- paste(getwd(),"/5 - ordinal1", sep = "")
setwd(wdname)

Ordinal1_noeffect <- read.csv("Practical_ordinal1_noeffect.csv")
Ordinal1_imp_both <- read.csv("Practical_ordinal1_imp_both.csv")
Ordinal1_imp_lib_worse_mort <- read.csv("Practical_ordinal1_imp_lib_worse_mort.csv")
Ordinal1_imp_mort_worse_lib <- read.csv("Practical_ordinal1_imp_mort_worse_lib.csv")
Ordinal1_worsen_both <- read.csv("Practical_ordinal1_worsen_both.csv")
Ordinal1_imp_lib_no_mort <- read.csv("Practical_ordinal1_imp_lib_no_mort.csv")
Ordinal1_imp_mort_no_lib <- read.csv("Practical_ordinal1_imp_mort_no_lib.csv")

#-------calculating lambda that controls type one error----------

Ordinal1_noeffect <- na.omit(Ordinal1_noeffect) %>% 
  filter(V1 != "Error") 
Ordinal1_noeffect <- Ordinal1_noeffect %>%
  mutate(X = seq(1, nrow(Ordinal1_noeffect), 1)) %>%
  rename(post_p = V1,
         trt_effect = V2,
         sd = V3) 

Ordinal1_noeffect <- Ordinal1_noeffect %>%
  arrange(desc(post_p))

#0.9605856 
(lambda_Ordinal1 <- quantile(as.numeric(Ordinal1_noeffect$post_p), probs = c(0.95)))

lambda_Ordinal1 <- 0.955
(Ordinal1_post_p_noeffect <- nrow(Ordinal1_noeffect[Ordinal1_noeffect$post_p > lambda_Ordinal1,])/nrow(Ordinal1_noeffect))


#--------based on the lambda, calculating "power" for each scenario----------
Ordinal1_worsen_both <- na.omit(Ordinal1_worsen_both) %>% 
  filter(V1 != "Error") 
Ordinal1_worsen_both <- Ordinal1_worsen_both %>% 
  mutate(X = seq(1, nrow(Ordinal1_worsen_both), 1)) %>%
  rename(post_p = V1,
         trt_effect = V2,
         sd = V3)
(Ordinal1_post_p_worsen_both <- nrow(Ordinal1_worsen_both[Ordinal1_worsen_both$post_p > lambda_Ordinal1,])/nrow(Ordinal1_worsen_both))

Ordinal1_imp_both <- na.omit(Ordinal1_imp_both) %>% 
  filter(V1 != "Error") 
Ordinal1_imp_both <- Ordinal1_imp_both %>% 
  mutate(X = seq(1, nrow(Ordinal1_imp_both), 1)) %>%
  rename(post_p = V1,
         trt_effect = V2,
         sd = V3)
(Ordinal1_post_p_imp_both <- nrow(Ordinal1_imp_both[Ordinal1_imp_both$post_p > lambda_Ordinal1,])/nrow(Ordinal1_imp_both))

Ordinal1_imp_lib_worse_mort <- na.omit(Ordinal1_imp_lib_worse_mort) %>% 
  filter(V1 != "Error") 
Ordinal1_imp_lib_worse_mort <- Ordinal1_imp_lib_worse_mort %>% 
  mutate(X = seq(1, nrow(Ordinal1_imp_lib_worse_mort), 1)) %>%
  rename(post_p = V1,
         trt_effect = V2,
         sd = V3)
(Ordinal1_post_imp_lib_worse_mort <- nrow(Ordinal1_imp_lib_worse_mort[Ordinal1_imp_lib_worse_mort$post_p > lambda_Ordinal1,])/nrow(Ordinal1_imp_lib_worse_mort))

Ordinal1_imp_mort_worse_lib <- na.omit(Ordinal1_imp_mort_worse_lib) %>% 
  filter(V1 != "Error") 
Ordinal1_imp_mort_worse_lib <- Ordinal1_imp_mort_worse_lib %>% 
  mutate(X = seq(1, nrow(Ordinal1_imp_mort_worse_lib), 1)) %>%
  rename(post_p = V1,
         trt_effect = V2,
         sd = V3)
(Ordinal1_post_imp_mort_worse_lib <- nrow(Ordinal1_imp_mort_worse_lib[Ordinal1_imp_mort_worse_lib$post_p > lambda_Ordinal1,])/nrow(Ordinal1_imp_mort_worse_lib))

Ordinal1_imp_lib_no_mort <- na.omit(Ordinal1_imp_lib_no_mort) %>% 
  filter(V1 != "Error") 
Ordinal1_imp_lib_no_mort <- Ordinal1_imp_lib_no_mort %>% 
  mutate(X = seq(1, nrow(Ordinal1_imp_lib_no_mort), 1)) %>%
  rename(post_p = V1,
         trt_effect = V2,
         sd = V3)
(Ordinal1_post_imp_lib_no_mort <- nrow(Ordinal1_imp_lib_no_mort[Ordinal1_imp_lib_no_mort$post_p > lambda_Ordinal1,])/nrow(Ordinal1_imp_lib_no_mort))

Ordinal1_imp_mort_no_lib <- na.omit(Ordinal1_imp_mort_no_lib) %>% 
  filter(V1 != "Error") 
Ordinal1_imp_mort_no_lib <- Ordinal1_imp_mort_no_lib %>% 
  mutate(X = seq(1, nrow(Ordinal1_imp_mort_no_lib), 1)) %>%
  rename(post_p = V1,
         trt_effect = V2,
         sd = V3)
(Ordinal1_post_imp_mort_no_lib <- nrow(Ordinal1_imp_mort_no_lib[Ordinal1_imp_mort_no_lib$post_p > lambda_Ordinal1,])/nrow(Ordinal1_imp_mort_no_lib))


Ordinal1_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "worsen_lib_imp_mort", "imp_lib_no_mort", "imp_mort_no_lib"),
                         post_p = c(Ordinal1_post_p_noeffect, Ordinal1_post_p_worsen_both, Ordinal1_post_p_imp_both, Ordinal1_post_imp_lib_worse_mort, Ordinal1_post_imp_mort_worse_lib, Ordinal1_post_imp_lib_no_mort, Ordinal1_post_imp_mort_no_lib),
                         trt_est_mean = c(mean(as.numeric(Ordinal1_noeffect$trt_effect)), mean(as.numeric(Ordinal1_worsen_both$trt_effect)), mean(as.numeric(Ordinal1_imp_both$trt_effect)), mean(as.numeric(Ordinal1_imp_lib_worse_mort$trt_effect)), mean(as.numeric(Ordinal1_imp_mort_worse_lib$trt_effect)), mean(as.numeric(Ordinal1_imp_lib_worse_mort$trt_effect)), mean(as.numeric(Ordinal1_imp_mort_worse_lib$trt_effect))),
                         sd_est_mean = c(mean(as.numeric(Ordinal1_noeffect$sd)), mean(as.numeric(Ordinal1_worsen_both$sd)), mean(as.numeric(Ordinal1_imp_both$sd)), mean(as.numeric(Ordinal1_imp_lib_worse_mort$sd)), mean(as.numeric(Ordinal1_imp_mort_worse_lib$sd)), mean(as.numeric(Ordinal1_imp_lib_worse_mort$sd)), mean(as.numeric(Ordinal1_imp_mort_worse_lib$sd))))

sqrt(Ordinal1_1$post_p*(1 - Ordinal1_1$post_p)/nrow(Ordinal1_imp_lib_no_mort))[c(1,3,2,5,4,7,6)]
