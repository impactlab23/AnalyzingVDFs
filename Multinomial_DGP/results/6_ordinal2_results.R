library(tidyverse)
library(dplyr)


setwd("/GitHub/AnalyzingVDFs/Multinomial_DGP")
wdname <- paste(getwd(),"/6 - ordinal2", sep = "")
setwd(wdname)

Ordinal3_noeffect <- read.csv("Ordinal3_noeffect.csv")
Ordinal3_imp_both <- read.csv("Ordinal3_imp_both.csv")
Ordinal3_imp_lib_worse_mort <- read.csv("Ordinal3_imp_lib_worse_mort.csv")
Ordinal3_imp_mort_worse_lib <- read.csv("Ordinal3_imp_mort_worse_lib.csv")
Ordinal3_worsen_both <- read.csv("Ordinal3_worsen_both.csv")
Ordinal3_imp_lib_no_mort <- read.csv("ordinal3_imp_lib_no_mort.csv")
Ordinal3_imp_mort_no_lib <- read.csv("ordinal3_imp_mort_no_lib.csv")

#-------calculating lambda that controls type one error----------
#hist(Ordinal3_imp_both$post_p)
#var(Ordinal3_imp_both$post_p)
Ordinal3_noeffect <- na.omit(Ordinal3_noeffect) %>% 
  mutate(X = seq(1, 2000, 1)) %>% filter(V1 != "Error") %>%
  mutate(post_p = as.numeric(V1),
         trt_effect = as.numeric(V2),
         sd = as.numeric(V3))

Ordinal3_noeffect <- Ordinal3_noeffect %>%
  arrange(desc(post_p))

#0.9468489 
(lambda_Ordinal3 <- quantile(Ordinal3_noeffect$post_p, probs = c(0.95)))

lambda_Ordinal3 <- 0.965
(Ordinal3_post_p_noeffect <- nrow(Ordinal3_noeffect[Ordinal3_noeffect$post_p >= lambda_Ordinal3,])/nrow(Ordinal3_noeffect))


#--------based on the lambda, calculating "power" for each scenario----------
Ordinal3_worsen_both <- na.omit(Ordinal3_worsen_both) %>% 
  mutate(X = seq(1, 2000, 1)) %>% filter(V1 != "Error") %>%
  mutate(post_p = as.numeric(V1),
         trt_effect = as.numeric(V2),
         sd = as.numeric(V3))
(Ordinal3_post_p_worsen_both <- nrow(Ordinal3_worsen_both[Ordinal3_worsen_both$post_p > lambda_Ordinal3,])/nrow(Ordinal3_worsen_both))

Ordinal3_imp_both <- na.omit(Ordinal3_imp_both) %>% 
  mutate(X = seq(1, 2000, 1)) %>% filter(V1 != "Error") %>%
  mutate(post_p = as.numeric(V1),
         trt_effect = as.numeric(V2),
         sd = as.numeric(V3))
(Ordinal3_post_p_imp_both <- nrow(Ordinal3_imp_both[Ordinal3_imp_both$post_p > lambda_Ordinal3,])/nrow(Ordinal3_imp_both))

Ordinal3_imp_lib_worse_mort <- na.omit(Ordinal3_imp_lib_worse_mort) %>% 
  mutate(X = seq(1, 2000, 1)) %>% filter(V1 != "Error") %>%
  mutate(post_p = as.numeric(V1),
         trt_effect = as.numeric(V2),
         sd = as.numeric(V3))
(Ordinal3_post_imp_lib_worse_mort <- nrow(Ordinal3_imp_lib_worse_mort[Ordinal3_imp_lib_worse_mort$post_p > lambda_Ordinal3,])/nrow(Ordinal3_imp_lib_worse_mort))

Ordinal3_imp_mort_worse_lib <- na.omit(Ordinal3_imp_mort_worse_lib) %>% 
  mutate(X = seq(1, 2000, 1)) %>% filter(V1 != "Error") %>%
  mutate(post_p = as.numeric(V1),
         trt_effect = as.numeric(V2),
         sd = as.numeric(V3))
(Ordinal3_post_imp_mort_worse_lib <- nrow(Ordinal3_imp_mort_worse_lib[Ordinal3_imp_mort_worse_lib$post_p > lambda_Ordinal3,])/nrow(Ordinal3_imp_mort_worse_lib))

Ordinal3_imp_lib_no_mort <- na.omit(Ordinal3_imp_lib_no_mort) %>% 
  filter(V1 != "Error") 
Ordinal3_imp_lib_no_mort <- Ordinal3_imp_lib_no_mort %>% 
  mutate(X = seq(1, nrow(Ordinal3_imp_lib_no_mort), 1)) %>%
  rename(post_p = V1,
         trt_effect = V2,
         sd = V3)
(Ordinal3_post_imp_lib_no_mort <- nrow(Ordinal3_imp_lib_no_mort[Ordinal3_imp_lib_no_mort$post_p > lambda_Ordinal3,])/nrow(Ordinal3_imp_lib_no_mort))

Ordinal3_imp_mort_no_lib <- na.omit(Ordinal3_imp_mort_no_lib) %>% 
  filter(V1 != "Error") 
Ordinal3_imp_mort_no_lib <- Ordinal3_imp_mort_no_lib %>% 
  mutate(X = seq(1, nrow(Ordinal3_imp_mort_no_lib), 1)) %>%
  rename(post_p = V1,
         trt_effect = V2,
         sd = V3)
(Ordinal3_post_imp_mort_no_lib <- nrow(Ordinal3_imp_mort_no_lib[Ordinal3_imp_mort_no_lib$post_p > lambda_Ordinal3,])/nrow(Ordinal3_imp_mort_no_lib))


Ordinal3_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "worsen_lib_imp_mort", "imp_lib_no_mort", "imp_mort_no_lib"),
                         post_p = c(Ordinal3_post_p_noeffect, Ordinal3_post_p_worsen_both, Ordinal3_post_p_imp_both, Ordinal3_post_imp_lib_worse_mort, Ordinal3_post_imp_mort_worse_lib, Ordinal3_post_imp_lib_no_mort, Ordinal3_post_imp_mort_no_lib),
                         trt_est_mean = c(mean(as.numeric(Ordinal3_noeffect$trt_effect)), mean(as.numeric(Ordinal3_worsen_both$trt_effect)), mean(as.numeric(Ordinal3_imp_both$trt_effect)), mean(as.numeric(Ordinal3_imp_lib_worse_mort$trt_effect)), mean(as.numeric(Ordinal3_imp_mort_worse_lib$trt_effect)), mean(as.numeric(Ordinal3_imp_lib_worse_mort$trt_effect)), mean(as.numeric(Ordinal3_imp_mort_worse_lib$trt_effect))),
                         sd_est_mean = c(mean(as.numeric(Ordinal3_noeffect$sd)), mean(as.numeric(Ordinal3_worsen_both$sd)), mean(as.numeric(Ordinal3_imp_both$sd)), mean(as.numeric(Ordinal3_imp_lib_worse_mort$sd)), mean(as.numeric(Ordinal3_imp_mort_worse_lib$sd)), mean(as.numeric(Ordinal3_imp_lib_worse_mort$sd)), mean(as.numeric(Ordinal3_imp_mort_worse_lib$sd))))

sqrt(Ordinal3_1$post_p*(1 - Ordinal3_1$post_p)/nrow(Ordinal3_imp_mort_no_lib))[c(1,3,2,5,4,7,6)]
