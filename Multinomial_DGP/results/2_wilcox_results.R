library(tidyverse)
library(dplyr)

setwd("/GitHub/AnalyzingVDFs/Multinomial_DGP")
wdname <- paste(getwd(),"/2 - wilcox", sep = "")
setwd(wdname)

wilcox_noeffect <- read.csv("wilcox_noeffect.csv")
wilcox_imp_both <- read.csv("wilcox_imp_both.csv")
wilcox_imp_lib_worse_mort <- read.csv("wilcox_imp_lib_worse_mort.csv")
wilcox_imp_mort_worse_lib <- read.csv("wilcox_imp_mort_worse_lib.csv")
wilcox_worsen_both <- read.csv("wilcox_worsen_both.csv")
wilcox_imp_lib_no_mort <- read.csv("wilcox_imp_lib_no_mort.csv")
wilcox_imp_mort_no_lib <- read.csv("wilcox_imp_mort_no_lib.csv")

#-------calculating lambda that controls type one error----------

wilcox_noeffect <- na.omit(wilcox_noeffect) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2)

#type 1 error

(wilcox_type1_error <- nrow(wilcox_noeffect[wilcox_noeffect$positive_difference == "Yes",])/nrow(wilcox_noeffect))


#--------based on the lambda, calculating "power" for each scenario----------
wilcox_worsen_both <- na.omit(wilcox_worsen_both) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2
  )
(wilcox_power_worsen_both <- nrow(wilcox_worsen_both[wilcox_worsen_both$positive_difference == "Yes",])/nrow(wilcox_worsen_both))

wilcox_imp_both <- na.omit(wilcox_imp_both) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2
  )
(wilcox_power_imp_both <- nrow(wilcox_imp_both[wilcox_imp_both$positive_difference == "Yes",])/nrow(wilcox_imp_both))

wilcox_imp_lib_worse_mort <- na.omit(wilcox_imp_lib_worse_mort) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2
  )
(wilcox_post_imp_lib_worse_mort <- nrow(wilcox_imp_lib_worse_mort[wilcox_imp_lib_worse_mort$positive_difference == "Yes",])/nrow(wilcox_imp_lib_worse_mort))

wilcox_imp_mort_worse_lib <- na.omit(wilcox_imp_mort_worse_lib) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2
  )
(wilcox_post_imp_mort_worse_lib <- nrow(wilcox_imp_mort_worse_lib[wilcox_imp_mort_worse_lib$positive_difference == "Yes",])/nrow(wilcox_imp_mort_worse_lib))

wilcox_imp_lib_no_mort <- na.omit(wilcox_imp_lib_no_mort) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2
  )
(wilcox_post_imp_lib_no_mort <- nrow(wilcox_imp_lib_no_mort[wilcox_imp_lib_no_mort$positive_difference == "Yes",])/nrow(wilcox_imp_lib_no_mort))

wilcox_imp_mort_no_lib <- na.omit(wilcox_imp_mort_no_lib) %>% 
  mutate(X = seq(1, 2000, 1)) %>%
  rename(positive_difference = V1,
         p_val = V2
  )
(wilcox_post_imp_mort_no_lib <- nrow(wilcox_imp_mort_no_lib[wilcox_imp_mort_no_lib$positive_difference == "Yes",])/nrow(wilcox_imp_mort_no_lib))


wilcox_1 <- data.frame(scenario = c("no_eff", "worsen_both", "imp_both", "worsen_mort_imp_lib", "worsen_lib_imp_mort", "imp_lib_no_mort", "imp_mort_no_lib"),
                       power = c(wilcox_type1_error, wilcox_power_worsen_both, wilcox_power_imp_both, wilcox_post_imp_lib_worse_mort, wilcox_post_imp_mort_worse_lib, wilcox_post_imp_lib_no_mort, wilcox_post_imp_mort_no_lib))
sqrt(wilcox_1$power*(1 - wilcox_1$power)/nrow(wilcox_imp_mort_no_lib))[c(1,3,2,5,4,7,6)]
