library(boot)
library(dplyr)
library(copula)

multinomial_data_gen_diff_effect <- function(OR_mort, OR_lib){
  VFDs <- c(seq(-1, 28, 1))[-30]
  N <- c(115, 14, 2, 2, 7, 4, 2, 3, 2, 1, 2, 2, 3, 5, 3, 5, 5, 6, 6, 6, 5, 10, 11, 10, 8, 13, 25, 21, 17, 15)[-30]
  
  #control group baseline proportion in each category
  ctl_prop <- N/sum(N)
  ctl_prop_cumsum <- cumsum(ctl_prop)
  ctl_odds <- ctl_prop_cumsum/(1-ctl_prop_cumsum)
  
  #treatment no effect
  trt_odds <- OR_lib*ctl_odds
  trt_prop_cumsum <- c((trt_odds/(1+trt_odds))[-29], 1)
  trt_prop <- c(trt_prop_cumsum[1],diff(trt_prop_cumsum))
  trt_prop <- trt_prop[-1] / sum(trt_prop[-1])
  
  ##############################################
  #### trt imp lib worse mort       ############(follow improving both outcomes)
  ##############################################
  
  trt_diff_effect_odds1 <- OR_mort * ctl_odds[1]
  trt_diff_effect_prop1 <- trt_diff_effect_odds1/(1+trt_diff_effect_odds1)
  trt_diff_effect_prop <- c(trt_diff_effect_prop1, trt_prop * (1 - trt_diff_effect_prop1))
  
  return(trt_diff_effect_prop)
}

multinomial_data_gen_same_effect <- function(OR){
  VFDs <- c(seq(-1, 28, 1))[-30]
  N <- c(115, 14, 2, 2, 7, 4, 2, 3, 2, 1, 2, 2, 3, 5, 3, 5, 5, 6,
         6, 6, 5, 10, 11, 10, 8, 13, 25, 21, 17, 15)[-30]
  
  #control group baseline proportion in each category
  ctl_prop <- N/sum(N)
  ctl_prop_cumsum <- cumsum(ctl_prop)
  ctl_odds <- ctl_prop_cumsum/(1-ctl_prop_cumsum)
  
  #treatment no effect
  #OR <- 1
  trt_odds <- OR*ctl_odds
  trt_prop_cumsum <- c((trt_odds/(1+trt_odds))[-29], 1)
  trt_prop <- c(trt_prop_cumsum[1],diff(trt_prop_cumsum))

  return(trt_prop)
}

## NULL
null_d_multi <- multinomial_data_gen_same_effect(1)[1]
null_l_multi <- multinomial_data_gen_same_effect(1)[-1]%*%28:1 / 
  sum(multinomial_data_gen_same_effect(1)[-1])
null_l_sd_multi <- sqrt((multinomial_data_gen_same_effect(1)[-1]%*%(28:1)^2 / 
  sum(multinomial_data_gen_same_effect(1)[-1])) - null_l_multi^2)
## IMPROVE BOTH
imp_d_multi <- multinomial_data_gen_same_effect(0.75)[1]
imp_l_multi <- multinomial_data_gen_same_effect(0.75)[-1]%*%28:1 / 
  sum(multinomial_data_gen_same_effect(0.75)[-1])
imp_l_sd_multi <- sqrt((multinomial_data_gen_same_effect(0.75)[-1]%*%(28:1)^2 / 
                           sum(multinomial_data_gen_same_effect(0.75)[-1])) - 
                         imp_l_multi^2)
## WORSEN BOTH
wor_d_multi <- multinomial_data_gen_same_effect(1.25)[1]
wor_l_multi <- multinomial_data_gen_same_effect(1.25)[-1]%*%28:1 / 
  sum(multinomial_data_gen_same_effect(1.25)[-1])
wor_l_sd_multi <- sqrt((multinomial_data_gen_same_effect(1.25)[-1]%*%(28:1)^2 / 
                          sum(multinomial_data_gen_same_effect(1.25)[-1])) - 
                         wor_l_multi^2)
## IMPROVE MORT; WORSEN LIB
imp_d_wor_l_d_multi <- multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1.25)[1]
imp_d_wor_l_l_multi <- multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1.25)[-1]%*%28:1 / 
  sum(multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1.25)[-1])
imp_d_wor_l_l_multi_sd <- sqrt((multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1.25)[-1]%*%(28:1)^2 / 
                          sum(multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1.25)[-1])) - 
                            imp_d_wor_l_l_multi^2)
## IMPROVE LIB; WORSEN MORT
imp_l_wor_d_d_multi <- multinomial_data_gen_diff_effect(OR_mort = 1.25, OR_lib = 0.75)[1]
imp_l_wor_d_l_multi <- multinomial_data_gen_diff_effect(OR_mort = 1.25, OR_lib = 0.75)[-1]%*%28:1 / 
  sum(multinomial_data_gen_diff_effect(OR_mort = 1.25, OR_lib = 0.75)[-1])
imp_l_wor_d_l_multi_sd <- sqrt((multinomial_data_gen_diff_effect(OR_mort = 1.25, OR_lib = 0.75)[-1]%*%(28:1)^2 / 
                                  sum(multinomial_data_gen_diff_effect(OR_mort = 1.25, OR_lib = 0.75)[-1])) - 
                                 imp_l_wor_d_l_multi^2)

## IMPROVE MORT; NO LIB
imp_d_no_l_d_multi <- multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1)[1]
imp_d_no_l_l_multi <- multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1)[-1]%*%28:1 / 
  sum(multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1)[-1])
imp_d_no_l_l_multi_sd <- sqrt((multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1)[-1]%*%(28:1)^2 / 
                                  sum(multinomial_data_gen_diff_effect(OR_mort = 0.75, OR_lib = 1)[-1])) - 
                                imp_d_no_l_l_multi^2)

## IMPROVE LIB; NO MORT
imp_l_no_d_d_multi <- multinomial_data_gen_diff_effect(OR_mort = 1, OR_lib = 0.75)[1]
imp_l_no_d_l_multi <- multinomial_data_gen_diff_effect(OR_mort = 1, OR_lib = 0.75)[-1]%*%28:1 / 
  sum(multinomial_data_gen_diff_effect(OR_mort = 1, OR_lib = 0.75)[-1])
imp_l_no_d_l_multi_sd <- sqrt((multinomial_data_gen_diff_effect(OR_mort = 1, OR_lib = 0.75)[-1]%*%(28:1)^2 / 
                                  sum(multinomial_data_gen_diff_effect(OR_mort = 1, OR_lib = 0.75)[-1])) - 
                                imp_l_no_d_l_multi^2)

cbind(
  c(null_d_multi, imp_d_multi, wor_d_multi, imp_d_wor_l_d_multi, 
    imp_l_wor_d_d_multi, imp_d_no_l_d_multi, imp_l_no_d_d_multi),
  c(null_l_multi, imp_l_multi, wor_l_multi, imp_d_wor_l_l_multi, 
    imp_l_wor_d_l_multi, imp_d_no_l_l_multi, imp_l_no_d_l_multi),
  c(null_l_sd_multi, imp_l_sd_multi, wor_l_sd_multi, imp_d_wor_l_l_multi_sd, 
    imp_l_wor_d_l_multi_sd, imp_d_no_l_l_multi_sd, imp_l_no_d_l_multi_sd))


data_gen_copula <- function(n_sample,        #trial sample size
                            F_90_d_trt,    #cumulative hazard at day 90
                            F_90_d_ctl,    ##cumulative hazard at day 90
                            F_28_v_trt,    #cumulative liberation rate from ventilator at day 28
                            F_28_v_ctl,    #cumulative liberation rate from ventilator at day 28
                            alpha_d,         #weibull distribution shape parameter for death
                            alpha_v         #weibull distribution shape parameter for liberation
){
  #generating patients -- making sure that treatment and control group are balanced
  treatment <- c("trt", "ctl")
  trt <- sample(treatment, n_sample, prob = c(0.5, 0.5), replace = TRUE)
  patient_id <- seq(1: n_sample)                                                    
  data <- data.frame(id = vector("integer", n_sample),
                     treatment = vector("character", n_sample))
  data <- data %>%
    mutate(id = patient_id,
           treatment = trt) 
  #calculating parameter values based on pre-set 90 day mortality
  #assuming the shape of the distribution is the same for trt and ctl group
  (beta_d_trt <- 90/((-log(1 - F_90_d_trt))^(1/alpha_d)))
  (beta_d_ctl <- 90/((-log(1 - F_90_d_ctl))^(1/alpha_d)))
  
  #calculating parameter values based on pre-set 28 day liberation
  #assuming the shape of the distribution is the same for trt and ctl group
  #assuming the increase in hazard rate for extubation is faster than mortality
  (beta_v_trt <- 28/((-log(1 - F_28_v_trt))^(1/alpha_v)))
  (beta_v_ctl <- 28/((-log(1 - F_28_v_ctl))^(1/alpha_v)))
  
  #generating correlated variables using copula
  cop_trt <- normalCopula(-0.4, dim = 2, dispstr = "ex")
  mvd_trt <- mvdc(copula = cop_trt, 
                  margins = c("weibull", "weibull"), 
                  paramMargins = list(list(shape = alpha_d, scale = beta_d_trt),
                                      list(shape = alpha_v, scale = beta_v_trt)))
  
  cop_ctl <- normalCopula(-0.4, dim = 2, dispstr = "ex")
  mvd_ctl <- mvdc(copula = cop_ctl, 
                  margins = c("weibull", "weibull"), 
                  paramMargins = list(list(shape = alpha_d, scale = beta_d_ctl),
                                      list(shape = alpha_v, scale = beta_v_ctl)))
  #initiating binary variables, time to event variables and parameter values for death and extubation
  death_bin <- NULL
  VF_bin <- NULL
  time_to_death <- NULL
  time_to_VF <- NULL
  
  
  for(i in 1: n_sample){
    if(data$treatment[i] == "trt"){
      var_trt <- rMvdc(1, mvd_trt)
      time_to_death[i] <- ceiling(var_trt[1])
      time_to_VF[i] <- ceiling(var_trt[2])
    }
    else if(data$treatment[i] == "ctl"){
      var_ctl <- rMvdc(1, mvd_ctl)
      time_to_death[i] <- ceiling(var_ctl[1])
      time_to_VF[i] <- ceiling(var_ctl[2])
    }
  }
  
  data <- data %>%
    mutate(time_to_death = time_to_death,
           time_to_VF = time_to_VF)
  return(data)
}


### NULL
dat_to_try <- data_gen_copula(50000, 0.37, 0.37, 0.759, 0.89, 0.9, 1.1)

null_d_wie <- mean(dat_to_try$time_to_death <= 90)
dat_to_try$censor_VF <- 1
dat_to_try$censor_VF[dat_to_try$time_to_VF >= 28] <- 0
dat_to_try$time_to_VF[dat_to_try$time_to_VF >= 28] <- 28
null_l_wie <- mean(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])
null_l_wie_sd <- sd(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])

##PLOTTING
dat_to_try

### IMPROVE BOTH
dat_to_try <- data_gen_copula(50000, 0.32, 0.32, 0.93, 0.93, 0.9, 1.1)

imp_d_wie <- mean(dat_to_try$time_to_death <= 90)
dat_to_try$censor_VF <- 1
dat_to_try$censor_VF[dat_to_try$time_to_VF >= 28] <- 0
dat_to_try$time_to_VF[dat_to_try$time_to_VF >= 28] <- 28
imp_l_wie <- mean(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])
imp_l_wie_sd <- sd(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])

## WORSEN BOTH
dat_to_try <- data_gen_copula(50000, 0.42, 0.42, 0.86, 0.86, 0.9, 1.1)

wor_d_wie <- mean(dat_to_try$time_to_death <= 90)
dat_to_try$censor_VF <- 1
dat_to_try$censor_VF[dat_to_try$time_to_VF >= 28] <- 0
dat_to_try$time_to_VF[dat_to_try$time_to_VF >= 28] <- 28
wor_l_wie <- mean(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])
wor_l_wie_sd <- sd(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])

### IMPROVE MORT; WORSEN LIB
dat_to_try <- data_gen_copula(50000, 0.32, 0.32, 0.86, 0.86, 0.9, 1.1)

imp_d_wor_l_d_wie <- mean(dat_to_try$time_to_death <= 90)
dat_to_try$censor_VF <- 1
dat_to_try$censor_VF[dat_to_try$time_to_VF >= 28] <- 0
dat_to_try$time_to_VF[dat_to_try$time_to_VF >= 28] <- 28
imp_d_wor_l_l_wie <- mean(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])
imp_d_wor_l_l_wie_sd <- sd(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])

### IMPROVE LIB; WORSEN MORT
dat_to_try <- data_gen_copula(50000, 0.42, 0.42, 0.93, 0.93, 0.9, 1.1)

imp_l_wor_d_d_wie <- mean(dat_to_try$time_to_death <= 90)
dat_to_try$censor_VF <- 1
dat_to_try$censor_VF[dat_to_try$time_to_VF >= 28] <- 0
dat_to_try$time_to_VF[dat_to_try$time_to_VF >= 28] <- 28
imp_l_wor_d_l_wie <- mean(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])
imp_l_wor_d_l_wie_sd <- sd(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])

### IMPROVE MORT; NO LIB
dat_to_try <- data_gen_copula(50000, 0.32, 0.32, 0.89, 0.89, 0.9, 1.1)

imp_d_no_l_d_wie <- mean(dat_to_try$time_to_death <= 90)
dat_to_try$censor_VF <- 1
dat_to_try$censor_VF[dat_to_try$time_to_VF >= 28] <- 0
dat_to_try$time_to_VF[dat_to_try$time_to_VF >= 28] <- 28
imp_d_no_l_l_wie <- mean(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])
imp_d_no_l_l_wie_sd <- sd(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])

### IMPROVE LIB; WORSEN MORT
dat_to_try <- data_gen_copula(50000, 0.37, 0.37, 0.93, 0.93, 0.9, 1.1)

imp_l_no_d_d_wie <- mean(dat_to_try$time_to_death <= 90)
dat_to_try$censor_VF <- 1
dat_to_try$censor_VF[dat_to_try$time_to_VF >= 28] <- 0
dat_to_try$time_to_VF[dat_to_try$time_to_VF >= 28] <- 28
imp_l_no_d_l_wie <- mean(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])
imp_l_no_d_l_wie_sd <- sd(dat_to_try[dat_to_try$time_to_death > 90,"time_to_VF"])

round(cbind(
  c(null_d_wie, imp_d_wie, wor_d_wie, imp_d_wor_l_d_wie, 
    imp_l_wor_d_d_wie, imp_d_no_l_d_wie, imp_l_no_d_d_wie),
  c(null_l_wie, imp_l_wie, wor_l_wie, imp_d_wor_l_l_wie, 
    imp_l_wor_d_l_wie, imp_d_no_l_l_wie, imp_l_no_d_l_wie),
  c(null_l_wie_sd, imp_l_wie_sd, wor_l_wie_sd, imp_d_wor_l_l_wie_sd, 
    imp_l_wor_d_l_wie_sd, imp_d_no_l_l_wie_sd, imp_l_no_d_l_wie_sd)
),2)
