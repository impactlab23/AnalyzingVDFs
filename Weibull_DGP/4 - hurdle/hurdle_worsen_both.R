################################################################################################
##########                    Simulations for PRACTICAL                   ######################
########                     Zero-inflated Poisson hurdle                            ###########
################################################################################################
library(boot)
library(dplyr)
library(doParallel)
library(foreach)
library(tidyverse)
library(brms)
library(copula)

#####################################################
#####    Data Generating Function      ##############
#####################################################
#--------------- simulating correlated data ------------------------
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
  
  
  for(i in 1: 1500){
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

#----------------binary + survival-------------------------
data_gen_weibull_binsurv <- function(n_sample,        #trial sample size
                                     F_90_d_trt,    #cumulative hazard at day 90
                                     F_90_d_ctl,    ##cumulative hazard at day 90
                                     F_28_v_trt,    #cumulative liberation rate from ventilator at day 28
                                     F_28_v_ctl,    #cumulative liberation rate from ventilator at day 28
                                     alpha_d,         #weibull distribution shape parameter for death
                                     alpha_v         #weibull distribution shape parameter for liberation
){
  data <- data_gen_copula(n_sample, F_90_d_trt, F_90_d_ctl, F_28_v_trt, F_28_v_ctl, alpha_d, alpha_v)
  
  data_CR_E <- data.frame(id = vector("integer", n_sample),
                          treatment = vector("character", n_sample),
                          censor = vector("integer", n_sample),
                          ftime = vector("integer", n_sample),
                          time_to_death = vector("integer", n_sample),
                          time_to_VF = vector("integer", n_sample),
                          death_status = vector("integer", n_sample))
  #ventilator status is not monitored after day 28
  ftime_CR_E <- NULL
  trt_CR_E <- NULL
  censoring_CR_E <- NULL
  ID_CR_E <- NULL
  time_to_death_E <- NULL
  time_to_VF_E <- NULL
  death_status_E <- NULL
  
  for(i in 1: n_sample){
    ID_CR_E[i] <- data$id[i]
    trt_CR_E[i] <- data$treatment[i]
    time_to_death_E[i] <- data$time_to_death[i]
    time_to_VF_E[i] <- data$time_to_VF[i]
    
    if(data$time_to_VF[i] <= 28){
      if(data$time_to_death[i] <= data$time_to_VF[i]){
        ftime_CR_E[i] <- data$time_to_death[i]
        censoring_CR_E[i] <- 1
        death_status_E[i] <- 1
      }
      else if(data$time_to_death[i] > data$time_to_VF[i] & data$time_to_death[i] < 90){
        ftime_CR_E[i] <- data$time_to_VF[i]
        censoring_CR_E[i] <- 0
        death_status_E[i] <- 1
      }
      else if(data$time_to_death[i] > data$time_to_VF[i] & data$time_to_death[i] >= 90){
        ftime_CR_E[i] <- data$time_to_VF[i]
        censoring_CR_E[i] <- 0
        death_status_E[i] <- 0
      }
    }
    else if(data$time_to_VF[i] > 28){
      if(data$time_to_death[i] <= 28){
        ftime_CR_E[i] <- data$time_to_death[i]
        censoring_CR_E[i] <- 1
        death_status_E[i] <- 1
      }
      else if(data$time_to_death[i] > 28 & data$time_to_death[i] < 90){
        ftime_CR_E[i] <- 28
        censoring_CR_E[i] <- 1
        death_status_E[i] <- 1
      }
      else if(data$time_to_death[i] >= 90){
        ftime_CR_E[i] <- 28
        censoring_CR_E[i] <- 1
        death_status_E[i] <- 0
      }
    }
  }
  #this can be used for the survival analysis later
  data_CR_E <- data_CR_E %>%
    mutate(id = ID_CR_E,
           treatment = trt_CR_E,
           censor = censoring_CR_E,
           ftime = ftime_CR_E,
           time_to_death = time_to_death_E,
           time_to_VF = time_to_VF_E,
           death_status = death_status_E)
  
  #using r-inla package, the censoring status should be 1 for event and 0 for censored
  data_binsurv <- data_CR_E %>%
    mutate(censor = ifelse(censor == 1, 0, 1))
  return(data_binsurv = data_binsurv)
}

#----------------Ordinal Outcome 2-------------------------
data_gen_weibull_ordinal_2 <- function(n_sample,        #trial sample size
                                       F_90_d_trt,    #cumulative hazard at day 90
                                       F_90_d_ctl,    ##cumulative hazard at day 90
                                       F_28_v_trt,    #cumulative liberation rate from ventilator at day 28
                                       F_28_v_ctl,    #cumulative liberation rate from ventilator at day 28
                                       alpha_d,         #weibull distribution shape parameter for death
                                       alpha_v         #weibull distribution shape parameter for liberation
){
  data_binsurv <- data_gen_weibull_binsurv(n_sample, F_90_d_trt, F_90_d_ctl, F_28_v_trt, F_28_v_ctl, alpha_d, alpha_v)
  
  #------------------VFDs based on SAP definition--------------------
  VFD_OO_2 <- NULL
  for(i in 1: n_sample){
    if(data_binsurv$death_status[i] == 1){
      VFD_OO_2[i] <- -1 + 2
    }
    else if(data_binsurv$death_status[i] == 0 & data_binsurv$time_to_VF[i] >= 28){
      VFD_OO_2[i] <- 0 + 2
    }
    else if(data_binsurv$death_status[i] == 0 & data_binsurv$time_to_VF[i] < 28){
      VFD_OO_2[i] <- 28 - data_binsurv$time_to_VF[i] + 2
    }
  }
  data_ordinal_2 <- data_binsurv %>% mutate(VFD_OO_2 = VFD_OO_2)
  return(data_ordinal_2 = data_ordinal_2)
}



PoissonHurdle_single_sim <- function(n_sample,        #trial sample size
                                     F_90_d_trt,    #cumulative hazard at day 90
                                     F_90_d_ctl,    ##cumulative hazard at day 90
                                     F_28_v_trt,    #cumulative liberation rate from ventilator at day 28
                                     F_28_v_ctl,    #cumulative liberation rate from ventilator at day 28
                                     alpha_d,         #weibull distribution shape parameter for death
                                     alpha_v         #weibull distribution shape parameter for liberation
){
  data <- data_gen_weibull_ordinal_2(n_sample, F_90_d_trt, F_90_d_ctl, F_28_v_trt, F_28_v_ctl, alpha_d, alpha_v)
  data <- data %>% 
    mutate(VFD_OO_2 = VFD_OO_2 - 1)
  #####################################################
  #####    Analysis               #####################
  #####################################################
  tryCatch({
    #----------------------------------------------------------------------
    #-------------  logit-poisson hurdle model   --------------------------
    #----------------------------------------------------------------------
    
    #-------------third definition of ordinal outcome using INLA------------
    hurdle <- brm(bf(VFD_OO_2 ~ treatment, 
                     hu ~ treatment), 
                  data = data,
                  set_prior("normal(0, 10)", class = "b", coef = "treatmenttrt"),
                  family = "hurdle_poisson")
    a <- summary(hurdle)
    ###################################################
    ###        posterior probability           ########
    ###################################################
    draw_stan_hurdle <- as.data.frame(hurdle)
    post_prob_hurdle_ber <- length(draw_stan_hurdle$b_hu_treatmenttrt[draw_stan_hurdle$b_hu_treatmenttrt < 0])/nrow(draw_stan_hurdle)
    post_prob_hurdle_poi <- length(draw_stan_hurdle$b_treatmenttrt[draw_stan_hurdle$b_treatmenttrt > 0])/nrow(draw_stan_hurdle)
    ###################################################
    ###        treatment effect estimates      ########
    ###################################################
    trt_hurdle_ber  <- summary(hurdle)$fixed[4,1]
    trt_hurdle_poi  <- summary(hurdle)$fixed[3,1]
    ###################################################
    ###        sd estimates                    ########
    ###################################################
    sd_hurdle_ber  <- summary(hurdle)$fixed[4,2]
    sd_hurdle_poi  <- summary(hurdle)$fixed[3,2]
    
    return(c(post_prob_hurdle_ber,
             post_prob_hurdle_poi,
             trt_hurdle_ber,
             trt_hurdle_poi,
             sd_hurdle_ber,
             sd_hurdle_poi
    ))
  },
  error = function(e){
    print("Error in INLA")
    post_prob_hurdle_ber <- "Error"
    post_prob_hurdle_poi <- "Error"
    ###################################################
    ###        treatment effect estimates      ########
    ###################################################
    trt_hurdle_ber  <- "Error"
    trt_hurdle_poi  <- "Error"
    ###################################################
    ###        sd estimates                    ########
    ###################################################
    sd_hurdle_ber  <- "Error"
    sd_hurdle_poi  <- "Error"
    
    return(c(post_prob_hurdle_ber,
             post_prob_hurdle_poi,
             trt_hurdle_ber,
             trt_hurdle_poi,
             sd_hurdle_ber,
             sd_hurdle_poi
    ))
  }
  )
}


#running the simulation for 1000 times to calculate the lambda to control for type one error 
n_sim <- 2000
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)
Sim_Res <- foreach(k = 1:n_sim, .combine = rbind,
                   .packages = c("boot", "dplyr", "tidyverse", "brms", "copula")
) %dopar% { 
  set.seed(k + 5125)
  PoissonHurdle_single_sim(n_sample = 1500,        #trial sample size
                           F_90_d_trt = 0.42,    #cumulative hazard at day 90
                           F_90_d_ctl = 0.37,    ##cumulative hazard at day 90
                           F_28_v_trt = 0.86,    #cumulative liberation rate from ventilator at day 28
                           F_28_v_ctl = 0.89,    #cumulative liberation rate from ventilator at day 28
                           alpha_d = 0.9,         #weibull distribution shape parameter for death
                           alpha_v = 1.1         #weibull distribution shape parameter for liberation
  )     
}
stopCluster(cl)

results <- as.data.frame(Sim_Res)

#write.csv(results, file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.csv")

#save.image(file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.Rdata")

write.csv(results, file = "C:/Users/Anna Heath/OneDrive - SickKids/Admin Anna/Hiring and Staff/Ziming Chen 2022/PRACTICAL/SimulationMar20/Weibull_AH/hurdle/hurdle_worsen_both.csv")


