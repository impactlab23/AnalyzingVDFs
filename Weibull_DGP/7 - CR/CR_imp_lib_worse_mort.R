################################################################################################
##########                    Simulations for PRACTICAL                   ######################
########                        Competing Risk  no effect                            ###########
################################################################################################
library(boot)
library(dplyr)
library(INLA)
library(doParallel)
library(foreach)
library(tidyverse)
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

#------------competing risks-----------------------
data_gen_weibull_CR <- function(n_sample,        #trial sample size
                                F_90_d_trt,    #cumulative hazard at day 90
                                F_90_d_ctl,    ##cumulative hazard at day 90
                                F_28_v_trt,    #cumulative liberation rate from ventilator at day 28
                                F_28_v_ctl,    #cumulative liberation rate from ventilator at day 28
                                alpha_d,         #weibull distribution shape parameter for death
                                alpha_v         #weibull distribution shape parameter for liberation
){
  data <- data_gen_copula(n_sample, F_90_d_trt, F_90_d_ctl, F_28_v_trt, F_28_v_ctl, alpha_d, alpha_v)
  data_CR_D <- data.frame(id = vector("integer", n_sample),
                          treatment = vector("character", n_sample),
                          censor = vector("integer", n_sample),
                          event = vector("character", n_sample),
                          ftime = vector("integer", n_sample),
                          time_to_death = vector("integer", n_sample),
                          time_to_VF = vector("integer", n_sample))
  
  data_CR_E <- data.frame(id = vector("integer", n_sample),
                          treatment = vector("character", n_sample),
                          censor = vector("integer", n_sample),
                          event = vector("character", n_sample),
                          ftime = vector("integer", n_sample),
                          time_to_death = vector("integer", n_sample),
                          time_to_VF = vector("integer", n_sample))
  #ventilator status is not monitored after day 28
  ftime_CR_D <- NULL
  trt_CR_D <- NULL
  event_CR_D <- NULL
  censoring_CR_D <- NULL
  ID_CR_D <- NULL
  time_to_death_D <- NULL
  time_to_VF_D <- NULL
  
  ftime_CR_E <- NULL
  trt_CR_E <- NULL
  event_CR_E <- NULL
  censoring_CR_E <- NULL
  ID_CR_E <- NULL
  time_to_death_E <- NULL
  time_to_VF_E <- NULL
  
  for(i in 1: n_sample){
    ID_CR_D[i] <- data$id[i]
    ID_CR_E[i] <- data$id[i]
    trt_CR_D[i] <- data$treatment[i]
    trt_CR_E[i] <- data$treatment[i]
    event_CR_D[i] <- "Death"
    event_CR_E[i] <- "Extubation"
    time_to_death_D[i] <- data$time_to_death[i]
    time_to_VF_D[i] <- data$time_to_VF[i]
    time_to_death_E[i] <- data$time_to_death[i]
    time_to_VF_E[i] <- data$time_to_VF[i]
    
    if(data$time_to_VF[i] <= 28){
      if(data$time_to_death[i] <= data$time_to_VF[i]){
        ftime_CR_D[i] <- data$time_to_death[i]
        ftime_CR_E[i] <- data$time_to_death[i]
        censoring_CR_D[i] <- 0
        censoring_CR_E[i] <- 1
      }
      else if(data$time_to_death[i] > data$time_to_VF[i] & data$time_to_death[i] < 90){
        ftime_CR_D[i] <- data$time_to_death[i]
        ftime_CR_E[i] <- data$time_to_VF[i]
        censoring_CR_D[i] <- 0
        censoring_CR_E[i] <- 0
      }
      else if(data$time_to_death[i] > data$time_to_VF[i] & data$time_to_death[i] >= 90){
        ftime_CR_D[i] <- 90
        ftime_CR_E[i] <- data$time_to_VF[i]
        censoring_CR_D[i] <- 1
        censoring_CR_E[i] <- 0
      }
    }
    else if(data$time_to_VF[i] > 28){
      if(data$time_to_death[i] <= 28){
        ftime_CR_D[i] <- data$time_to_death[i]
        ftime_CR_E[i] <- data$time_to_death[i]
        censoring_CR_D[i] <- 0
        censoring_CR_E[i] <- 1
      }
      else if(data$time_to_death[i] > 28 & data$time_to_death[i] < 90){
        ftime_CR_D[i] <- data$time_to_death[i]
        ftime_CR_E[i] <- 28
        censoring_CR_D[i] <- 0
        censoring_CR_E[i] <- 1
      }
      else if(data$time_to_death[i] >= 90){
        ftime_CR_D[i] <- 90
        ftime_CR_E[i] <- 28
        censoring_CR_D[i] <- 1
        censoring_CR_E[i] <- 1
      }
    }
  }
  data_CR_D <- data_CR_D %>%
    mutate(id = ID_CR_D,
           treatment = trt_CR_D,
           censor = censoring_CR_D,
           event = event_CR_D,
           ftime = ftime_CR_D,
           time_to_death = time_to_death_D,
           time_to_VF = time_to_VF_D)
  
  #this can be used for the survival analysis later
  data_CR_E <- data_CR_E %>%
    mutate(id = ID_CR_E,
           treatment = trt_CR_E,
           censor = censoring_CR_E,
           event = event_CR_E,
           ftime = ftime_CR_E,
           time_to_death = time_to_death_E,
           time_to_VF = time_to_VF_E)
  
  data_CR <- rbind(data_CR_D, data_CR_E) 
  data_CR <- data_CR[order(data_CR$id),]
  #using r-inla package, the censoring status should be 1 for event and 0 for censored
  data_CR <- data_CR %>%
    mutate(status_CR_inla = ifelse(censor == 1, 0, 1))
  return(data_CR = data_CR)
}

CR_single_sim <- function(n_sample,        #trial sample size
                          F_90_d_trt,    #cumulative hazard at day 90
                          F_90_d_ctl,    ##cumulative hazard at day 90
                          F_28_v_trt,    #cumulative liberation rate from ventilator at day 28
                          F_28_v_ctl,    #cumulative liberation rate from ventilator at day 28
                          alpha_d,         #weibull distribution shape parameter for death
                          alpha_v         #weibull distribution shape parameter for liberation
){
  data_CR <- data_gen_weibull_CR(n_sample, F_90_d_trt, F_90_d_ctl, F_28_v_trt, F_28_v_ctl, alpha_d, alpha_v)
  #####################################################
  #####    Analysis               #####################
  #####################################################
  tryCatch({
    lc <- inla.make.lincombs('treatmenttrt:eventExtubation' = c(0, 1), treatmenttrt = c(1, 1), eventExtubation = c(0, 0))
    CR_inla <- inla(inla.surv(ftime, status_CR_inla) ~ treatment * event, data = data_CR, family = "coxph",
                    control.fixed = list(mean= 0, prec= 0.1),
                    lincomb = lc,
                    control.compute = list(dic = T, cpo = T, config = T))
    summary(CR_inla)
    
    ###################################################
    ###        posterior probabilities         ########
    ###################################################
    post_prob_death_CR <- inla.pmarginal(0, CR_inla$marginals.lincomb.derived$lc1)
    post_prob_liberation_CR <- 1 - inla.pmarginal(0, CR_inla$marginals.lincomb.derived$lc2)
    
    #trt_CR_inla is towards mortality improvement
    #trt_CR_inla + extubation_CR_inla + interaction_trt_ext_CR_inla is towards liberation rate improvement
    
    ###################################################
    ###        treatment effect estimates      ########
    ###################################################
    trt_death_CR <- summary(CR_inla)$fixed[2,1] 
    trt_liberation_CR <- summary(CR_inla)$fixed[3,1]
    trt_interaction_trt_ext_CR_inla <- summary(CR_inla)$fixed[4,1]
    
    ###################################################
    ###        sd estimates                    ########
    ###################################################
    sd_death_CR <- summary(CR_inla)$fixed[2,2] 
    sd_liberation_CR <- summary(CR_inla)$fixed[3,2]
    sd_interaction_trt_ext_CR_inla <- summary(CR_inla)$fixed[4,2]
    
    return(c(post_prob_death_CR = post_prob_death_CR, 
             post_prob_liberation_CR = post_prob_liberation_CR,
             trt_death_CR = trt_death_CR,
             trt_liberation_CR = trt_liberation_CR,
             trt_interaction_trt_ext_CR_inla = trt_interaction_trt_ext_CR_inla,
             sd_death_CR = sd_death_CR,
             sd_liberation_CR = sd_liberation_CR,
             sd_interaction_trt_ext_CR_inla = sd_interaction_trt_ext_CR_inla))
    
  },
  error = function(e){
    print("Error in INLA")
    ###################################################
    ###        posterior probabilities         ########
    ###################################################
    post_prob_death_CR <- "Error"
    post_prob_liberation_CR <- "Error"
    
    ###################################################
    ###        treatment effect estimates      ########
    ###################################################
    trt_death_CR <- "Error"
    trt_liberation_CR <- "Error"
    trt_interaction_trt_ext_CR_inla <- "Error"
    ###################################################
    ###        sd estimates                    ########
    ###################################################
    sd_death_CR <- "Error"
    sd_liberation_CR <- "Error"
    sd_interaction_trt_ext_CR_inla <- "Error"
    
    return(c(post_prob_death_CR = post_prob_death_CR, 
             post_prob_liberation_CR = post_prob_liberation_CR,
             trt_death_CR = trt_death_CR,
             trt_liberation_CR = trt_liberation_CR,
             trt_interaction_trt_ext_CR_inla = trt_interaction_trt_ext_CR_inla,
             sd_death_CR = sd_death_CR,
             sd_liberation_CR = sd_liberation_CR,
             sd_interaction_trt_ext_CR_inla = sd_interaction_trt_ext_CR_inla))
  }
  )
}


#running the simulation for 1000 times to calculate the lambda to control for type one error 
inla.setOption(num.threads = 8)
n_sim <- 2000
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)
Sim_Res <- foreach(k = 1:n_sim, .combine = rbind,
                   .packages = c("boot", "INLA", "dplyr", "tidyverse", "copula")
) %dopar% { 
  set.seed(k + 5125)
  CR_single_sim(n_sample = 1500,        #trial sample size
                F_90_d_trt = 0.42,    #cumulative hazard at day 90
                F_90_d_ctl = 0.37,    ##cumulative hazard at day 90
                F_28_v_trt = 0.93,    #cumulative liberation rate from ventilator at day 28
                F_28_v_ctl = 0.86,    #cumulative liberation rate from ventilator at day 28
                alpha_d = 0.9,         #weibull distribution shape parameter for death
                alpha_v = 2       #weibull distribution shape parameter for liberation
  )     
}
stopCluster(cl)

results <- as.data.frame(Sim_Res)

#write.csv(results, file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/CR_noeffect.csv")

#save.image(file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/CR_noeffect.Rdata")

write.csv(results, file = "C:/Users/Anna Heath/OneDrive - SickKids/Admin Anna/Hiring and Staff/Ziming Chen 2022/PRACTICAL/SimulationMar20/Weibull_AH/CR/CR_imp_lib_worse_mort.csv")
