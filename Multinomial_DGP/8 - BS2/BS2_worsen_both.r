################################################################################################
##########                    Simulations for PRACTICAL                   ######################
########                    Linear Regression+KG&J                                   ###########
################################################################################################
library(boot)
library(dplyr)
library(doParallel)
library(foreach)
library(tidyverse)
library(INLA)

#this is generating multinomial probabilities for this data generation mechanism
multinomial_data_gen_same_effect <- function(OR){
  VFDs <- c(seq(-1, 28, 1))[-30]
  N <- c(115, 14, 2, 2, 7, 4, 2, 3, 2, 1, 2, 2, 3, 5, 3, 5, 5, 6, 6, 6, 5, 10, 11, 10, 8, 13, 25, 21, 17, 15)[-30]
  
  #control group baseline proportion in each category
  ctl_prop <- N/sum(N)
  ctl_prop_cumsum <- cumsum(ctl_prop)
  ctl_odds <- ctl_prop_cumsum/(1-ctl_prop_cumsum)
  
  #treatment no effect
  #OR <- 1
  trt_odds <- OR*ctl_odds
  trt_prop_cumsum <- c((trt_odds/(1+trt_odds))[-29], 1)
  trt_prop <- c(trt_prop_cumsum[1],diff(trt_prop_cumsum))
  
  #######################################################################
  #####   Generating the patient matrix based on above ##################
  #######################################################################
  treatment <- c("trt", "ctl")
  
  data <- as.data.frame(matrix(nrow = 1500, ncol = 4))
  colnames(data) <- c("ID", "treatment", "VFDs", "death_time")
  for(i in 1 : nrow(data)){
    data[i, 1] <- i
    data[i, 4] <- ceiling(runif(1, min = 0, max = 90))
    data[i, 2] <- sample(treatment, 1, prob = c(0.5, 0.5), replace = TRUE)
    if(data[i, 2] == "trt"){
      data[i, 3] <- sample(VFDs, 1, prob = trt_prop)
    }
    else{
      data[i, 3] <- sample(VFDs, 1, prob = ctl_prop)
    }
  }
  data <- data %>%
    mutate(death = ifelse(VFDs == -1, 1, 0))
  return(data)
}


BS2_single_sim <- function(OR){
  data <- multinomial_data_gen_same_effect(OR)
  sub_data <- data %>% 
    filter(death != 1) %>%
    mutate(censor = ifelse(VFDs == 0, 0, 1),
           time = ifelse(VFDs == 0, 28, 28 - VFDs))
  #####################################################
  #####    Analysis               #####################
  #####################################################
  tryCatch({
    #-------------third definition of ordinal outcome using INLA------------
    #################      modeling binary death           #######################
    result_BM_death <- inla(death ~ treatment, 
                            data = data,
                            family = "binomial",
                            control.fixed = list(mean= 0, prec= 0.1),
                            control.compute = list(dic = T, cpo = T, config = T),
                            verbose = F)
    
    ##########               modeling liberation using survival analysis            #########################
    result_surv_extubation <- inla(inla.surv(time, censor) ~ treatment, 
                                   data = sub_data,
                                   family = "coxph",
                                   control.fixed = list(mean= 0, prec= 0.1),
                                   control.compute = list(dic = T, cpo = T, config = T),
                                   verbose = F)
    ###################################################
    ###        posterior probabilities         ########
    post_prob_death_BM <- inla.pmarginal(0, result_BM_death$marginals.fixed$treatmenttrt)
    post_prob_extubation_surv <- 1 - inla.pmarginal(0, result_surv_extubation$marginals.fixed$treatmenttrt)
    
    ###################################################
    ###        treatment effect estimates      ########
    ###################################################
    trt_death_BM  <- summary(result_BM_death)$fixed[2,1]
    trt_extubation_surv  <- summary(result_surv_extubation)$fixed[2,1]
    
    ###################################################
    ###        sd estimates                    ########
    ###################################################
    sd_death_BM  <- summary(result_BM_death)$fixed[2,2]
    sd_extubation_surv  <- summary(result_surv_extubation)$fixed[2,2]
    
    return(c(post_prob_death_BM = post_prob_death_BM, 
             post_prob_extubation_surv = post_prob_extubation_surv,
             trt_death_BM = trt_death_BM,
             trt_extubation_surv = trt_extubation_surv,
             sd_death_BM = sd_death_BM,
             sd_extubation_surv = sd_extubation_surv))
    
  },
  error = function(e){
    print("Error in INLA")
    ###################################################
    ###        posterior probabilities         ########
    ###################################################
    post_prob_death_BM <- "Error"
    post_prob_extubation_surv <- "Error"
    
    ###################################################
    ###        treatment effect estimates      ########
    ###################################################
    trt_death_BM <- "Error"
    trt_extubation_surv <- "Error"
    
    ###################################################
    ###        sd estimates                    ########
    ###################################################
    sd_death_BM <- "Error"
    sd_extubation_surv <- "Error"
    
    return(c(post_prob_death_BM = post_prob_death_BM, 
             post_prob_extubation_surv = post_prob_extubation_surv,
             trt_death_BM = trt_death_BM,
             trt_extubation_surv = trt_extubation_surv,
             sd_death_BM = sd_death_BM,
             sd_extubation_surv = sd_extubation_surv))
  }
  )
}


#running the simulation for 1000 times to calculate the lambda to control for type one error 
n_sim <- 2000
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)
Sim_Res <- foreach(k = 1:n_sim, .combine = rbind,
                   .packages = c("boot", "dplyr", "tidyverse", "INLA")
) %dopar% { 
  set.seed(k + 5125)
  BS2_single_sim(OR = 1.3)     
}
stopCluster(cl)

results <- as.data.frame(Sim_Res)

#write.csv(results, file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.csv")

#save.image(file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.Rdata")

write.csv(results, file = "C:/Users/Anna Heath/OneDrive - SickKids/Admin Anna/Hiring and Staff/Ziming Chen 2022/PRACTICAL/SimulationMar20/Multinomial_AH/BS2/BS2_worsen_both.csv")

save.image(file = "C:/Users/Anna Heath/OneDrive - SickKids/Admin Anna/Hiring and Staff/Ziming Chen 2022/PRACTICAL/SimulationMar20/Multinomial_AH/BS2/BS2_worsen_both.Rdata")


