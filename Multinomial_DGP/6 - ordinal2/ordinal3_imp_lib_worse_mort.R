######################################################################
###########     multinomial ordinal 3 imp lib worsen mort        #####
######################################################################

library(boot)
library(dplyr)
library(doParallel)
library(foreach)
library(tidyverse)
library(INLA)

#this is generating multinomial probabilities for this data generation mechanism
multinomial_data_gen_diff_effect <- function(OR_mort, OR_lib){
  VFDs <- c(seq(-1, 28, 1))[-30]
  N <- c(115, 14, 2, 2, 7, 4, 2, 3, 2, 1, 2, 2, 3, 5, 3, 
         5, 5, 6, 6, 6, 5, 10, 11, 10, 8, 13, 25, 21, 17, 15)[-30]
  
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
  
  #######################################################################
  #####   Generating the patient matrix based on above ##################
  #######################################################################
  treatment <- c("trt", "ctl")
  
  data <- as.data.frame(matrix(nrow = 1500, ncol = 4))
  colnames(data) <- c("ID", "treatment", "VFDs","death_time")
  for(i in 1 : nrow(data)){
    data[i, 1] <- i
    data[i, 4] <- ceiling(runif(1, min = 0, max = 90))
    data[i, 2] <- sample(treatment, 1, prob = c(0.5, 0.5), replace = TRUE)
    if(data[i, 2] == "trt"){
      data[i, 3] <- sample(VFDs, 1, prob = trt_diff_effect_prop)
    }
    else{
      data[i, 3] <- sample(VFDs, 1, prob = ctl_prop)
    }
  }
  data <- data %>%
    mutate(VFDs2 = case_when(VFDs == -1 ~ 1,
                             VFDs == 0 ~ 2,
                             VFDs >= 1 & VFDs <= 5 ~ 3,
                             VFDs >= 6 & VFDs <= 12 ~ 4,
                             VFDs >= 13 & VFDs <= 17 ~ 5,
                             VFDs >= 18 & VFDs <= 21 ~ 6,
                             VFDs >= 22 & VFDs <= 24 ~ 7,
                             VFDs == 25 ~ 8,
                             VFDs == 26 ~ 9,
                             VFDs == 27 ~ 10))
  return(data)
}#OR = 1

ordinal3_single_sim <- function(OR_mort, OR_lib){
  data <- multinomial_data_gen_diff_effect(OR_mort, OR_lib)
  #####################################################
  #####    Analysis               #####################
  #####################################################
  tryCatch({
    #----------------------------------------------------------------------
    #-------------   Ordinal Outcome   ------------------------------------
    #----------------------------------------------------------------------
    
    #-------------third definition of ordinal outcome using INLA------------
    lc <- inla.make.lincombs(treatmenttrt = -1)
    result_OO3 <- inla(VFDs2 ~ treatment, 
                       data = data,
                       family = "pom",
                       lincomb = lc,
                       control.fixed = list(mean= 0, prec= 0.1),
                       control.family = list(hyper = list(theta1 = list(param = 10))),
                       control.compute = list(dic = T, cpo = T, config = T),
                       verbose = F)
    
    post_prob_OO3 <- inla.pmarginal(0, result_OO3$marginals.lincomb.derived$lc00000001)
    
    
    ###################################################
    ###        treatment effect estimates      ########
    ###################################################
    trt_OO3  <- summary(result_OO3)$fixed[2,1]
    
    ###################################################
    ###        sd estimates                    ########
    ###################################################
    sd_OO3  <- summary(result_OO3)$fixed[2,2]
    
    return(c(post_prob_OO3,
             trt_OO3,
             sd_OO3
    ))
  },
  error = function(e){
    print("Error in INLA")
    post_prob_OO3 <- "Error"
    
    trt_OO3 <- "Error"
    
    sd_OO3 <- "Error"
    
    return(c(post_prob_OO3,
             trt_OO3,
             sd_OO3))
  }
  )
}

#OR_mort = 1.2 
#OR_lib = 0.8

#running the simulation for 1000 times to calculate the lambda to control for type one error 
inla.setOption(num.threads = 8)
n_sim <- 2000
#ncores <- detectCores() - 1
#cl <- makeCluster(ncores)
#registerDoParallel(cl)
Sim_Res <- foreach(k = 1:n_sim, .combine = rbind,
                   .packages = c("boot", "dplyr", "tidyverse", "INLA")
) %do% { 
  set.seed(k + 5125)
  ordinal3_single_sim(OR_mort = 1.25,
                      OR_lib = 0.75)     
}
stopCluster(cl)

results <- as.data.frame(Sim_Res)

#write.csv(results, file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.csv")

#save.image(file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.Rdata")

write.csv(results, file = "C:/Users/Anna Heath/OneDrive/OneDrive - SickKids/Admin Anna/Hiring and Staff/Ziming Chen 2022/PRACTICAL/SimulationMar20/Multinomial_AH/ordinal3/ordinal3_imp_lib_worse_mort.csv")

save.image(file = "/home/aheath/Ziming_Sim/ordinal3/ordinal3_imp_lib_worse_mort.Rdata")