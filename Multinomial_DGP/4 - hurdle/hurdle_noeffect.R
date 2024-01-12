#############################################################
########### multinomial hurdle model no effect ##############
#############################################################
library(tidyverse)
library(dplyr)
library(boot)
library(doParallel)
library(foreach)
library(brms)

#this is generating multinomial probabilities for this data generation mechanism
multinomial_data_gen_same_effect <- function(OR){
  VFDs <- c(seq(-1, 28, 1))[-2]
  N <- c(115, 14, 2, 2, 7, 4, 2, 3, 2, 1, 2, 2, 3, 5, 3, 5, 5, 6, 6, 6, 5, 10, 11, 10, 8, 13, 25, 21, 17, 15)[-2]
  
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
  
  data <- as.data.frame(matrix(nrow = 1500, ncol = 3))
  colnames(data) <- c("ID", "treatment", "VFDs")
  for(i in 1 : nrow(data)){
    data[i, 1] <- i
    data[i, 2] <- sample(treatment, 1, prob = c(0.5, 0.5), replace = TRUE)
    if(data[i, 2] == "trt"){
      data[i, 3] <- sample(VFDs, 1, prob = trt_prop)
    }
    else{
      data[i, 3] <- sample(VFDs, 1, prob = ctl_prop)
    }
  }
  return(data)
}

#OR = 1
PoissonHurdle_single_sim <- function(OR){
  data <- multinomial_data_gen_same_effect(OR)
  data <- data %>% 
    mutate(VFDs = VFDs + 1)
  #####################################################
  #####    Analysis               #####################
  #####################################################
  tryCatch({
    #----------------------------------------------------------------------
    #-------------  logit-poisson hurdle model   --------------------------
    #----------------------------------------------------------------------
    
    #-------------third definition of ordinal outcome using INLA------------
    hurdle <- brm(bf(VFDs ~ treatment, 
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
ncores <- detectCores() - 2
cl <- makeCluster(ncores)
registerDoParallel(cl)
Sim_Res <- foreach(k = 1:n_sim, .combine = rbind,
                   .packages = c("boot", "dplyr", "tidyverse", "brms", "copula")
) %dopar% { 
  set.seed(k + 5125)
  PoissonHurdle_single_sim(OR = 1)     
}
stopCluster(cl)

results <- as.data.frame(Sim_Res)

#write.csv(results, file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.csv")

#save.image(file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.Rdata")

write.csv(results, file = "C:/Users/Anna Heath/OneDrive/OneDrive - SickKids/Admin Anna/Hiring and Staff/Ziming Chen 2022/PRACTICAL/SimulationMar20/Multinomial_AH/hurdle/hurdle_noeffect.csv")
