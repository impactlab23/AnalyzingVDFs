#############################################################
########### multinomial t test no effect       ##############
#############################################################
library(tidyverse)
library(dplyr)
library(boot)
library(doParallel)
library(foreach)

#this is generating multinomial probabilities for this data generation mechanism
multinomial_data_gen_same_effect <- function(OR){
  VFDs <- c(seq(0, 28, 1))[-30]
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
  data[data == -1] <- 0
  return(data)
}

#OR = 1
ttest_single_sim <- function(OR){
  data <- multinomial_data_gen_same_effect(OR)
  data <- data %>% 
    mutate(VFDs = VFDs + 1)
  data <- data %>% 
    mutate(treatment = as.factor(ifelse(treatment == "trt", 1, 0)))
  #####################################################
  #####    Analysis               #####################
  #####################################################
  tryCatch({
    #-------------third definition of ordinal outcome using INLA------------
    t_test <- t.test(VFDs ~ treatment, data = data, alternative = "less")
    
    ###################################################
    ###        p values                        ########
    ###################################################
    p_val <- t_test$p.value
    
    ###################################################
    ###      treatment effect mean difference  ########
    ###################################################
    
    #mean in treatment group - mean in control group
    trt_effect <- as.numeric(t_test$estimate[2] - t_test$estimate[1])
    
    ###################################################
    ###        treatment effect yes or no      ########
    ###################################################
    trt_effect_YN <- ifelse(p_val < 0.05, "Yes", "No")
    
    return(c(trt_effect_YN,
             p_val,
             trt_effect))
  },
  error = function(e){
    print("Error")
    trt_effect_YN <- "Error"
    p_val <- "Error"
    trt_effect <- "Error"
    
    return(c(trt_effect_YN,
             p_val,
             trt_effect))
  }
  )
}


#running the simulation for 1000 times to calculate the lambda to control for type one error 
n_sim <- 2000
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)
Sim_Res <- foreach(k = 1:n_sim, .combine = rbind,
                   .packages = c("boot", "dplyr", "tidyverse")
) %dopar% { 
  set.seed(k + 5125)
  ttest_single_sim(OR = 1)     
}
stopCluster(cl)

results <- as.data.frame(Sim_Res)

write.csv(results, file = "C:/Users/Anna Heath/OneDrive/OneDrive - SickKids/Admin Anna/Hiring and Staff/Ziming Chen 2022/PRACTICAL/SimulationMar20/Multinomial_AH/t test/ttest_noeffect.csv")
