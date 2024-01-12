################################################################################################
##########                    Simulations for PRACTICAL                   ######################
########                    Linear Regression+KG&J                                   ###########
################################################################################################
library(boot)
library(dplyr)
library(doParallel)
library(foreach)
library(tidyverse)
#install_github('aejensen/TruncComp')
library(TruncComp)

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




KJL_single_sim <- function(OR){
  data <- multinomial_data_gen_same_effect(OR)
    data <- data %>% 
    mutate(treatment = ifelse(treatment == "trt", 1, 0))
  #####################################################
  #####    Analysis               #####################
  #####################################################
  tryCatch({
    #-------------third definition of ordinal outcome using INLA------------
    model <- truncComp(VFDs ~ treatment, atom = 0, data = data, method="SPLRT")
    trt_effect_vector <- NULL
    for(b in 1: 1000){
      data2 <- data[sample(1:nrow(data), size = nrow(data), replace = TRUE), ]
      model_lm <- lm(VFDs ~ treatment, data = data2)
      trt_effect_vector[b] <- coef(model_lm)[2]
    }
    trt <- mean(na.omit(trt_effect_vector))
    sd <- sd(na.omit(trt_effect_vector))
    ###################################################
    ###        p values                        ########
    ###################################################
    p_val <- model$p
    
    ###################################################
    ###  positive treatment effect yes or no   ########
    ###################################################
    trt_effect_YN <- ifelse(p_val < 0.05 & trt > 0, "Yes", "No")
    
    return(c(trt_effect_YN,
             p_val,
             trt,
             sd))
  },
  error = function(e){
    print("Error")
    trt_effect_YN <- "Error"
    p_val <- "Error"
    trt <- "Error"
    sd <- "Error"
    
    return(c(trt_effect_YN,
             p_val,
             trt,
             sd))
  }
  )
}


#running the simulation for 1000 times to calculate the lambda to control for type one error 
n_sim <- 2000
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)
Sim_Res <- foreach(k = 1:n_sim, .combine = rbind,
                   .packages = c("boot", "dplyr", "tidyverse", "TruncComp")
) %dopar% { 
  set.seed(k + 5125)
  KJL_single_sim(OR = 1.25)     
}
stopCluster(cl)

results <- as.data.frame(Sim_Res)

#write.csv(results, file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.csv")

#save.image(file = "C:/Users/Ziming Chen/OneDrive - SickKids/PRACTICAL/SimulationFinal/ordinal1_noeffect.Rdata")

write.csv(results, file = "C:/Users/Anna Heath/OneDrive/OneDrive - SickKids/Admin Anna/Hiring and Staff/Ziming Chen 2022/PRACTICAL/SimulationMar20/Multinomial_AH/KJL/LR_worsen_both.csv")



