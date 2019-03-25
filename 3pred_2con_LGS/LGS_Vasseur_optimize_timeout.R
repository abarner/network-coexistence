library(deSolve)
library(ggplot2)
library(reshape)
library(cowplot)
library(tidyverse)
library(DEoptim)
library(R.utils)



optimize_fun <- function(params) {
  
  
  VassFox_Cvar_2P_3C <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      #M_C1_j = M_C1_0 * exp()
      #M_C2_j = M_C2_0 * exp()
      
      
      #variables to make math easier to read
      
      
      all_c_p1 <-
        ((O_P1_C1 * C1) +  (O_P1_C2 * C2) + (1 - (O_P1_C1 + O_P1_C2)) * C3)
      
      
      all_c_p2 <-
        ((O_P2_C1 * C1) +  (O_P2_C2 * C2) + (1 - (O_P2_C1 + O_P2_C2)) * C3)
      
      both_preds_eat_C1 <-
        ((O_P1_C1 * J_P1 * P1 * C1) / (all_c_p1 + all_c_p2 + C_0)) + ((O_P2_C1 * J_P2 * P2 * C1) / (all_c_p1 +
                                                                                                      all_c_p2 + C_0))
      
      both_preds_eat_C2 <-
        ((O_P1_C2 * J_P1 * P1 * C2) / (all_c_p1 + all_c_p2 + C_0)) + ((O_P2_C2 * J_P2 * P2 * C2) / (all_c_p1 +
                                                                                                      all_c_p2 + C_0))
      
      
      both_preds_eat_C3 <-
        (((1 - (
          O_P1_C1 + O_P1_C2
        )) * J_P1 * P1 * C3) / (all_c_p1 + all_c_p2 + C_0)) + (((1 - (
          O_P2_C1 + O_P2_C2
        )) * J_P2 * P2 * C3) / (all_c_p1 + all_c_p2 + C_0))
      
      
      
      dP_1 = -(M_P1 * P1) + (((J_P1 * P1) * (
        (O_P1_C1 * C1) + (O_P1_C2 * C2) + (1 - (O_P1_C1 + O_P1_C2)) * C3
      )) / (all_c_p1 + all_c_p2 + C_0))
      
      
      dP_2 = -(M_P2 * P2) + (((J_P2 * P2) * (
        (O_P2_C1 * C1) + (O_P2_C2 * C2) + (1 - (O_P2_C1 + O_P2_C2)) * C3
      )) / (all_c_p1 + all_c_p2 + C_0))
      
      
      #dP = - (M_P * P) + ( ( (J_P * P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
      #dC1 = - (M_C1 * C1) + ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_P1_C1 * J_P1 * P1 * C1) / ( (O_P1_C1 * C1) + ((1 - O_P1_C1) * C2) + C_0) )
      
      
      
      
      dC1 = -(M_C1 * C1) + ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - both_preds_eat_C1
      dC2 = -(M_C2 * C2) + ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2)) - both_preds_eat_C2
      dC3 = -(M_C3 * C3) + ((O_C3_R * J_C3 * C3 * R) / (R + R_0_3)) - both_preds_eat_C3
      
      dR = r * R * (1 - (R / K)) - ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2)) -  ((O_C3_R * J_C3 * C3 * R) / (R + R_0_3))
      
      return(list(c(dP_1, dP_2, dC1, dC2, dC3, dR)))
      
    })
    
    
  }
  
  
  #params 1,2, and 3 are the ro values, hence they can be negaitive (-1 to 1)
  #all other parameters are between 0 and 1
  
  #enviro flux parameters
  
  ro <- params[1]
  pred_ro <- params[2]
  res_ro <- params[3]
  
  
  sigma <- params[4]
  pred_sigma <- params[5]
  res_sigma <- params[6]
  
  
  #calc enivor flux on consumers
  z <- matrix(data = rnorm(
    n = 2 * 1000,
    mean = 0,
    sd = 1
  ), nrow = 2)
  cholesky <-
    matrix(data = c(sigma ^ 2, ro * (sigma ^ 2), ro * (sigma ^ 2), sigma ^
                      2),
           nrow = 2)
  g <- cholesky %*% z
  M_C1_temp <- 0.4 * exp(g[1, ])
  M_C2_temp <- 0.2 * exp(g[2, ])
  
  #redraw fluxs for C3, rather then make a 3x matrix, its easier to do this.
  z <- matrix(data = rnorm(
    n = 2 * 1000,
    mean = 0,
    sd = 1
  ), nrow = 2)
  
  g <- cholesky %*% z
  #we added consumer 3 so we estimate its growth rate
  M_C3_temp <- params[7] * exp(g[1, ])
  
  
  
  #enviro flux preds
  z_pred <-
    matrix(data = rnorm(
      n = 2 * 1000,
      mean = 0,
      sd = 1
    ), nrow = 2)
  
  pred_cholesky <-
    matrix(
      data = c(
        pred_sigma ^ 2,
        pred_ro * (pred_sigma ^ 2),
        pred_ro * (pred_sigma ^ 2),
        pred_sigma ^ 2
      ),
      nrow = 2
    )
  
  g_pred <- pred_cholesky %*% z_pred
  
  M_P1_temp <- 0.08 * exp(g_pred[1, ])
  #we added pred 2 so its growth needs to be estimated
  M_P2_temp <- params[8] * exp(g_pred[2, ])
  
  
  #enviro flux on res
  z_res <- rnorm(n = 1000, mean = 0, sd = 1)
  
  res_cholesky <-
    matrix(data = c(res_sigma ^ 2, res_ro * (res_sigma ^ 2)), ncol = 1)
  
  g_res <- res_cholesky %*% z_res
  
  r_temp <- 1 * exp(g_res[1, ])
  
  
  #define state variables
  State <- c(
    P1 = 1,
    P2 = 1,
    C1 = 1,
    C2 = 1,
    C3 = 1,
    R = 1
  )
  Time <- seq(0, 1000, by = 1)
  
  results <- matrix(data = NA,
                    nrow = 1000,
                    ncol = 6)
  results[1, ] <- State
  
  for (t in 2:(1000)) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    M_C3 <- M_C3_temp[t]
    
    M_P1 <- M_P1_temp[t]
    M_P2 <- M_P2_temp[t]
    r <- r_temp[t]
    # from Table 1
    pars <- c(
      # resource intrinsic rate of growth
      #r = 1.0,
      r = r,
      # resource carrying capacity
      K = 1.0,
      # consumer 1 ingestion rate
      J_C1 = 0.8036,
      # consumer 2 ingestion rate
      J_C2 = 0.7,
      # consumer 3 ingestion rate (we added this)
      #J_C3 = 0.923,
      J_C3 = params[9],
      # predator ingestion rate
      J_P1 = 0.4,
      # predator 2 ingestion rate
      J_P2 = params[10],
      # medial consumer 1 mortality rate
      M_C1 = M_C1,
      # medial consumer 2 mortality rate
      M_C2 = M_C2,
      # medial consumer 3 mortality rate
      M_C3 = M_C3,
      # predator mortality rate
      #M_P = 0.08,
      M_P1 = M_P1,
      M_P2 = M_P2,
      # half saturation constant
      R_0_1 = 0.16129,
      R_0_2 = 0.9,
      #need to extimate because we added another consumer
      R_0_3 = params[11],
      
      C_0 = 0.5,
      # preference coefficients (need to estimate all of these)
      # preference for consumer 3 is just 1-(pref1+pref2), so pref1+pref2 <=1
      
      
      O_P1_C1 = params[12],
      O_P2_C1 = params[13],
      
      O_P1_C2 = params[14],
      O_P2_C2 = params[15],
      
      
      O_C1_R = params[16],
      O_C2_R = params[17],
      O_C3_R = params[18],
      
      #don't need to pass these
      sigma = sigma,
      ro = ro
    )
    
    
    ## preference for consumer 3 is just 1-(pref1+pref2), so pref1+pref2 <=1
    if ((params[12] + params[14]) > 1) {
      return(Inf)
    }
    if ((params[13] + params[15]) > 1) {
      return(Inf)
    }
    
    
    
    State <-
      c(
        P1 = results[t - 1, 1],
        P2 = results[t - 1, 2],
        C1 = results[t - 1, 3],
        C2 = results[t - 1, 4],
        C3 = results[t - 1, 5],
        R = results[t - 1, 6]
      )
    
    an.error.occured <- FALSE

    
    #uf the eval produces errorrs or takes a long time don't bother just throw an error
    tryCatch({
      VF_out <- withTimeout( {
        as.data.frame(
          ode(
            func = VassFox_Cvar_2P_3C,
            y = State,
            parms = pars,
            times = seq(0, 1)
          ) ) },
          timeout = 2,
          events = list(func = eventfun) )
    }, 
    error = function(e) {
     an.error.occured <<- TRUE
    },
    TimeoutException = function(e) { an.error.occured <<- TRUE
    }
    
    )
    
    #if the call to the ode solver produced and error this is a bad parameter set
    if (an.error.occured) {
      return(Inf)
    }
    
    results[t, ] <-
      c(VF_out[2, 2], VF_out[2, 3], VF_out[2, 4], VF_out[2, 5], VF_out[2, 6], VF_out[2, 7])
    
    #d_check<- data.frame(pred1=results[,1],pred2=results[,2],con1=results[,3],con2=results[,4],con3=results[,5],res=results[,6])
    
    #suppressMessages({s<-melt(d_check)})
    
    
    #we don't want dynamics that crash so if any of the pops get below some lower theshold, we consider that the pop crashed.
    #if(any(na.omit(s$value) < 0.00001)){
    # return(Inf)
    #}
    
  }
  
  dat <-
    data.frame(
      pred1 = results[, 1],
      pred2 = results[, 2],
      con1 = results[, 3],
      con2 = results[, 4],
      con3 = results[, 5],
      res = results[, 6]
    )
  suppressMessages({
    d <- melt(dat)
  })
  
  d$time <- c(rep(seq(1, 1000), 6))
  
  if(any(is.na(d$value))){
    return(Inf)
  }
  
  d <- d[d$time > 500, ]
  #ggplot(d, aes(x=time,y=value,color = variable))+geom_line()+facet_wrap( ~ variable, scales="free")
  
  
 # if(any(d$value > 5)){
  #  return(Inf)
  #}
  #anytime a population crashes below this threshhold it get penalized in the scoring function
  Crash_penality <- sum(d$value < 0.05)
  

  return(Crash_penality)
}

#define lower bounds on params, the first 3 are the sig values for envir var they can be negative, the rest of the parameters (4-18) are bound between zero and 1
lower <- c(rep(-1, 3), rep(0, 15))
high <- c(rep(1, 18))
opt_out<-DEoptim(
  optimize_fun,
  lower = lower,
  upper = high,
  DEoptim.control(
    NP=2000,
    itermax=125,
    strategy=3,
    c=0.05,
    parallelType = 1,
    packages = c("deSolve", "reshape", "tidyverse", "R.utils")
  )
)

best_memb <- opt_out$member$bestmemit
best_val <-opt_out$member$bestvalit

all<-as.data.frame(cbind(best_memb,best_val))

results<-all[all$best_val == 0,]

write.csv(results, file="Estimated_stable_params.csv")
