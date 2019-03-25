library(deSolve)
library(ggplot2)
library(reshape)
library(cowplot)
library(tidyverse)
library(DEoptim)
library(R.utils)



draw_optimize_fun <- function(params, check) {
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
    n = 2 * 5000,
    mean = 0,
    sd = 1
  ), nrow = 2)
  cholesky <-
    matrix(data = c(sigma ^ 2, ro * (sigma ^ 2), ro * (sigma ^ 2), sigma ^
                      2),
           nrow = 2)
  g <- cholesky %*% z
  M_C1_temp <- 0.4 * exp(g[1,])
  M_C2_temp <- 0.2 * exp(g[2,])
  
  #redraw fluxs for C3, rather then make a 3x matrix, its easier to do this.
  z <- matrix(data = rnorm(
    n = 2 * 5000,
    mean = 0,
    sd = 1
  ), nrow = 2)
  
  g <- cholesky %*% z
  #we added consumer 3 so we estimate its growth rate
  M_C3_temp <- params[7] * exp(g[1,])
  
  
  
  #enviro flux preds
  z_pred <-
    matrix(data = rnorm(
      n = 2 * 5000,
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
  
  M_P1_temp <- 0.08 * exp(g_pred[1,])
  #we added pred 2 so its growth needs to be estimated
  M_P2_temp <- params[8] * exp(g_pred[2,])
  
  
  #enviro flux on res
  z_res <- rnorm(n = 5000, mean = 0, sd = 1)
  
  res_cholesky <-
    matrix(data = c(res_sigma ^ 2, res_ro * (res_sigma ^ 2)), ncol = 1)
  
  g_res <- res_cholesky %*% z_res
  
  r_temp <- 1 * exp(g_res[1,])
  
  
  #define state variables
  State <- c(
    P1 = 1,
    P2 = 1,
    C1 = 1,
    C2 = 1,
    C3 = 1,
    R = 1
  )
  Time <- seq(0, 5000, by = 1)
  
  results <- matrix(data = NA,
                    nrow = 5000,
                    ncol = 6)
  results[1,] <- State
  
  for (t in 2:(5000)) {
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
    
    
    #if the eval produces errorrs or takes a long time don't bother just throw an error
    tryCatch({
      VF_out <- withTimeout({
        as.data.frame(ode(
          func = VassFox_Cvar_2P_3C,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ))
      },
      timeout = 2,
      events = list(func = eventfun))
    },
    error = function(e) {
      an.error.occured <<- TRUE
    },
    TimeoutException = function(e) {
      an.error.occured <<- TRUE
    })
    
    #if the call to the ode solver produced and error this is a bad parameter set
    if (an.error.occured) {
      return(Inf)
    }
    
    results[t,] <-
      c(VF_out[2, 2], VF_out[2, 3], VF_out[2, 4], VF_out[2, 5], VF_out[2, 6], VF_out[2, 7])
    
    
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
  
  d$time <- c(rep(seq(1, 5000), 6))
  
  #if just checking to see which work return text
  if (check == TRUE) {
    if (any(d$value < 0.05)) {
      return("Crash")
    }
    else{
      return("Good")
    }
  }
  
  #if check is false return pictures
  if (check == FALSE) {
    g <-
      ggplot(d, aes(x = time, y = value, color = variable)) + geom_line() + facet_wrap(~ variable, scales =
                                                                                         "free") + ylim(0, max(d$value))
    
    
    return(g)
    
  }
}


#read in parameter values that might work
opted_parm <-
  read.csv("./eco_working_group/network-coexistence/Estimated_stable_params.csv",
           header = TRUE)


opt_remove_index <- opted_parm[, -1]
opt_remove_index <- opt_remove_index[, -19]

data <- unique(opt_remove_index)

#check that none of the pops drop below a cut off
no_crash <- NULL
for (i in 1:length(data$par1)) {
  g <- draw_optimize_fun(as.numeric(data[i, ]), TRUE)
  if (g == "Good") {
    no_crash <- c(no_crash, i)
  }
}

#this is the list of good parameter values
print(no_crash)

#working param sets: 2(good),15(pop migh crash),17(good),36(uncompled, pretty flat dynamics),46(neat larger dyanmics),68 (not stable looking)

#what parameters are
#ro <- params[1]
#pred_ro <- params[2]
#res_ro <- params[3]
#sigma <- params[4]
#pred_sigma <- params[5]
#res_sigma <- params[6]
#M_C3_temp <- params[7]
#M_P2_temp <- params[8]
#J_C3 = params[9]
#J_P2 = params[10]
#R_0_3 = params[11]
#O_P1_C1 = params[12]
#O_P2_C1 = params[13]
#O_P1_C2 = params[14]
#O_P2_C2 = params[15]
#O_C1_R = params[16]
#O_C2_R = params[17]
#O_C3_R = params[18]


#these are the 5 parameter sets that seem to have resonable dynamics. 

data[2, ]
#par1       par2        par3       par4      par5      par6    par7       par8      par9     par10     par11     par12      par13     par14     par15     par16     par17     par18
#-0.9948819 -0.1463557 -0.08771538 0.09810206 0.4959867 0.6458228 0.13142 0.08108551 0.8224542 0.2919573 0.9526727 0.7185076 0.08430302 0.1647937 0.5594141 0.8047696 0.9846117 0.5688471

data[17, ]
#par1       par2       par3      par4      par5      par6       par7      par8      par9     par10   par11     par12     par13     par14     par15     par16     par17     par18
#0.6204637 0.05813391 -0.7418263 0.2215021 0.4562604 0.4393107 0.07764343 0.0765969 0.4086876 0.2356045 0.67062 0.6451473 0.0298092 0.2679468 0.7543071 0.8062275 0.9879175 0.5240379

data[36, ]
#par1       par2       par3      par4      par5      par6       par7      par8      par9     par10     par11     par12      par13     par14     par15     par16     par17    par18
#0.4733374 -0.3684155 -0.5820099 0.1528808 0.3657811 0.3046128 0.06546211 0.0785416 0.5845653 0.2682388 0.9016482 0.7754392 0.02568254 0.2145904 0.8412819 0.7867506 0.9803509 0.322561

data[46, ]
#par1      par2       par3      par4      par5      par6       par7       par8      par9     par10     par11     par12      par13     par14     par15     par16     par17     par18
#-0.1047472 0.2080639 -0.6394717 0.1878601 0.3905516 0.5217772 0.04218716 0.08408349 0.6186747 0.2919329 0.9273697 0.6955827 0.01753484 0.2386237 0.8902442 0.7653685 0.9834837 0.2176153

data[68, ]
#par1       par2      par3      par4      par5      par6       par7      par8      par9     par10     par11     par12      par13      par14    par15     par16     par17    par18
#0.2447766 0.03046945 0.1919667 0.1613392 0.1381634 0.3399561 0.07790641 0.1778277 0.4890899 0.6040426 0.7010799 0.8418758 0.08309623 0.09946049 0.721384 0.8183416 0.9833202 0.451942

g1 <- draw_optimize_fun(as.numeric(data[2, ]), FALSE)


g2 <- draw_optimize_fun(as.numeric(data[17, ]), FALSE)

g3 <- draw_optimize_fun(as.numeric(data[36, ]), FALSE)

g4 <- draw_optimize_fun(as.numeric(data[46, ]), FALSE)

g5 <- draw_optimize_fun(as.numeric(data[68, ]), FALSE)


save_plot("param1.pdf", g1, base_aspect_ratio = 3)
save_plot("param2.pdf", g2, base_aspect_ratio = 3)
save_plot("param3.pdf", g3, base_aspect_ratio = 3)
save_plot("param4.pdf", g4, base_aspect_ratio = 3)
save_plot("param5.pdf", g5, base_aspect_ratio = 3)
