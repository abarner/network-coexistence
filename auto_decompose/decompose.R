# Implementation of the Vasseur and Fox Model with temporal variation in mortality rates
# Compare positive, no, and negative autocorrelation
# For tree species mortality rates through time (e.g. from Vasseur and Fox appendix)
# decompose into mechanisms, a rehash of Lauren's code. This gives the same results that she got so I'm pretty sure it works
rm(list = ls())
library(deSolve)
library(here)
library(foreach)
library(doSNOW)
# ----------------------------------------------------------------------------------------------------
# Run model with both species to get overall dynamics

VassFox_Cvar <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dP = -(M_P * P) + (((J_P * P) * ((O_P_C1 * C1) + ((1 - O_P_C1) * C2))) / ((O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0))
    dC1 = -(M_C1 * C1) + ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - ((O_P_C1 * J_P * P * C1) / ((O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0))
    dC2 = -(M_C2 * C2) + ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2)) - (((1 -
                                                                         O_P_C1) * J_P * P * C2) / ((O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0))
    dR = r * R * (1 - (R / K)) - ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2))
    return(list(c(dP, dC1, dC2, dR)))
    
  })
}

  # parameters
  # resource intrinsic rate of growth
  r <- 1.0
  # resource carrying capacity
  K <- 1.0
  # consumer 1 ingestion rate
  J_C1 <- 0.8036
  # consumer 2 ingestion rate
  J_C2 <- 0.7
  # predator ingestion rate
  J_P <- 0.4
  # predator mortality rate
  M_P <- 0.08
  # half saturation constant
  R_0_1 <- 0.16129
  R_0_2 <- 0.9
  C_0 <- 0.5
  # preference coefficient
  O_P_C1 <- 0.92
  O_C1_R <- 1.0
  O_C2_R <- 0.98
  time <- 5000
  
runtoequal <- function(State_int, M_C1_temp, M_C2_temp) {
  
  # Run to equilibrium
  mat <- matrix(data = NA,
                nrow = time,
                ncol = 4)
  mat[1, ] <- State_int
  
  for (t in 2:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    
    pars <-
      c(
        r = r,
        K = K,
        J_C1 = J_C1,
        J_C2 = J_C2,
        J_P = J_P,
        M_C1 = M_C1,
        M_C2 = M_C2,
        M_P = M_P,
        R_0_1 = R_0_1,
        R_0_2 = R_0_2,
        C_0 = C_0,
        O_P_C1 = O_P_C1,
        O_C1_R = O_C1_R,
        O_C2_R = O_C2_R
      )
    
    # Udate state variables to output from last timestep
    State <- c(
      P = mat[t - 1, 1],
      C1 = mat[t - 1, 2],
      C2 = mat[t - 1, 3],
      R = mat[t - 1, 4]
    )
    
    # run ODE solver
    VF <-
      as.data.frame(
        ode(
          func = VassFox_Cvar,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ),
        events = list(func = eventfun)
      )
    
    # Update results matrix
    mat[t, ] <- c(VF[2, 2], VF[2, 3], VF[2, 4], VF[2, 5])
    
    
  }
  return(mat)
  
}



runtoequal_set <- function(State_int, M_C1_temp, M_C2_temp, P) {
#sets pred to mean value
  mat <- matrix(data = NA,
                nrow = time,
                ncol = 4)
  mat[1, ] <- State_int
  
  for (t in 2:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    pars <-
      c(
        r = r,
        K = K,
        J_C1 = J_C1,
        J_C2 = J_C2,
        J_P = J_P,
        M_C1 = M_C1,
        M_C2 = M_C2,
        M_P = M_P,
        R_0_1 = R_0_1,
        R_0_2 = R_0_2,
        C_0 = C_0,
        O_P_C1 = O_P_C1,
        O_C1_R = O_C1_R,
        O_C2_R = O_C2_R
      )
    
    # Udate state variables to output from last timestep
    State <- c(
      P = P,
      C1 = mat[t - 1, 2],
      C2 = mat[t - 1, 3],
      R = mat[t - 1, 4]
    )
    
    # run ODE solver
    VF <-
      as.data.frame(
        ode(
          func = VassFox_Cvar,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ),
        events = list(func = eventfun)
      )
    
    # Update results matrix
    mat[t, ] <- c(P, VF[2, 3], VF[2, 4], VF[2, 5])
    
  }
  return(mat)
  
}


invade_C1 <- function(mat, M_C1_temp, M_C2_temp) {
#calculates what happens when C1 invades
  invade_start_time <- time / 2
  
  # now invade C1
  C1_ldgr <- matrix(
    data = NA,
    nrow = (time - invade_start_time),
    ncol = 1
  )
  C2_resident <-
    matrix(
      data = NA,
      nrow = (time - invade_start_time),
      ncol = 1
    )
  counter <- 1
  
  for (t in invade_start_time:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    
    pars <-
      c(
        r = r,
        K = K,
        J_C1 = J_C1,
        J_C2 = J_C2,
        J_P = J_P,
        M_C1 = M_C1,
        M_C2 = M_C2,
        M_P = M_P,
        R_0_1 = R_0_1,
        R_0_2 = R_0_2,
        C_0 = C_0,
        O_P_C1 = O_P_C1,
        O_C1_R = O_C1_R,
        O_C2_R = O_C2_R
      )
    
    # Udate state variables to output from last timestep
    State <- c(
      P = mat[t - 1, 1],
      C1 = 0.001,
      C2 = mat[t - 1, 3],
      R = mat[t - 1, 4]
    )
    
    # run ODE solver
    VF_out <-
      as.data.frame(
        ode(
          func = VassFox_Cvar,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ),
        events = list(func = eventfun)
      )
    C1_ldgr[counter] <- log (VF_out[2, 3] / .001)
    C2_resident[counter] <- log (VF_out[2, 4] / VF_out[1, 4])
    
    counter <- counter + 1
  }
  
  
  
  result <- data.frame(C1_ldgr = C1_ldgr, C2_resident = C2_resident)
  return(result)
  
}

invade_C2 <- function(mat, M_C1_temp, M_C2_temp) {
#calculate what happens when C2 invades
  invade_start_time <- time / 2
  
  # now invade C1
  C2_ldgr <- matrix(
    data = NA,
    nrow = (time - invade_start_time),
    ncol = 1
  )
  C1_resident <-
    matrix(
      data = NA,
      nrow = (time - invade_start_time),
      ncol = 1
    )
  counter <- 1
  
  for (t in invade_start_time:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    
    pars <-
      c(
        r = r,
        K = K,
        J_C1 = J_C1,
        J_C2 = J_C2,
        J_P = J_P,
        M_C1 = M_C1,
        M_C2 = M_C2,
        M_P = M_P,
        R_0_1 = R_0_1,
        R_0_2 = R_0_2,
        C_0 = C_0,
        O_P_C1 = O_P_C1,
        O_C1_R = O_C1_R,
        O_C2_R = O_C2_R
      )
    
    # run ODE solver
    State <- c(
      P = mat[t - 1, 1],
      C1 = mat[t - 1, 2],
      C2 = .001,
      R = mat[t - 1, 4]
    )
    
    # run ODE solver
    VF_out <-
      as.data.frame(
        ode(
          func = VassFox_Cvar,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ),
        events = list(func = eventfun)
      )
    C2_ldgr[counter] <- log (VF_out[2, 4] / .001)
    C1_resident[counter] <- log (VF_out[2, 3] / VF_out[1, 3])
    
    counter <- counter + 1
  }
  
  result <- data.frame(C2_ldgr = C2_ldgr, C1_resident = C1_resident)
  return(result)
  
}



runall <- function(sigma, rho, time) {
  #evaluate all combinations of invaders, wrapped up in a function so that it can be called in parallel
  # ----------------------------------------------------------------------------------------------------
  # starting conditions
  State <- c(P = 1,
             C1 = 1,
             C2 = 1,
             R = 1) # starting parameters
  
  # Calculate temporal variation in mortality timeseries
  z <- matrix(data = rnorm(
    n = 2 * time,
    mean = 0,
    sd = 1
  ), nrow = 2)
  
  cholesky <-
    matrix(data = c(sigma ^ 2, rho * (sigma ^ 2), rho * (sigma ^ 2), sigma ^
                      2),
           nrow = 2)
  
  g <- cholesky %*% z
  M_C1_temp <- 0.4 * exp(g[1, ])
  M_C2_temp <- 0.2 * exp(g[2, ])
  
  
  #results <- runtoequal(State, M_C1_temp, M_C2_temp)
  invade_start_time <- time / 2
  
  State <- c(P = 1,
             C1 = 0,
             C2 = 1,
             R = 1) # starting parameters
  
  C1_invade_C2_resident <- runtoequal(State, M_C1_temp, M_C2_temp)
  # now invade C1
  invadeC1 <- invade_C1(C1_invade_C2_resident, M_C1_temp, M_C2_temp)
  C1_ldgr <- invadeC1$C1_ldgr
  C2_resident <- invadeC1$C2_resident
  
  # C2 as the invader
  State <- c(P = 1,
             C1 = 1,
             C2 = 0,
             R = 1) # starting parameters
  
  C2_invade_C1_resident <- runtoequal(State, M_C1_temp, M_C2_temp)
  # now invade C2
  invadeC2 <- invade_C2(C2_invade_C1_resident, M_C1_temp, M_C2_temp)
  C2_ldgr <- invadeC2$C2_ldgr
  C1_resident <- invadeC2$C1_resident
  
  # calculate r_bar
  C1_r_bar <- mean(C1_ldgr) - mean(C2_resident)
  C2_r_bar <- mean(C2_ldgr) - mean(C1_resident)
  
  # ----------------------------------------------------------------------------------------------------
  # Partitioning coexistence mechanisms
  # set up non-fluctuating conditions
  avg_M_C1 <- mean(M_C1_temp[invade_start_time:time])
  avg_M_C2 <- mean(M_C2_temp[invade_start_time:time])
  
  avg_predator_C1_invade_C2_resident <-
    mean(C1_invade_C2_resident[invade_start_time:time, 1])
  avg_predator_C2_invade_C1_resident <-
    mean(C2_invade_C1_resident[invade_start_time:time, 1])
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_0
  
  # C1 as the invader
  State <-
    c(P = avg_predator_C1_invade_C2_resident,
      C1 = 0,
      C2 = 1,
      R = 1) # starting parameters
  
  
  epsilon_0_C1_invade_C2_resident <-
    runtoequal_set(
      State,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      avg_predator_C1_invade_C2_resident
    )
  
  # now invade C1
  invadeC1 <-
    invade_C1(epsilon_0_C1_invade_C2_resident,
              rep(avg_M_C1, length(M_C1_temp)),
              rep(avg_M_C2, length(M_C2_temp)))
  
  
  C1_epsilon_0 <- invadeC1$C1_ldgr
  C2_resident_epsilon_0 <- invadeC1$C2_resident
  
  
  # C2 as the invader
  State <-
    c(P = avg_predator_C2_invade_C1_resident,
      C1 = 1,
      C2 = 0,
      R = 1) # starting parameters
  
  epsilon_0_C2_invade_C1_resident <-
    runtoequal_set(
      State,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      avg_predator_C2_invade_C1_resident
    )
  
  # now invade C2
  invadeC2 <-
    invade_C2(epsilon_0_C2_invade_C1_resident,
              rep(avg_M_C1, length(M_C1_temp)),
              rep(avg_M_C2, length(M_C2_temp)))
  
  C2_epsilon_0 <- invadeC2$C2_ldgr
  C1_resident_epsilon_0 <- invadeC2$C1_resident
  
  C1_delta_0 <- mean(C1_epsilon_0) - mean(C2_resident_epsilon_0)
  C2_delta_0 <- mean(C2_epsilon_0) - mean(C1_resident_epsilon_0)
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_E
  
  # C1 as the invader
  State <-
    c(P = avg_predator_C1_invade_C2_resident,
      C1 = 0,
      C2 = 1,
      R = 1) # starting parameters
  
  
  epsilon_E_C1_invade_C2_resident <-
    runtoequal_set(State,
                   M_C1_temp,
                   M_C2_temp,
                   avg_predator_C1_invade_C2_resident)
  
  # now invade C1
  invadeC1 <-
    invade_C1(epsilon_E_C1_invade_C2_resident, M_C1_temp, M_C2_temp)
  
  
  C1_epsilon_E <- invadeC1$C1_ldgr
  
  C2_resident_epsilon_E <- invadeC1$C2_resident
  
  
  # C2 as the invader
  State <-
    c(P = avg_predator_C2_invade_C1_resident,
      C1 = 1,
      C2 = 0,
      R = 1) 
  
  epsilon_E_C2_invade_C1_resident <-
    runtoequal_set(State,
                   M_C1_temp,
                   M_C2_temp,
                   avg_predator_C2_invade_C1_resident)
  # now invade C2
  invadeC2 <-
    invade_C2(epsilon_E_C2_invade_C1_resident, M_C1_temp, M_C2_temp)
  C2_epsilon_E <- invadeC2$C2_ldgr
  
  C1_resident_epsilon_E <- invadeC2$C1_resident
  
  
  C1_delta_E <-
    mean(C1_epsilon_E) - mean(C2_resident_epsilon_E) - C1_delta_0
  C2_delta_E <-
    mean(C2_epsilon_E) - mean(C1_resident_epsilon_E) - C2_delta_0
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_P
  # predator population size varies, but mortality remains constant
  
  # C1 as the invader
  State <- c(P = 1,
             C1 = 0,
             C2 = 1,
             R = 1) 
  
  epsilon_P_C1_invade_C2_resident <-
    runtoequal(State, rep(avg_M_C1, length(M_C1_temp)), rep(avg_M_C2, length(M_C2_temp)))
  
  
  # now invade C1
  invadeC1 <-
    invade_C1(epsilon_P_C1_invade_C2_resident,
              rep(avg_M_C1, length(M_C1_temp)),
              rep(avg_M_C2, length(M_C2_temp)))
  
  C1_epsilon_P <- invadeC1$C1_ldgr
  
  C2_resident_epsilon_P <- invadeC1$C2_resident
  
  # C2 as the invader
  State <- c(P = 1,
             C1 = 1,
             C2 = 0,
             R = 1) # starting parameters
  
  epsilon_P_C2_invade_C1_resident <-
    runtoequal(State, rep(avg_M_C1, length(M_C1_temp)), rep(avg_M_C2, length(M_C2_temp)))
  
  # now invade C2
  invadeC2 <-
    invade_C2(epsilon_P_C2_invade_C1_resident ,
              rep(avg_M_C1, length(M_C1_temp)),
              rep(avg_M_C2, length(M_C2_temp)))
  

  C2_epsilon_P <- invadeC2$C2_ldgr
  C1_resident_epsilon_P <- invadeC2$C1_resident
  
  C1_delta_P <-
    mean(C1_epsilon_P) - mean(C2_resident_epsilon_P) - C1_delta_0
  C2_delta_P <-
    mean(C2_epsilon_P) - mean(C1_resident_epsilon_P) - C2_delta_0
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_EP
  C1_delta_EP <- C1_r_bar - (C1_delta_0 + C1_delta_P + C1_delta_E)
  C2_delta_EP <- C2_r_bar - (C2_delta_0 + C2_delta_P + C2_delta_E)
  
  # ----------------------------------------------------------------------------------------------------
  # Update final results vector
  C1_results <-
    c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)
  C2_results <-
    c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)
  
  done <- NULL
  #record results
  done$C1_final_mechanisms <- C1_results
  done$C2_final_mechanisms <- C2_results
  
  return(done)
  
}

# strength of env. on mortality rate
sigma <- 0.55

# cross-correlation of C1 and C2
rho <- .75

#set up for parallization, change the number based on how many cores your machine has. Most personal laptops have 4 or 8. If you don't know how many you have its probably 4. If you have a fancy cluster crank this up to like 100.  
cluster = makeCluster(8, type = "SOCK")
registerDoSNOW(cluster)
# looping over multiple runs
runs <- 100

#call the function that will do all of the calculations and store the results runs with do par
results <-
  foreach (
    run_loop = 1:runs,
    .packages = c("deSolve") ,
    .combine = 'rbind'
  ) %dopar% {
    d = runall(sigma, rho, time)
  }


#format things nicely and save
C1_final_mechanisms <- t(as.data.frame(results[, 1]))
colnames(C1_final_mechanisms) <-
  c("C1_r_bar",
    "C1_delta_0",
    "C1_delta_P",
    "C1_delta_E",
    "C1_delta_EP")

C2_final_mechanisms <-  t(as.data.frame(results[, 2]))
colnames(C2_final_mechanisms) <-
  c("C2_r_bar",
    "C2_delta_0",
    "C2_delta_P",
    "C2_delta_E",
    "C2_delta_EP")

write.csv(
  C1_final_mechanisms,
  file = here("C1_final_mechanisms_positive_cross.csv"),
  row.names = FALSE
)
write.csv(
  C2_final_mechanisms,
  file = here("C2_final_mechanisms_positive_cross.csv"),
  row.names = FALSE
)

#release the cores we used, this is important! Always remember to run this step if you have created a cluster, or else your cores won't be released until you restart the R session. 
stopCluster(cluster)

# ----------------------------------------------------------------------------------------------------
