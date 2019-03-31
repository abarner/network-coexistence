# Implementation of the Vasseur and Fox Model with temporal variation in mortality rates
# Compare positive, no, and negative autocorrelation
# For tree species mortality rates through time (e.g. from Vasseur and Fox appendix)

library(here)
library(deSolve)
# ----------------------------------------------------------------------------------------------------
# Run model with both species to get overall dynamics

VassFox_Cvar <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dP = - (M_P * P) + ( ( (J_P * P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
    dC1 = - (M_C1 * C1) + ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P * C1) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dC2 = - (M_C2 * C2) + ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P * C2) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dR = r * R * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) )
    
    return(list(c(dP, dC1, dC2, dR)))
    
  })
}

# ----------------------------------------------------------------------------------------------------
# parameters
# resource intrinsic rate of growth
r = 1.0
# resource carrying capacity
K = 1.0
# consumer 1 ingestion rate
J_C1 = 0.8036
# consumer 2 ingestion rate
J_C2 = 0.7
# predator ingestion rate
J_P = 0.4
# medial consumer 1 mortality rate
M_C1 = M_C1
# medial consumer 2 mortality rate
M_C2 = M_C2
# predator mortality rate
M_P = 0.08
# half saturation constant
R_0_1 = 0.16129
R_0_2 = 0.9
C_0 = 0.5
# preference coefficient
#O_P_C1 = 0.92
O_C1_R = 1.0
O_C2_R = 0.98

# strength of env. on mortality rate
#sigma=0.55

# cross-correlation of C1 and C2
rho=0

time  <- 5000 # number of timesteps to run the model 

# ----------------------------------------------------------------------------------------------------
# looping over parameter space

sigma_seq <- seq(0,.75,.05)
omega_seq <- seq(.5, 1, .05)
tracker <- 1

coexistence_space <- matrix(data=NA, nrow=length(omega_seq), ncol=length(sigma_seq))
rownames(coexistence_space) <- omega_seq
colnames(coexistence_space) <- sigma_seq
C1_space <- coexistence_space
C2_space <- coexistence_space


for (seq1 in 1:length(sigma_seq)) {
  for (seq2 in 1:length(omega_seq)) {
    
    sigma <- sigma_seq[seq1]
    O_P_C1 = omega_seq[seq2]
    # ----------------------------------------------------------------------------------------------------
    # starting conditions
    State <- c(P = 1, C1 = 1, C2 = 1, R = 1) # starting parameters
    
    # Calculate temporal variation in mortality timeseries
    z <- matrix(data = rnorm(n = 2*time, mean = 0, sd = 1), nrow = 2)
    cholesky <- matrix(data = c(sigma^2, rho * (sigma ^2), rho * (sigma ^2), sigma ^2), 
                       nrow = 2)
    g <- cholesky %*% z
    M_C1_temp <- 0.4 * exp(g[1,])
    M_C2_temp <- 0.2 * exp(g[2,])
    
    results <- matrix(data=NA, nrow=time, ncol=4)
    results[1,] <- State
    
    for (t in 2:time) {
      M_C1 <- M_C1_temp[t]
      M_C2 <- M_C2_temp[t]

      pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)

      # Udate state variables to output from last timestep
      State <- c(P = results[t-1,1], C1 = results[t-1,2], C2 = results[t-1,3], R = results[t-1,4])

      # run ODE solver
      VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)),
                              events = list(func = eventfun))

      # Update results matrix
      results[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])

    }
    

    # ----------------------------------------------------------------------------------------------------
    # Low density growth rate calculations
    
    # C1 as the invader
    State <- c(P = 1, C1 = 0, C2 = 1, R = 1) # starting parameters
    
    # invasion time
    invade_start_time <- time/2
    
    # Run to equilibrium
    C1_invade_C2_resident <- matrix(data=NA, nrow=time, ncol=4)
    C1_invade_C2_resident[1,] <- State
    
    for (t in 2:time) {
      M_C1 <- M_C1_temp[t]
      M_C2 <- M_C2_temp[t]
      
      pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
      
      # Udate state variables to output from last timestep
      State <- c(P = C1_invade_C2_resident[t-1,1], C1 = C1_invade_C2_resident[t-1,2], 
                 C2 = C1_invade_C2_resident[t-1,3], R = C1_invade_C2_resident[t-1,4])
      
      # run ODE solver
      VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                              events = list(func = eventfun))
      
      # Update results matrix
      C1_invade_C2_resident[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
      
    }
    
    # now invade C1
    C1_ldgr <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    C2_resident <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    counter <- 1
    
    for (t in invade_start_time:time) {
      
      M_C1 <- M_C1_temp[t]
      M_C2 <- M_C2_temp[t]
      
      pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
      
      # Udate state variables to output from last timestep
      State <- c(P = C1_invade_C2_resident[t-1,1], C1 = 0.001, 
                 C2 = C1_invade_C2_resident[t-1,3], R = C1_invade_C2_resident[t-1,4])
      
      # run ODE solver
      VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                              events = list(func = eventfun))
      C1_ldgr[counter] <- log (VF_out[2,3]/.001)
      C2_resident[counter] <- log (VF_out[2,4]/VF_out[1,4])
      
      counter <- counter + 1
    }
    
    # C2 as the invader
    State <- c(P = 1, C1 = 1, C2 = 0, R = 1) # starting parameters
    
    # invasion time
    invade_start_time <- time/2
    
    # Run to equilibrium
    C2_invade_C1_resident <- matrix(data=NA, nrow=time, ncol=4)
    C2_invade_C1_resident[1,] <- State
    
    for (t in 2:time) {
      M_C1 <- M_C1_temp[t]
      M_C2 <- M_C2_temp[t]
      
      pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
      
      # Udate state variables to output from last timestep
      State <- c(P = C2_invade_C1_resident[t-1,1], C1 = C2_invade_C1_resident[t-1,2], 
                 C2 = C2_invade_C1_resident[t-1,3], R = C2_invade_C1_resident[t-1,4])
      
      # run ODE solver
      VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                              events = list(func = eventfun))
      
      # Update results matrix
      C2_invade_C1_resident[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
      
    }
    
    # now invade C2
    C2_ldgr <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    C1_resident <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    counter <- 1
    
    for (t in invade_start_time:time) {
      
      M_C1 <- M_C1_temp[t]
      M_C2 <- M_C2_temp[t]
      
      pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
      
      # Udate state variables to output from last timestep
      State <- c(P = C2_invade_C1_resident[t-1,1], C1 = C2_invade_C1_resident[t-1,2], 
                 C2 = .001, R = C2_invade_C1_resident[t-1,4])
      
      # run ODE solver
      VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                              events = list(func = eventfun))
      C2_ldgr[counter] <- log (VF_out[2,4]/.001)
      C1_resident[counter] <- log (VF_out[2,3]/VF_out[1,3])
      
      counter <- counter + 1
    }
    
    # calculate r_bar
    C1_r_bar <- mean(C1_ldgr)-mean(C2_resident)
    C2_r_bar <- mean(C2_ldgr)-mean(C1_resident)
    
    # if (C1_r_bar > 0 & C2_r_bar > 0) {
    #   # ----------------------------------------------------------------------------------------------------
    #   # Partitioning coexistence mechanisms
    #   # set up non-fluctuating conditions
    #   avg_M_C1 <- mean(M_C1_temp[invade_start_time:time])
    #   avg_M_C2 <- mean(M_C2_temp[invade_start_time:time])
    #   
    #   avg_predator_C1_invade_C2_resident <- mean(C1_invade_C2_resident[invade_start_time:time,1])
    #   avg_predator_C2_invade_C1_resident <- mean(C2_invade_C1_resident[invade_start_time:time,1])
    #   
    #   # ----------------------------------------------------------------------------------------------------
    #   # calculate delta_0 
    #   
    #   # C1 as the invader
    #   State <- c(P = avg_predator_C1_invade_C2_resident, C1 = 0, C2 = 1, R = 1) # starting parameters
    #   
    #   # Run to equilibrium
    #   epsilon_0_C1_invade_C2_resident <- matrix(data=NA, nrow=time, ncol=4)
    #   epsilon_0_C1_invade_C2_resident[1,] <- State
    #   
    #   M_C1 <- avg_M_C1 
    #   M_C2 <- avg_M_C2 
    #   pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #   
    #   for (t in 2:time) {
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = avg_predator_C1_invade_C2_resident, C1 = epsilon_0_C1_invade_C2_resident[t-1,2], 
    #                C2 = epsilon_0_C1_invade_C2_resident[t-1,3], R = epsilon_0_C1_invade_C2_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     
    #     # Update results matrix
    #     epsilon_0_C1_invade_C2_resident[t,] <- c(avg_predator_C1_invade_C2_resident, VF_out[2,3], VF_out[2,4], VF_out[2,5])
    #     
    #   }
    #   
    #   
    #   # now invade C1
    #   C1_epsilon_0 <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   C2_resident_epsilon_0 <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   counter <- 1
    #   
    #   M_C1 <- avg_M_C1 
    #   M_C2 <- avg_M_C2 
    #   pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #   
    #   for (t in invade_start_time:time) {
    #     
    #     pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = epsilon_0_C1_invade_C2_resident[t-1,1], C1 = 0.001, 
    #                C2 = epsilon_0_C1_invade_C2_resident[t-1,3], R = epsilon_0_C1_invade_C2_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     C1_epsilon_0[counter] <- log (VF_out[2,3]/.001)
    #     C2_resident_epsilon_0[counter] <- log (VF_out[2,4]/VF_out[1,4])
    #     
    #     counter <- counter + 1
    #   }
    #   
    #   # C2 as the invader
    #   State <- c(P = avg_predator_C2_invade_C1_resident, C1 = 1, C2 = 0, R = 1) # starting parameters
    #   
    #   # Run to equilibrium
    #   epsilon_0_C2_invade_C1_resident <- matrix(data=NA, nrow=time, ncol=4)
    #   epsilon_0_C2_invade_C1_resident[1,] <- State
    #   
    #   M_C1 <- avg_M_C1 
    #   M_C2 <- avg_M_C2 
    #   pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #   
    #   for (t in 2:time) {
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = avg_predator_C2_invade_C1_resident, C1 = epsilon_0_C2_invade_C1_resident[t-1,2], 
    #                C2 = epsilon_0_C2_invade_C1_resident[t-1,3], R = epsilon_0_C2_invade_C1_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     
    #     # Update results matrix
    #     epsilon_0_C2_invade_C1_resident[t,] <- c(avg_predator_C2_invade_C1_resident, VF_out[2,3], VF_out[2,4], VF_out[2,5])
    #     
    #   }
    #   
    #   
    #   # now invade C2
    #   C2_epsilon_0 <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   C1_resident_epsilon_0 <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   counter <- 1
    #   
    #   M_C1 <- avg_M_C1 
    #   M_C2 <- avg_M_C2 
    #   pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #   
    #   for (t in invade_start_time:time) {
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = epsilon_0_C2_invade_C1_resident[t-1,1], C1 = epsilon_0_C2_invade_C1_resident[t-1,2], 
    #                C2 = 0.001, R = epsilon_0_C2_invade_C1_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     C2_epsilon_0[counter] <- log (VF_out[2,4]/.001)
    #     C1_resident_epsilon_0[counter] <- log (VF_out[2,3]/VF_out[1,3])
    #     
    #     counter <- counter + 1
    #   }
    #   
    #   C1_delta_0 <- mean(C1_epsilon_0)-mean(C2_resident_epsilon_0)
    #   C2_delta_0 <- mean(C2_epsilon_0)-mean(C1_resident_epsilon_0)
    #   
    #   # ----------------------------------------------------------------------------------------------------
    #   # calculate delta_E 
    #   
    #   # C1 as the invader
    #   State <- c(P = avg_predator_C1_invade_C2_resident, C1 = 0, C2 = 1, R = 1) # starting parameters
    #   
    #   # Run to equilibrium
    #   epsilon_E_C1_invade_C2_resident <- matrix(data=NA, nrow=time, ncol=4)
    #   epsilon_E_C1_invade_C2_resident[1,] <- State
    #   
    #   for (t in 2:time) {
    #     M_C1 <- M_C1_temp[t]
    #     M_C2 <- M_C2_temp[t]
    #     
    #     pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = avg_predator_C1_invade_C2_resident, C1 = epsilon_E_C1_invade_C2_resident[t-1,2], 
    #                C2 = epsilon_E_C1_invade_C2_resident[t-1,3], R = epsilon_E_C1_invade_C2_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     
    #     # Update results matrix
    #     epsilon_E_C1_invade_C2_resident[t,] <- c(avg_predator_C1_invade_C2_resident, VF_out[2,3], VF_out[2,4], VF_out[2,5])
    #     
    #   }
    #   
    #   
    #   # now invade C1
    #   C1_epsilon_E <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   C2_resident_epsilon_E <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   counter <- 1
    #   
    #   for (t in invade_start_time:time) {
    #     
    #     M_C1 <- M_C1_temp[t]
    #     M_C2 <- M_C2_temp[t]
    #     
    #     pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = epsilon_E_C1_invade_C2_resident[t-1,1], C1 = 0.001, 
    #                C2 = epsilon_E_C1_invade_C2_resident[t-1,3], R = epsilon_E_C1_invade_C2_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     C1_epsilon_E[counter] <- log (VF_out[2,3]/.001)
    #     C2_resident_epsilon_E[counter] <- log (VF_out[2,4]/VF_out[1,4])
    #     
    #     counter <- counter + 1
    #   }
    #   
    #   # C2 as the invader
    #   State <- c(P = avg_predator_C2_invade_C1_resident, C1 = 1, C2 = 0, R = 1) # starting parameters
    #   
    #   # Run to equilibrium
    #   epsilon_E_C2_invade_C1_resident <- matrix(data=NA, nrow=time, ncol=4)
    #   epsilon_E_C2_invade_C1_resident[1,] <- State
    #   
    #   for (t in 2:time) {
    #     
    #     M_C1 <- M_C1_temp[t]
    #     M_C2 <- M_C2_temp[t]
    #     
    #     pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = avg_predator_C2_invade_C1_resident, C1 = epsilon_E_C2_invade_C1_resident[t-1,2], 
    #                C2 = epsilon_E_C2_invade_C1_resident[t-1,3], R = epsilon_E_C2_invade_C1_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     
    #     # Update results matrix
    #     epsilon_E_C2_invade_C1_resident[t,] <- c(avg_predator_C2_invade_C1_resident, VF_out[2,3], VF_out[2,4], VF_out[2,5])
    #     
    #   }
    #   
    #   
    #   # now invade C2
    #   C2_epsilon_E <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   C1_resident_epsilon_E <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   counter <- 1
    #   
    #   for (t in invade_start_time:time) {
    #     
    #     M_C1 <- M_C1_temp[t]
    #     M_C2 <- M_C2_temp[t]
    #     
    #     pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #     
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = epsilon_E_C2_invade_C1_resident[t-1,1], C1 = epsilon_E_C2_invade_C1_resident[t-1,2], 
    #                C2 = 0.001, R = epsilon_E_C2_invade_C1_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     C2_epsilon_E[counter] <- log (VF_out[2,4]/.001)
    #     C1_resident_epsilon_E[counter] <- log (VF_out[2,3]/VF_out[1,3])
    #     
    #     counter <- counter + 1
    #   }
    #   
    #   C1_delta_E <- mean(C1_epsilon_E)-mean(C2_resident_epsilon_E) - C1_delta_0
    #   C2_delta_E <- mean(C2_epsilon_E)-mean(C1_resident_epsilon_E) - C2_delta_0
    #   
    #   # ----------------------------------------------------------------------------------------------------
    #   # calculate delta_P 
    #   # predator population size varies, but mortality remains constant
    #   
    #   # C1 as the invader
    #   State <- c(P = 1, C1 = 0, C2 = 1, R = 1) # starting parameters
    #   
    #   # Run to equilibrium
    #   epsilon_P_C1_invade_C2_resident <- matrix(data=NA, nrow=time, ncol=4)
    #   epsilon_P_C1_invade_C2_resident[1,] <- State
    #   
    #   M_C1 <- avg_M_C1 
    #   M_C2 <- avg_M_C2 
    #   pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #   
    #   for (t in 2:time) {
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = epsilon_P_C1_invade_C2_resident[t-1,1], C1 = epsilon_P_C1_invade_C2_resident[t-1,2], 
    #                C2 = epsilon_P_C1_invade_C2_resident[t-1,3], R = epsilon_P_C1_invade_C2_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     
    #     # Update results matrix
    #     epsilon_P_C1_invade_C2_resident[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
    #     
    #   }
    #   
    #   
    #   # now invade C1
    #   C1_epsilon_P <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   C2_resident_epsilon_P <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   counter <- 1
    #   
    #   M_C1 <- avg_M_C1 
    #   M_C2 <- avg_M_C2 
    #   pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #   
    #   for (t in invade_start_time:time) {
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = epsilon_P_C1_invade_C2_resident[t-1,1], C1 = 0.001, 
    #                C2 = epsilon_P_C1_invade_C2_resident[t-1,3], R = epsilon_P_C1_invade_C2_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     C1_epsilon_P[counter] <- log (VF_out[2,3]/.001)
    #     C2_resident_epsilon_P[counter] <- log (VF_out[2,4]/VF_out[1,4])
    #     
    #     counter <- counter + 1
    #   }
    #   
    #   # C2 as the invader
    #   State <- c(P = 1, C1 = 1, C2 = 0, R = 1) # starting parameters
    #   
    #   # Run to equilibrium
    #   epsilon_P_C2_invade_C1_resident <- matrix(data=NA, nrow=time, ncol=4)
    #   epsilon_P_C2_invade_C1_resident[1,] <- State
    #   
    #   M_C1 <- avg_M_C1 
    #   M_C2 <- avg_M_C2 
    #   pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #   
    #   for (t in 2:time) {
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = epsilon_P_C2_invade_C1_resident[t-1,1], C1 = epsilon_P_C2_invade_C1_resident[t-1,2], 
    #                C2 = epsilon_P_C2_invade_C1_resident[t-1,3], R = epsilon_P_C2_invade_C1_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     
    #     # Update results matrix
    #     epsilon_P_C2_invade_C1_resident[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
    #     
    #   }
    #   
    #   
    #   # now invade C2
    #   C2_epsilon_P <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   C1_resident_epsilon_P <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
    #   counter <- 1
    #   
    #   M_C1 <- avg_M_C1 
    #   M_C2 <- avg_M_C2 
    #   pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    #   
    #   for (t in invade_start_time:time) {
    #     
    #     # Udate state variables to output from last timestep
    #     State <- c(P = epsilon_P_C2_invade_C1_resident[t-1,1], C1 = epsilon_P_C2_invade_C1_resident[t-1,2], 
    #                C2 = 0.001, R = epsilon_P_C2_invade_C1_resident[t-1,4])
    #     
    #     # run ODE solver
    #     VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
    #                             events = list(func = eventfun))
    #     C2_epsilon_P[counter] <- log (VF_out[2,4]/.001)
    #     C1_resident_epsilon_P[counter] <- log (VF_out[2,3]/VF_out[1,3])
    #     
    #     counter <- counter + 1
    #   }
    #   
    #   C1_delta_P <- mean(C1_epsilon_P)-mean(C2_resident_epsilon_P) - C1_delta_0
    #   C2_delta_P <- mean(C2_epsilon_P)-mean(C1_resident_epsilon_P) - C2_delta_0
    #   
    #   # ----------------------------------------------------------------------------------------------------
    #   # calculate delta_EP 
    #   C1_delta_EP <- C1_r_bar - (C1_delta_0 + C1_delta_P + C1_delta_E)
    #   C2_delta_EP <- C2_r_bar - (C2_delta_0 + C2_delta_P + C2_delta_E)
    #   
    # } else{
    #   C1_delta_0 <- NA
    #   C1_delta_P <- NA 
    #   C1_delta_E <- NA 
    #   C1_delta_EP <- NA
    #   C2_delta_0 <- NA 
    #   C2_delta_P <- NA 
    #   C2_delta_E <- NA 
    #   C2_delta_EP <- NA
    # }
    # 
    # ----------------------------------------------------------------------------------------------------
    # Update final results vector
    coexist <- NA
    if (C1_r_bar > 0 & C2_r_bar > 0)
    {
      coexist <- 3
    } else if (C1_r_bar > 0 & C2_r_bar <= 0) {
      coexist <- 1
    } else if (C1_r_bar <= 0 & C2_r_bar > 0) {
      coexist <- 2
    } else {
      coexist <- 1
    }
    # 
    # C1_results <- c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)
    # C2_results <- c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)
    # 

    # update matrix of final results  
    coexistence_space[seq2, seq1] <- coexist
    C1_space[seq2, seq1] <- C1_r_bar
    C2_space[seq2, seq1] <- C2_r_bar

    
    tracker <- tracker + 1
    # print how far we've come
    if(tracker%%25==0) {
      print(tracker)
    }
    
  }
}


write.csv(coexistence_space, file = here("foodwebs", "cr_variation_examplecode_files", 
                                           "mechanism_figs", "parameter_sweep.csv"), row.names=FALSE)
write.csv(C1_space, file = here("foodwebs", "cr_variation_examplecode_files", 
                                         "mechanism_figs", "C1_parameter_sweet.csv"), row.names=FALSE)
write.csv(C2_space, file = here("foodwebs", "cr_variation_examplecode_files", 
                                "mechanism_figs", "C2_parameter_sweet.csv"), row.names=FALSE)

heatmap(coexistence_space)
# ----------------------------------------------------------------------------------------------------
# to plot run plotting_mechanisms.R
library(plot3D)

quartz(width=6, height=6)
par(mar=c(2,2,2,2), oma=c(2,2,.5, .5), mfrow=c(2,2))
image2D(coexistence_space, x=omega_seq, y=sigma_seq, border="black", 
        col=c("white", "grey", "black"), yaxt="n", xaxt="n", colkey=FALSE, xlab="", ylab="",
        main= "Coexistence")
axis(side=1, tick=seq(.5, 1, .1))
mtext(side=2, line=2.5, "Environmental Variation")
axis(side=2, at=seq(0, .75, .25), tick=TRUE)

plot.new()
legend(x=0.25, y=0.75, c("Species 1 Only", "Species 2 Only", "Coexistence"), 
       col=c("white", "grey", "black"), pch=15, xpd=TRUE, bty="n", pt.cex=2)

legend(x=0.25, y=0.75, c("", "", ""), 
       col="black", pch=22, xpd=TRUE, bty="n", pt.cex=2)

library(RColorBrewer)
#display.brewer.all(colorblindFriendly = T)

#quartz(width=5, height = 5)
#par(mfrow=c(2,1), oma=c(1,1,1,1), mar=c(2,2,2,2))
min_val <- -.221
max_val <- .161

col_pos=brewer.pal(9,"YlOrRd")
col_neg=brewer.pal(9,"Blues")

col1_pos_temp=round(9*(max(C1_space)/max_val))
col2_pos_temp=round(9*(max(C2_space)/max_val))
col1_neg_temp=round(9*(min(C1_space)/min_val))
col2_neg_temp=round(9*(min(C2_space)/min_val))

cols1_neg <- col_neg[1:col1_neg_temp]
cols1_neg <- rev(cols1_neg)
cols2_neg <- col_neg[1:col2_neg_temp]
cols2_neg <- rev(cols2_neg)

cols1_pos<- col_pos[1:col1_pos_temp]
cols2_pos<- col_pos[1:col2_pos_temp]

cols1 <- c(cols1_neg, cols1_pos)
cols2 <- c(cols2_neg, cols2_pos)

image2D(C1_space, x=omega_seq, y=sigma_seq, border="black", col=cols1,
         yaxt="n", xaxt="n", xlab="", ylab="", contour=FALSE, main="Species 1")
axis(side=1, at=seq(.5, 1, .25), tick=TRUE)
axis(side=2, at=seq(0, .75, .25), tick=TRUE)
mtext(side=1, line=2.5, "Predation Preference")
mtext(side=2, line=2.5, "Environmental Variation")

image2D(C2_space, x=omega_seq, y=sigma_seq, border="black", col=cols2,
        yaxt="n", xaxt="n", xlab="", ylab="", contour=FALSE, main="Species 2", colkey=FALSE)
axis(side=1, tick=seq(.5, 1, .25))
axis(side=2, tick=seq(0, .75, .25))
mtext(side=1, line=2.5, "Predation Preference")
