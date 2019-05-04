# Implementation of the Vasseur and Fox Model with temporal variation in mortality rates
# Compare positive, no, and negative autocorrelation
# For tree species mortality rates through time (e.g. from Vasseur and Fox appendix)

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
#M_C1 = M_C1  # set later
# medial consumer 2 mortality rate
#M_C2 = M_C2  # set later
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

# cross-correlation of C1 and C2
rho=0

time  <- 5000 # number of timesteps to run the model 

# ----------------------------------------------------------------------------------------------------
# looping over parameter space

sigma_seq <- seq(0,.75,.05)  # set to 0.05
omega_seq <- seq(.5, 1, .05) # set to 0.05
runs <- 10
tracker <- 1

coexistence_space <- array(data=NA, c(length(omega_seq), length(sigma_seq), runs))

C1_space <- coexistence_space
C2_space <- coexistence_space

for (run_num in 1:runs) {
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
      
      # results <- matrix(data=NA, nrow=time, ncol=4)
      # results[1,] <- State
      # 
      # for (t in 2:time) {
      #   M_C1 <- M_C1_temp[t]
      #   M_C2 <- M_C2_temp[t]
      #   
      #   pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
      #   
      #   # Udate state variables to output from last timestep
      #   State <- c(P = results[t-1,1], C1 = results[t-1,2], C2 = results[t-1,3], R = results[t-1,4])
      #   
      #   # run ODE solver
      #   VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)),
      #                           events = list(func = eventfun))
      #   
      #   # Update results matrix
      #   results[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
      #   
      # }
      
      
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
      
      # update matrix of final results  
      coexistence_space[seq2, seq1, run_num] <- coexist
      C1_space[seq2, seq1, run_num] <- C1_r_bar
      C2_space[seq2, seq1, run_num] <- C2_r_bar
      
    }
  }

}

save(coexistence_space, C1_space, C2_space, sigma_seq, omega_seq, file="MSIresults.RData")
