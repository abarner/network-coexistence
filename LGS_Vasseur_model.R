
VassFox_Cvar <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    #M_C1_j = M_C1_0 * exp()
    #M_C2_j = M_C2_0 * exp()
    
    dP = - (M_P * P) + ( ( (J_P * P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
    dC1 = - (M_C1 * C1) + ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P * C1) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dC2 = - (M_C2 * C2) + ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P * C2) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dR = r * R * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) )
    
    return(list(c(dP, dC1, dC2, dR)))
    
  })
}


sigma <- 0.55
ro <- -0.75

z <- matrix(data = rnorm(n = 2*1000, mean = 0, sd = 1), nrow = 2)
cholesky <- matrix(data = c(sigma^2, ro * (sigma ^2), ro * (sigma ^2), sigma ^2), 
                   nrow = 2)
g <- cholesky %*% z
M_C1_temp <- 0.4 * exp(g[1,])
M_C2_temp <- 0.2 * exp(g[2,])

State <- c(P = 1, C1 = 1, C2 = 1, R = 1)
Time <- seq(0, 1000, by = 1)

results <- matrix(data=NA, nrow=1000, ncol=4)
results[1,] <- State

for (t in 2:(1000)) {
  M_C1 <- M_C1_temp[t]
  M_C2 <- M_C2_temp[t]
 # from Table 1
pars <- c(
  # resource intrinsic rate of growth
  r = 1.0,
  # resource carrying capacity
  K = 1.0,
  # consumer 1 ingestion rate
  J_C1 = 0.8036,
  # consumer 2 ingestion rate
  J_C2 = 0.7,
  # predator ingestion rate
  J_P = 0.4,
  # medial consumer 1 mortality rate
  M_C1 = M_C1,
  # medial consumer 2 mortality rate
  M_C2 = M_C2,
  # predator mortality rate
  M_P = 0.08,
  # half saturation constant
  R_0_1 = 0.16129,
  R_0_2 = 0.9,
  C_0 = 0.5,
  # preference coefficient
  O_P_C1 = 0.92,
  O_C1_R = 1.0,
  O_C2_R = 0.98,
  sigma=sigma,
  ro=ro
)

State <- c(P = results[t-1,1], C1 = results[t-1,2], C2 = results[t-1,3], R = results[t-1,4])
VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                        events = list(func = eventfun))

results[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])

}

quartz(width=5, height=5)
plot(results[,1], type="l", col="black")
lines(results[,2], col="blue")
lines(results[,3], col="purple")
lines(results[,4], col="green")
legend("topright", c("predator", "prey 1", "prey 2", "resource"), col=c("black", "blue", "purple", "green"), lty=1, lwd=2)
