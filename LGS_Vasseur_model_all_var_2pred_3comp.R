library(deSolve)
library(ggplot2)
library(reshape)
library(cowplot)
library(tidyverse)

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





VassFox_Cvar_2P_3C <-function(Time, State, Pars){
  
  with(as.list(c(State, Pars)), {
    #M_C1_j = M_C1_0 * exp()
    #M_C2_j = M_C2_0 * exp()
    
    
    #variables to make math easier to read
    
    
    all_c_p1 <- ( (O_P1_C1 * C1) +  (O_P1_C2 * C2) + (1-(O_P1_C1 + O_P1_C2)) * C3) 
    
    
    all_c_p2 <- ( (O_P2_C1 * C1) +  (O_P2_C2 * C2) + (1-(O_P2_C1 + O_P2_C2)) * C3) 
    
    both_preds_eat_C1 <- ( (O_P1_C1 * J_P1 * P1 * C1) / ( all_c_p1+all_c_p2+C_0 )) + ( (O_P2_C1 * J_P2 * P2 * C1) / ( all_c_p1+all_c_p2+C_0 ))
    
    both_preds_eat_C2<- ( ( O_P1_C2 * J_P1 * P1 * C2) / ( all_c_p1+all_c_p2+C_0 )) + ( (O_P2_C2 * J_P2 * P2 * C2) / ( all_c_p1+all_c_p2+C_0 ))
    
    
    both_preds_eat_C3<- ( ( (1-(O_P1_C1 + O_P1_C2  )) * J_P1 * P1 * C3) / ( all_c_p1+all_c_p2+C_0 )) + ( ((1-(O_P2_C1 + O_P2_C2  ))* J_P2 * P2 * C3) / ( all_c_p1+all_c_p2+C_0 ))
    
    
    
    dP_1 = - (M_P1 * P1) + ( ( (J_P1 * P1) * ( (O_P1_C1 * C1) + ( O_P1_C2 * C2) + (1-(O_P1_C1 + O_P1_C2  )) * C3  ) ) / (all_c_p1 +all_c_p2+C_0))
  
  
    dP_2 = - (M_P2 * P2) + ( ( (J_P2 * P2) * ( (O_P2_C1 * C1) + (O_P2_C2 * C2) + (1-(O_P2_C1 + O_P2_C2  )) * C3 ) ) / (all_c_p1 +all_c_p2+C_0))
    
    
    #dP = - (M_P * P) + ( ( (J_P * P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
    #dC1 = - (M_C1 * C1) + ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_P1_C1 * J_P1 * P1 * C1) / ( (O_P1_C1 * C1) + ((1 - O_P1_C1) * C2) + C_0) )

    
    
    
    dC1 = - (M_C1 * C1) + ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - both_preds_eat_C1
    dC2 = - (M_C2 * C2) + ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) ) - both_preds_eat_C2
    dC3 = - (M_C3 * C3) + ( (O_C3_R * J_C3 * C3 * R) / (R + R_0_3) ) - both_preds_eat_C3
    
    dR = r * R * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) ) -  ( (O_C3_R * J_C3 * C3 * R) / (R + R_0_3) )
    
    return(list(c(dP_1, dP_2, dC1, dC2,dC3, dR)))
    
  })


}


#enviro flux consumers
sigma <- 0.8
ro <- -0.2

z <- matrix(data = rnorm(n = 2*5000, mean = 0, sd = 1), nrow = 2)
cholesky <- matrix(data = c(sigma^2, ro * (sigma ^2), ro * (sigma ^2), sigma ^2), 
                   nrow = 2)
g <- cholesky %*% z
M_C1_temp <- 0.4 * exp(g[1,])
M_C2_temp <- 0.2 * exp(g[2,])

#redraw fluxs for C3, rather then make a 3x matrix, its easier to do this.
z <- matrix(data = rnorm(n = 2*5000, mean = 0, sd = 1), nrow = 2)
g <- cholesky %*% z

M_C3_temp <- 0.25 * exp(g[1,])



#enviro flux preds
pred_sigma <-.1
pred_ro <- -0.2

z_pred <- matrix(data = rnorm(n = 2*5000, mean = 0, sd = 1), nrow = 2)

pred_cholesky <- matrix(data = c(pred_sigma^2, pred_ro * (pred_sigma ^2), pred_ro * (pred_sigma ^2),pred_sigma^2 ),nrow=2)

g_pred <- pred_cholesky %*% z_pred

M_P1_temp <-0.09*exp(g_pred[1,])
M_P2_temp <-0.12*exp(g_pred[2,])



#enviro flux resource
res_sigma <-0.9
res_ro <- -0.1

z_res<-rnorm(n=5000,mean=0,sd=1)

res_cholesky<-matrix(data = c(res_sigma^2, res_ro*(res_sigma ^2)),ncol=1)

g_res <-res_cholesky %*% z_res

r_temp <- 1*exp(g_res[1,])


State <- c(P1 = 1, P2=1, C1 = 1, C2 = 1,C3 =1, R = 1)
Time <- seq(0, 5000, by = 1)

results <- matrix(data=NA, nrow=5000, ncol=6)
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
    r=r,
    # resource carrying capacity
    K = 1.0,
    # consumer 1 ingestion rate
    J_C1 = 0.8036,
    # consumer 2 ingestion rate
    J_C2 = 0.7036,
    # consumer 3 ingestion rate
    J_C3 = 0.923,
    # predator ingestion rate
    J_P1 = 0.85,
    J_P2 = 0.4,
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
    R_0_3 = 0.66,
    
    C_0 = 0.5,
    # preference coefficient
    O_P1_C1 = 0.45,
    O_P2_C1 = 0.1,
   
    O_P1_C2= 0.1,
    O_P2_C2 = 0.42,
    
    
    O_C1_R = .91,
    O_C2_R = 0.98,
    O_C3_R = 0.9,
    sigma=sigma,
    ro=ro
  )
  
  State <- c(P1 = results[t-1,1], P2 = results[t-1,2],C1 = results[t-1,3], C2 = results[t-1,4], C3 = results[t-1,5],R = results[t-1,6])
  VF_out <- as.data.frame(ode(func = VassFox_Cvar_2P_3C, y = State, parms = pars, times = seq(0,1)), 
                          events = list(func = eventfun))
  
  results[t,] <- c(VF_out[2,2],VF_out[2,3], VF_out[2,4], VF_out[2,5], VF_out[2,6],VF_out[2,7])
  
}

dat<- data.frame(pred1=results[,1],pred2=results[,2],con1=results[,3],con2=results[,4],con3=results[,5],res=results[,6])
d<-melt(dat)
d$time<-c(rep(seq(1,5000),6))


d<-d[d$time>2500,]
#ggplot(d, aes(x=time,y=value,color = variable))+geom_line()+facet_wrap( ~ variable, scales="free")


scale<-1
#quartz(width=5, height=5)
#plot(results[2500:5000,1]*scale, type="l", col="black")
#lines(results[2500:5000,2]*scale, col="red")
#lines(results[2500:5000,3]*scale, col="blue")
#lines(results[2500:5000,4]*scale, col="purple")
#lines(results[2500:5000,5]*scale, col="green")
#legend("topright", c("predator 1", "predator 2", "prey 1", "prey 2", "resource"), col=c("black", "red", "blue", "purple", "green"), lty=1, lwd=2)



ggplot(d, aes(x=time,y=value,color = variable))+geom_line()+ylim(0,1)+facet_wrap( ~ variable, scales="free")
