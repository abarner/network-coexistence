## updated version of the forde & doak 2004 model
## functions only

library(tidyverse)

#### functions internal to simulation ####

do.free.space.calculation <- function(total, B, size.B, C, size.C, L, size.L) {
  # Same as F&D model
  
  # total = total area
  # B is the population size of Balanus glandula
  # size.B = average size of adult B. glandula
  # C is population size of Chthamalus fissus/dalli
  # size.C = average size of C. fissus/dalli
  # L is population size of limpets
  # size.L = average size of an adult limpet
  free <- total - (B*size.B + C*size.C + L*size.L)
  return(free)
}

# for barnacles and limpets
do.potential.recruitment <- function(free, size.x, size.recruit.x, larvae.x) {
  # Same as F&D model
  
  # F is the amount of free space on the rock
  # size.x is the average size (e.g. size.B, size.C, or size.L)
  # size.recruit.x is the size of a recruit
  # larvae.x is the number of settling larvae of species x
  L <- (free/size.x)*(1-exp(-(size.recruit.x*larvae.x)/free))
  return(L)
}

do.actual.recruitment <- function(free, L, C) {
  # Updated from F&D: now based on iwasa & roughgarden equation 1
  
  # C = c_i of iwasa & roughgarden, larval settlement rate/recruit survival of target species
  # free = amount of free space in the system
  # L = potential number of recruits of target species coming from the larval pool
  return(C * L * free)
}

do.population.size.barnacles <- function(S, p.whelk, W.prev, X.prev, S.r, R, p.star, P.prev) {
  # Same as F&D model
  
  # S is survivorship (before whelk predation)
  # p.whelk is the encounter rate of barnacles by whelks
  # W.prev is the whelk population size in the previous month
  # X.prev is the barnacle population size in the previous month
  # S.r is the survivorship of recruits
  # R is the number of recruits
  # p.star is the encounter rate of barnacles by sea stars
  # P.prev is the whelk population size in the previous month
  
  # rewritten from previous form to make more intuitive sense
  X_before_predation <- S*X.prev
  
  if (X_before_predation < 0) { X_before_predation <- 0 } 
  
  X <- X_before_predation + S.r*R - p.whelk*W.prev*X_before_predation - p.star*P.prev*X_before_predation
  
  # for some reason, X might return as NA - add in function to see when
  if (is.na(X)) {
    print(c(X_before_predation, X.prev, R, W.prev, P.prev))
  }
  
  if (X < 0) {
    X <- 0
    return(X)
  } else {
    return (X)
    }
  
}

do.population.size.limpets <- function(S, L.prev, S.r, R, delta) {
  # Same as F&D model
  
  # S is survivorship
  # L.prev is the limpet population size in the previous month
  # S.r is the survivorship of recruits
  # R is the number of recruits
  # delta is density-dependence 
  L <- S*L.prev + S.r*R*exp(delta*L.prev)
  return(L)
}


# one of the main problems that causes the system to be "unbalanced" 
# (e.g., major barnacle mortality) is a lack of logistic growth in the
# whelk population. the time steps are too large to result in logistic growth
# and whelk recruitment is >> 90 (which is the set max population size) 
# in f & d. however, if we set maximum whelk population size at a flat number,
# that results in a constant removal of barnacle population size, which 
# doesn't quite represent predator dynamics. instead, 
do.population.size.whelks <- function(W.prev, 
                                      #W.feeding, 
                                      S, 
                                      R) {
  # Updated from F&D: Now follows simple similar form as other predator equation (sea star)
  
  # W.prev is whelk population size in the previous month
  # p is per capita predation rate
  # B.prev is the B. glandula population
  # C.prev is the C. fissus/dalli population
  # R is the recruitment rate
  # S is the survival rate
  # W.feeding is food intake
  
  W <- (W.prev * S) + R
  
  # note that the functional form of this population model is unclear
  # from the manuscript. Alternate model forms could be:
  # W <- W.prev * W. feeding + R
  # W <- W.prev*W.feeding*S + R
  
  return(W)
}

do.whelk.recruitment <- function(avg.C, avg.B, p.B, p.C, Y, W.prev, 
                                 B.prev, C.prev, S, 
                                 r = .3) {
  # Updated from F&D to explicitly incorporate density dependence
  # (Which was previously set at a hard limit of 90)
  
  # avg.C is the average number of C. fissus/dalli from April through June
  # avg. B is the average number of B. glandula from April through June
  # p is the per capita predation rate
  
  R <- (avg.C + avg.B)*3*Y*W.prev*S*((p.B*B.prev) + (p.C*C.prev))
  
  R_logistic <- (R / 90) * exp(r * (1 - (R / 90)))
  
  return(R*R_logistic)
}

do.population.size.seastar <- function(S, P.prev, R, survival.recruit.P, r = 1) {
  # Updated from F&D to explicitly incorporate density dependence
  # (Which was previously set at a hard limit of 6)
  
  # S is seastar adult survival
  # P.prev is previous seastar population size
  # R is abundance of recruits
  # survival.recruit.P is seastar recruit survival
  # Prey.prev is total population size of prey at last time point
  # r is population growth rate 
  
  # note that there is a typo in equation 7 in Forde & Doak appendix 
  # where "P.prev" (= P at time t-1) is incorrectly included as 
  # P at time t+1
  
  # add equation for seastar recruits?
  # sea stars grow very slowly, and should only enter into the "adult"
  # population once a year                                                                                                           
  
  P <- S*P.prev + survival.recruit.P*R
  
  P_logistic <- (P / 6) * exp(r * (1 - (P / 6)))
  
  return(P_logistic * P)
}

#### code to run simulation ####

do.intertidal.simulation <- function(
  
  years_set = 50, # number of years to run simulation
  
  # in connolly & roughgarden 1999, table 1, settlement is 0.02/hour/m2 for barn
  # estimates for limpet settlement is set to same as to barnacle settlement
  # don't have an estimate otherwise
  # use Fig 5a in Gilman 2006 Ecography (5-10 recruits max / 625 cm2)
  # vs Figures from Menge et al. 2011 JEMBE (50-100 recruits / 100 cm2)
  settlement.B = .002 * 30 * 24,
  settlement.C = .002 * 30 * 24,
  settlement.L = .002 * 30 * 24,
  
  # all means set to "high" recruitment scenarios 
  # (from forde & doak table 1)
  B.mean, # = 90000,
  B.stdev, # = 67823,
  C.mean, # = 70000,
  C.stdev, # = 52440,
  L.mean, # = 3000,
  L.stdev, # = 1949,
  P.mean, # = 6873,
  P.stdev, # = 5495,
  
  var_P = NULL, # if want constant value for pisaster recruitment, set value here
  var_B = NULL, # for balanus recruitment
  var_C = NULL, # for chthamalus recruitment
  var_L = NULL, # for limpet recruitment
  
  P_avg = NULL, # if want constant value for predator abundance, set value here
  W_avg = NULL,
  
  B_1 = 4100, # defaults to starting conditions given by forde & doak
  C_1 = 4100, # default was = 11000, but set to same as balanus
  L_1 = 239,
  #L_1 = 400, # after gilman 2006
  W_1 = 93,
  P_1 = 1,
  total_1 = 1
  
) {
  
  # ----------------------------------------------------------------------------------------------------
  # Model parameters
  
  size.B <- .000098
  size.C <- .000032
  size.L <- .00008
  # size.L <- .0001 # after gilman 2006 ecography
  
  size.recruit.B <- .000003
  size.recruit.C <- .000003
  size.recruit.L <- .000003
  
  survival.B <- .7 # asymmetry based on connell 1961 .75
  survival.C <- .7
  survival.L <- .97
  survival.W <- .94
  survival.P <- .992
  
  # connell 1961 suggests that survival of balanus > chthamalus
  survival.recruit.B <- .7
  survival.recruit.C <- .7
  survival.recruit.L <- .88
  survival.recruit.W <- .88
  survival.recruit.P <- 0.001
  # from menge 1975:
  # average annual survival of spawned gametes to postmaturity longevity = 
  # 1.46 x 10-9/m2/year, and annual mortality of gametes is 0.999
  
  delta <- -.02 # density dependence for limpets
  Y <- .01 # whelk conversion rate
  # connell 1961 suggests that whelk predation rate on balanus > chthamalus
  p.whelk.b <- 0.02 # per capita whelk predation rate on balanus 
  p.whelk.c <- 0.02 # per capita whelk predation rate on chth.
  # navarrete (year?) suggests sea star predation rate on balanus > chthamalus
  p.seastar.b <- 0.02 # per capita sea star predation rate
  p.seastar.c <- 0.02 # per capita sea star predation rate on chth.
  
  total <- total_1
  
  # ----------------------------------------------------------------------------------------------------
  # Initial conditions
  
  years <- years_set
  timesteps <- years*12
  month <- rep(c("Apr", "May", "Jun", "Jul", "Aug", 
                 "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar"), years)
  summer <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep")
  winter <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar")
  
  B <- C <- L <- W <- P <- free <- rep(NA, timesteps)
  B[1] <- B_1
  C[1] <- C_1
  L[1] <- L_1
  W[1] <- W_1
  P[1] <- P_1
  
  free[1] <- do.free.space.calculation(total=total, B=B[1], size.B=size.B, C=C[1], 
                                    size.C=size.C, L=L[1], size.L=size.L)
  
  # B.mean <- 90000
  # B.stdev <- sqrt(4.6*10^9)
  location.B <- log(B.mean^2 / sqrt(B.stdev^2 + B.mean^2))
  shape.B <- sqrt(log(1 + (B.stdev^2 / B.mean^2)))
  larvae.B <- rlnorm(n=timesteps, location.B, shape.B)
  
  # C.mean <- 70000
  # C.stdev <- sqrt(2.75*10^9)
  location.C <- log(C.mean^2 / sqrt(C.stdev^2 + C.mean^2))
  shape.C <- sqrt(log(1 + (C.stdev^2 / C.mean^2)))
  larvae.C <- rlnorm(n=timesteps, location.C, shape.C)
  
  # L.mean <- 3000
  # L.stdev <- sqrt(3.8*10^6)
  location.L <- log(L.mean^2 / sqrt(L.stdev^2 + L.mean^2))
  shape.L <- sqrt(log(1 + (L.stdev^2 / L.mean^2)))
  larvae.L <- rlnorm(n=timesteps, location.L, shape.L)
  
  # P.mean <- 727
  # P.stdev <- sqrt(3.4*10^5)
  location.P <- log(P.mean^2 / sqrt(P.stdev^2 + P.mean^2))
  shape.P <- sqrt(log(1 + (P.stdev^2 / P.mean^2)))
  larvae.P <- rlnorm(n=timesteps, location.P, shape.P)
  
  
  # to try to understand dynamics, want to track each variable at each time step
  results <- data.frame(
    timesteps = rep(NA, timesteps),
    B.potential.recruits = rep(NA, timesteps),
    C.potential.recruits = rep(NA, timesteps),
    L.potential.recruits = rep(NA, timesteps),
    B.recruits = rep(NA, timesteps),
    C.recruits = rep(NA, timesteps),
    L.recruits = rep(NA, timesteps),
    P.recruits = rep(NA, timesteps),
    B = rep(NA, timesteps),
    C = rep(NA, timesteps),
    L = rep(NA, timesteps),
    W.recruits = rep(NA, timesteps),
    W = rep(NA, timesteps),
    P = rep(NA, timesteps),
    free = rep(NA, timesteps),
    larvae.B = larvae.B,
    larvae.C = larvae.C,
    larvae.L = larvae.L,
    larvae.P = larvae.P
  )
  
  for (t in 2:timesteps) {
    B.potential.recruits <- do.potential.recruitment(free=free[t-1], size.x=size.B, size.recruit.x=size.recruit.B, 
                                                     larvae.x= ifelse(is.null(var_B),
                                                                      larvae.B[t],
                                                                      var_B))
    C.potential.recruits <- do.potential.recruitment(free=free[t-1], size.x=size.C, size.recruit.x=size.recruit.C, 
                                                     larvae.x=ifelse(is.null(var_C),
                                                                     larvae.C[t],
                                                                     var_C))
    L.potential.recruits <- do.potential.recruitment(free=free[t-1], size.x=size.L, size.recruit.x=size.recruit.L, 
                                                     larvae.x=ifelse(is.null(var_L),
                                                                     larvae.L[t],
                                                                     var_L))
    
    # note that there are more recruits than potential recruits :(
    B.recruits <- do.actual.recruitment(free=free[t-1], L= B.potential.recruits, 
                                        C = settlement.B)
    C.recruits <- do.actual.recruitment(free=free[t-1], L= C.potential.recruits, 
                                        C = settlement.C)
    L.recruits <- do.actual.recruitment(free=free[t-1], L= L.potential.recruits, 
                                        C = settlement.L)
    P.recruits <- ifelse(is.null(var_P),
                         larvae.P[t],
                         var_P)
    
    
    B[t] <- do.population.size.barnacles(S=survival.B, p.whelk=p.whelk.b, W.prev = W[t-1], X.prev=B[t-1],
                                         S.r = survival.recruit.B, R=B.recruits, p.star=p.seastar.b, P.prev=P[t-1])
    C[t] <- do.population.size.barnacles(S=survival.C, p.whelk=p.whelk.c, W.prev = W[t-1], X.prev=C[t-1],
                                         S.r = survival.recruit.C, R=C.recruits, p.star=p.seastar.c, P.prev=P[t-1])
    L[t] <- do.population.size.limpets(S=survival.L, L.prev=L[t-1], S.r=survival.recruit.L, R=L.recruits, 
                                       delta=delta)
    
    # do whelk recruitment only in june
    if(month[t] == "Jun") {
      W.recruits <- do.whelk.recruitment(avg.C=mean(c(C[t], C[t-1], C[t-2])), 
                                         avg.B=mean(c(B[t], B[t-1], B[t-2])),
                                         S = survival.recruit.W,
                                         p.B=p.whelk.b, p.C = p.whelk.c,
                                         Y=Y, 
                                         W.prev=W[t-1], B.prev=B[t], C.prev=C[t])   
    } else {
      W.recruits <- 0
    }
    
    # if holding predation constant...
    if (is.null(W_avg)) {
      W[t] <- do.population.size.whelks(W.prev=W[t-1],
                                        S = survival.W,
                                        R=W.recruits)
    } else {
      W[t] <- W_avg
    }
    
    if (is.null(P_avg)) {
      P[t] <- do.population.size.seastar(S=survival.P, P.prev=P[t-1], 
                                         survival.recruit.P = survival.recruit.P,
                                         R=P.recruits)
    } else {
      P[t] <- P_avg
    }
    
    # to do coexistence invasion, need ability to set starting (and total) population size
    # to 0
    
    if (B_1 == 0) B[t] <- 0
    if (C_1 == 0) C[t] <- 0
    if (L_1 == 0) L[t] <- 0
    if (W_1 == 0) W[t] <- 0
    if (P_1 == 0) P[t] <- 0
    
    free[t] <- do.free.space.calculation(total=total, B=B[t], size.B=size.B, C=C[t], 
                                      size.C=size.C, L=L[t], size.L=size.L)
    
    # track all results
    results$timesteps[t] <- t
    results$B.potential.recruits[t] <- B.potential.recruits
    results$C.potential.recruits[t] <- C.potential.recruits
    results$L.potential.recruits[t] <- L.potential.recruits
    results$B.recruits[t] <- B.recruits
    results$C.recruits[t] <- C.recruits
    results$L.recruits[t] <- L.recruits
    results$P.recruits[t] <- P.recruits
    results$B[t] <- B[t]
    results$C[t] <- C[t]
    results$L[t] <- L[t]
    results$W.recruits[t] <- W.recruits
    results$W[t] <- W[t]
    results$P[t] <- P[t]
    results$free[t] <- free[t]
  }
  
  fd_results <- tibble(time = seq(1:length(B)),
                       balanus_glandula = B,
                       chthamalus_dalli = C,
                       limpets = L,
                       whelks = W,
                       pisaster_ochraceus = P,
                       free_space = free,
                       larvae.B = larvae.B,
                       larvae.C = larvae.C,
                       larvae.L = larvae.L,
                       larvae.P = larvae.P)
  return(fd_results)
  
}

do.growth.rates <- function(results,
                            col_nums # numeric vector of which columns to calculate growth rates for
                            ) {
  tmp <- as.matrix(results[, c(col_nums)])
  tmp_out <- matrix(NA, nrow = nrow(tmp), ncol= ncol(tmp))
  colnames(tmp_out) <- colnames(tmp)
  
  for(i in 2:nrow(results)) {
    tmp_out[i, ] <- tmp[i, ] / tmp[i-1, ]
  }
  
  # sometimes get "inf" responses when dividor = 0
  tmp_out[is.infinite(tmp_out)] <- 0
  
  return(tmp_out)
  
}

