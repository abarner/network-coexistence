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
  if (X_before_predation < 0) {X_before_predation <- 0}
  
  X <- X_before_predation + S.r*R - p.whelk*W.prev*X_before_predation - p.star*P.prev*X_before_predation
  if (X < 0) {X <- 0}
  
  return (X)
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
  
  W <- W.prev * S + R
  
  # note that the functional form of this population model is unclear
  # from the manuscript. Alternate model forms could be:
  # W <- W.prev * W. feeding + R
  # W <- W.prev*W.feeding*S + R
  
  return(W)
}

do.whelk.recruitment <- function(avg.C, avg.B, p, Y, W.prev, B.prev, C.prev, S, 
                                 r = 1) {
  # Updated from F&D to explicitly incorporate density dependence
  # (Which was previously set at a hard limit of 90)
  
  # avg.C is the average number of C. fissus/dalli from April through June
  # avg. B is the average number of B. glandula from April through June
  # p is the per capita predation rate
  
  R <- (avg.C + avg.B)*3*p*Y*W.prev*S*(B.prev + C.prev)
  
  R_logistic <- (R / 90) * exp(r * (1 - (R / 90)))
  
  return(R_logistic * R)
}

do.population.size.seastar <- function(S, P.prev, R, r = 1) {
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
  settlement.B = .002 * 30 * 24,
  settlement.C = .002 * 30 * 24,
  settlement.L = .0002 * 30 * 24,
  B.mean = 90000,
  B.stdev = sqrt(4.6*10^9),
  C.mean = 70000,
  C.stdev = sqrt(2.75*10^9),
  L.mean = 3000,
  L.stdev = sqrt(3.8*10^6),
  P.mean = 727,
  P.stdev = sqrt(3.4*10^5)
) {
  
  # ----------------------------------------------------------------------------------------------------
  # Model parameters
  
  size.B <- .000098
  size.C <- .000032
  size.L <- .00008
  
  size.recruit.B <- .000003
  size.recruit.C <- .000003
  size.recruit.L <- .000003
  
  survival.B <- .7
  survival.C <- .7
  survival.L <- .97
  survival.W <- .94
  survival.P <- .992
  
  # in connolly & roughgarden 1999, table 1, settlement is 0.02/hour/m2 for barn
    # settlement.B <- .002 * 30 * 24
    # settlement.C <- .002 * 30 * 24
    # settlement.L <- .0002 * 30 * 24
  # estimates for limpet settlement is relative to barnacle settlement
  # don't have an estimate otherwise
  # use Fig 5a in Gilman 2006 Ecography (5-10 recruits max / 625 cm2)
  # vs Figures from Menge et al. 2011 JEMBE (50-100 recruits / 100 cm2)
  # so should be about ~ an order of magnitude lower than barnacle settlement?
  
  survival.recruit.B <- .7
  survival.recruit.C <- .7
  survival.recruit.L <- .88
  survival.recruit.W <- .88
  survival.recruit.P <- 0.001
  # from menge 1975:
  # average annual survival of spawned gametes to postmaturity longevity = 
  # 1.46 x 10-9/m2/year, and annual mortality of gametes is 0.999
  
  delta <- -.02 # density dependence for limpets
  p.whelk <- .001 # per capita whelk predation rate
  Y <- .001 # whelk conversion rate
  p.seastar <- .007 # per capita whelk predation rate
  
  total <- 1
  
  # ----------------------------------------------------------------------------------------------------
  # Initial conditions
  
  years <- 50
  timesteps <- years*12
  month <- rep(c("Apr", "May", "Jun", "Jul", "Aug", 
                 "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar"), years)
  summer <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep")
  winter <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar")
  
  B <- C <- L <- W <- P <- free <- rep(NA, timesteps)
  B[1] <- 4100
  C[1] <- 11000
  #B[1] <- 0
  #C[1] <- 0
  L[1] <- 239
  W[1] <- 93
  P[1] <- 1
  
  free[1] <- do.free.space.calculation(total=total, B=B[1], size.B=size.B, C=C[1], 
                                    size.C=size.C, L=L[1], size.L=size.L)
  
  # B.mean <- 90000
  # B.stdev <- sqrt(4.6*10^9)
  location.B <- log(B.mean^2 / sqrt(B.stdev^2 + B.mean^2))
  shape.B <- sqrt(log(1 + (B.stdev^2 / B.mean^2)))
  
  # C.mean <- 70000
  # C.stdev <- sqrt(2.75*10^9)
  location.C <- log(C.mean^2 / sqrt(C.stdev^2 + C.mean^2))
  shape.C <- sqrt(log(1 + (C.stdev^2 / C.mean^2)))
  
  # L.mean <- 3000
  # L.stdev <- sqrt(3.8*10^6)
  location.L <- log(L.mean^2 / sqrt(L.stdev^2 + L.mean^2))
  shape.L <- sqrt(log(1 + (L.stdev^2 / L.mean^2)))
  
  # P.mean <- 727
  # P.stdev <- sqrt(3.4*10^5)
  location.P <- log(P.mean^2 / sqrt(P.stdev^2 + P.mean^2))
  shape.P <- sqrt(log(1 + (P.stdev^2 / P.mean^2)))
  
  
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
    free = rep(NA, timesteps)
  )
  
  for (t in 2:timesteps) {
    B.potential.recruits <- do.potential.recruitment(free=free[t-1], size.x=size.B, size.recruit.x=size.recruit.B, 
                                                     larvae.x=rlnorm(n=1, location.B, shape.B))
    C.potential.recruits <- do.potential.recruitment(free=free[t-1], size.x=size.C, size.recruit.x=size.recruit.C, 
                                                     larvae.x=rlnorm(n=1, location.C, shape.C))
    L.potential.recruits <- do.potential.recruitment(free=free[t-1], size.x=size.L, size.recruit.x=size.recruit.L, 
                                                     larvae.x=rlnorm(n=1, location.L, shape.L))
    
    # note that there are more recruits than potential recruits :(
    B.recruits <- do.actual.recruitment(free=free[t-1], L= B.potential.recruits, 
                                        C = settlement.B)
    C.recruits <- do.actual.recruitment(free=free[t-1], L= C.potential.recruits, 
                                        C = settlement.C)
    L.recruits <- do.actual.recruitment(free=free[t-1], L= L.potential.recruits, 
                                        C = settlement.L)
    P.recruits <- rlnorm(n=1, location.P, shape.P)
    
    
    B[t] <- do.population.size.barnacles(S=survival.B, p.whelk=p.whelk, W.prev = W[t-1], X.prev=B[t-1],
                                         S.r = survival.recruit.B, R=B.recruits, p.star=p.seastar, P.prev=P[t-1])
    C[t] <- do.population.size.barnacles(S=survival.C, p.whelk=p.whelk, W.prev = W[t-1], X.prev=C[t-1],
                                         S.r = survival.recruit.C, R=C.recruits, p.star=p.seastar, P.prev=P[t-1])
    L[t] <- do.population.size.limpets(S=survival.L, L.prev=L[t-1], S.r=survival.recruit.L, R=L.recruits, delta=delta)
    
    if(month[t] == "Jun") {
      W.recruits <- do.whelk.recruitment(avg.C=mean(c(C[t], C[t-1], C[t-2])), 
                                         avg.B=mean(c(B[t], B[t-1], B[t-2])),
                                         S = survival.recruit.W,
                                         p=p.whelk, Y=Y, W.prev=W[t-1], B.prev=B[t], C.prev=C[t])   
    } else {
      W.recruits <- 0
    }
    
    ## Include this if using form of whelk population model that depends 
    ## on barnacle population size
    
    # if(month[t] %in% summer) {
    #   W.feeding <- p.whelk*(B[t-1]+C[t-1]) / (1+ (p.whelk*(B[t-1]+C[t-1])))
    # } else {
    #   W.feeding <- p.whelk*(B[t-6]+C[t-6]) / (1+ (p.whelk*(B[t-6]+C[t-6])))
    # }
    
    W[t] <- do.population.size.whelks(W.prev=W[t-1],
                                      #W.feeding = W.feeding, 
                                      S = survival.W,
                                      R=W.recruits)
    
    P[t] <- do.population.size.seastar(S=survival.P, P.prev=P[t-1], R=P.recruits)
    
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
                       free_space = free)
  return(fd_results)
  
}


# quartz(width=8, height=4)
# par(mfrow=c(1,3))
# plot(B, xlab="Time", ylab="B. glandula")
# plot(C, xlab="Time", ylab="C. fissus")
# plot(L, xlab="Time", ylab="Limpets")
# 
# quartz(width=6, height=4)
# par(mfrow=c(1,3))
# plot(W, xlab="Time", ylab="Whelk")
# plot(P, xlab="Time", ylab="Seastars")
# plot(F, xlab = "Time", ylab = "Free Space")

## plot ggplot results
# fd_results <- tibble(time = seq(1:length(B)),
#                      balanus_glandula = B,
#                      chthamalus_dalli = C,
#                      limpets = L,
#                      whelks = W,
#                      pisaster_ochraceus = P,
#                      free_space = F)
# plot.new()
#png(filename = "forde_and_doak_new_dynamics.png", width = 10, height = 7, units = "in", res = 300)
# fd_results %>%
#   gather(balanus_glandula:free_space, key = "species", 
#          value = "abundance") %>%
#   mutate(neg_value = ifelse(abundance < 0, "yes", "no")) %>%
#   ggplot(aes(x = time, y = abundance)) +
#   geom_point(size = 2, aes(col = neg_value)) +
#   geom_line() + 
#   facet_wrap(~species, scales = "free") +
#   theme(legend.position = "none")
#dev.off()

