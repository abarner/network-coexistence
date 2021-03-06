## updated version of the forde & doak 2004 model
## changes primarily to the whelk and seastar predation models

do.free.space.calculation <- function(T, B, size.B, C, size.C, L, size.L) {
  # T = total area
  # B is the population size of Balanus glandula
  # size.B = average size of adult B. glandula
  # C is population size of Chthamalus fissus/dalli
  # size.C = average size of C. fissus/dalli
  # L is population size of limpets
  # size.L = average size of an adult limpet
  F <- T - (B*size.B + C*size.C + L*size.L)
  return(F)
}

# for barnacles and limpets
do.potential.recruitment <- function(F, size.x, size.recruit.x, larvae.x) {
  # F is the amount of free space on the rock
  # size.x is the average size (e.g. size.B, size.C, or size.L)
  # size.recruit.x is the size of a recruit
  # larvae.x is the number of settling larvae of species x
  L <- (F/size.x)*(1-exp(-(size.recruit.x*larvae.x)/F))
  return(L)
}

do.actual.recruitment <- function(F, L, S, L.B, size.B, L.C, size.C, L.L, size.L) {
  # L is the potential recruits of target species
  # S is the survival of recruits of target species
  # L.B is potential recruits of B. glandula
  # L.C is potential recruits of C. fissus/dalli
  # L.L is potential recruits of limpets
  # size.B = average size of adult B. glandula
  # size.C = average size of C. fissus/dalli
  # size.L = average size of an adult limpet
  
  ## note: below is eqn 2 in appendix, causes populations to crash
  ## because the fraction is >>1
  R <- ((F*S) / (L.B*size.B + L.C*size.C + L.L*size.L)) * L
  
  ## rewrite logically as fraction of occupied space
  #R <- ((L.B*size.B + L.C*size.C + L.L*size.L) / F) * S * L
  return(R)
}

do.population.size.barnacles <- function(S, p.whelk, W.prev, X.prev, S.r, R, p.star, P.prev) {
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
  X <- X_before_predation + S.r*R - p.whelk*W.prev*X_before_predation - 
    p.star*P.prev*X_before_predation
  return (X)
}

do.population.size.limpets <- function(S, L.prev, S.r, R, delta) {
  # S is survivorship
  # L.prev is the limpet population size in the previous month
  # S.r is the survivorship of recruits
  # R is the number of recruits
  # delta is density-dependence 
  L <- S*L.prev + S.r*R*exp(delta*L.prev)
  return(L)
}

do.population.size.whelks <- function(W.prev, 
                                      #W.feeding, 
                                      S, 
                                      R) {
  # W.prev is whelk population size in the previous month
  # p is per capita predation rate
  # B.prev is the B. glandula population
  # C.prev is the C. fissus/dalli population
  # R is the recruitment rate
  # S is the survival rate
  # W.feeding is food intake
  
  # note that the functional form of this population model is unclear
  # from the manuscript. Alternate model forms could be:
  W <- W.prev * S + R
  # W <- W.prev * W. feeding + R
  # W <- W.prev*W.feeding*S + R
  
  return(W)
}

do.whelk.recruitment <- function(avg.C, avg.B, p, Y, W.prev, B.prev, C.prev) {
  # avg.C is the average number of C. fissus/dalli from April through June
  # avg. B is the average number of B. glandula from April through June
  # p is the per capita predation rate
  R <- (avg.C + avg.B)*3*p*Y*W.prev*(B.prev+C.prev)
  if (R >= 90) {
    # note that if the max. recruitment rate is not set, then
    # whelk recruitment becomes very large very quickly
    R <- 90
  }
  return(R)
}

do.population.size.seastar <- function(S, P.prev, R) {
  # S is seastar adult survival
  # P.prev is previous seastar population size
  
  # note that there is a typo in equation 7 in Forde & Doak appendix 
  # where "P.prev" (= P at time t-1) is incorrectly included as 
  # P at time t+1
  
  # adult population size should be a function of mortality & prey consumption
  # how to set up predator density dependence? type 2 functional response
  # 'attack rate' = 'encounter rate'
  #type_ii_barn <- (attack * X) / (1 + (attack * handling * X))
  # or use estimates from novak et al. 2017 (calculate myself)
  #attack <- .
  
  #P <- p.po.balanus*B*P.prev + p.po.chtham*C*P.prev + survival.recruit.P*R - (1 -S)*P.prev 
  
  P <- S*P.prev + survival.recruit.P*R
  
  if (P > 6) {
    P <- 6
  }
  if(P<0) {
    P <- 0
  }
  return(P)
}

# ----------------------------------------------------------------------------------------------------
# model parameters
size.B <- .98
size.C <- .32
size.L <- .8

size.recruit.B <- .03
size.recruit.C <- .03
size.recruit.L <- .03

survival.B <- .7
survival.C <- .7
survival.L <- .97
survival.W <- .94
survival.P <- .992


survival.recruit.B <- .7
survival.recruit.C <- .7
survival.recruit.L <- .88
survival.recruit.P <- .998

delta <- -.02 # density dependence for limpets
p.whelk <- .01 # per capita whelk predation rate
Y <- .01 # whelk conversion rate
p.seastar <- .007 # per capita whelk predation rate

T <- 100

# ----------------------------------------------------------------------------------------------------
# Initial conditions

years <- 50
timesteps <- years*12
month <- rep(c("Apr", "May", "Jun", "Jul", "Aug", 
               "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar"), years)
summer <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep")
winter <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar")

B <- C <- L <- W <- P <- F <- rep(NA, timesteps)
B[1] <- 4100
C[1] <- 11000
#B[1] <- 0
#C[1] <- 0
L[1] <- 239
W[1] <- 93
P[1] <- 1

F[1] <- do.free.space.calculation(T=T, B=B[1], size.B=size.B, C=C[1], 
                                  size.C=size.C, L=L[1], size.L=size.L)

B.mean <- 50000
B.stdev <-sqrt(3.24*10^10)
location.B <- log(B.mean^2 / sqrt(B.stdev^2 + B.mean^2))
shape.B <- sqrt(log(1 + (B.stdev^2 / B.mean^2)))

C.mean <- 30000
C.stdev <- sqrt(3.6*10^9)
location.C <- log(C.mean^2 / sqrt(C.stdev^2 + C.mean^2))
shape.C <- sqrt(log(1 + (C.stdev^2 / C.mean^2)))

L.mean <- 2400
L.stdev <- sqrt(1.7*10^7)
location.L <- log(L.mean^2 / sqrt(L.stdev^2 + L.mean^2))
shape.L <- sqrt(log(1 + (L.stdev^2 / L.mean^2)))

P.mean <- 3800
P.stdev <- sqrt(3.7*10^7)
location.P <- log(P.mean^2 / sqrt(P.stdev^2 + P.mean^2))
shape.P <- sqrt(log(1 + (P.stdev^2 / P.mean^2)))
  
 
  
  
for (t in 2:timesteps) {
  B.potential.recruits <- do.potential.recruitment(F=F[t-1], size.x=size.B, size.recruit.x=size.recruit.B, 
                                                   larvae.x=rlnorm(n=1, location.B, shape.B))
  C.potential.recruits <- do.potential.recruitment(F=F[t-1], size.x=size.C, size.recruit.x=size.recruit.C, 
                                                   larvae.x=rlnorm(n=1, location.C, shape.C))
  L.potential.recruits <- do.potential.recruitment(F=F[t-1], size.x=size.L, size.recruit.x=size.recruit.L, 
                                                   larvae.x=rlnorm(n=1, location.L, shape.L))
  
  # note that there are more recruits than potential recruits :(

  B.recruits<-0
  C.recruits<-0
  L.recruits<-0
  P.recruits<-0
  
  if(any(month[t] == c("Apr", "May", "Jun", "Jul", "Aug", "Sep"))) {
  
  B.recruits <- do.actual.recruitment(F=F[t-1], L= B.potential.recruits, 
                                      S = survival.B,
                                      L.B=B.potential.recruits, 
                                      size.B=size.B, L.C=C.potential.recruits, size.C=size.C,
                                      L.L=L.potential.recruits, size.L=size.L)
  
  }
  
  C.recruits <- do.actual.recruitment(F=F[t-1], L= C.potential.recruits, 
                                      S = survival.C,
                                      L.B=B.potential.recruits, 
                                      size.B=size.B, L.C=C.potential.recruits, size.C=size.C,
                                      L.L=L.potential.recruits, size.L=size.L)
  
  
  if(any(month[t] == c("Jul", "Aug", "Sep"))) {
  
  
  L.recruits <- do.actual.recruitment(F=F[t-1], L= L.potential.recruits, 
                                      S = survival.L,
                                      L.B=B.potential.recruits, 
                                      size.B=size.B, L.C=C.potential.recruits, size.C=size.C,
                                      L.L=L.potential.recruits, size.L=size.L)
  }
  
  if(any(month[t] == c("Jun", "Jul", "Aug", "Sep"))) {
  
  P.recruits <- rlnorm(n=1, location.P, shape.P)
  
  }
  
  
  B[t] <- do.population.size.barnacles(S=survival.B, p.whelk=p.whelk, W.prev = W[t-1], X.prev=B[t-1],
                                       S.r = survival.recruit.B, R=B.recruits, p.star=p.seastar, P.prev=P[t-1])
  C[t] <- do.population.size.barnacles(S=survival.C, p.whelk=p.whelk, W.prev = W[t-1], X.prev=C[t-1],
                                       S.r = survival.recruit.C, R=C.recruits, p.star=p.seastar, P.prev=P[t-1])
  L[t] <- do.population.size.limpets(S=survival.L, L.prev=L[t-1], S.r=survival.recruit.L, R=L.recruits, delta=delta)
  
  if(month[t] == "Jun") {
    W.recruits <- do.whelk.recruitment(avg.C=mean(c(C[t], C[t-1], C[t-2])), 
                                       avg.B=mean(c(B[t], B[t-1], B[t-2])),
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
  
  F[t] <- do.free.space.calculation(T=T, B=B[t], size.B=size.B, C=C[t], 
                                    size.C=size.C, L=L[t], size.L=size.L)
  
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
library(tidyverse)
library(ggplot2)
library(tidyr)
fd_results <- tibble(time = seq(1:length(B)),
                     balanus_glandula = B,
                     chthamalus_dalli = C,
                     limpets = L,
                     whelks = W,
                     pisaster_ochraceus = P,
                     free_space = F)

fd_results<-fd_results[fd_results$time>100,]

plot.new()
fd_results %>%
  gather(balanus_glandula:free_space, key = "species", 
         value = "abundance") %>%
  ggplot(aes(x = time, y = abundance)) +
  geom_point(size = 0.5) +
  geom_line() + 
  facet_wrap(~species, scales = "free")
