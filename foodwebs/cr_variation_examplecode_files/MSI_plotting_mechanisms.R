# Create graph of mechanistic decomposition
#source("LGS_Vasseur_model.R") # if wanting to run a full set of simulations
library(here)

load(file = here("foodwebs", "cr_variation_examplecode_files", 
                 "MSI_results", "MSI_sigma1_omega90.RData"))

C1_r_bar <- mean(C1_final_mechanisms[,1])
C1_delta_0 <- mean(C1_final_mechanisms[,2])
C1_delta_P <- mean(C1_final_mechanisms[,3])
C1_delta_E <- mean(C1_final_mechanisms[,4])
C1_delta_EP <- mean(C1_final_mechanisms[,5])

C2_r_bar <- mean(C2_final_mechanisms[,1])
C2_delta_0 <- mean(C2_final_mechanisms[,2])
C2_delta_P <- mean(C2_final_mechanisms[,3])
C2_delta_E <- mean(C2_final_mechanisms[,4])
C2_delta_EP <- mean(C2_final_mechanisms[,5])

C1_results_sigma1 <- c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)
C2_results_sigma1 <- c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)

# and calculate standard deviation
sdC1_r_bar <- sd(C1_final_mechanisms[,1])
sdC1_delta_0 <- sd(C1_final_mechanisms[,2])
sdC1_delta_P <- sd(C1_final_mechanisms[,3])
sdC1_delta_E <- sd(C1_final_mechanisms[,4])
sdC1_delta_EP <- sd(C1_final_mechanisms[,5])

sdC2_r_bar <- sd(C2_final_mechanisms[,1])
sdC2_delta_0 <- sd(C2_final_mechanisms[,2])
sdC2_delta_P <- sd(C2_final_mechanisms[,3])
sdC2_delta_E <- sd(C2_final_mechanisms[,4])
sdC2_delta_EP <- sd(C2_final_mechanisms[,5])

sdC1_results_sigma1 <- c(sdC1_r_bar, sdC1_delta_0, sdC1_delta_P, sdC1_delta_E, sdC1_delta_EP)
sdC2_results_sigma1 <- c(sdC2_r_bar, sdC2_delta_0, sdC2_delta_P, sdC2_delta_E, sdC2_delta_EP)

# ----------------------------------------------------------------------------------------------------
# Plot results
#quartz(width=4.5, height=6)
pdf(file=here("foodwebs", "cr_variation_examplecode_files", 
              "MSI_results", "sigma_variation_results.pdf"), width=4.5, height=7)
par(mfrow=c(3,2), oma=c(4,2, 1.5, 1), mar=c(.5,2,.5,0))
x <- barplot(C1_results_sigma1, ylim=c(-.1, .17), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")
abline(h=0)
axis(side=2, at=c(-.1, 0, .1), labels=TRUE, tick=TRUE)
#axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
#                                                 "b" = expression(bar(Delta)[i]^0),
#                                                 "c" = expression(bar(Delta)[i]^P),
#                                                 "d" = expression(bar(Delta)[i]^E),
#                                                 "e" = expression(bar(Delta)[i]^{E*P})))


box(which = "plot", lty = "solid")
mtext(expression("Competitor 1"), side=3, outer=FALSE, adj=0.5)
text(x=5.9, y=.155, "A)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=C1_results_sigma1-sdC1_results_sigma1, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=C1_results_sigma1+sdC1_results_sigma1, length=.05,
       angle=90, col=c("black"), code=3)

barplot(C2_results_sigma1, ylim=c(-.1, .17), xlab="", ylab=c(""),
        col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")
abline(h=0)
#axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
#                                                 "b" = expression(bar(Delta)[i]^0),
#                                                 "c" = expression(bar(Delta)[i]^P),
#                                                 "d" = expression(bar(Delta)[i]^E),
#                                                 "e" = expression(bar(Delta)[i]^{E*P})))

axis(side=2, at=c(-.1, .0, .1), lab=c("", "", ""))
box(which = "plot", lty = "solid")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=C2_results_sigma1-sdC2_results_sigma1, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=C2_results_sigma1+sdC2_results_sigma1, length=.05,
       angle=90, col=c("black"), code=3)
text(x=5.9, y=.155, "B)")
#mtext("Mechanistic Partitioning", side=1, outer=TRUE, adj=0.5, line=2.25)
#mtext("Growth Rate When Rare", side=2, outer=TRUE, adj=0.5, line=.5)
mtext(expression("Competitor 2"), side=3, outer=FALSE, adj=0.5)


# ------------------------------------------------------------------------------------------------------------
load(file = here("foodwebs", "cr_variation_examplecode_files", 
                 "MSI_results", "MSI_sigma4_omega90.RData"))

C1_r_bar <- mean(C1_final_mechanisms[,1])
C1_delta_0 <- mean(C1_final_mechanisms[,2])
C1_delta_P <- mean(C1_final_mechanisms[,3])
C1_delta_E <- mean(C1_final_mechanisms[,4])
C1_delta_EP <- mean(C1_final_mechanisms[,5])

C2_r_bar <- mean(C2_final_mechanisms[,1])
C2_delta_0 <- mean(C2_final_mechanisms[,2])
C2_delta_P <- mean(C2_final_mechanisms[,3])
C2_delta_E <- mean(C2_final_mechanisms[,4])
C2_delta_EP <- mean(C2_final_mechanisms[,5])

C1_results_sigma4 <- c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)
C2_results_sigma4 <- c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)

# and calculate standard deviation
sdC1_r_bar <- sd(C1_final_mechanisms[,1])
sdC1_delta_0 <- sd(C1_final_mechanisms[,2])
sdC1_delta_P <- sd(C1_final_mechanisms[,3])
sdC1_delta_E <- sd(C1_final_mechanisms[,4])
sdC1_delta_EP <- sd(C1_final_mechanisms[,5])

sdC2_r_bar <- sd(C2_final_mechanisms[,1])
sdC2_delta_0 <- sd(C2_final_mechanisms[,2])
sdC2_delta_P <- sd(C2_final_mechanisms[,3])
sdC2_delta_E <- sd(C2_final_mechanisms[,4])
sdC2_delta_EP <- sd(C2_final_mechanisms[,5])

sdC1_results_sigma4 <- c(sdC1_r_bar, sdC1_delta_0, sdC1_delta_P, sdC1_delta_E, sdC1_delta_EP)
sdC2_results_sigma4 <- c(sdC2_r_bar, sdC2_delta_0, sdC2_delta_P, sdC2_delta_E, sdC2_delta_EP)

# ----------------------------------------------------------------------------------------------------
# Plot results
#quartz(width=6, height=4)
#par(mfrow=c(1,2), oma=c(4,2, 1.5, 1), mar=c(0,2,0,0))
x <- barplot(C1_results_sigma4, ylim=c(-.1, .17), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")
abline(h=0)
#axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
#                                                 "b" = expression(bar(Delta)[i]^0),
#                                                 "c" = expression(bar(Delta)[i]^P),
#                                                 "d" = expression(bar(Delta)[i]^E),
#                                                 "e" = expression(bar(Delta)[i]^{E*P})))

axis(side=2, at=c(-.1, 0, .1), labels=TRUE, tick=TRUE)
box(which = "plot", lty = "solid")
#mtext(expression("Competitor 1"), side=3, outer=FALSE, adj=0.5)
text(x=5.9, y=.155, "C)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=C1_results_sigma4-sdC1_results_sigma4, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=C1_results_sigma4+sdC1_results_sigma4, length=.05,
       angle=90, col=c("black"), code=3)

barplot(C2_results_sigma4, ylim=c(-.1, .17), xlab="", ylab=c(""),
        col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")
abline(h=0)
#axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
#                                                 "b" = expression(bar(Delta)[i]^0),
#                                                 "c" = expression(bar(Delta)[i]^P),
#                                                 "d" = expression(bar(Delta)[i]^E),
#                                                 "e" = expression(bar(Delta)[i]^{E*P})))

axis(side=2, at=c(-.1, 0, .1), labels=FALSE)
box(which = "plot", lty = "solid")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=C2_results_sigma4-sdC2_results_sigma4, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=C2_results_sigma4+sdC2_results_sigma4, length=.05,
       angle=90, col=c("black"), code=3)
text(x=5.9, y=.155, "D)")
#mtext("Mechanistic Partitioning", side=1, outer=TRUE, adj=0.5, line=2.25)
mtext("Growth Rate When Rare", side=2, outer=TRUE, adj=0.5, line=.5)
#mtext(expression("Competitor 2"), side=3, outer=FALSE, adj=0.5)


# ------------------------------------------------------------------------------------------------------------
load(file = here("foodwebs", "cr_variation_examplecode_files", 
                 "MSI_results", "MSI_sigma7_omega90.RData"))

C1_r_bar <- mean(C1_final_mechanisms[,1])
C1_delta_0 <- mean(C1_final_mechanisms[,2])
C1_delta_P <- mean(C1_final_mechanisms[,3])
C1_delta_E <- mean(C1_final_mechanisms[,4])
C1_delta_EP <- mean(C1_final_mechanisms[,5])

C2_r_bar <- mean(C2_final_mechanisms[,1])
C2_delta_0 <- mean(C2_final_mechanisms[,2])
C2_delta_P <- mean(C2_final_mechanisms[,3])
C2_delta_E <- mean(C2_final_mechanisms[,4])
C2_delta_EP <- mean(C2_final_mechanisms[,5])

C1_results_sigma7 <- c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)
C2_results_sigma7 <- c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)

# and calculate standard deviation
sdC1_r_bar <- sd(C1_final_mechanisms[,1])
sdC1_delta_0 <- sd(C1_final_mechanisms[,2])
sdC1_delta_P <- sd(C1_final_mechanisms[,3])
sdC1_delta_E <- sd(C1_final_mechanisms[,4])
sdC1_delta_EP <- sd(C1_final_mechanisms[,5])

sdC2_r_bar <- sd(C2_final_mechanisms[,1])
sdC2_delta_0 <- sd(C2_final_mechanisms[,2])
sdC2_delta_P <- sd(C2_final_mechanisms[,3])
sdC2_delta_E <- sd(C2_final_mechanisms[,4])
sdC2_delta_EP <- sd(C2_final_mechanisms[,5])

sdC1_results_sigma7 <- c(sdC1_r_bar, sdC1_delta_0, sdC1_delta_P, sdC1_delta_E, sdC1_delta_EP)
sdC2_results_sigma7 <- c(sdC2_r_bar, sdC2_delta_0, sdC2_delta_P, sdC2_delta_E, sdC2_delta_EP)

# ----------------------------------------------------------------------------------------------------
# Plot results
#quartz(width=6, height=4)
#par(mfrow=c(1,2), oma=c(4,2, 1.5, 1), mar=c(0,2,0,0))
x <- barplot(C1_results_sigma7, ylim=c(-.1, .17), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")
abline(h=0)
axis(side=2, at=c(-.1, 0., .1), labels=TRUE, tick=TRUE)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^P),
                                                 "d" = expression(bar(Delta)[i]^E),
                                                 "e" = expression(bar(Delta)[i]^{E*P})))


box(which = "plot", lty = "solid")
#mtext(expression("Competitor 1"), side=3, outer=FALSE, adj=0.5)
text(x=5.9, y=.155, "E)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=C1_results_sigma7-sdC1_results_sigma7, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=C1_results_sigma7+sdC1_results_sigma7, length=.05,
       angle=90, col=c("black"), code=3)

barplot(C2_results_sigma7, ylim=c(-.1, .17), xlab="", ylab=c(""),
        col=c("grey40", "grey70", "grey70", "grey70", "grey70"), yaxt="n")
abline(h=0)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^P),
                                                 "d" = expression(bar(Delta)[i]^E),
                                                 "e" = expression(bar(Delta)[i]^{E*P})))

axis(side=2, at=c(-.1, 0, .1), labels=FALSE)
box(which = "plot", lty = "solid")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=C2_results_sigma7-sdC2_results_sigma7, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=C2_results_sigma7+sdC2_results_sigma7, length=.05,
       angle=90, col=c("black"), code=3)
text(x=5.9, y=.155, "F)")
mtext("Mechanistic Partitioning", side=1, outer=TRUE, adj=0.5, line=2.25)
#mtext("Growth Rate When Rare", side=2, outer=TRUE, adj=0.5, line=.5)
#mtext(expression("Competitor 2"), side=3, outer=FALSE, adj=0.5)
dev.off()
