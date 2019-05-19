# Create graph of mechanistic decomposition
#source("LGS_Vasseur_model.R") # if wanting to run a full set of simulations
library(here)

# or if wanting to read them in
C1_final_mechanisms <- read.csv(file = here("foodwebs", "cr_variation_examplecode_files", 
                                           "mechanism_figs", "C1_final_mechanisms_positive_cross.csv"))

C2_final_mechanisms <- read.csv(file = here("foodwebs", "cr_variation_examplecode_files", 
                                            "mechanism_figs", "C2_final_mechanisms_positive_cross.csv"))

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

C1_results <- c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)
C2_results <- c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)

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

sdC1_results <- c(sdC1_r_bar, sdC1_delta_0, sdC1_delta_P, sdC1_delta_E, sdC1_delta_EP)
sdC2_results <- c(sdC2_r_bar, sdC2_delta_0, sdC2_delta_P, sdC2_delta_E, sdC2_delta_EP)

# ----------------------------------------------------------------------------------------------------
# Plot results
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

quartz(width=6, height=4)
par(mfrow=c(1,2), oma=c(4,2, 1.5, 1), mar=c(0,2,0,0))
x <- barplot(C1_results, ylim=c(-.06, .17), xlab="", ylab=c("Growth Rate When Rare"),
             col=c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"))

abline(h=0)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^P),
                                                 "d" = expression(bar(Delta)[i]^E),
                                                 "e" = expression(bar(Delta)[i]^{E*P})))


box(which = "plot", lty = "solid")
mtext(expression("Competitor 1"), side=3, outer=FALSE, adj=0.5)
text(x=.3, y=.16, "A)")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=C1_results-sdC1_results, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=C1_results+sdC1_results, length=.05,
       angle=90, col=c("black"), code=3)

barplot(C2_results, ylim=c(-.06, .17), xlab="", ylab=c(""),
        col=c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"), yaxt="n")
abline(h=0)
axis(side=1, at=c(.7, 1.9, 3.1, 4.3, 5.5), lab=c("a" = expression(bar("r")[i]-bar("r")[r]) ,
                                                 "b" = expression(bar(Delta)[i]^0),
                                                 "c" = expression(bar(Delta)[i]^P),
                                                 "d" = expression(bar(Delta)[i]^E),
                                                 "e" = expression(bar(Delta)[i]^{E*P})))

axis(side=2, at=seq(-.05, .15, .05), lab=c("", "", "", "", ""))
box(which = "plot", lty = "solid")
arrows(x0=c(.7, 1.9, 3.1, 4.3, 5.5), y0=C2_results-sdC2_results, 
       x1=c(.7, 1.9, 3.1, 4.3, 5.5), y1=C2_results+sdC2_results, length=.05,
       angle=90, col=c("black"), code=3)
text(x=.3, y=.16, "B)")
mtext("Mechanistic Partitioning", side=1, outer=TRUE, adj=0.5, line=2.25)
mtext("Growth Rate When Rare", side=2, outer=TRUE, adj=0.5, line=.5)
mtext(expression("Competitor 2"), side=3, outer=FALSE, adj=0.5)

quartz(5,5)
plot.new()
legend("topleft", c("Low-density growth rate", "Growth rate under avg. conditions", 
                    "Effect of variation in predation rates", "Effect of environmental fluctuations", 
                    "Interactive effect of varying predation \n and environmental fluctuations"), lwd=4,
       col=c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"), bty="n", cex=.75)
