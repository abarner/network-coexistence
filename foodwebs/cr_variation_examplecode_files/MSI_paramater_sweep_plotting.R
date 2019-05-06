# to plot run plotting_mechanisms.R
library(here)

read.csv(file = here("foodwebs", "cr_variation_examplecode_files", 
                     "MSI_results", "MSIresults.RData"))

coexistence_space <- read.csv(file = here("foodwebs", "cr_variation_examplecode_files", 
                                          "mechanism_figs", "parameter_sweep.csv"))

coexistence_avg <- matrix(NA, nrow=length(omega_seq), ncol=length(sigma_seq))
for (o in 1:length(omega_seq)) {
  for (s in 1:length(sigma_seq)) {
    coexistence_avg[o,s] <- mean(coexistence_space[o,s,])
  }
}

C1_avg <- matrix(NA, nrow=length(omega_seq), ncol=length(sigma_seq))
C1_stdev <- matrix(NA, nrow=length(omega_seq), ncol=length(sigma_seq))
C2_avg <- matrix(NA, nrow=length(omega_seq), ncol=length(sigma_seq))
C2_stdev <- matrix(NA, nrow=length(omega_seq), ncol=length(sigma_seq))

for (o in 1:length(omega_seq)) {
  for (s in 1:length(sigma_seq)) {
    C1_avg[o,s] <- mean(C1_space[o,s,])
    C1_stdev[o,s] <- sd(C1_space[o,s,])
    
    C2_avg[o,s] <- mean(C2_space[o,s,])
    C2_stdev[o,s] <- sd(C2_space[o,s,])
  }
}

library(plyr)
C1_round <- round_any(C1_avg, .05)
C2_round <- round_any(C2_avg, .05)

coexistence_avg[9,16] <- 5
coexistence_avg[5,6] <- 4
coexistence_avg[5,7] <- 4
coexistence_avg[5,8] <- 4
coexistence_avg[5,9] <- 4
coexistence_avg[6,12] <- 4
coexistence_avg[6,13] <- 4
coexistence_avg[7,13] <- 4
coexistence_avg[7,14] <- 4
coexistence_avg[7,15] <- 4
coexistence_avg[8,15] <- 4

library(viridis)
temp <- viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1,
            option = "D")

library(plot3D)

quartz(width=6, height=6)
par(mar=c(2,2,2,2), oma=c(2,2,.5, .5), mfrow=c(2,2))
image2D(coexistence_avg, x=omega_seq, y=sigma_seq, border="black", col=c("#FDE725FF", "#33638DFF", "#440154FF",
                                                                         "#55C667FF", "black"),
         yaxt="n", xaxt="n", colkey=FALSE, xlab="", ylab="",
        main= "Coexistence")
axis(side=1, tick=seq(.5, 1, .1))
mtext(side=2, line=2.5, "Environmental Variation")
axis(side=2, at=seq(0, .75, .25), tick=TRUE)

plot.new()
legend("topleft", c("Species 1 Wins", "Species 2 Wins", "Coexistence", "Species 1 or 2 Wins", "Coexistence or Species 1 Wins"), 
       col=c("#FDE725FF", "#33638DFF", "#440154FF",
             "#55C667FF", "black"), pch=15, xpd=TRUE, bty="n", pt.cex=2)

legend("topleft", c("", "", "", "", ""), 
       col="black", pch=22, xpd=TRUE, bty="n", pt.cex=2)

library(RColorBrewer)
#display.brewer.all(colorblindFriendly = T)

#quartz(width=5, height = 5)
#par(mfrow=c(2,1), oma=c(1,1,1,1), mar=c(2,2,2,2))
min_val <- -.221
max_val <- .159

col_pos=brewer.pal(4,"YlOrRd")
col_neg=brewer.pal(6,"Blues")
col_neg <- rev(col_neg)


#quartz()
#plot.new()
#legend("topleft", c("", "", "", "", "", "", "", "", "", ""),  col=c(col_neg, col_pos), cex=2, pch=15)

cols1 <- c(col_neg[c(1,2,4,5,6)], col_pos[2:4])
cols2 <- c(col_neg[c(5,6)], col_pos[2:3])

image2D(C1_round, x=omega_seq, y=sigma_seq, border="black", col=cols1,
        yaxt="n", xaxt="n", xlab="", ylab="", contour=FALSE, main="Species 1")
axis(side=1, at=seq(.5, 1, .25), tick=TRUE)
axis(side=2, at=seq(0, .75, .25), tick=TRUE)
mtext(side=1, line=2.5, "Predation Preference")
mtext(side=2, line=2.5, "Environmental Variation")

image2D(C2_round, x=omega_seq, y=sigma_seq, border="black", col=cols2,
        yaxt="n", xaxt="n", xlab="", ylab="", contour=FALSE, main="Species 2", colkey=FALSE)
axis(side=1, at=seq(.5, 1, .25), tick=TRUE)
axis(side=2, at=seq(0, .75, .25), tick=TRUE, label=FALSE)
mtext(side=1, line=2.5, "Predation Preference")
