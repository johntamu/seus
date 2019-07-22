####################################################
## This is analysis for ALV3808-5, copied from    ##
## the main script file, but also includes edits. ##
####################################################

################################
## Do not expect code here    ##
## to be exactly the same     ##
## as in the original script  ##
## file.                      ##
################################

# Load the necessary libraries
library(Rssa)
library(zoo)
library(signal)
library(spectral.methods)
library(ggplot2)
library(lattice)
library(TSA)
library(pracma)
library(readxl)
library(xts)
library(astsa)
library(forecast)

########################################
## Bamboo coral 3808-5, DISK 1        ##
##                                    ##
########################################

# remember that analyses start from the outside, so need to reverse
bam5 <- read_excel("3808 #5 D1 TopC1 graph.xls")
bam5 <- bam5[-c(1397),]
colnames(bam5)[colnames(bam5)=='Sr/Ca#5 Disk 1 TopC1 T1'] <- "srca"
colnames(bam5)[colnames(bam5)=='Distance from outer edge 1'] <- "distance"
colnames(bam5)[colnames(bam5)=='Cd (ppm) ALV-3808 #5 Disk 1 TopC1 T1'] <- "cd"
bam5.rev <- bam5 %>% arrange(desc(distance)) # Reorder rows, so that the inner portion (oldest) is first
bam5.rev <- bam5.rev %>% dplyr::select(distance, cd, srca)

plot(bam5.rev$cd ~ bam5.rev$distance, type = "l")
plot(bam5.rev$anom ~ bam5.rev$distance, type = "l")

########################
## Calculate anomaly  ##
########################

anom1 <- (bam5.rev$srca - mean(bam5.rev$srca))/sd(bam5.rev$srca)
bam5.rev$anom <- anom1

par(mfrow = c(2,1))
plot(detrend(bam5.rev$cd) ~ bam5.rev$distance, type = "l", col = alpha("darkgrey", 0.4),
     xlab = NA, ylab = "Cadmium", bty = "n")
lines(forecast::ma(detrend(bam5.rev$cd), order = 12) ~ bam5.rev$distance, 
      type = "l", lwd = 2)
abline(h=0, lty = 5)
plot(bam5.rev$anom ~ bam5.rev$distance, type = "l", col = alpha("darkgrey", 0.4),
     ylim = c(-1,1),
     xlab = "Distance (mm)", ylab = "Sr/Ca Anomaly", bty = "n")
lines(forecast::ma(bam5.rev$anom, order = 12) ~ bam5.rev$distance, 
      type = "l", lwd = 2)
abline(h=0, lty = 5)

###########################
## SSA with Anomaly data ##
###########################

ssa5 <- ssa(bam5.rev$anom, L = 120)
plot(ssa5)
plot(ssa5, type = "vector")
