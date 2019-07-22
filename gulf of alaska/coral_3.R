####################################################
## This is analysis for ALV3808-3, copied from    ##
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
library(dplyr)
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
## Bamboo coral 3808-3, DISK 2        ##
## Sr/Ca                              ##
########################################

#' Estimated age for this bamboo coral is 75 +/- 5 yrs
#' So, max = 80, min = 70, median = 75
#' From Sr/Ca cycle counting, the age is approximately 84 years.
#' from the Roark et al (2005) paper
#' Need to doublecheck this data with Brendan

bam3.d2 <- read.csv("3808-3.csv")
bam3.d2 <- bam3.d2[-c(nrow(bam3.d2)),-c(3:ncol(bam3.d2))] # Get rid of unnecessary columns
# bam3.d2 <- bam3.d2 %>% dplyr::filter(Dist.from.outside > 0.150)
colnames(bam3.d2)[colnames(bam3.d2)=='Sr.Ca.3808.3D2.Alpha.T1'] <- "srca"
colnames(bam3.d2)[colnames(bam3.d2)=='Dist.from.outside'] <- "distance"
# bam3.d2$scaled <- scale(bam3.d2$srca)
bam3.d2.rev <- bam3.d2 %>% arrange(desc(distance))
bam3.d2.rev <- bam3.d2.rev %>% dplyr::select(distance, srca)
growthrate <- 0.120 # mm per year, +/- 0.10 mm
yr.ad <- 2002.5 - (bam3.d2$distance)/(growthrate)
bam3.d2.rev$yr.ad <- yr.ad

####################################
## Find Peaks in the Sr/Ca data.  ##
## Still determining best way     ##
## to do this.                    ##
####################################






########################
## Calculate anomaly  ##
########################

# bam3.anomaly <- bam3.d2.rev$srca - mean(bam3.d2.rev$srca)
bam3.anomaly <- (bam3.d2.rev$srca - mean(bam3.d2.rev$srca))/sd(bam3.d2.rev$srca)
bam3.d2.rev$anomaly <- forecast::ma(bam3.anomaly, order = 3, centre = TRUE)

##########################
## Create a time series ##
##########################

# bam3.d2.ts <- ts(bam3.d2.rev$srca, 
#                  start=c(1927), end=c(2001), frequency = 12) # Convert to time series
# # Above, going by 84 years from the text (not the Table, for now), collected in 2002
# plot(bam3.d2.ts)
# plot(diff(bam3.d2.ts))

bam3.d2.ssa <- ssa((bam3.d2.rev$srca), L = 80) # Research to determine best window length for your data
plot(bam3.d2.ssa) 
plot(bam3.d2.ssa, "series", groups=1:20)
plot(bam3.d2.ssa, type="vectors")
plot(bam3.d2.ssa, type="paired")
plot(wcor(bam3.d2.ssa))
recon1 <- reconstruct(bam3.d2.ssa, groups=list(1, c(2,3), c(4,5), c(6,7)), drop.attributes = TRUE)
plot(bam3.d2.ssa)
plot(recon1$F1)

##########################################
## Basic moving averages with the data	##
##########################################

sma.bam3 <- smooth::sma(bam3.d2.rev$srca, h = 20)
sma.bam3 <- forecast::ma(bam3.d2.rev$srca, n = 12)
bam3.d2.rev$ma <- sma.bam3

##############
## Plotting ##
##############

plot(bam3.d2.rev$ma ~ bam3.d2.rev$distance, type="l",
     col = "blue",
     lwd = 1.5)
lines(bam3.d2.rev$srca ~ bam3.d2.rev$distance,
      col = alpha("black", 0.4))

par(mfrow=c(2,1))
astsa::tsplot(bam3ts2, ylab='Sr/Ca Anomaly', col=alpha(4, 0.5), main='ALV 3808-3', cex.main=1.5,
              xlim=c(1925,2002),
              ylim=c(-2.5, 2.5))
lines(forecast::ma(bam3ts2, order = 12, centre = TRUE), lwd = 2)
abline(h=0, lty=5)
astsa::tsplot(pdo.ts, ylab="PDO Index", col=ifelse(pdo.ts < 0, "blue","red"),
              xlim=c(1925,2002), type = "h")
lines(forecast::ma(pdo.ts, order = 24, centre = TRUE))
abline(h=0, lty=5)

######################################
## Plots combining Sr/Ca anomalies  ##
## from Coral 3 and Coral 4         ##
######################################

par(mfrow=c(3,1))
astsa::tsplot(bam3ts2, ylab='Sr/Ca Ano., Coral 3', col=alpha(4, 0.3), cex.main=1.5,
              xlim=c(1925,2002),
              ylim=c(-2.5, 2.5))
lines(forecast::ma(bam3ts2, n = 12), lwd = 2)
abline(h=0, lty=5)
astsa::tsplot(sma.ts2, ylab='Sr/Ca Ano., Coral 4', col=alpha(4, 0.3), cex.main=1.5,
              xlim=c(1925,2002))
lines(forecast::ma(sma.ts2, n = 12), lwd = 2)
abline(h=0, lty = "dashed")
# astsa::tsplot(pdo.ts, ylab="PDO Index", col=6,
#               xlim=c(1925,2002))
# lines(forecast::ma(pdo.ts, n = 24))
astsa::tsplot(norm.ts, xlim = c(1925, 2002), col = alpha(6, 0.3))
lines(forecast::ma(norm.ts, n = 12), lwd = 2)
abline(h=0, lty=5)

##########################
## Create a time series ##
##########################

##########################
## Worth noting: Corals ##
## were collected in    ##
## June 2002.           ##
##########################

# bam3ts1 <- ts(bam3.d2.rev$ma, end = c(2002, 6), frequency = 12) # Using end date rather than beginning date
bam3ts2 <- ts(bam3.d2.rev$anomaly, end = c(2002, 6), frequency = 12)

###########################
## SSA with Anomaly data ##
###########################

data <- bam3ts2
# data <- bam3.d2.rev$anomaly

ssa3 <- ssa(data, L = 120)
plot(ssa3)
plot(ssa3, type = "series", groups = 1:20)
plot(ssa3, type = "vector")
plot(ssa3, type = "paired")
plot(ssa3, type = "wcor")

recon <- reconstruct(ssa3, groups = list(c(1),c(2),c(5,6),c(4,5,6,7),c(8,9,10)))
plot(recon$F3, xlim =c(1940, 1950))
abline(h=0, lty = 5)

########################
## UNUSED CODE BELOW  ##
########################

# bam3ts1 <- ts(bam3.d2.rev$ma, start=c(1922), frequency = 12) # 80 years old
# bam3ts2 <- ts(bam3.d2.rev$ma, start=c(1927), frequency = 12) # 75 years old
# bam3ts3 <- ts(bam3.d2.rev$ma, start=c(1932), frequency = 13) # 70 years old
# bam3ts4 <- ts(bam3.d2.rev$ma, start=c(1918), frequency = 11) # 84 years old

