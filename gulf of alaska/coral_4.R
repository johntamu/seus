####################################################
## This is analysis for ALV3808-4, copied from    ##
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
## Bamboo coral 3808-4, DISK 1        ##
##                                    ##
########################################

# remember that analyses start from the outside, so need to reverse
bam4.d1 <- read_excel("3808 #4 disk 1 C2 T1 graph.xls")
bam4.d1 <- bam4.d1[-c(698),]
colnames(bam4.d1)[colnames(bam4.d1)=='Sr/Ca ALV-3808 #4 disk 1 C2 T1'] <- "srca"
colnames(bam4.d1)[colnames(bam4.d1)=='Distance from outer edge 1'] <- "distance"
colnames(bam4.d1)[colnames(bam4.d1)=='Cd (ppm)#4 disk 1 C2 T1'] <- "cd"
bam4.d1.rev <- bam4.d1 %>% arrange(desc(distance)) # Reorder rows, so that the inner portion (oldest) is first
bam4.d1.rev <- bam4.d1.rev %>% dplyr::select(distance, cd, srca)
growthrate <- 0.09 # mm per year, +/- 0.10 mm
yr.ad <- 2002.5 - (bam4.d1$distance)/(growthrate)
bam4.d1.rev$yr.ad<- yr.ad

########################################
## Bamboo coral 3808-4, DISK 2        ##
##                                    ##
########################################

bam4.d2 <- read_excel("ALV-3808#4 D2 C3 working.xls")
colnames(bam4.d2)[colnames(bam4.d2)=='Distance from outer edge 1'] <- "distance"
colnames(bam4.d2)[colnames(bam4.d2)=='Cd (ppm) ALV-3808 #4 disk 2 C3 T1'] <- "cd"
bam4.d2.rev <- bam4.d2 %>% arrange(desc(distance))
bam4.d2.rev <- bam4.d2.rev %>% dplyr::select(distance, cd)
growthrate <- 0.09 # mm per year, +/- 0.10 mm
yr.ad <- 2002.5 - (bam4.d2$distance)/(growthrate)
bam4.d2.rev$yr.ad <- yr.ad

########################################
## Bamboo coral 3808-4, DISK 3        ##
##                                    ##
########################################

bam4.d3 <- read_excel("ALV-3808#4 D3 C3 working.xls")
colnames(bam4.d3)[colnames(bam4.d3)=='Distance from outer edge 1'] <- "distance"
colnames(bam4.d3)[colnames(bam4.d3)=='Cd (ppm) ALV-3808 #4 Disk 3 C3 T1'] <- "cd"
bam4.d3.rev <- bam4.d3 %>% arrange(desc(distance))
bam4.d3.rev <- bam4.d3.rev %>% dplyr::select(distance, cd)

View(bam4.d1.rev)
View(bam4.d2.rev)
View(bam4.d3.rev)

#' -------------------------------------------------------------------------
#' Based on the data from Roark et al (2005), we know the following things:
#' 1. Corals were collected in 2002 (end date)
#' Based on radiocarbon, Coral 4 is 64 +/- 4 years old.
#' So, max = 68, min = 60, median = 64
#' We do not know the estimated age based on Sr/Ca cycles, however.
#' But, we can use R to find peaks in the Sr/Ca record.
#' Did Roark et al (2005) do cycle count after 3pt smoothing? Why?
#' -------------------------------------------------------------------------

##########################################
## Basic moving averages with the data	##
##########################################

sma.bam4 <- smooth::sma(bam4.d1.rev$srca, h = 20)
sma.bam4 <- TTR::SMA(bam4.d1.rev$srca, n = 12)
sma.bam4 <- forecast::ma(bam4.d1.rev$srca, order = 12, centre = TRUE)
bam4.d1.rev$ma <- sma.bam4

########################
## Calculate anomaly  ##
########################

# bam.anomaly <- bam4.d1.rev$srca - mean(bam4.d1.rev$srca)
bam.anomaly <- (bam4.d1.rev$srca - mean(bam4.d1.rev$srca))/sd(bam4.d1.rev$srca) # Disks 2 and 3 do not have Sr/Ca
# Roark (2005) just subtracted but other sources say to subtract the mean
# and then divide by the standard deviation

plot(bam.anomaly, type = "l")
plot(bam4.d1.rev$srca, type = "l")
bam4.d1.rev$anomaly <- forecast::ma(bam.anomaly, order = 1)

##############
## Plotting ##
##############

plot(sma.bam4, type = "l")
plot(bam4.d1.rev$ma ~ bam4.d1.rev$distance, type="l",
     col = "blue",
     lwd = 1.5)
lines(bam4.d1.rev$srca ~ bam4.d1.rev$distance,
      col = alpha("black", 0.4))
# Add Cadmium to this plot, so: Sr/Ca, Cadmium, PDO Index

##########################
## Create a time series ##
##########################

sma.ts1 <- ts(bam4.d1.rev$ma, end=c(2002, 6), frequency = 12)
sma.ts2 <- ts(bam4.d1.rev$anomaly, end=c(2002, 6), frequency = 12)

# Cd time series
cdts1 <- ts(bam4.d1.rev$cd, end = c(2002, 6), frequency = 12)
cdts2 <- ts(bam4.d2.rev$cd, end = c(2002, 6), frequency = 12)
cdts3 <- ts(bam4.d3.rev$cd, end = c(2002, 6), frequency = 12)

detrend1 <- ts(pracma::detrend(bam4.d1.rev$cd), end = c(2002, 6), frequency = 12)
detrend2 <- ts(pracma::detrend(bam4.d2.rev$cd), end = c(2002, 6), frequency = 12)
detrend3 <- ts(pracma::detrend(bam4.d3.rev$cd), end = c(2002, 6), frequency = 12)

# Create function to normalize data
normalize <- function(data_vector) {
  data_vector <- ((data_vector - mean(data_vector)))
  return(data_vector)
}

################
## Find Peaks ##
################

##############
## Plotting ##
##############

par(mfrow=c(2,1))
astsa::tsplot(sma.ts2, ylab='Sr/Ca Anomaly', col=4, main='ALV 3808-4', cex.main=1.5,
              xlim=c(1940,2002))
abline(h=0, lty = "dashed")
astsa::tsplot(pdo.ts, ylab="PDO Index", col=6,
              xlim=c(1940,2002))
abline(h=0, lty=5)

par(mfrow=c(4,1))
astsa::tsplot(cdts1, ylab='Cd Disk 1', col=4, cex.main=1.5,
              xlim=c(1940,2002),
              xlab = NULL)
astsa::tsplot(cdts2, ylab='Cd Disk 2', col=4, cex.main=1.5,
              xlim=c(1940,2002),
              xlab = NULL)
astsa::tsplot(cdts3, ylab='Cd Disk 3', col=4,  cex.main=1.5,
              xlim=c(1940,2002),
              xlab = NULL)
astsa::tsplot(pdo.ts, ylab="PDO Index", col="darkgrey",
              xlim=c(1940,2002), type = "h",
              xlab = "Calendar Year (C.E)")
lines(forecast::ma(pdo.ts, order = 36))
abline(h=0, lty=5)


######################################
## Plotting Cd against distance in  ##
## each disk.                       ##
######################################
cd1 <- pracma::detrend(bam4.d1$cd)
cd2 <- pracma::detrend(bam4.d2$cd)
cd3 <- pracma::detrend(bam4.d3$cd)

par(mfrow=c(3,1))
plot(cd1 ~ bam4.d1$distance, type="l", bty = "l",
     xlab = NA, ylab = "Cd, Disk 1")
lines(forecast::ma(cd1, order = 12, centre = TRUE))
abline(h=0, lty = 5)
plot(cd2 ~ bam4.d2$distance, type="l", bty = "l",
     xlab = NA, ylab = "Cd, Disk 2")
abline(h=0, lty = 5)
plot(cd3 ~ bam4.d3$distance, type="l", bty = "l",
     xlab = "Distance (mm)", ylab = "Cd, Disk 3")
abline(h=0, lty = 5)


######################################
## Plotting with the detrended data ##
######################################

par(mfrow=c(4,1))
astsa::tsplot(detrend2, ylab='Cd Disk 2, Detrended', col=alpha("darkgrey", 0.4), cex.main=1.5,
              xlim=c(1940,2002),
              xlab = NULL, bty = "n")
lines(forecast::ma(detrend2, order = 12, centre = TRUE), lwd = 2) # Moving average for Cd
abline(v=1972, col = "blue", lwd = 1.5)
abline(v=1982, col = "blue", lwd = 1.5)
abline(v=1997, col = "blue", lwd = 1.5)
abline(h=0, lty=5)
# astsa::tsplot(detrend2, ylab='Cd Disk 2', col=4, cex.main=1.5,
#               xlim=c(1940,2002),
#               xlab = NULL)
# abline(h=0, lty=5)
# astsa::tsplot(detrend3, ylab='Cd Disk 3', col=4,  cex.main=1.5,
#               xlim=c(1940,2002),
#               xlab = NULL)
# abline(h=0, lty=5)
astsa::tsplot(sma.ts2, ylab='Sr/Ca Anomaly, Disk 1', col=alpha("darkgrey", 0.4), cex.main=1.5,
              xlim=c(1940,2002),
              xlab = NULL, type = "h", bty = "n",
              ylim = c(-1, 1))
lines(forecast::ma(sma.ts2, order = 12, centre = TRUE), lwd = 2) # Moving average for Sr/Ca anomaly, Coral 4
abline(v=1972, col = "blue", lwd = 1.5)
abline(v=1982, col = "blue", lwd = 1.5)
abline(v=1997, col = "blue", lwd = 1.5)
abline(h=0, lty = 5)
astsa::tsplot(bam3ts2, ylab='Sr/Ca Ano., Coral 3', col=alpha("darkgrey", 0.4), cex.main=1.5,
              xlim=c(1940,2002),
              ylim=c(-1, 1), type = "h", bty = "n")
lines(forecast::ma(bam3ts2, order = 12, centre = TRUE), lwd = 2) # Moving average for Sr/Ca anomaly, Coral 3
abline(v=1972, col = "blue", lwd = 1.5)
abline(v=1982, col = "blue", lwd = 1.5)
abline(v=1997, col = "blue", lwd = 1.5)
abline(h=0, lty=5)

astsa::tsplot(pdo.ts, ylab="PDO Index", col=ifelse(pdo.ts < 0, "blue","red"),
              xlim=c(1940,2002), type = "h",
              xlab = "Calendar Year (C.E)",
              bty = "n")
lines(forecast::ma(pdo.ts, order = 12, centre = TRUE), lwd = 2)
abline(h=0, lty=5)

# astsa::tsplot(forecast::ma(pdo.ts, 12), ylab="PDO Index", col=6,
#               xlim=c(1935,2005))

############################
## SSA with the Data but  ## 
## without moving average ##
## or anomaly             ##
############################

bam4.d1.ssa <- ssa(bam4.d1.rev$srca, L = 80, neig = 20)
plot(bam4.d1.ssa)
plot(bam4.d1.ssa, "series", groups=1:15)
plot(bam4.d1.ssa, type="vectors")
plot(bam4.d1.ssa, type="paired")
plot(wcor(bam4.d1.ssa))
recon1 <- reconstruct(bam4.d1.ssa, groups=list(2))
bam4.ts <- ts(recon1$F1, start=c(1938), end=c(2003), frequency=12)
plot(recon1)
# resid1 <- residuals(recon1)
# spec.pgram(resid1, detrend = FALSE, log = "no")
# plot(resid1)

# t.trend <- recon1$F2
# plot(t.trend)
# abline(v=c(1927, 1977), col="red", lty=5)
# abline(v=1947, col="blue")
# abline(h=0.003096223, lty=5)

# bam4.d1$testnorm <- (srca - min(srca))/(max(srca) - min(srca))
# bam4.d1$testnorm <- (srca - mean(srca))*1000
# bam4.d1.sr <- bam4.d1$`Sr/Ca ALV-3808 #4 disk 1 C2 T1`[1:697]
# bam4.d1.ts <- ts(bam4.d1.rev$`Sr/Ca ALV-3808 #4 disk 1 C2 T1`, start=c(1938), end=c(2002), frequency=11) # Convert to time series

# SSA with anomaly data
output <- ssa(bam4.d1$anomaly)
plot(output)
plot(output, "series", groups = 1:20)

########################################
## SSA with Sr/Ca after converting to ## 
## a time series object in R.         ##
########################################



########################################
## Bamboo coral 3808-4, DISK 1        ##
## Cadmium                            ##
########################################

bam4.d1.cd <- bam4.d1$`Cd (ppm)#4 disk 1 C2 T1`[1:697]
bam4.d1.cd.rev <- rev(bam4.d1.cd)
detrend.cd <- detrend(bam4.d1.cd.rev, tt="linear")
plot(detrend.cd)
plot(diff(detrend.cd))
diff.cd <- diff(bam4.d1.cd.rev)
# What I'm doing here is detrending the data because of that overall increase
# Not sure why that trend is there?
# b4.cd.ts <- ts(b4.cd.diff, start=c(1938), end=c(2002), frequency=10)
bam4.cd.ts <- ts(detrend.cd, start=c(1938), end=c(2002), frequency=12) # Time series of the detrended data
plot(bam4.cd.ts1, lty=1, xlim=c(1935,2005))
abline(v=c(1927, 1977), col="red", lty=5)
abline(v=1947, col="blue")

# SSA on the Cadmium
bam4.cd.ssa <- ssa(bam4.cd.ts, L = 60)
plot(bam4.cd.ssa)
plot(bam4.cd.ssa, type="vectors")
plot(bam4.cd.ssa, type="paired")
plot(wcor(bam4.cd.ssa))
t.recon1 <- reconstruct(bam4.cd.ssa, groups=list(1,2,3,4,5,6))
plot(bam4.cd.ssa, "series", groups=1:10)
plot(t.recon1$F1)
abline(v=c(1927, 1977), col="red", lty=5)
abline(v=1947, col="blue")


########################################
## Bamboo coral 3808-4, DISK 2        ##
## Cadmium                            ##
########################################

bam4.d2.cd <- read_excel("ALV-3808#4 D2 C3 working.xls")
bam4.d2.cd <- bam4.d2.cd$`Cd (ppm) ALV-3808 #4 disk 2 C3 T1`[1:564]
bam4.d2.cd.rev <- rev(bam4.d2.cd)
bam4.d2.rev.det <- detrend(bam4.d2.cd.rev, tt="linear")
plot(bam4.d2.rev.det)
bam4.d2.ts <- ts(bam4.d2.rev.det, start=c(1938), end=c(2002), frequency=12)
plot(bam4.d2.ts, lty=1, bty="n")
title(main="Bamboo 3808-4, Disk 2, Cadmium, Detrended")

bam4.d2.cd.ssa <- ssa(bam4.d2.ts, L = 60) # Step: Do SSA
plot(bam4.d2.cd.ssa)
plot(bam4.d2.cd.ssa, type="vectors")
plot(bam4.d2.cd.ssa, type="paired")
plot(wcor(bam4.d2.cd.ssa))
t.recon2 <- reconstruct(bam4.d2.cd.ssa, groups=list(1:2,3,4,5,6))
# plot(bam4.d2.cd.ssa, "series", groups=1:10)
# plot(reconstruct(bam4.d2.cd.ssa, groups=list(1:2)))
plot(t.recon2$F1)
title(main="Bamboo 3808-4, Disk 2, SSA-Cadmium")

###########################
## SSA with Anomaly data ##
###########################

bam.anomaly <- (bam4.d1.rev$srca - mean(bam4.d1.rev$srca))/sd(bam4.d1.rev$srca)
bam4.d1.rev$anomaly <- forecast::ma(bam.anomaly, order = 1)

##################
## UNUSED CODE  ##
##################

# bam4.d1.rev$ma.normalized --> I was going to use a column or normalized data?
# Below, I have used findPeaks and set that to p; n is the vector. But findPeaks is not a great function
# See: https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
# plot(n, col = ifelse(n %in% n[p], "red", "black"), type = "o", xlim = c(0,50)) <- find peaks and highlight

# sma.ts <- ts(bam4.d1.rev$ma, start=c(1934), frequency = 10)
# sma.ts2 <- ts(bam4.d1.rev$ma, start=c(1938), frequency = 11)
# sma.ts3 <- ts(bam4.d1.rev$ma, start=c(1930), frequency = 10)

# bam4.d1$scaled <- scale(bam4.d1$srca)

# bam4.d1.ts <- ts(bam4.d1.rev$srca, start=c(1938), end=c(2003), frequency=11) # Convert to time series
# Key here: start and end years are based on the ages reported in Roark et al (2005)
# plot(bam4.d1.ts)
# plot(diff(bam4.d1.ts))
## Research to determine best window length for your data, probably n/2