# --------------------
# John Schiff
# November 7, 2018
# --------------------

################################################################
## Code in this script is for using SSA and other time series ##
## analyses to look at the trace metal measurements in        ##
## bamboo corals                                              ##
################################################################

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

# ---------------------------------------------------------------
# Load the climate data you will compare to
# ---------------------------------------------------------------
npi <- read.csv("npindex_monthly.csv") # North Pacific Index, 1899 - March 2018 (3rd month of the year)
ali <- read.csv("ALPI_1900_2015_EN.csv") # Aleutian Low, annual
nam <- read.delim("SV-NAM-data.txt", sep="") # Northern Annual Mode, monthly
pdo <- read.csv("pdo_data_edit.csv", header = TRUE) # Pacific Decadal Oscillation, from https://www.ncdc.noaa.gov/teleconnections/pdo/

################################################################
# Edit the data depending on the format the file is in         #
################################################################

# Convert North Pacific Index to a time series object, do SSA, and plot a rolling average (moving average)
# ---------------------
# North Pacific Index
# ---------------------

row.names(npi) <- (npi$date)
npi <- npi[1:1428,2]
npi <- npi[-c(552)]
npi.norm <- (npi - mean(npi))/sd(npi)
View(npi)

npi.ts <- ts(npi, start=c(1899), end=c(2018), frequency=12)
norm.ts <- ts(npi.norm, start = c(1899), end = c(2018), frequency = 12)
View(npi.ts)
time(npi.ts)
plot(norm.ts)

npi.ssa <- ssa(npi.ts, L=60)
plot(npi.ssa, "series", groups=c(1:10))
plot(npi.ssa, type="paired")
plot(reconstruct(npi.ssa, groups=list(1)))
npi.recon <- reconstruct(npi.ssa, groups=list(1))
plot(npi.recon$F1)

npi.ma <- ma(npi.ts, order=20, centre=TRUE)
plot(npi.ma)

# ---------------------
# Aleutian Low
# ---------------------

ali <- ali[,2]
ali.ts <- ts(ali, start=c(1900), end=c(2015))
View(ali.ts)
plot(ali.ts)

ali.ma <- ma(ali.ts, order=10, centre=TRUE)
plot(ali.ma)

# ---------------------
# Northern Annual Mode
# ---------------------

nam <- read.delim("SV-NAM-data.txt", sep="")
nam <- nam[,2:13]
nam.ts <- ts(as.vector(t(as.matrix(nam))), 
             start=c(1948,1), end=c(2016,12), frequency=12)

nam.ssa <- ssa(nam.ts, L=60)
plot(nam.ssa, "series", groups=c(1:10))
plot(nam.ssa, type="paired")
plot(reconstruct(nam.ssa, groups=list(1)))
nam.recon <- reconstruct(nam.ssa, groups=list(1))
plot(nam.recon$F1)

plot(nam.ts)

# ----------------------------
# Pacific Decadal Oscillation
# ----------------------------

pdo <- pdo
pdo$Date <- lubridate::ymd(pdo$Date)
pdo.ts <- ts(pdo$Value, start=c(1854), end=c(2019), frequency=12)
pdo.xts <- xts(pdo$Value, pdo$Date)
# pdo.ts <- as.ts(pdo.xts)
plot(pdo.xts)

pdo.ssa <- ssa(pdo.ts, L=984)
plot(pdo.ssa , "series", groups=c(1:10))
plot(pdo.ssa , type="paired")
plot(reconstruct(pdo.ssa , groups=list(1,2,3,4,5)))
pdo.recon <- reconstruct(pdo.ssa, groups=list(2:3))
plot(pdo.recon$F1 ~ pdo$Date, type = "n")
pdo.xts2 <- xts(pdo.recon$F1, pdo$Date)
xyplot(pdo.recon$F1 ~ pdo$Date, type = "l")

###################################################################
## Now we have the climate data converted into time series with  ##
## appropriate spectral analyses and moving averages applied.    ##
## This is useful for when we want to compare them to the actual ##
## bamboo coral trace metal data                                 ##
###################################################################

########################################
## Bamboo coral 3808-4, DISK 1        ##
## Sr/Ca                              ##
########################################

# remember that analyses start from the outside, so need to reverse
bam4.d1 <- read_excel("3808 #4 disk 1 C2 T1 graph.xls")
bam4.d1 <- bam4.d1[-c(698),]
# bam4.d1$scaled <- scale(bam4.d1$srca)
colnames(bam4.d1)[colnames(bam4.d1)=='Sr/Ca ALV-3808 #4 disk 1 C2 T1'] <- "srca"
colnames(bam4.d1)[colnames(bam4.d1)=='Distance from outer edge 1'] <- "distance"
bam4.d1.rev <- rev(bam4.d1)
bam4.d1.rev <- bam4.d1.rev %>% dplyr::select(distance, srca)

# Write vectors so they can be analyzed in SSA-MTM
write.table()
##########################################
## Basic moving averages with the data	##
##########################################

sma.bam4 <- fore(bam4.d1.rev$srca*1000, h = 20)
bam4.d1.rev$ma <- sma.bam4
bam4.d1.rev$ma.normalized

plot(sma.bam4, type = "l")
plot(bam4.d1.rev$ma ~ bam4.d1.rev$distance, type="l",
     col = "blue",
     lwd = 1.5,
     ylim=c(2.85, 3.25))
lines(bam4.d1.rev$srca*1000 ~ bam4.d1.rev$distance,
      col = alpha("black", 0.4))
sma.ts <- ts(bam4.d1.rev$ma, start=c(1938), end = c(2003), frequency = 12)
par(mfrow=c(2,1))
astsa::tsplot(sma.ts, ylab='Sr/Ca 12pt MA', col=4, main='ALV 3808-4', cex.main=1.5,
              xlim=c(1935,2005))
astsa::tsplot(pdo.ts, ylab="PDO Index", col=6,
              xlim=c(1935,2005))
abline(h=0, lty=5)

# astsa::tsplot(forecast::ma(pdo.ts, 12), ylab="PDO Index", col=6,
#               xlim=c(1935,2005))

##########################
## SSA with the Data	##
##########################

bam4.d1.ssa <- ssa(bam4.d1.rev$srca*1000, L = 336, neig = 20)
plot(bam4.d1.ssa, "series", groups=1:15)
plot(bam4.d1.ssa, type="vectors")
plot(bam4.d1.ssa, type="paired")
plot(wcor(bam4.d1.ssa))
recon1 <- reconstruct(bam4.d1.ssa, groups=list(2:3))
bam4.ts <- ts(recon1$F1, start=c(1938), end=c(2003), frequency=12)
plot(bam4.ts)

t.trend <- recon1$F2
plot(t.trend)
abline(v=c(1927, 1977), col="red", lty=5)
abline(v=1947, col="blue")
abline(h=0.003096223, lty=5)

#####################################
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

######################################## EDIT THIS!
## Bamboo coral 3808-4, DISK 3        ##
## Cadmium                            ##
########################################

bam4.d3.cd <- read_excel("ALV-3808#4 D3 C3 working.xls")
bam4.d3.cd <- bam4.d3.cd$`Cd (ppm) ALV-3808 #4 Disk 3 C3 T1`[1:552]
bam4.d3.cd.rev <- rev(bam4.d3.cd)
bam4.d3.rev.det <- detrend(bam4.d3.cd.rev, tt="linear")
plot(bam4.d3.rev.det)
bam4.d3.ts <- ts(bam4.d3.rev.det, start=c(1942), end=c(2002), frequency=12)
plot(bam4.d3.ts, lty=1, bty="n")

bam4.d3.cd.ssa <- ssa(bam4.d3.ts, L = 60)
plot(bam4.d3.cd.ssa)
plot(bam4.d3.cd.ssa, type="vectors")
plot(bam4.d3.cd.ssa, type="paired")
plot(wcor(bam4.d3.cd.ssa))
t.recon3 <- reconstruct(bam4.d3.cd.ssa, groups=list(1:2,3,4,5,6))
# plot(bam4.d3.cd.ssa, "series", groups=1:10)
# plot(reconstruct(bam4.d3.cd.ssa, groups=list(1:2)))
plot(t.recon3$F1)

print(ts.plot(t.recon1$F1,
              t.recon2$F1,
              t.recon3$F1, 
              col=c("#253494", "#2c7fb8", "#a1dab4"), 
              ylim=c(-0.65, 0.6), 
              lty=c(1,2,3), 
              lwd=c(0.5, 0.5, 3)))

abline(v=c(1927, 1977), col="red")
abline(v=1947, col="blue")
abline(h=0, lty=5)

xyplot(t.recon1$F1 + # Disk 1 cadmium
         t.recon2$F1 + # Disk 2 cadmium
         t.recon3$F1, # Disk 3 cadmium
       main="corals")


########################################
## Bamboo coral 3808-3, DISK 2        ##
## Sr/Ca                              ##
########################################

#' Estimated age for this bamboo coral is 75 +/- 5 yrs
#' from the Roark et al (2005) paper
#' Need to doublecheck this data with Brendan

bam3.d2 <- read.csv("3808-3.csv")
bam3.d2 <- bam3.d2[-c(nrow(bam3.d2)),-c(4:ncol(bam3.d2))]
bam3.d2 <- bam3.d2 %>% dplyr::filter(Dist.from.outside > 0.150)
colnames(bam3.d2)[colnames(bam3.d2)=='Sr.Ca.3808.3D2.Alpha.T1'] <- "srca"
# bam3.d2$scaled <- scale(bam3.d2$srca)
bam3.d2$anomaly <- bam3.d2$srca - mean(bam3.d2$srca)
bam3.d2.rev <- rev(bam3.d2)
bam3.d2.rev <- bam3.d2.rev %>% dplyr::select(Dist.from.outside, srca)
bam3.d2.ts <- ts(bam3.d2.rev$srca, 
                 start=c(1927), end=c(2001), frequency = 12) # Convert to time series
# Above, going by 84 years from the text (not the Table, for now), collected in 2002
plot(bam3.d2.ts)
plot(diff(bam3.d2.ts))

bam3.d2.ssa <- ssa((bam3.d2.rev$srca*1000), L = 456) # Research to determine best window length for your data
plot(bam3.d2.ssa) 
plot(bam3.d2.ssa, "series", groups=1:20)
plot(bam3.d2.ssa, type="vectors")
plot(bam3.d2.ssa, type="paired")
plot(wcor(bam3.d2.ssa))
recon1 <- reconstruct(bam3.d2.ssa, groups=list(1, c(2,3), c(4,5), c(6,7)), drop.attributes = TRUE)
plot(bam3.d2.ssa)
plot(recon1$F1)



###############
# UNUSED CODE #
###############
# bam3.d2.ssa <- ssa(bam3.d2.ts, L = 36)
# plot(bam3.d2.ssa, "series", groups=1:20)
# plot(bam3.d2.ssa, type="vectors")
# plot(bam3.d2.ssa, type="paired")
# plot(wcor(bam3.d2.ssa))
# recon1 <- reconstruct(bam3.d2.ssa, groups=list(1))
# plot(recon1$F1)

#' Different way of going about this
#' Link: https://www.r-bloggers.com/wheres-the-magic-emd-and-ssa-in-r/
#' Trying out EMD (empirical decomposition mode) to look at this data
# library(hht)
# ee <- EEMD(c(bam3.d2.ts), time(bam3.d2.ts), 250, 6, "trials")

# resid1 <- residuals(recon1)
# spec.pgram(resid1, detrend = FALSE, log = "no")
# plot(resid1)


# bam4.d1.ts <- ts(bam4.d1.rev$srca, start=c(1938), end=c(2003), frequency=11) # Convert to time series
# Key here: start and end years are based on the ages reported in Roark et al (2005)
# plot(bam4.d1.ts)
# plot(diff(bam4.d1.ts))
# bam4.d1.ssa <- ssa(bam4.d1.rev$srca*1000, L = 336, neig = 20) # Research to determine best window length for your data, probably n/2

# bam4.d1$testnorm <- (srca - min(srca))/(max(srca) - min(srca))
# bam4.d1$testnorm <- (srca - mean(srca))*1000
# bam4.d1.sr <- bam4.d1$`Sr/Ca ALV-3808 #4 disk 1 C2 T1`[1:697]
# bam4.d1.ts <- ts(bam4.d1.rev$`Sr/Ca ALV-3808 #4 disk 1 C2 T1`, start=c(1938), end=c(2002), frequency=11) # Convert to time series

