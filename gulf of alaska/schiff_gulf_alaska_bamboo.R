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

################################
##                            ##
##                            ##
##                            ##
## Sr/Ca in each Bamboo coral ##
##                            ##
##                            ##
##                            ##
################################
#' Note, I am no longer going to try to use
#' a package such as XTS or ts() to convert
#' these to timeseries objects. I can still analyze them
#' with just a vector (and this is what you would do anyway
#' with SSA-MTM or kSpectra)

########################################
## Bamboo coral 3808-4, DISK 1        ##
## Sr/Ca                              ##
########################################

bam4.d1 <- read_excel("3808 #4 disk 1 C2 T1 graph.xls")
bam4.d1 <- bam4.d1[-c(698),]
# bam4.d1$scaled <- scale(bam4.d1$srca)
colnames(bam4.d1)[colnames(bam4.d1)=='Sr/Ca ALV-3808 #4 disk 1 C2 T1'] <- "srca"
colnames(bam4.d1)[colnames(bam4.d1)=='Distance from outer edge 1'] <- "distance"
bam4.d1.rev <- bam4.d1.rev %>% dplyr::select(distance, srca)

########################################
## Bamboo coral 3808-3, DISK 2        ##
## Sr/Ca                              ##
########################################

# Write vectors so they can be analyzed in SSA-MTM or kSpectra

##########################################
## Basic moving averages with the data	##
##########################################

########################################
## Bamboo coral 3808-5, DISK 1        ##
## Sr/Ca                              ##
########################################

##################################
##                              ##
##                              ##
##                              ##
## Cadmium in each Bamboo coral ##
##                              ##
##                              ##
##                              ##
##################################

########################################
## Bamboo coral 3808-4, DISK 1        ##
## Cadmium                            ##
########################################

########################################
## Bamboo coral 3808-4, DISK 2        ##
## Cadmium                            ##
########################################
