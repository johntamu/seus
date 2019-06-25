## ---
## Topic: Iodine spectral and data analysis in deep-sea coral samples
## John Schiff
## Date: January 18, 2019
# ---

############################################################################
## This code will be used to analyze any Iodine data from deep-sea corals ##
## that we may have. This will be for spectral and statisical analysis.   ##
############################################################################

###########################################################################
## Note: We are not converting the Iodine to a time series object.       ##
## We do not know the age of the coral or any temporal context. We are   ##
## USING the Iodine data to determine a temporal context and chronology. ##
###########################################################################

################
## Libraries  ##
################

# Link: https://stackoverflow.com/questions/29321696/what-is-a-spectrogram-and-how-do-i-set-its-parameters
# Link: https://stackoverflow.com/questions/16341717/detecting-cycle-maxima-peaks-in-noisy-time-series-in-

library(lattice)
library(Rssa) # Singular spectrum analysis
library(zoo)
library(signal)
# library(spectral.methods) # Also does singular spectrum analysis
library(ggplot2)
library(data.table)
library(lattice)
library(pracma)
library(TTR)
# library(hht) # Empirical mode decomposition

# First, import the Iodine data
# Stetson-4904 BC1 Part 1
link1 <- '~/Google Drive/projects/rproj/seus/data/STET_JSL05_4904_BC1_Part 1.csv'
link2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/STET_JSL05_4904_BC1_Part 1.csv'

stetson <- read.csv(link2)
stetson <- read.csv(link1)

stetdata <- as.data.table(stetson) # To use data.table() functions, but could also use dplyr instead
stetdata <- stetdata[, .(Distance, Total.Iodine, Total.Br, BSE.1)]
stetdata <- na.omit(stetdata)
stetdata$i.ma <- forecast::ma(stetdata$Total.Iodine, order = 6, centre = TRUE)
stetdata$b.ma <- forecast::ma(stetdata$Total.Br, order = 6, centre = TRUE)
names(stetdata) <- c("distance", "iodine", "bromine", "BSE", "i.ma", "b.ma") # rename columns for easier coding

s <- stetdata
plot(iodine ~ distance, data = s, type = "l", xlim = c(0, 7000),
     ylab = "Counts", col = alpha("black", 0.4))
lines(bromine ~ distance, data = s, col = "blue")
# From plotting the Iodine and looking at figure 7 in Prouty (2018), Iodine @ 1500 could be a good cutoff

################################
## Spectrum Analysis          ##
################################

################################
## Will use Rssa package here ##
################################

#' A note about SSA and filtering (such as FIR)
#' SSA breaks up a time series into its own
#' different components, similar to PCA. An FIR
#' is just a low-pass filter to get rid of noise
#' 
#' A potential alternative to SSA to exlplore
#' would be the multitaper method.

############################################
## Assign ages to the Stetson record for  ##
## Iodine and Bromine                     ##
############################################

#' This will be useful for putting years bp
#' (or AD) to the Iodine and Bromine trends.
#' Use radiocarbon data (John analysis) and
#' Mitty age estimates to do comparisons.
#' 

#' A note from the Rssa paper by Golyandina:
#' The first component in ssa() is going to be the trend
#' To extract seasonality, you need to remove this trend
#' first. See pg. 24 in the paper for details
#' 
#' 
s$year.ad <- 2005 - (max(s$distance) - s$distance)/growthrate2 # Growth rate from schiff_linear_age_models.R
s$year.ad2 <- 2005 - (max(s$distance) - s$distance)/gr.stet2 # Growth rate using the 4-data point radiocarbon run


#' Note to self: Isolate Iodine pattern after accounting for recon$F1 trend
#' and see how those peaks look, or how the spectrogram looks after that
#' 
#' 
#' 
#' Creating time series of Bromine and Iodine.
#' This is likely unnecessary, and I do the same
#' SSA analysis but without them as time series in a section
#' after this one.
#' 
t.iodine <- ts(s$iodine, end = c(2005, 4), frequency = 8)
t.bromine <- ts(s$bromine, end = c(2005, 4), frequency = 8)
t.bse <- ts(s$BSE, end = c(2005, 4), frequency = 8)

plot(t.iodine)
plot(t.bromine)

t.ssa <- ssa(t.bromine, L = 5000)
plot(t.ssa, groups = 1:100)
plot(t.ssa, "wcor", groups = 1:20)
plot(t.ssa, "series", groups = 30:50)
t.recon1 <- reconstruct(t.ssa, groups = list(1, 2:50))
plot(t.recon1, plot.method = "xyplot", type = "cumsum")
plot(t.bromine, col = alpha("black", 0.4))
lines(t.recon1$F2, col = "blue")

res.trend <- residuals(t.recon1)
spec.pgram(res.trend, detrend = FALSE, log = "no")
t.ssa <- ssa(res.trend)
plot(res.trend, xlim = c(1300, 1500))
plot(t.ssa, "series", groups = 1:20)
data <- as.data.frame(res.trend)
data <- data %>% filter(x > 250)
ts <- ts(data$x)
plot(ts)
#' 02-27-2019
#' Note: It looks like there are peaks in the
#' bromine record. However, I cannot find a clear reason why
#' these peaks exist, and therefore I can't think of an avenue 
#' to explore in regards to this.

res.ssa <- ssa(t.recon1$residuals)

t.ssa <- ssa(t.iodine, L = 256)
plot(t.ssa)
plot(t.ssa, "wcor", groups = 1:20)
plot(t.ssa, "vectors")
plot(t.ssa, "series", groups = 1:20)
t.recon2 <- reconstruct(t.ssa, groups = list(1,2,3,c(4,5,6)))
plot(t.recon2$F1, col = alpha("black", 0.4))
lines(t.recon$F2, col = "blue")

t.ssa <- ssa(t.bse, L = 256)

plot(t.bromine, type = "l", col = alpha("black", 0.3),
     xlim = c(500, 2000))
lines(t.recon$F1, col = "blue", lwd = 2)

plot(t.recon2$F1, type = "l", xlim = c(500, 2000), ylim = c(100, 1000))
lines(t.recon1$F1, col = "blue")

#'
#'
#' Separate section with Bromine and Iodine
#' But I do not convert them to time series
#' objects in R. I keep them as vectors.
#'
#'




################################
## Finding peaks in a series  ##
################################

##############################
## Looking at three things: ## 
## Iodine, Bromine, and BSE ##
##############################
par(mfrow=c(2,1))
plot(iodine ~ distance, data = s, type = "l",
     ylab = "Counts", xlim = c(0, 2000))
plot(bromine ~ year.ad, data = s, type = "l", col = "blue",
     xlim = c(450, 2000))
abline(h=750, lty = "dashed")

plot(BSE ~ distance, data = s, type = "l",
     ylab = "Counts", xlim = c(300,400))

# Isolate the Iodine peaks according to Prouty et al (2018)

# Consider seewave or tuneR packages? Research this

#' Research question:
#' What are the controls on Iodine and Bromine in the 
#' black coral? Assuming the growth ages are true, what
#' can they tell us about climatic events, if anything?
#' Consider the peak number age and radiocarbon age
#' 
#' Also make sure to compare BSE estimates with the Iodine estimates
#' 

##################
##################
##              ##
##              ##
## UNUSED CODE  ##
##              ##
##              ##
##################
##################

# library(pracma)
# p <- findpeaks(iodine, minpeakheight = 1000, minpeakdistance = 2, threshold = 1000)
# length(p)
# plot(p, type = "l")

#   
# xyplot(Total.Br ~ Total.Iodine,
#        stetson)
# xyplot(Total.Iodine ~ Distance,
#        stetson,
#        type="l")
# xyplot(Total.Iodine ~ Distance, data = stetdata, type = "l", xlim = c(3500, 3600))

# s <- ssa(iodine, L = (length(iodine)/4))
# spec <- signal::specgram(iodine,
#                          n = 1024,
#                          Fs = 1000,
#                          window = 1000) # Need some sort of sampling rate, such as how manyu data points per second when scanning
# 
# # discard phase information
# P = abs(spec$S)
# 
# # normalize
# P = P/max(P)
# 
# # convert to dB
# P = 10*log10(P)
# 
# # config time axis
# t = spec$t
# 
# # plot spectrogram
# imagep(x = t,
#        y = spec$f,
#        z = t(P),
#        col = oce.colorsViridis,
#        ylab = 'Frequency [Hz]',
#        xlab = 'Time [s]',
#        drawPalette = T,
#        decimate = F)
# 

#' ############
#' ## Iodine ##
#' ############
#' 
#' 
#' ssai <- ssa(iodine, L = (256))
#' plot(ssai)
#' plot(ssai, "series", groups = 1:20)
#' plot(ssai, type = "vector")
#' plot(ssai, type = "paired")
#' plot(ssai, type = "wcor")
#' # recon <- reconstruct(s, groups = list(c(2), c(3), c(4,5), c(6,7,8,9))) # based on plot(s)
#' reconi <- reconstruct(ssai, groups = list(1, 2, c(3), c(4,5), 6, c(7:10)))
#' plot(wcor(ssai, groups = 1:20))
#' plot(reconi)
#' plot(reconi$F1 ~ s$distance, type = "l", xlim = c(0,12000)) # Compare to estimated age years
#' plot(reconi$F4, type = "l", xlim = c(6000,7000))
#' # Define parameters to construct spectrogram
#' # Link: https://stackoverflow.com/questions/29321696/what-is-a-spectrogram-and-how-do-i-set-its-parameters
#' # Link: https://stackoverflow.com/questions/16341717/detecting-cycle-maxima-peaks-in-noisy-time-series-in-r
#' 
#' #' Below will be a table comparing the Iodine and Radiocarbon ages
#' #' between Mitty's thesis and my work
#' #' Her radioarbon ages are different
#' #' 
#' bc <- c("STET4904-BC1", "JACK4686-BC1")
#' iodine.table <- data.frame()
#' 
#' ##############
#' ## Bromine  ##
#' ##############
#' 
#' bromine <- stetdata$Total.Br
#' b <- ssa(bromine, L = (length(bromine)/3))
#' plot(b)
#' plot(b, "series", groups = 1:50)
#' plot(b, type = "vector")
#' plot(b, type = "paired")
#' plot(b, type = "wcor")
#' reconb <- reconstruct(b, groups = list(c(1), c(2), c(3,4,5), c(6,7,8,9), c(10:15), c(16:21), c(22:24)))
#' plot(reconb, plot.method = "xyplot", add.residuals = FALSE, superpose = TRUE,
#'      auto.key = list(columns = 2))
#' rate <- 8 # estimated growth rate in microns per year
#' dis <- stetdata$Distance
#' yrsbp <- stetdata$Distance/rate
#' yearad <- 2005 - yrsbp
#' plot(bromine-reconb$F2, type = "l")
#' plot(reconb$F4 ~ dis, type = "l", xlim = c(2000,2100))
#' 
#' par(mfrow = c(2,1))
#' plot(reconb$F4 ~ dis, type = "l", xlim = c(2000,2100))
#' plot(bromine ~ dis, type = "l", xlim = c(2000,2100))
#' abline(h = 50)