## ---
## Topic: Climate data from North Atlantic (e.g., AMO Index) to compare with deep-sea coral isotope data
## John Schiff
## Date: January 28, 2019
## ---

##################################################################
## This code will be used to do some data analysis              ##          
## on AMO index data (since it is SST) and see how it compares  ##
## at all to the SEUS deep-sea coral data                       ##
##################################################################

paths.win <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/amo index from noaa esrl.txt'
# paths.mac 

amo <- read.table(paths.win,fill = TRUE)
amo <- amo[-c(1,165:168),] # remove superfluous metadata
# amo <- amo[,-c(1)]
# colnames(amo) <- c("Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
months <- format(seq.Date(as.Date("2013-01-01"), as.Date("2013-12-01"), 
                          by = "month"), format = "%b")
colnames(amo) <- c("Year", months)
# Issue here is that the data is in a monthly format, which is good, but we need to consider that when analyzing
# amo.ts <- ts(amo, start=c(1856), frequency = 12)

######################################
## Example data below to set up a   ## 
## dataframe with monthly data over ## 
## 100 years                        ##
######################################
# Link: https://stackoverflow.com/questions/16250724/converting-a-data-frame-to-monthly-time-series

# set.seed(12)
# dummy.df <- as.data.frame(matrix(round(rnorm(1200),digits=2),nrow=100,ncol=12))
# months <- format(seq.Date(as.Date("2013-01-01"), as.Date("2013-12-01"), 
#                           by = "month"), format = "%b")
# colnames(dummy.df) <- months
# dummy.df$Year <- seq(1901, 2000)

amo <- melt(amo, id.vars = "Year")

amo$Date <- as.Date(paste(amo$Year, amo$variable, "01", sep = "-"),
                         format = ("%Y-%b-%d"))
amo <- amo[-c(nrow),]
amo <- head(amo, -1) # Remove last row

amo.ts <- ts(amo$value, start=c(1856,1), end=c(2018,12), frequency=12)

library(forecast)
plot(ma(amo.ts, 240))
aw
