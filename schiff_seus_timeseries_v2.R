#' ------------
#' John Schiff
#' 02/26/2019
#' ------------

##################################         
## Importing and plotting       ##
## bulk time series, including  ##
## calculating moving average   ## 
## and SSA reconstruction.      ##
##################################

####################################
## Before using this script file, ##
## make sure that you have        ##
## established age models and     ##
## they are in your global        ##
## environment. You will attach   ##
## them to the bulk records here. ##
####################################

################
## Libraries  ##
################
library(dplyr) # easy data manipulation
library(Rssa) # for SSA reconstruction
library(lattice)
library(ggplot2)
library(forecast) # for moving averages
library(data.table) # alternative for data manipulation

# Establish a pathfile for the data
pathmac1 <- '~/Google Drive/projects/rproj/seus/data/schiff cn age models 02-05-2019.csv'
pathwin1 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff cn age models 02-05-2019.csv'

# Load the data
seus.bulk <- read.csv(pathmac1) # If on MacOS, don't alter this dataframe
seus.bulk <- read.csv(pathwin1) # If on Windows, don't alter this dataframe

#' Note on filtering/selecting rows and columns
#' It is useful to know both dplyr and data.table() 
#' packages. I will use dplyr for most of my code
#' but cheat sheets for data.table() are available
#' if you prefer that.

stetson <- seus.bulk %>% filter(coral.id == 'stet-4904-bc1-d2')
jack <- seus.bulk %>% filter(coral.id == 'jack-4907-bc1-d3')
jack2 <- seus.bulk %>% filter(coral.id == 'jack-4907-bc1-d1')
jack4684 <- seus.bulk %>% filter(coral.id == 'jack-4684-bc-unk')
jack4686 <- seus.bulk %>% filter(coral.id == 'jack-4686-bc-d1-t1')
sav <- seus.bulk %>% filter(coral.id == 'sav-4902-bc1-unk')

####################################
## Add linear growth rate ages to ##
## bulk records.                  ## 
####################################

jack$linear.ad <- binded$yrs # Add the combined columns for two growth rates you calculated in the age model script file
stetson$linear.ad <- stet.linear.ad1 # Add the calculated years based on the linear growth rate
sav$linear.ad <- linear.ad2
jack2$linear.ad <- binded2$yrs
jack4684$linear.ad <- ad5
jack4686$linear.ad <- NA

jack$linear.ad2 <- yrs # From age models script. This is based on the overall growth rate for the coral, approximately 
jack2$linear.ad2 <- NA
stetson$linear.ad2 <- stet.linear.ad2 # This is to keep the number of columns the same so we can bind them later
sav$linear.ad2 <- NA
jack4684$linear.ad2 <- NA
jack4686$linear.ad2 <- NA

##########################################
## Add secondary (i.e., min and max)    ##
## growth rate ages to the bulk records ##
##########################################



##########################################
## Calculate moving averages using the  ##
## forecast package                     ##
##########################################

# 3 pt moving average
jack$d15n.3pt <- as.numeric(forecast::ma(jack$d15n.vs.air, order = 3, centre = TRUE)) # Forecast::ma converts them to ts objects, which can't be binded to a dataframe
jack2$d15n.3pt <- as.numeric(forecast::ma(jack2$d15n.vs.air, order = 3, centre = TRUE))
jack4684$d15n.3pt <- as.numeric(forecast::ma(jack4684$d15n.vs.air, order = 3, centre = TRUE))
stetson$d15n.3pt <- as.numeric(forecast::ma(stetson$d15n.vs.air, order = 3, centre = TRUE))
sav$d15n.3pt <- as.numeric(forecast::ma(sav$d15n.vs.air, order = 3, centre = TRUE))
jack4686$d15n.3pt <- NA

jack$d13c.3pt <- as.numeric(forecast::ma(jack$d13c.vs.vpdb, order = 3, centre = TRUE))
jack2$d13c.3pt <- as.numeric(forecast::ma(jack2$d13c.vs.vpdb, order = 3, centre = TRUE))
jack4684$d13c.3pt <- as.numeric(forecast::ma(jack4684$d13c.vs.vpdb, order = 3, centre = TRUE))
stetson$d13c.3pt <- as.numeric(forecast::ma(stetson$d13c.vs.vpdb, order = 3, centre = TRUE))
sav$d13c.3pt <- as.numeric(forecast::ma(sav$d13c.vs.vpdb, order = 3, centre = TRUE))
jack4686$d13c.3pt <- NA

# 5 point moving average
jack$d15n.5pt <- as.numeric(forecast::ma(jack$d15n.vs.air, order = 5, centre = TRUE))
jack2$d15n.5pt <- as.numeric(forecast::ma(jack2$d15n.vs.air, order = 5, centre = TRUE))
jack4684$d15n.5pt <- NA
stetson$d15n.5pt <- as.numeric(forecast::ma(stetson$d15n.vs.air, order = 5, centre = TRUE))
sav$d15n.5pt <- as.numeric(forecast::ma(sav$d15n.vs.air, order = 5, centre = TRUE))
jack4686$d15n.5pt <- NA

jack$d13c.5pt <- as.numeric(forecast::ma(jack$d13c.vs.vpdb, order = 5, centre = TRUE))
jack2$d13c.5pt <- as.numeric(forecast::ma(jack2$d13c.vs.vpdb, order = 5, centre = TRUE))
jack4684$d13c.5pt <- NA
stetson$d13c.5pt <- as.numeric(forecast::ma(stetson$d13c.vs.vpdb, order = 5, centre = TRUE))
sav$d13c.5pt <- as.numeric(forecast::ma(sav$d13c.vs.vpdb, order = 5, centre = TRUE))
jack4686$d13c.5pt <- NA

# 12 point moving average
jack$d15n.12pt <- as.numeric(forecast::ma(jack$d15n.vs.air, order = 12, centre = TRUE))
jack2$d15n.12pt <- as.numeric(forecast::ma(jack2$d15n.vs.air, order = 12, centre = TRUE))
jack4684$d15n.12pt <- NA
stetson$d15n.12pt <- as.numeric(forecast::ma(stetson$d15n.vs.air, order = 12, centre = TRUE))
sav$d15n.12pt <- as.numeric(forecast::ma(sav$d15n.vs.air, order = 12, centre = TRUE))
jack4686$d15n.12pt <- NA

jack$d13c.12pt <- as.numeric(forecast::ma(jack$d13c.vs.vpdb, order = 12, centre = TRUE))
jack2$d13c.12pt <- as.numeric(forecast::ma(jack2$d13c.vs.vpdb, order = 12, centre = TRUE))
jack4684$d13c.12pt <- NA
stetson$d13c.12pt <- as.numeric(forecast::ma(stetson$d13c.vs.vpdb, order = 12, centre = TRUE))
sav$d13c.12pt <- as.numeric(forecast::ma(sav$d13c.vs.vpdb, order = 12, centre = TRUE))
jack4686$d13c.12pt <- NA

######################################################################
## Bind the data together to create a new data frame with the years ##
######################################################################

df.bulk <- rbind(jack, jack2, jack4684, jack4686, stetson, sav)
write.csv(df.bulk, '~/Google Drive/projects/rproj/seus/data/schiff_bulk_years_06-18-2019.csv')
write.csv(df.bulk, 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_bulk_years_06-18-2019.csv')

#' ----------------------------------------------
#' Add clam linear regression or smooth spline
#' ages to bulk records if you did those.
#' ----------------------------------------------

################################
## PART 2: SSA reconstruction ##
################################
library(Rssa)

#' Procedure:
#' 1. Use ssa() function and plot(ssa) to identify pairs
#' and lone points to use in reconstruction
#' 2. recon$F1 will generally be the trend; isolate this
#' by then analyzing the residuals without recon$F1
#' 
#' 

##########################
## Do SSA on each coral ##
##########################

######################
## Nitrogen records ##
######################

ssa.jack <- ssa(df.jack$d15n, L = 30)
ssa.stet <- ssa(df.stet$d15n, L = 30) # Need to remove two values from linear.ad column because the isotope record has two NAs
ssa.sav <- ssa(df.sav$d15n, L = 30)
ssa.jack2 <- ssa(jack4684$d15n.vs.air, L = 30)

######################
## Carbon records   ##
######################

cssa.jack <- ssa(df.jack$d13c, L = 30)
cssa.stet <- ssa(df.stet$d13c, L = 30) # Need to remove two values from linear.ad column because the isotope record has two NAs
cssa.sav <- ssa(df.sav$d13c, L = 20)

################################
## Looking at SSA components  ##
################################

################
# Jacksonville #
################

plot(ssa.jack)
plot(ssa.jack, "paired")
plot(ssa.jack, "series", groups = 1:20)
jack.recon <- reconstruct(ssa.jack, groups = list(1,c(2,3),4:7))
plot(jack.recon, plot.method = "xyplot", type = "cumsum")

# Test plot
plot(d15n ~ linear.ad, df.jack, type = "l", col = alpha("black", 0.3))
lines(d15n.5pt ~ linear.ad, df.jack, col = "blue")
lines(jack.recon$F1 + jack.recon$F2 ~ df.jack$linear.ad, col = "red")

############
# Savannah #
############

plot(ssa.sav)
plot(ssa.sav, "paired")
plot(ssa.sav, "series", groups = 1:20)
sav.recon <- reconstruct(ssa.sav, groups = list(1,2,3,4:7))
plot(sav.recon, plot.method = "xyplot", type = "cumsum")

# Test plot
plot(sav.recon$F1 + sav.recon$F2 + sav.recon$F3 + sav.recon$F4 ~ linear.ad, sav, type = "l", ylim = c(7, 11))
lines(jack.recon$F1 + jack.recon$F2 + jack.recon$F3 ~ df.jack$linear.ad, col = "red")


##############################
## Suess effect correction  ##
##############################

#' Also relevant to bulk d13C data
#' McMahon (2015) used this method to correct
#' for the Suess effect in their data.
#' Two references: Francey et al, 1999; Quay et al, 2013
#' 
#' 0.16permil per decade since 1960
#' 0.05permil per decade between 1860 - 1960
#' 
#' One method: Make a for-loop to write a vector with corrected bulk d13C values
#' or just do it in Excel and re-import the corrected d13C values as a vector
#'
df.stet %>%
  select(linear.ad2, d13c) %>%
  write.csv(., file = "stet_bulk_c2.csv", row.names = FALSE)

new_df <- read.csv("~/Google Drive/projects/rproj/seus/stet_bulk_c_corrected.csv")
df.stet$d13c.corrected <- new_df$d13c.suess.corrected


########################################
## Statistical analyses on bulk data  ##
########################################

#' First separate the different specimens
#' 

df.bulk %>%
  filter(coral.id == "jack-4907-bc1-d3"
         | coral.id == "sav-4902-bc1-unk"
         | coral.id == "stet-4904-bc1-d2") -> new.bulk

n.aov <- aov(d15n ~ coral.id, new.bulk)
c.aov <- aov(d13c ~ coral.id, new.bulk)
summary(n.aov)
summary(c.aov)

TukeyHSD(n.aov)
TukeyHSD(c.aov)

#' Tests with pre- and post-human activity d15N datasets
#' in Stetson-4904 BC1
#' 

df.stet %>%
  filter(linear.ad2 > 1850) -> stet.post
stet.post$V1 <- 'post'
df.stet %>%
  filter(linear.ad2 < 1850 & linear.ad2 > 1700) -> stet.pre
stet.pre$V1 <- 'pre'
V1 <- rbind(stet.post, stet.pre)

t.aov <- aov(d15n ~ V1, V1)
summary(t.aov)

summary(lm(d15n ~ linear.ad2, stet.post))

t.test(stet.post$d15n, stet.pre$d15n)

df.stet %>%
  filter(linear.ad2 > 800 & linear.ad2 < 950) -> spike
summary(lm(d15n ~ d13c, spike))
