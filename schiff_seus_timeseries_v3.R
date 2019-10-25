#' ------------
#' John Schiff
#' 08/04/2019
#' ------------
#' 
#' Note: Timeseries v2 was getting pretty hefty.
#' I am taking the code to load the data from the
#' previous file and importing it here.
#' I am then writing code to simply put in the Year AD 
#' of choice rather than all of them. This hopefully reduces the dependence
#' of running the age models.v2 file and the timeseries script file
#' together.
#' 



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
path1 <- '~/Documents/GitHub/data/schiff cn age models 09-04-2019.csv'
# path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff cn age models 09-04-2019.csv'
# path3 <- '/home/Desktop/data' # Linux computer link

load <- path1

seus.bulk <- read.csv(load, comment.char = '#')

# Data is now loaded

stetson <- seus.bulk %>% dplyr::filter(coral.id == 'stet-4904-bc1-d2')
jack <- seus.bulk %>% dplyr::filter(coral.id == 'jack-4907-bc1-d3')
jack2 <- seus.bulk %>% dplyr::filter(coral.id == 'jack-4907-bc1-d1')
jack4684 <- seus.bulk %>% dplyr::filter(coral.id == 'jack-4684-bc-unk')
jack4686 <- seus.bulk %>% dplyr::filter(coral.id == 'jack-4686-bc-d1-t1')
sav <- seus.bulk %>% dplyr::filter(coral.id == 'sav-4902-bc1-unk')

####################################
## Add year AD chronologies to    ##
## bulk records.                  ## 
####################################

# jack4907.ad <- binded$yrs NOTE: These linear ages need to be cleaned up in the age models v3 script file
# jack4907.usgs.ad <- binded2$yrs
# jack4686.ad <- NA
# jack4684.ad <- ad2
# stet4904.ad <- stet.linear.ad3 # Stetson age model with split chronologies
# sav4902.ad <- sav.predict.ad

# Comment out above and select below to include only clam ages
jack4907.ad <- 1950 - jackdepths$best
jack4907.usgs.ad <- NA
jack4686.ad <- NA
jack4684.ad <- 1950 - jack4684depths$best
jack4684.ad <- ad2 # 23 microns per year growth rate
# stet4904.ad <- 1950 - stetdepths$best
stet4904.ad <- 1950 - stet2depths$best # from the coarser resolution Stetson radiocarbon record
sav4902.ad <- 1950 - savdepths$best

# Add to data frames
jack$linear.ad <- jack4907.ad
stetson$linear.ad <- stet4904.ad
sav$linear.ad <- sav4902.ad
jack2$linear.ad <- jack4907.usgs.ad
jack4684$linear.ad <- jack4684.ad
jack4686$linear.ad <- jack4686.ad

string <- '~/Documents/GitHub/data/schiff_bulk_years_09-04-2019.csv'
  
df.bulk <- rbind(jack, jack2, jack4684, jack4686, stetson, sav)
write.csv(df.bulk, string, row.names = FALSE)

#---------------------
# Creating two separate age records
# One linear, the other uses multiple
# linear growth rates
# Good comparison
#---------------------

# Only after making sure you are importing the linear growth rate model
jack.agemodel.linear <- jack
jack.agemodel.linear$yr.ad <- jack4907.ad
jack.agemodel.linear$yr.bp <- 1950 - jack.agemodel.linear$yr.ad

jack$yr.bp <- 1950 - jack$linear.ad

jack.smoothspline <- jack
jack.smoothspline$yr.ad <- jack4907.ad
jack.smoothspline$yr.bp <- 1950 - jack.smoothspline$yr.ad