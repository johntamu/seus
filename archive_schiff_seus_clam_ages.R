# ----------------------------------------------
# John Schiff
# October 26, 2018
# https://cran.r-project.org/web/packages/clam/index.html
# ----------------------------------------------

library(clam)

#################################################
# This script shows how to use the clam package #
# in R to do robust radiocarbon age models      #
# beyond the linear age calculations in the     #
# other script. E.g., smoothing spline          #
# across the data points                        #
#################################################

# Assign a folder where you will store your age model output

core.dir1 <- "~/Google Drive/projects/rproj/seus/data/clam/Cores"
core.dir2 <- "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores"

# ----------------------------------
# Age model for Jacksonville-4907
# ----------------------------------

clam(core="4907bc1john", type = 4, smooth = 0.6, prob = 0.95, its = 1000,
     coredir = core.dir1, cc = 2, ignore = c(), BCAD = FALSE, depth = "mm", proxies = TRUE)
# Sherwood et al (2014) used 0.6 smoothing, though did not explain why

# when proxies = TRUE, clam will look for a .csv file with the isotope values and their respective depths
# these depths are in centimeters (cm) by default

# ----------------------------------
# Age model for Savannah-4902
# ----------------------------------
sav.depth <- bulk.sav$distance..mm.
write.table(sav.depth, file = "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/4902usgs/4902usgs_depths.txt",
            row.names = FALSE, col.names = FALSE)

clam(core="4902usgs", type = 2, smooth = 0.6, prob = 0.95, its = 1000,
     coredir = core.dir2, cc = 2, ignore = c(), BCAD = FALSE, depth = "mm",
     depths.file = TRUE)
clam(core="4902usgs", type = 2, smooth = 0.6, prob = 0.95, its = 1000,
     coredir = core.dir1, cc = 2, ignore = c(), BCAD = FALSE, depth = "mm",
     depths.file = TRUE)
file <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/4902usgs/4902usgs_polyn_regr_ages.txt'
sav.depths <- read.delim(file, header = TRUE, sep = "\t", dec = ".")


# ----------------------------------
# Age model for Jacksonville-4684
# ----------------------------------
jack4684.depth <- bulk.jack4684$distance..mm.
write.table(youngjack.depth, file = "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/4684usgs/4684usgs_depths.txt",
            row.names = FALSE, col.names = FALSE)
clam(core="4684usgs", type = 2, smooth = 0.5, prob = 0.95, its = 1000,
     coredir = core.dir2, cc = 2, ignore = c(2:6), BCAD = FALSE, depth = "mm",
     depths.file = TRUE)
file <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/4684usgs/4684usgs_polyn_regr_ages.txt'
file <- '~/Google Drive/projects/rproj/seus/data/clam/Cores/4684usgs/4684usgs_polyn_regr_ages.txt'
jack4684.depths <- read.delim(file, header = TRUE, sep = "\t", dec = ".")

#' Doing something different here.
#' I'm going to re-import the radiocarbon data without
#' the Calib calibrations and redo them entirely using
#' the clam package. 
#' 
#' First, import the radiocarbon data and remove calibrated ages. 
#' Make sure the resulting data fram matches what clam wants to
#' work with. 

# Borrowing from Rscript for linear age models
# path1 <- '~/Google Drive/projects/rproj/seus/data/09-05-2018-schiff radiocarbon.csv'
path1 <- '~/Google Drive/projects/rproj/seus/data/schiff radiocarbon 02-05-2019.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff radiocarbon 02-05-2019.csv'
radiocarbon <- read.csv(path2)
radiocarbon <- read.csv(path1)
# radiocarbon <- radiocarbon[,c(1:13)]

dir1 <- "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/"
dir2 <- "~/Google Drive/Google Drive/projects/rproj/seus/data/clam/Cores/"

# Create new data frame fitting clam() syntax
df <- data.frame(lab_ID = radiocarbon$ID.2,
                 C14_age = radiocarbon$X14C.Age,
                 cal_age = NA,
                 error = radiocarbon$error.1,
                 reservoir = NA,
                 depth = radiocarbon$Distance..um./1000) # Convert micron to mm

write.csv(df, file = "schiff radiocarbon data for clam.csv", row.names = FALSE)
df <- read.csv("schiff radiocarbon data for clam.csv")

## --------------------------------------- ##
##
## Clam ages for Steson-4904 BC1
##
## --------------------------------------- ##

stetson <- df[c(63:90),]
stetson$thickness <- 0.035 # Estimated thickness of each peel

write.csv(stetson, file = "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson.csv", row.names = FALSE)
write.csv(stetson, file = "~/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson.csv", row.names = FALSE)

stetson.depth <- bulk.stet$distance..mm.

write.table(stetson.depth, file = "C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson_depths.txt",
            row.names = FALSE, col.names = FALSE)
write.table(stetson.depth, file = "~/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson_depths.txt",
            row.names = FALSE, col.names = FALSE)

clam(core="stetson", type = 2, smooth = 0.9, prob = 0.95, its = 1000,
     coredir = core.dir1, cc = 2, ignore = c(2:5), BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = FALSE)
file <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson_polyn_regr_ages.txt'
file <- '~/Google Drive/projects/rproj/seus/data/clam/Cores/stetson/stetson_polyn_regr_ages.txt'
stet.depths <- read.delim(file, header = TRUE, sep = "\t", dec = ".")

##' ---------------------------------------
##' Clam ages for Steson-4904 BC1
##' But this time with the 4 data point
##' radiocarbon dataset.
##' ---------------------------------------

clam(core="stetson2", type = 2, smooth = 0.9, prob = 0.95, its = 1000,
     coredir = core.dir1, cc = 2, ignore = c(), BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE)

write.table(stetson.depth, file = "~/Google Drive/projects/rproj/seus/data/clam/Cores/stetson2/stetson2_depths.txt",
            row.names = FALSE, col.names = FALSE)
file <- '~/Google Drive/projects/rproj/seus/data/clam/Cores/stetson2/stetson2_polyn_regr_ages.txt'
stet2.depths <- read.delim(file, header = TRUE, sep = "\t", dec = ".")
# Now we just bind the depths vector to the stetson bulk isotope dataframe!

#' Note to self: Investigate further the radiocarbon data around the Bond event
#' There could be some indication that the Bond event influenced the radiocarbon age
#' by changing circulation of older/younger water in the Gulf Stream
#' around this time. A bit of a stretch, but consider it for later.
#' 
#' Could be an addition to the Iodine chapter?
#' 





