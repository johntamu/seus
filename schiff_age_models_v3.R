#' ------------------------------------------------------
#' Cleaned up age model code goes here when finalized
#' 08/14/2019
#' ------------------------------------------------------
#' 

#' --------------------------------------------
#' Load packages and the data
#' --------------------------------------------

library(clam)
library(dplyr)

path1 <- '~/Documents/GitHub/data/schiff radiocarbon 02-05-2019.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff radiocarbon 02-05-2019.csv'
path3 <- '/home/john/Desktop/data/schiff radiocarbon 02-05-2019.csv'

df <- read.csv(path1, header = TRUE)

radiocarbon <- df

##################################
## Create separate dataframes.  ##
##################################
r.jack <- radiocarbon %>% filter(Coral == 'jack-4907-bc1-d1')
r.jack2 <- radiocarbon %>% filter(Coral == 'jack-4684-bc1')
r.stet <- radiocarbon %>% filter(Coral == 'stet-4904-bc1-d5')
r.sav <- radiocarbon %>% filter(Coral == 'sav-4902-bc1')
r.stet2 <- radiocarbon %>% filter(Coral == 'stet-4904-bc1')

####################################
## Attach yrmin, yrmax estimates  ## 
## from clam CALIBRATION output   ##
####################################

# Jacksvonille-4907 BC1

# Savannah-4902 BC1

# Stetson-4904 BC1

# Jacksonville-4684 BC1



#' --------------------------------------------
#' PART 1: Linear Age Models
#' 
#' 
#' --------------------------------------------

# **********************
# Jacksonville-4907 BC1 
# **********************

# **********************
# Stetson-4904 BC1
# **********************

# **********************
# Savannah-4902 BC1
# **********************

# **********************
# Jacksonville-4684 BC1
# **********************

# Note: This one is particularly important because we can als use the bomb spike

################
# Growth rates #
################

gr1 <- 270/(2004-1975)
gr2 <- (682-270)/(1975-1957) # 23 microns per year
gr3 <- 682/(2004-1957) # 15 microns per year
gr4 <- as.numeric(1/lm.whole$coefficients[2])

print(gr1)
print(gr2)
print(gr3)
print(gr4)

m <- mean(c(gr1, gr2, gr3))
sd(c(gr1, gr2, gr3))

ad1 <- 2004 - ((jack4684$distance..mm.*1000)/gr3) # 15 microns per year, from rounding up average overall growth rate
ad2 <- 2004 - ((jack4684$distance..mm.*1000)/23) # 23 microns per year, this also fits with overall linear trend
ad3 <- 2004 - ((jack4684$distance..mm.*1000)/gr1) # 9 microns per year
ad4 <- 2004 - ((jack4684$distance..mm.*1000)/20) # Test, based on sheet from Nancy
ad5 <- 2004 - ((jack4684$distance..mm.*1000)/16) # from 'm'

t.df <- cbind(ad1, ad2, ad3, ad4, ad5)
t.df <- as.data.frame(t.df)


#' --------------------------------------------
#' PART 2: clam Age Models (with clam R package)
#' 
#' 
#' --------------------------------------------
#' 

# **********************
# Stetson-4904 BC1
# ********************** 
stetson %>%
  select(distance..mm.) -> depths
write.table(depths, '~/Documents/GitHub/data/clam/Cores/stetson/stetson_depths.txt', 
            sep = '\t', row.names = FALSE, col.names = FALSE)

clam(core="stetson", type = 1, smooth = 0.8, prob = 0.99, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 2)

stetdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/stetson/stetson_smooth_spline_ages.txt')

# **********************
# Stetson-4904 BC1 PART 2
# ********************** 

stetson %>%
  select(distance..mm.) -> depths
write.table(depths, '~/Documents/GitHub/data/clam/Cores/stetson2/stetson2_depths.txt', 
            sep = '\t', row.names = FALSE, col.names = FALSE)

clam(core="stetson2", type = 4, smooth = 0.2, prob = 0.99, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 1, bty = "o")

clam(core="stetson2", type = 1, prob = 0.99, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 1, bty = "o", cmyr = TRUE)

# stet2depths <- read.delim('~/Documents/GitHub/data/clam/Cores/stetson2/stetson2_smooth_spline_ages.txt')
stet2depths <- read.delim('~/Documents/GitHub/data/clam/Cores/stetson2/stetson2_interpolated_ages.txt')
stet2depths$rate <- stet2depths$acc.rate*1000

# **********************
# Jacksonville-4907 BC1 
# **********************
jack %>%
  select(distance..mm.) -> depths
write.table(depths, '~/Documents/GitHub/data/clam/Cores/jack4907/jack4907_depths.txt', 
            sep = '\t', row.names = FALSE, col.names = FALSE)

clam(core="jack4907", type = 4, smooth = 0.5, prob = 0.99, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 5, cmyr = TRUE)

clam(core="jack4907", type = 2, prob = 0.99, its = 1000, # Use linear regression with a hiatus
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, thickness = 0.1, est = 1, cmyr = TRUE, hiatus = c(0.95, 3.5), # hiatus argument is depth in cm ***
     bty = "o") 

# jackdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4907/jack4907_smooth_spline_ages.txt')
jackdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4907/jack4907_polyn_regr_ages.txt')
jackdepths$rate <- jackdepths$acc.rate*1000

# **********************
# Savannah-4907 BC1 
# **********************
sav %>%
  select(distance..mm.) -> depths
write.table(depths, '~/Documents/GitHub/data/clam/Cores/sav4902/sav4902_depths.txt', 
            sep = '\t', row.names = FALSE, col.names = FALSE)

clam(core="sav4902", type = 4, smooth = 0.2, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 1, youngest = c(1))

clam(core="sav4902", type = 4, smooth = 0.5, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 1, youngest = c(1), cmyr = TRUE)

savdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/sav4902/sav4902_smooth_spline_ages.txt')
savdepths$rate = savdepths$acc.rate*1000

# **********************
# Jacksonville-4684 BC1 
# **********************

jack4684 %>%
  select(distance..mm.) -> depths
write.table(depths, '~/Documents/GitHub/data/clam/Cores/jack4684/jack4684_depths.txt', 
            sep = '\t', row.names = FALSE, col.names = FALSE)

clam(core="jack4684", type = 2, smooth = 0.4, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 2, ignore = c(2:7), youngest = c(1))

# jack4684depths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4684/jack4684_smooth_spline_ages.txt')
jack4684depths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4684/jack4684_polyn_regr_ages.txt')

#' --------------------------------------------
#' PART 3: Bacon Age Models (with Bacon R package)
#' 
#' 
#' --------------------------------------------