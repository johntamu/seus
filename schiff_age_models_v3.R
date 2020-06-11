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
core.dir <- '~/Documents/GitHub/data/clam/Cores'

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

# Note: This one is particularly important because we can also use the bomb spike

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

clam(core="stetson", type = 1, prob = 0.95, its = 1000, # Linear regression
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, thickness = 0.1, est=1, outlier=c(2:11,13:15,17:24,26:27), cmyr = TRUE, youngest = -55, bty = "o")

clam(core="stetson", type = 4, smooth = 0.8, prob = 0.95, its = 1000, # Spline model
    coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
    depths.file = TRUE, thickness = 0.1, est = 2)

clam(core="stetson", type = 1, smooth = 0.3, prob = 0.95, its = 1000, # Linear interpolation
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 2)

stetdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/stetson/stetson_polyn_regr_ages.txt')

# **********************
# Stetson-4904 BC1 PART 2
# ********************** 

stetson %>%
  select(distance..mm.) -> depths
write.table(depths, '~/Documents/GitHub/data/clam/Cores/stetson2/stetson2_depths.txt', 
            sep = '\t', row.names = FALSE, col.names = FALSE)

par(pty = "s", mfrow = c(1,2))
clam(core="stetson2", type = 2, smooth = 0.3, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, thickness = 0.1, est = 1, bty = "o", cmyr = TRUE, youngest = (1950-2005))

clam(core="stetson2", type = 1, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, thickness = 0.1, est = 1, bty = "o", cmyr = TRUE, youngest = (1950-2005))

mix.calibrationcurves(proportion=0.5, cc1 = "IntCal13.14C", cc2 = "Marine13.14C", name = "mixed.14C", dir = ".", offset = c(2,20), sep = "\t")

# stet2depths <- read.delim('~/Documents/GitHub/data/clam/Cores/stetson2/stetson2_smooth_spline_ages.txt')
stet2depths <- read.delim('~/Documents/GitHub/data/clam/Cores/stetson2/stetson2_interpolated_ages.txt')
stet2depths$rate <- stet2depths$acc.rate*1000

#-------------------------
# Stetson-4904 model
# comparison
#-------------------------

pryr.stetmodel1 %<a-% {jack.agemodel.linear %>%
    plot(d15n.vs.air ~ yr.bp, .,
         type = "l",
         bty = "n",
         col = alpha("black", 0.3),
         xlab = "Years BP",
         ylab = n,
         # main = "Age Model 1",
         xaxt = "n",
         xlim = c(500, 3500),
         ylim = c(7,10))
  jack.agemodel.linear %>%
    lines(rollmean(d15n.vs.air, 3, na.pad = TRUE) ~ yr.bp, ., col = "black")
}


# **********************
# Jacksonville-4907 BC1 
# **********************
jack %>%
  select(distance..mm.) -> depths
write.table(depths, '~/Documents/GitHub/data/clam/Cores/jack4907/jack4907_depths.txt', 
            sep = '\t', row.names = FALSE, col.names = FALSE)

par(pty = "s", mfrow=c(1,3))
clam(core="jack4907", type = 4, smooth = 0.5, prob = 0.99, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 5, cmyr = TRUE, bty = "o", mix.calibrationcurves(offset = c(-37,5)))


clam(core="jack4907", type = 1, prob = 0.95, its = 1000, # Use linear regression with a hiatus
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, outlier = c(4:9,11:16), thickness = 0.05, est = 1, cmyr = TRUE, 
     bty = "o") 

clam(core="jack4907", type = 2, prob = 0.99, its = 1000, # Use linear regression with a hiatus
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, thickness = 0.05, est = 1, cmyr = TRUE, hiatus = c(0.935, 3.5), # hiatus argument is depth in cm ***
     bty = "o", youngest = c(1)) 

clam(core="jackusgs", type = 2, prob = 0.99, its = 1000, # Use linear regression with a hiatus
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, thickness = 0.05, est = 1, cmyr = TRUE, hiatus = c(0.6, 2.8), # hiatus argument is depth in cm ***
     bty = "o",  youngest = c(1)) 

clam(core="jack4907", type = 2, prob = 0.95, its = 1000, # Use linear regression with a hiatus
     coredir = core.dir, cc = 4, ccdir=".", BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, thickness = 0.1, est = 7, cmyr = TRUE, hiatus = c(1.35), # 1.1 works; hiatus argument is depth in cm ***
     bty = "o", wghts = 0)

clam(core="jack4907", type = 2, prob = 0.95, its = 1000, # Use linear regression with a hiatus
     coredir = core.dir, cc = 4, ccdir=".", BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.05, est = 1, cmyr = TRUE, hiatus = c(1.1), # hiatus argument is depth in cm ***
     bty = "o", youngest = c(1))

mix.calibrationcurves(proportion=0.0, cc1 = "IntCal13.14C", cc2 = "Marine13.14C", name = "mixed.14C", dir = ".", offset = c(2,20), sep = "\t")
calibrate(1010, 35, cc=2, reservoir=c(2,20), prob = 0.95)
# One linear regression

clam(core="jack4907", type = 2, prob = 0.99, its = 1000, # Use linear regression with a hiatus
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 1, cmyr = TRUE,
     bty = "o", youngest = 500) 

# Smoothing spline with clam package

clam(core="jack4907", type = 4, smooth = 0.5, prob = 0.99, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 5, cmyr = TRUE)

# jackdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4907/jack4907_smooth_spline_ages.txt')
jackdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4907/jack4907_polyn_regr_ages.txt')
# jackdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4907/jack4907_interpolated_ages.txt')
jackdepths$rate <- jackdepths$acc.rate*1000

# --------------------------------
# Jacksonville-4907 
# model comparison
# --------------------------------

pryr.agemodel1 %<a-% {jack.agemodel.linear %>%
    plot(d15n.vs.air ~ yr.bp, .,
         type = "l",
         bty = "n",
         col = alpha("black", 0.3),
         xlab = "Years BP",
         ylab = n,
         # main = "Age Model 1",
         xaxt = "n",
         xlim = c(500, 3500),
         ylim = c(7,10))
  jack.agemodel.linear %>%
    lines(rollmean(d15n.vs.air, 3, na.pad = TRUE) ~ yr.bp, ., col = "black")
}
pryr.agemodel2 %<a-% {jack %>%
    plot(d15n.vs.air ~ yr.bp, .,
         type = "l",
         bty = "l",
         col = alpha("black", 0.3),
         xlab = "Years BP",
         ylab = n,
         # main = "Age Model 2",
         # xaxt = "n",
         xlim = c(500, 3500),
         ylim = c(7,10))
  axis(side = 2)
  jack %>%
    lines(rollmean(d15n.vs.air, 3, na.pad = TRUE) ~ yr.bp, ., col = "black")
}
pryr.agemodel3 %<a-% {jack.smoothspline %>%
    plot(d15n.vs.air ~ yr.bp, .,
         type = "l",
         bty = "l",
         col = alpha("black", 0.3),
         xlab = "Years BP",
         ylab = n,
         # main = "Age Model 1",
         # xaxt = "n",
         xlim = c(500, 3500),
         ylim = c(7,10))
  jack.smoothspline %>%
    lines(rollmean(d15n.vs.air, 3, na.pad = TRUE) ~ yr.bp, ., col = "black")
}

par.set <- par(pty = "s", mfrow=c(3,1), oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
               mar = c(1, 1, 0, 0)+0.3) # space for one row of text at ticks and to separate plots

agemodel1 <- ggplot() + # With ggplot2
  geom_line(data=jack.agemodel.linear, aes(x = yr.bp, y = d15n.vs.air), color = "gray", alpha = 0.6, size = 0.5) +
  geom_line(data=jack.agemodel.linear, aes(x = yr.bp,
                              y = rollmean(d15n.vs.air, 3, na.pad = TRUE)), color = "black", alpha = 0.99, size = 0.5) +
  ylab(n) +
  theme_classic() +
  # xlab(x) +
  xlim(500, 3500) +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(500, 3500)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

agemodel2 <- ggplot() + # With ggplot2
  geom_line(data=jack, aes(x = yr.bp, y = d15n.vs.air), color = "gray", alpha = 0.6, size = 0.5) +
  geom_line(data=jack, aes(x = yr.bp, y = rollmean(d15n.vs.air, 3, na.pad = TRUE)), color = "black", alpha = 0.99, size = 0.5) +
  ylab(n) +
  theme_classic() +
  xlab("Years BP") +
  # xlim(500, 3500) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(500, 3500)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(agemodel1), ggplotGrob(agemodel2), size = "last"))

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
par(pty = "s")
clam(core="sav4902", type = 1, smooth = 0.5, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, thickness = 0.1, outlier=c(2,5), cmyr = TRUE, bty ="o")

clam(core="sav4902", type = 2, smooth = 0.5, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 1, youngest = c(1), cmyr = TRUE, hiatus = c(20.5))

# savdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/sav4902/sav4902_polyn_regr_ages.txt')
savdepths <- read.delim('~/Documents/GitHub/data/clam/Cores/sav4902/sav4902_interpolated_ages.txt')
savdepths$rate = savdepths$acc.rate*1000

# **********************
# Jacksonville-4684 BC1 
# **********************

jack4684 %>%
  select(distance..mm.) -> depths
write.table(depths, '~/Documents/GitHub/data/clam/Cores/jack4684/jack4684_depths.txt', 
            sep = '\t', row.names = FALSE, col.names = FALSE)

clam(core="jack4684", type = 2, smooth = 0.4, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = TRUE,
     depths.file = TRUE, thickness = 0.1, est = 2, ignore = c(2:7), cmyr = TRUE, youngest = c(1), bty = "o")

# jack4684depths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4684/jack4684_smooth_spline_ages.txt')
jack4684depths <- read.delim('~/Documents/GitHub/data/clam/Cores/jack4684/jack4684_polyn_regr_ages.txt')

# *********************************
#
# Create some figures for Results 
# section of dissertation
#
# *********************************

# Linear regression figure

par(pty = "s", mfrow = c(2,2))
clam(core="stetson", type = 2, smooth = 0.6, prob = 0.99, its = 1000, # Linear regression
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 2, bty = "o")

clam(core="jack4907", type = 2, prob = 0.99, its = 1000, # Use linear regression with a hiatus
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 1, cmyr = TRUE,
     bty = "o") 

clam(core="sav4902", type = 2, smooth = 0.5, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 1, youngest = c(1), cmyr = TRUE, bty = "o")

clam(core="jack4684", type = 2, smooth = 0.4, prob = 0.95, its = 1000,
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.1, est = 2, ignore = c(2:7), youngest = c(1), bty = "o")

# Breaking into multiple regressions 

par(pty = "s", mfrow = c(1,1))
clam(core="jack4907", type = 2, prob = 0.99, its = 1000, # Use linear regression with a hiatus
     coredir = core.dir, cc = 2, BCAD = FALSE, depth = "mm", plotname = FALSE,
     depths.file = TRUE, thickness = 0.05, est = 1, cmyr = TRUE, hiatus = c(0.935, 3.5), # hiatus argument is depth in cm ***
     bty = "o") 




#' --------------------------------------------
#' PART 3: Bacon Age Models (with Bacon R package)
#' 
#' 
#' --------------------------------------------