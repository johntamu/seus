#' ------------------------------------------
#' Note: Make sure to have data loaded from
#' figures script file v2.
#' ------------------------------------------
#' 08/03/2019
#' 

# Required packages
library(maps)
library(ggplot2)
library(marmap)
library(oceanmap)
library(dplyr)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(reshape2)
library(Hmisc)
library(grid)
library(zoo) 
library(forecast)

##################################
## Load df.bulk from timeseries ##
## script file                  ##
##################################
path1 <- '~/Documents/GitHub/data/schiff_bulk_years_08-05-2019.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_bulk_years_08-05-2019.csv'
path3 <- '/home/john/Desktop/data/schiff_bulk_years_08-05-2019.csv'

path <- path1

df.bulk <- read.csv(path)
View(df.bulk)

colnames(df.bulk)[names(df.bulk) == "distance..mm."] <- "distance" # Rename some columns for easier coding
colnames(df.bulk)[names(df.bulk) == "d15n.vs.air"] <- "d15n"
colnames(df.bulk)[names(df.bulk) == "d13c.vs.vpdb"] <- "d13c"

df.jack <- df.bulk %>% filter(coral.id == 'jack-4907-bc1-d3')
df.sav <- df.bulk %>% filter(coral.id == 'sav-4902-bc1-unk')
df.stet <- df.bulk %>% filter(coral.id == 'stet-4904-bc1-d2')
df.jack2 <- df.bulk %>% filter(coral.id == 'jack-4907-bc1-d1')
df.jack4684 <- df.bulk %>% filter(coral.id == 'jack-4684-bc-unk')
df.jack4686 <- df.bulk %>% filter(coral.id == 'jack-4686-bc-d1-t1')

x <- 'Year CE'
phe <- expression({delta}^15*"N"[" Phe"]*" (\u2030)")
n <- expression(delta^{15}*"N (\u2030)")
c <- expression(delta^{13}*"C (\u2030)")
eaa.neaa <- c('Phe', 'Thr', 'Ile', 'Leu', 'Val', 'Asx', 'Glx', 'Pro', 'Ala', 'Ser', 'Gly') # For Essential/Non-Essential ordering
tr.srcaa <- c('Glu', 'Asp', 'Ala', 'Ile', 'Leu', 'Pro', 'Val', 'Gly', 'Ser', 'Lys', 'Tyr', 'Phe', 'Thr') # For Trophic/Source AA ordering

# Below is a general function to generate a continuous color palette

#' Figure: 
#' --------------------------------------------------------------------------------------
#' Bathymetric map of study area
#' --------------------------------------------------------------------------------------

# get bathymetry data
b = getNOAA.bathy(lon1=-95, lon2=-65,
                  lat1=15, lat2=40, 
                  resolution=1)

## Querying NOAA database ...
## This may take seconds to minutes, depending on grid size
## Building bathy matrix ...

# make a simple track line
# lin = data.frame(
#   lon = c(-65.17536, -65.37423, -65.64541, -66.06122, -66.15161),  
#   lat = c(43.30837, 42.94679, 42.87448, 42.92871, 42.72985)
# )

# make a few points
pts = data.frame(
  lon = c(-079.6616667, -079.6416667, -079.1271667, -077.609956, -079.6515000),
  lat = c(30.51333333, 30.8012600, 31.7040000, 31.8441917, 30.5021667)
)

# # build a polygon (in this case the 'Roseway Basin Area To Be Avoided')
# ply = data.frame(
#   lon = c(-64.916667,-64.983333,-65.516667, -66.083333),
#   lat = c(43.266667,  42.783333, 42.65, 42.866667)
# )

# Import and load libraries
library(oce)
library(ocedata)
data("coastlineWorldFine")

# convert bathymetry
bathyLon = as.numeric(rownames(b))
bathyLat = as.numeric(colnames(b))
bathyZ = as.numeric(b)
dim(bathyZ) = dim(b)

# define plotting region
mlon = mean(pts$lon)
mlat = mean(pts$lat)
span = 1200
lonlim = c(-90, -70)
latlim = c(20, 35)

# plot coastline (with projection)
plot(coastlineWorldFine, clon = mlon, clat = mlat, span = span, 
     projection="+proj=merc", col = 'lightgrey')

# plot bathymetry
mapContour(bathyLon,bathyLat,bathyZ,
           levels = c(-500, -1000, -1500, -2000, -2500, -3000, -3500, -4000, -4500, -5000),
           # lwd = c(1, 1, 2, 2, 3),
           # lty = c(3, 1, 3, 1, 3),
           col = 'darkgray')

# # add depth legend
# legend("bottomright", seg.len = 3, cex = 0.7,
#        lwd = c(1, 1, 2, 2, 3),
#        lty = c(3, 1, 3, 1, 3),
#        legend = c("50", "100", "150", "200", "250"),
#        col = 'darkgray', title = "Depth [m]", bg = "white")

# add map data
mapPoints(longitude = pts$lon, latitude = pts$lat, pch = 21, col = 'black', bg = 'red')
mapLines(longitude = lin$lon, latitude = lin$lat, col = 'blue')

#' Figure: 
#' --------------------------------------------------------------------------------------
#' Chlorophyll map with gridded data
#' --------------------------------------------------------------------------------------

# I actually have the code to do this in Python with NetCDF data, so I do it there.



#' Figure: 
#' --------------------------------------------------------------------------------------
#' Linear age models (radiocarbon)
#' --------------------------------------------------------------------------------------

par(pty = "s", mfrow = c(2,2), mar=c(4,5,1,1))
plot(X14C.Age ~ Distance..um., r.jack,
     pch = 22, bg = '#ffeda0', col = "black", 
     ylab = expression({Delta}^14*"C Age (CRA)"), xlab = 'Distance from edge (um)')
# abline(mod1)
plot(X14C.Age ~ Distance..um., r.sav,
     pch = 23, bg = '#feb24c', col = "black", 
     ylab = expression({Delta}^14*"C Age (CRA)"), xlab = 'Distance from edge (um)')
# abline(mod2)
plot(X14C.Age ~ Distance..um., r.stet,
     pch = 21, bg = '#f03b20', col = "black",
     ylab = expression({Delta}^14*"C Age (CRA)"), xlab = 'Distance from edge (um)')
# abline(mod3)
r.jack$mm <- r.jack$Distance..um./1000
r.jack2$mm <- r.jack2$Distance..um./1000
r.stet$mm <- r.stet$Distance..um./1000
r.sav$mm <- r.sav$Distance..um./1000

par(pty = "s", mfrow = c(2,2), mar=c(4,5,1,1))
plot(mean ~ mm, r.jack,
     type = "o", bg = 'gray', col = "black", 
     ylab = "Years BP", xlab = 'Distance from edge (um)')
abline(jlm1.mm, lty = "dashed")
abline(jlm2.mm, lty = "dashed")
plot(X14C.Age ~ mm, whole,
     type = "o", bg = 'gray', col = "black", 
     ylab = "Years BP", xlab = 'Distance from edge (um)')
abline(lm.whole.mm, lty = "dashed")
plot(mean ~ mm, r.sav,
     type = "o", bg = 'gray', col = "black", 
     ylab = "Years BP", xlab = 'Distance from edge (um)')
abline(lm.sav.mm, lty = "dashed")
plot(mean ~ mm, r.stet,
     type = "o", bg = 'gray', col = "black",
     ylab = "Years BP", xlab = 'Distance from edge (um)')
abline(s1.mm, lty = "dashed")
abline(s2.mm, lty = "dashed")

#' Figure: 
#' --------------------------------------------------------------------------------------
#' Bomb spike age model for Jack-4684 BC1
#' --------------------------------------------------------------------------------------



#' Figure: 
#' --------------------------------------------------------------------------------------
#' Bulk data - overall
#' --------------------------------------------------------------------------------------

mod1 <- lm(d15n ~ d13c, df.jack)
mod2 <- lm(d15n ~ d13c, df.sav)
mod3 <- lm(d15n ~ d13c, df.stet)

summary(mod)
summary(mod2)
summary(mod3)

par(mfrow = c(1,3))
plot(d15n ~ d13c, df.jack,
     pch = 22, bg = '#ffeda0', col = "black", 
     ylab = n, xlab = c)
abline(mod1)
plot(d15n ~ d13c, df.sav,
     pch = 23, bg = '#feb24c', col = "black", 
     ylab = n, xlab = c)
abline(mod2)
plot(d15n ~ d13c, df.stet,
     pch = 21, bg = '#f03b20', col = "black",
     ylab = n, xlab = c)
abline(mod3)

#' Figure: 
#' --------------------------------------------------------------------------------------
#' Bulk data - testing reproducibility
#' --------------------------------------------------------------------------------------

# Jack-4907 bulk data vs distance from edge

par(pty = 's', mfrow = c(1,2), mar=c(5,6,4,1)+.05)
plot(d15n ~ distance, df.jack,
     xlab = 'Distance from edge (mm)',
     ylab = expression(delta^{15}*"N"),
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3))
lines(d15n ~ distance, df.jack2,
      col = alpha("black", 0.75), lwd = 1.5)

plot(d13c ~ distance, df.jack,
     xlab = 'Distance from edge (mm)',
     ylab = expression(delta^{13}*"C"),
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3))
lines(d13c ~ distance, df.jack2,
      col = alpha("black", 0.75), lwd = 1.5)

par(pty = "s")
plot(d15n ~ linear.ad, df.jack,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.0))
lines(forecast::ma(df.jack$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = alpha("black", 0.4), lwd = 1.5)
lines(forecast::ma(df.jack2$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack2,
      col = "#1f78b4", lwd = 1.5)
points(forecast::ma(df.jack2$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack2,
      col = "#1f78b4", lwd = 1.5)

par(pty = "s")
plot(d13c ~ linear.ad, df.jack,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.0))
lines(forecast::ma(df.jack$d13c, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = alpha("black", 0.4), lwd = 1.5)
lines(forecast::ma(df.jack2$d13c, order = 2, centre = TRUE) ~ linear.ad, df.jack2,
      col = "#1f78b4", lwd = 1.5)
points(forecast::ma(df.jack2$d13c, order = 2, centre = TRUE) ~ linear.ad, df.jack2,
      col = "#1f78b4", lwd = 1.5)

# For now, use the figure already generated and copy paste that into the file

#' Figure: Recent history
#' --------------------------------------------------------------------------------------
#' Bulk data - Stet4904 and Jack4684 d15N
#' --------------------------------------------------------------------------------------

par(pty = "s")
plot(d15n ~ linear.ad, df.stet,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     xlim = c(1500,2005),
     col = alpha("black", 0.0))
points(forecast::ma(df.stet$d15n, order = 1, centre = TRUE) ~ linear.ad, df.stet,
      col = "#33a02c", lwd = 1.5)
lines(forecast::ma(df.jack4684$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack4684,
      col = "#1f78b4", lwd = 1.5)

#' Figure: Recent history, bulk
#' --------------------------------------------------------------------------------------
#' Bulk data - Stet-4904 entire record
#' --------------------------------------------------------------------------------------

# Nitrogen
par(pty = "s", mfrow = c(2,2), mar=c(4,5,1,1))
plot(d15n ~ distance, df.stet,
     xlab = 'Distance from edge (mm)',
     ylab = n,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend(7,2, box.lty = 0, legend = "Stetson-4904-BC1 Disk 1", bg = NULL)
plot(d15n ~ distance, df.jack,
     xlab = 'Distance from edge (mm)',
     ylab = n,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Jacksonville-4907-BC1 Disk 3", bg = NULL)
plot(d15n ~ distance, df.jack4684,
     xlab = 'Distance from edge (mm)',
     ylab = n,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Jacksonville-4684-BC1 Disk 1", bg = NULL)
plot(d15n ~ distance, df.sav,
     xlab = 'Distance from edge (mm)',
     ylab = n,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Savannah Banks-BC1 Base 1", bg = NULL)

# Carbon
par(pty = "s", mfrow = c(2,2), mar=c(4,5,1,1))
plot(d13c ~ distance, df.stet,
     xlab = 'Distance from edge (mm)',
     ylab = c,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend(7,2, box.lty = 0, legend = "Stetson-4904-BC1 Disk 1", bg = NULL)
plot(d13c ~ distance, df.jack,
     xlab = 'Distance from edge (mm)',
     ylab = c,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Jacksonville-4907-BC1 Disk 3", bg = NULL)
df.jack4684 %>%
  filter(d13c )
plot(d13c ~ distance, df.jack4684,
     xlab = 'Distance from edge (mm)',
     ylab = c,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Jacksonville-4684-BC1 Disk 1", bg = NULL)
plot(d13c ~ distance, df.sav,
     xlab = 'Distance from edge (mm)',
     ylab = c,
     type = "l",
     cex = 0.75,
     col = alpha("#0072B2", 0.99))
# text(2,7, labels = "A", cex = 2.25)
# legend("bottomleft", box.lty = 0, legend = "Savannah Banks-BC1 Base 1", bg = NULL)

plot(d13c ~ linear.ad, df.stet,
     xlab = x,
     ylab = c,
     type = "l",
     cex = 0.5,
     # xlim = c(1500,2005),
     col = alpha("black", 0.3))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
lines(forecast::ma(df.stet$d13c, order = 3, centre = TRUE) ~ linear.ad, df.stet,
      col = "#33a02c", lwd = 1.5)

plot(d15n ~ linear.ad, df.stet,
     xlab = x,
     ylab = c,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
lines(forecast::ma(df.stet$d15n, order = 3, centre = TRUE) ~ linear.ad, df.stet,
      col = "#33a02c", lwd = 1.5)

obj2 <- xyplot(forecast::ma(df.stet$d15n, order = 3, centre = TRUE) ~ linear.ad, data = df.stet, type = "l", lwd = 1.5,
               ylab = n, xlab = x, col = '#081d58', ylim = c(5, 10))
obj1 <- xyplot(forecast::ma(df.stet$d13c, order = 3, centre = TRUE) ~ linear.ad, data = df.stet, type = "l", lwd = 1.5,
               xlab = x, ylab = c, col = '#225ea8', ylim = c(-18, -10))
doubleYScale(obj2, obj1, add.ylab2 = TRUE, style1= NULL, style2 = NULL)

splot1 <- df.stet %>%
  dplyr::select(linear.ad, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.65) +
  geom_line(aes(y=rollmean(d15n, 3, na.pad = TRUE)), size = 0.65, color = '#225ea8') +
  
  geom_vline(xintercept = 950, lty = 'dotted') +
  geom_vline(xintercept = 1250, lty = 'dotted') +
  geom_vline(xintercept = 1550, lty = 'dotted') +
  geom_vline(xintercept = 1850, lty = 'dotted') +
  
  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 3) +
  annotate("text", x = 1700, y = 7, label = "Little Ice Age", size = 3) +
  
  ylab(n) +
  xlab(NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

splot2 <- df.stet %>%
  dplyr::select(linear.ad, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.65) +
  geom_line(aes(y=rollmean(d13c, 3, na.pad = TRUE)), size = 0.65, color = '#081d58') +
  
  geom_vline(xintercept = 950, lty = 'dotted') +
  geom_vline(xintercept = 1250, lty = 'dotted') +
  geom_vline(xintercept = 1550, lty = 'dotted') +
  geom_vline(xintercept = 1850, lty = 'dotted') +
  
  ylab(c) +
  theme_classic() +
  xlab(x) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(splot1), ggplotGrob(splot2), size = "last"))

#' Figure: Recent history, bulk
#' --------------------------------------------------------------------------------------
#' Stet4904 d13C vs d15N XY plot during certain periods
#' --------------------------------------------------------------------------------------

df.stet %>%
  filter(linear.ad < 1100) -> mwp

plot(d13c ~ d15n, mwp)
abline(lm(d13c ~ d15n, mwp))
summary(lm(d13c ~ d15n, mwp))

#' Figure: Recent history, bulk
#' --------------------------------------------------------------------------------------
#' Bulk data - Jack4686 vs distance
#' --------------------------------------------------------------------------------------

plot(d15n ~ distance, df.jack4686,
     xlab = "Distance from edge (mm)",
     ylab = n,
     type = "o",
     cex = 0.5)
# xlim = c(1500,2005),
# col = alpha("black", 0.3))

#' Figure: Ancient bulk values
#' --------------------------------------------------------------------------------------
#' Savannah-4902 by itself
#' --------------------------------------------------------------------------------------

savplot1 <- df.sav %>%
  select(linear.ad, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(d15n, 3, na.pad = TRUE)), size = 0.75, color = '#02818a') +
  # geom_point(color = '#31a354', shape = 21) +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -900, lty = 'longdash') +
  geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +
  
  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 2.5) +
  annotate("text", x = -600, y = 7, label = "Iron Age Cold Epoch", size = 4) +
  annotate("text", x = 75, y = 7, label = "Roman Warm Period", size = 3) +
  
  ylab(n) +
  xlab(NULL) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(-300,1300)) +
  theme_classic() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

savplot2 <- df.sav %>%
  select(linear.ad, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(d13c, 3, na.pad = TRUE)), size = 0.75, color = '#3690c0') +
  # geom_point(color = '#addd8e', shape = 21) +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -900, lty = 'longdash') +
  geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +
  
  ylab(c) +
  theme_classic() +
  xlab(x) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(-300,1300)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(savplot1), ggplotGrob(savplot2), size = "last"))




#' Figure: Ancient bulk values
#' --------------------------------------------------------------------------------------
#' Jack-4907 by itself
#' --------------------------------------------------------------------------------------

plot(d13c ~ linear.ad, df.jack,
     xlab = x,
     ylab = c,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3))
abline(v = -250, col = "black", lty = "dashed")
abline(v = 400, col = "black", lty = "dashed")
abline(v = -900, col = alpha("black", 0.75), lty = "longdash")
abline(v = -300, col = alpha("black", 0.75), lty = "longdash")
lines(forecast::ma(df.jack$d13c, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5)

plot(d15n ~ linear.ad, df.jack,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.3))
abline(v = -250, col = "black", lty = "dashed")
abline(v = 400, col = "black", lty = "dashed")
abline(v = -900, col = alpha("black", 0.75), lty = "longdash")
abline(v = -300, col = alpha("black", 0.75), lty = "longdash")
lines(forecast::ma(df.jack$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5)

obj2 <- xyplot(forecast::ma(df.jack$d15n, order = 3, centre = TRUE) ~ linear.ad, data = df.jack, type = "l", lwd = 1.5,
               ylab = n, xlab = x, col = 'red', ylim = c(5, 10))
obj1 <- xyplot(forecast::ma(df.jack$d13c, order = 3, centre = TRUE) ~ linear.ad, data = df.jack, type = "l", lwd = 1.5,
               xlab = x, ylab = c, col = 'blue', ylim = c(-18, -10))
doubleYScale(obj2, obj1, add.ylab2 = TRUE, style1= NULL, style2 = NULL)

jplot1 <- df.jack %>%
  select(linear.ad, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(d15n, 3, na.pad = TRUE)), size = 0.85, color = '#cc4c02') +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -900, lty = 'longdash') +
  geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +
  
  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 2) +
  annotate("text", x = -600, y = 7, label = "Iron Age Cold Epoch", size = 4) +
  annotate("text", x = 75, y = 7, label = "Roman Warm Period", size = 2) +
  
  ylab(n) +
  xlab(NULL) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(-1200,0)) +
  theme_classic() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

jplot2 <- df.jack %>%
  select(linear.ad, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(d13c, 3, na.pad = TRUE)), size = 0.85, color = '#fe9929') +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = -900, lty = 'longdash') +
  geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +

  ylab(c) +
  theme_classic() +
  xlab(x) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(-1200,0)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(jplot1), ggplotGrob(jplot2), size = "last"))

#' Figure: Ancient bulk values
#' --------------------------------------------------------------------------------------
#' Jack-4907 vs Sav-4902
#' --------------------------------------------------------------------------------------

plot(d15n ~ linear.ad, df.sav, # With Base R
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     # xlim = c(0,1500),
     col = alpha("black", 0.0))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
abline(v = 400, col = alpha("black", 0.75), lty = "longdash")
lines(forecast::ma(df.jack$d15n, order = 2, centre = TRUE) ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5)
lines(forecast::ma(df.sav$d15n, order = 3, centre = TRUE) ~ linear.ad, df.sav,
      col = "#f46d43", lwd = 1.5)

p1 <- ggplot() + # With ggplot2
  geom_line(data=df.jack, aes(x = linear.ad, y = d15n), color = "gray", alpha = 0.0, size = 0.5) +
  geom_line(data=df.jack, aes(x = linear.ad,
                              y = rollmean(d15n, 3, na.pad = TRUE)), color = "#4393c3", alpha = 0.85, size = 0.75) +
  geom_line(data=df.sav, aes(x = linear.ad, y = d15n), color = "gray", alpha = 0.0, size = 0.5) +
  geom_line(data=df.sav, aes(x = linear.ad,
                             y = rollmean(d15n, 3, na.pad = TRUE)), color = "#d6604d", alpha = 0.85, size = 0.75) +
  
  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 3) +
  # annotate("text", x = -600, y = 7, label = "Iron Age Cold Epoch", size = 4) +
  annotate("text", x = 75, y = 7, label = "Roman Warm Period", size = 3) +
  annotate("pointrange", y = 7.5, x = 650, xmin = (650-200), xmax = (650+200), color = "black") +
  
  geom_vline(xintercept = 950, lty = 'dotted') +
  geom_vline(xintercept = 1250, lty = 'dotted') +
  # geom_vline(xintercept = -900, lty = 'longdash') +
  # geom_vline(xintercept = -300, lty = 'longdash') +
  geom_vline(xintercept = -250, lty = "dotted") +
  geom_vline(xintercept = 400, lty = "dotted") +
  
  ylab(c) +
  theme_classic() +
  xlab(x) +
  # xlim(-1200, 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(-300,1300)) +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_text(size=10, color = "black"),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p1

plot(d13c ~ linear.ad, df.sav,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     ylim = c(-19, -14),
     # xlim = c(0,1500),
     col = alpha("black", 0.0))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
abline(v = 400, col = alpha("black", 0.75), lty = "longdash")
lines(forecast::ma(df.jack$d13c, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5)
lines(forecast::ma(df.sav$d13c, order = 3, centre = TRUE) ~ linear.ad, df.sav,
      col = "#f46d43", lwd = 1.5)

plot(d15n ~ linear.ad, df.sav,
     xlab = x,
     ylab = n,
     type = "o",
     cex = 0.5,
     xlim = c(500,1500),
     col = alpha("black", 0.0))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
lines(d15n ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5, type = "o")
lines(d15n ~ linear.ad, df.sav,
      col = "#f46d43", lwd = 1.5, type = "l")

#' -----------------------------------------------------
#' Compound-specific figures below
#' 
#' 
#' 
#' 
#' 
#' -----------------------------------------------------

################################
##                            ##
##                            ##
## N-CSIAA Data visualization ##
##                            ##
##                            ##
################################
ndata <- read.csv("~/Google Drive/projects/rproj/seus/cleaned_ndata.csv") # Contains SEUS black coral data, GOM black corals, and some POM data from elsewhere

######################
## Overall AA plots ##
######################
# Helpful link: https://stackoverflow.com/questions/21982987/mean-per-group-in-a-data-frame

bwtheme <- standard.theme("pdf", color=FALSE)

myColours <- brewer.pal(9,'Greys')

my.settings <- list(
  superpose.symbol=list(fill=myColours[2:5],col = "black", border="transparent", pch=c(21,23,22,24,25,8)),
  strip.background=list(col=myColours[6]),
  strip.border=list(col='black'))

ndata %>%
  melt(., "Type") %>%
  filter(variable %in% tr.srcaa) %>% # Success! 03-04-2019
  xyplot(value ~ variable,
         .,
         panel = function(x,  y, ...){
           panel.xyplot(x, y, ...)
           panel.abline(h = 9, lty = 1)
           panel.abline(h = 9.25, lty = 2)
           panel.abline(h = 8.75, lty = 2)
         },
         cex=1.5,
         ylim=c(-20, 40), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
         group = Type,
         xlab = NULL,
         ylab = expression({delta}^15*"N (\u2030)"),
         auto.key=list(columns=2, cex = 0.75),
         par.settings = my.settings)


ndata.melt <- melt(ndata[1:45,], "Sample.ID2")
ndata.melt$variable1 <- factor(ndata.melt$variable, 
                               levels=tr.srcaa)

bwtheme <- standard.theme("pdf", color = FALSE) # Figure for black and white overall plot, using myColors above
xyplot(value ~ variable1,
       ndata.melt,
       panel = function(x,  y, ...){
         panel.xyplot(x, y, ...)
         panel.abline(h = 9, lty = 1)
         panel.abline(h = 9.25, lty = 2)
         panel.abline(h = 8.75, lty = 2)
       },
       cex=1.5,
       ylim=c(-20, 40), # ylim and xlim are "first class" parameters and don't need to be in scales=list()
       group = Sample.ID2,
       xlab = NULL,
       ylab = expression({delta}^15*"N (\u2030)"),
       auto.key=list(columns=2, cex = 0.75),
       par.settings = my.settings)

######################
## Trophic Dynamics ## 
######################
ndata %>%
  filter(Type == "Black Coral - GOM" | Type == "Black Coral - SEUS") -> bc.ndata

p <- ggplot(bc.ndata, aes(x=Sum.V, y=TP, shape = Region, fill = Region), color = "black")
p + geom_point(size = 3) +
  facet_wrap(~Region, ncol = 1) +
  scale_shape_manual(values = c(22, 23)) +
  scale_color_manual(values = c("#D55E00", "#E69F00")) +
  ylab("Trophic Position (Glu - Phe)") +
  xlab(expression(paste(Sigma,"V"))) +
  theme_bw() +
  theme(axis.text=element_text(size=11, color = "black"))


########################
## SumV through time  ##
########################
ndata %>%
  filter(Region == "SEUS") %>%
  droplevels(.) %>%
  xyplot(Sum.V ~ Year.CE,
         data = .,
         group = Sample.ID2,
         cex = 1.5,
         # xlim = c(1500,2010),
         par.settings = my.settings,
         auto.key = list(columns=c(2), cex = 0.95),
         xlab = x,
         ylab = expression(paste(Sigma,"V")))

########################
## TP through time    ##
########################
ndata %>%
  # filter(Region == "SEUS") %>%
  # droplevels(.) %>%
  xyplot(TP ~ Year.CE,
         data = .,
         group = Sample.ID2,
         cex = 1.5,
         # xlim = c(750, 1500),
         par.settings = my.settings,
         auto.key = list(columns=c(2), cex = 0.95),
         xlab = x,
         ylab = "Trophic Position (Glu - Phe)")

######################
## Phe through time ##
######################

ndata %>% # Showing Phe through time, but only for select specimens (modify code as needed)
  filter(Phe > 1) %>%
  filter(Phe < 13) %>%
  filter(Sample.ID2 == "Jacksonville-4907" | Sample.ID2 == "Savannah Banks-4902" | Sample.ID2 == "Jacksonville-4684") -> t.ndata

par(mfrow=c(1,2))
xyplot(Phe ~ Year.CE,
       data = t.ndata,
       groups = Sample.ID2,
       # xlim =c(1500,1900),
       xlab = x,
       # xlim = c(-250,1500),
       ylab = phe,
       cex = 3,
       ylim = c(4, 14))


################################
##                            ##
##                            ##
## C-CSIAA Data visualization ##
##                            ##
##                            ##
################################
seus_carbon <- read.csv("C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson_cleaned.csv")
seus_carbon <- read.csv("~/Google Drive/projects/rproj/seus/data/schiff c-csiaa stetson_cleaned.csv")

######################
## Overall AA Plots ##
######################
# Helpful link: https://stackoverflow.com/questions/21982987/mean-per-group-in-a-data-frame