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
path1 <- '~/Documents/GitHub/data/schiff_bulk_years_08-04-2019.csv'
path2 <- 'C:/Users/jschiff.GEOSAD/Google Drive/projects/rproj/seus/data/schiff_bulk_years_08-04-2019.csv'
path3 <- '/home/john/Desktop/data/schiff_bulk_years_08-04-2019.csv'

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

par(mfrow = c(1,3))
plot(X14C.Age ~ Distance..um., r.jack,
     pch = 22, bg = '#ffeda0', col = "black", 
     ylab = expression({Delta}^14*"C Age (CRA)", xlab = 'Distance from edge (um)'))
# abline(mod1)
plot(X14C.Age ~ Distance..um., r.sav,
     pch = 23, bg = '#feb24c', col = "black", 
     ylab = expression({Delta}^14*"C Age (CRA)", xlab = 'Distance from edge (um)'))
# abline(mod2)
plot(X14C.Age ~ Distance..um., r.stet,
     pch = 21, bg = '#f03b20', col = "black",
     ylab = expression({Delta}^14*"C Age (CRA)", xlab = 'Distance from edge (um)'))
# abline(mod3)

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

par(pty = "s")
plot(d13c ~ linear.ad, df.jack,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     col = alpha("black", 0.0))
lines(forecast::ma(df.jack$d13c, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = alpha("black", 0.4), lwd = 1.5)
lines(forecast::ma(df.jack2$d13c, order = 3, centre = TRUE) ~ linear.ad, df.jack2,
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
lines(forecast::ma(df.stet$d15n, order = 3, centre = TRUE) ~ linear.ad, df.stet,
      col = "#33a02c", lwd = 1.5)
lines(forecast::ma(df.jack4684$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack4684,
      col = "#1f78b4", lwd = 1.5)

#' Figure: Recent history, bulk
#' --------------------------------------------------------------------------------------
#' Bulk data - Stet-4904 entire record
#' --------------------------------------------------------------------------------------

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
  select(linear.ad, d15n) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d15n), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(d15n, 3, na.pad = TRUE)), size = 0.85, color = '#225ea8') +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
  geom_vline(xintercept = 1550, lty = 'dotted') +
  geom_vline(xintercept = 1850, lty = 'dotted') +
  
  annotate("text", x = 1100, y = 7, label = 'Medieval Warming', size = 3.75) +
  annotate("text", x = 1700, y = 7, label = "Little Ice Age", size = 4) +
  
  ylab(n) +
  xlab(NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(axis.text.y   = element_text(size=10, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

splot2 <- df.stet %>%
  select(linear.ad, d13c) %>%
  na.omit() %>%
  ggplot(aes(x = linear.ad, y = d13c), size = 0.5, alpha = 0.75) +
  geom_line(color = "gray", alpha = 0.4, size = 0.75) +
  geom_line(aes(y=rollmean(d13c, 3, na.pad = TRUE)), size = 0.85, color = '#081d58') +
  
  geom_vline(xintercept = 950, lty = 'dashed') +
  geom_vline(xintercept = 1250, lty = 'dashed') +
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
  annotate("text", x = -600, y = 7, label = "Iron Age Cold Epoch", size = 2) +
  annotate("text", x = 75, y = 7, label = "Roman Warm Period", size = 2) +
  
  ylab(n) +
  xlab(NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
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
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
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

plot(d15n ~ linear.ad, df.sav,
     xlab = x,
     ylab = n,
     type = "l",
     cex = 0.5,
     # xlim = c(0,1500),
     col = alpha("black", 0.0))
abline(v = 950, col = "black", lty = 'dashed')
abline(v = 1250, col = "black", lty = 'dashed')
abline(v = 400, col = alpha("black", 0.75), lty = "longdash")
lines(forecast::ma(df.jack$d15n, order = 3, centre = TRUE) ~ linear.ad, df.jack,
      col = "#4575b4", lwd = 1.5)
lines(forecast::ma(df.sav$d15n, order = 3, centre = TRUE) ~ linear.ad, df.sav,
      col = "#f46d43", lwd = 1.5)

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



