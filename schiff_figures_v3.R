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


# For now, use the figure already generated and copy paste that into the file

#' Figure: Recent history
#' --------------------------------------------------------------------------------------
#' Bulk data - Stet4904 and Jack4684 d15N
#' --------------------------------------------------------------------------------------

