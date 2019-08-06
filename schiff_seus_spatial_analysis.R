## ---
## Topic: Using Spatial Analysis Tools in R for SEUS Deep-sea Coral Research
## John Schiff
## Date: October 21, 2018
## ---

##################################################################
## The code in this script file is for plotting study sites     ##
## for a given project. This can be done with a variety of      ##
## packages dependening on what you want to visualize spatially ##
##################################################################

library(maps)
library(ggplot2)
library(marmap)
library(oceanmap) # A new (young) package to plot data sites, 
# but you can do this with other packages too

########################################################################
##                                                                    ##
## Plot the Black coral stations with plotmap() from oceanmap package ##
## Use figure() function from same package to save plots              ##
##                                                                    ##
########################################################################

names <- c("Jacksonville 4684BC1", "Jacksonville 4907BC1", "Savannah 4902BC1", 
           "Stetson Banks 4904BC1", "Jacksonville 4686BC1")
latitude <- c(30.51333333, 30.8012600, 31.7040000, 31.8441917, 30.5021667)
longitude <- c(-079.6616667, -079.6416667, -079.1271667, -077.609956, -079.6515000)
stations <- cbind(longitude, latitude)

lon1 <- c(-86,-75)
lat1 <- c(25,35)
plotmap(lon=lon1, lat=lat1, main="Black Coral Stations", grid=FALSE, col.bg="lightblue")

# text(longitude, latitude, labels=names, cex=0.45)
legend("bottomright", legend=names, pch=c(0,1,2,5,6), text.width=5, ncol=1, cex=0.55)

###########################################################################################
## Same plot as above but with three data points (Jacksonville, Stetson Banks, Savannah) ##
###########################################################################################

names2 <- c("Jacksonville", "Savannah", "Stetson")
latitude2 <- c(30.6, 31.7040000, 31.8441917)
longitude2 <- c(-079.65,-079.1271667, -077.609956)
plotmap(lon=lon1, lat=lat1, main="Black Coral Stations", grid=FALSE, col.bg="lightblue")
get.bathy(lon=lon1, lat=lat1, terrain=FALSE)
points(longitude2, latitude2, pch=c(21),col="black", bg="red",cex=1.75)
text(longitude2, latitude2, labels=names2, cex=0.65, pos=1)
# legend("topright", inset=c(0.55, 0.05), 
#        legend=names2, pch=c(21,22,23), col="black", text.width=3.5, ncol=1, cex=0.5)

map("worldHires", 
    "usa",
    xlim=c(-86.5,-76.5),
    ylim=c(25,34),
    col="gray90",
    fill=TRUE)
points(longitude,latitude,pch=1,col="red")
map.scale(-86,33.7,ratio=FALSE,relwidth=0.3,cex=1.0)

############################################################
##                                                        ##
## Bathymetric map of the study site using marmap package ##
##                                                        ##
############################################################
# Helpful link: https://remi-daigle.github.io/2017-CHONe-Data/shocknawe.nb.html (not by those who wrote the package)
# Helpful link: https://stackoverflow.com/questions/43543318/map-inset-in-marmap-including-a-marmap-file-in-ggplot2
seusbathy <- getNOAA.bathy(lon1=-85, lon2=-70,
                           lat1=24, lat2=39, 
                           resolution=1)
par(pty = "s") # Last ran on 3-2-2019. Do this one!
plot(seusbathy, image=TRUE, land = FALSE, col = c("black", "darkgrey"), n = 12, drawlabels = TRUE)
plot(seusbathy, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
box(which = "plot", lty = "solid")
points(longitude2,latitude2, pch=c(21,21,21), col="black", bg=c('#bd0026', '#f03b20', '#feb24c'), cex=1.75)
# text(longitude2, latitude2, labels=names2, cex=1, pos=1)
legend("topleft", 
       legend=names2, pch=c(21,21,21), col="black", pt.bg = c('#bd0026', '#f03b20', '#feb24c'), bty="n", cex = 0.95) # inset=c(0.55, 0.05),
# text.width=3.5, ncol=1, cex=0.5)
scaleBathy(seusbathy, deg=1.75, x="bottomleft", inset=10)

############################################################
## Transect data visualization using marmap               ## 
############################################################

transect <- get.transect(seusbathy, -82.5, 31.75, -76.5, 32.25,  
             distance=TRUE)
plotProfile(transect)



# ----

data(seus)
z <- as.bathy(seus)
library(lattice)
wireframe(unclass(z), 
          screen=list(x=-25, y=10, z=15), shade=TRUE, aspect=c(1,0.5), zlab="z", xlab="x", ylab="y")
# use screen = list() to rotate view

####################################
## Using code from a blog I found ##
####################################
# Link: https://www.molecularecologist.com/2015/07/marmap/ I have modified it for my study area
#  Fetch data on NOAA servers and write on disk
bat <- getNOAA.bathy(lon1=-95, lon2=-75,
                     lat1=15, lat2=35, 
                     resolution=1)

# Create nice looking color palettes
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))
# Plot
par(pty = "s")
plot(bat, image = TRUE, land = FALSE, lwd = 0.1, bpal = list(c(0, max(bat), greys), c(min(bat), 0, blues)))
plot(bat, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline
points(longitude2,latitude2, pch=c(21), col="black", bg="red", cex=1.5)
scaleBathy(seusbathy, deg=4.19, x="bottomleft", inset=5)

#' -------------------------------------------------------------
#' Using some different code I found
#' Link: https://hansenjohnson.org/post/bathymetric-maps-in-r/
#' Link: https://dankelley.github.io/r/2015/04/03/oce-proj.html
#'
#'
#'
#' -------------------------------------------------------------

# Get the data
library(marmap)

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
span2 = 300
lonlim = c(-90, -70)
latlim = c(20, 35)
# lonlim2 = c(-81, -77)
# latlim2 = c(27, 33)

# plot coastline (with projection)
# par(mfrow = c(1,2))
plot(coastlineWorldFine, clon = mlon, clat = mlat, span = span, 
     projection="+proj=merc", col = 'lightgrey', drawBox = TRUE)
# plot bathymetry
mapContour(bathyLon,bathyLat,bathyZ,
        levels = c(-500, -1000, -1500, -2000, -2500, -3000, -3500, -4000, -4500, -5000),
        # lwd = c(1, 1, 2, 2, 3),
        # lty = c(3, 1, 3, 1, 3),
        col = 'darkgray')

# # add depth legend
legend("topleft", 
       # seg.len = 3,
       pch = c(4,2,1,6,5),
       cex = 0.4,
       legend = c("Jacksonville-4907 BC1", "Jacksonville-4684 BC1", "Jacksonville-4686 BC1", "Savannah-4902 BC1", "Stetson-4904"),
       col = 'black', 
       # title = "Depth [m]",
       bg = "white")

# add map data
mapPoints(longitude = pts$lon, latitude = pts$lat, pch = c(4,2,1,6,5), lwd = 1.5, col = "black", cex = 1.25)
# mapLines(longitude = lin$lon, latitude = lin$lat, col = 'blue')
# mapPolygon(longitude = ply$lon, latitude = ply$lat, lty = 2)






##########################################################
## UNUSED CODE: bits and pieces not used in actual code ##
##########################################################

# get.bathy(lon=lon1, lat=lat1, cbpos='r', grid=FALSE, cby=c(-1.25,-1.25))
# points(longitude, latitude, pch=c(0,1,2,5,6),col="red")
# legend("topleft", 
#        legend=names, 
#        pch=c(0,1,2,5,6), 
#        text.width=3.5, 
#        ncol=1, 
#        cex=0.70, 
#        bg="white", 
#        inset=c(0.15, 0.03))



