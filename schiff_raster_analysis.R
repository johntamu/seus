library(raster)
r <- stack('~/Documents/GitHub/gridded-data/nesdisVHNSQchlaMonthly_47a5_7140_35b3.nc')
plot(r)
nlayers(r)

for(i in 1:nlayers(r)) {
  band <- r[[i]]
  # save raster in a separate file
  writeRaster(band,paste('band',i,'.tif', sep=''))
}
