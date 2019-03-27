library(rgdal)
library(raster)
library(ncdf4)

aquaRast = raster(ncvar_get(nc_open("Downloads/aqua_test.nc"), 'sst'))
aquaRast = flip(t(aquaRast), 'y')
extent(aquaRast) = c(18, 31, -98, -80)
crs(aquaRast) = "+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs"

newproj <- "+ellps=WGS72 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs"
projRast = projectRaster(aquaRast, crs = newproj)

newr <- raster(nrow=1314, ncol=1314, crs = newproj)
extent(newr) = c(18, 31, -98, -80)

projRast = resample(projRast, newr)

writeRaster(projRast, 'Downloads/test.tif', format = "GTiff", overwrite=TRUE, varname='sst')#, compression=7)
