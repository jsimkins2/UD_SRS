# converting AQUA data from netcdf to geotiff after a resampling
# the goal is to convert this to the ROFFS output
library(rgdal)
library(raster)
library(ncdf4)

f<- nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/Aqua1DayAggregate.nc")


tm<- ncvar_get(f, "time")
lon <- ncvar_get(f, "lon") ### loading the longitude from the THREDDS server
#lon<-ncvar_get(f, 'longitude')
lat <- ncvar_get(f, "lat") ### loading the latitude from the THREDDS server
#lat<-ncvar_get(f, 'latitude')
loncount = which.min(abs(lon - -80)) - which.min(abs(lon - -98))
latcount = which.min(abs(lat - 31)) - which.min(abs(lat - 18))
sst<- ncvar_get(f,"sst", start = c(which.min(abs(lon - -98)), which.min(abs(lat - 25)),length(tm)), count = c(loncount, latcount,1))

sstRast = raster(sst)
sstRast = flip(t(sstRast), 'y')
extent(sstRast) = c(-98, -80, 18, 31)
crs(sstRast) = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#sstRast = crop(sstRast,extent(-98, -80, 18, 31))

newproj <- "+proj=merc +lon_0=-89 +lat_ts=24.5 +x_0=0 +y_0=0 +ellps=WGS72 +towgs72=0,0,4.5,0,0,0.554,0.2263 +units=m +no_defs"
#newproj <- "+proj=longlat +ellps=WGS72 +towgs84=0,0,4.5,0,0,0.554,0.2263 +no_defs "

projRast = projectRaster(sstRast, crs = newproj)
rast = projRast
refTif = raster("Downloads/182971533TMO.wmex.tif")
extent(rast) = extent(refTif)
rast = resample(rast, refTif)


writeRaster(x = rast, filename = 'Downloads/test3.tif', format = "GTiff", overwrite=TRUE, varname='band1', 
            prj=TRUE, datatype = 'INT1U',options=c("COMPRESS=NONE"))



test = raster("Downloads/test.tif")

sst.aggregate = aggregate(projRast, fact=108910.9)



require(rgdal)
require(maptools)
require(raster)

data(wrld_simpl)
mollCRS = CRS('+proj=moll')
sst_moll <- projectRaster(sstRast, crs=mollCRS)
wrld <- spTransform(wrld_simpl, mollCRS)

data(wrld_simpl)
w <- crop(wrld, sst_m)

plot(sstRast)
plot(w, col='gray', add=TRUE)