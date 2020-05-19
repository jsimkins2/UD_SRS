# script for creating a bias mask for Sarah
library(R.matlab)
library(raster)
NYharbias = 0.397106319832402 #40.369 -73.703 91, 183
DBayBias =  0.351385446175637 #38.457 -74.702 196, 128
BarnBias = 0.065860630407911 #39.778 -73.769 123, 180



lons = readMat("Downloads/lon_bands.mat")
lats = readMat("Downloads/lat_bands.mat")

# create a raster 278x276 
# mins and maxes of lons/lats bands 
# throw in bias values above into array
# populate cells via bilinear and ngb
x <- matrix(data=NA,nrow=276,ncol=278)

# use this to find which values go where in the bias mask
which(abs(lons$lon.bands + 73.769 )==min(abs(lons$lon.bands + 73.769)))

#plug in the values
x[91, 183] = NYharbias
x[196, 128] = DBayBias 
x[123, 180] = BarnBias

rastX = raster(x)
extent(rastX) = c(min(lons$lon.bands), max(lons$lon.bands), min(lats$lat.bands), max(lats$lat.bands))
projection(rastX) <- CRS("+proj=longlat +datum=WGS84")
