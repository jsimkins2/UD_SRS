# this script is dedicated to repositioning GPT Aqua data to the APS Aqua data lat/lons
# I tried many different avenues but this is going to be the most accurate way to do it

# James Simkins 2019
library(ncdf4)
library(raster)
library(lubridate)

clim.nc = nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/aqua_clim_rolling8.nc")
rec_time  <- max(clim.nc$dim$time$vals)
rec_len   <- clim.nc$dim$time$len
sst.c <- ncvar_get(clim.nc, "sst",start = c(1,1,as.numeric(jdayGoes)), count = c(clim.nc$dim$lon$len,clim.nc$dim$lat$len,1))
sst.c = flip(t(raster(sst.c)),2)
proj4string(sst.c) = CRS("+init=epsg:3857")
extent(sst.c)=c(min(ncvar_get(clim.nc,'lon')), max(ncvar_get(clim.nc, 'lon')),min(ncvar_get(clim.nc,'lat')), max(ncvar_get(clim.nc, 'lat')))
writeRaster(sst.c, filename="Downloads/rolling8.nc", format="CDF", overwrite=TRUE,
            varname="sst", varunit="celsius", longname="aqua sst 8 day climatology",
            xname="longitude", yname="latitude")

writeRaster(sst.c, filename="Downloads/rolling8.tif", format="GTiff", overwrite=TRUE,
            varname="sst", varunit="celsius", longname="aqua sst 8 day climatology",
            xname="longitude", yname="latitude")

gdal_translate rolling8.nc rolling8.tif -of GTIFF -a_srs epsg:3857 -a_nodata -999
gdalwarp rolling8.tif rep_rolling8.tif -of GTIFF -tr 0.018000001395654 0.018000001395654 -te -99.9991505 15.9937824 -49.9951466 51.9946839 -s_srs '+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs' -t_srs '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
gdalwarp -t_srs '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0' rolling8.tif rep_rolling8.tif
gdalwarp rep_rolling8.tif trolling.tif -of GTIFF -tr 0.017993521913511 -0.018090900672387 
gdalwarp rolling8.tif rep_rolling8.tif -of GTIFF -ot Float32 -tr 0.017993521913511 -0.018090900672387 -te -99.9991505 15.9937824 -49.9951466 51.9946839 -srcnodata -999 -dstnodata -999 -r average -s_srs EPSG:3857 -t_srs EPSG:4326

rep = raster("Downloads/rolling8.tif")
plot(rep)

rastFile = nc_open("Downloads/GOES16_SST_8day.nc")
goesRast = raster(ncvar_get(rastFile,'sst') - 273.15)
jdayGoes = format(as.POSIXct(rastFile$dim$time$vals,origin = "1970-01-01",tz = "GMT"),format="%j")
extent(goesRast)=c(min(ncvar_get(rastFile,'latitude')), max(ncvar_get(rastFile, 'latitude')),min(ncvar_get(rastFile,'longitude')), max(ncvar_get(rastFile, 'longitude')))
proj4string(goesRast) = CRS("+init=epsg:4326")
goesRast= flip(t(goesRast),2)
writeRaster(goesRast, filename="Downloads/goesRast.tif", format="GTiff", overwrite=TRUE,
            varname="sst", varunit="celsius", longname="aqua sst 8 day climatology",
            xname="longitude", yname="latitude")
gdalwarp -t_srs EPSG:3857 goesRast.tif repGoes.tif

gdalwarp rolling8.tif rep_rolling8.tif -of GTIFF -tr 0.012952772295672 0.025131223668040 -te -99.990 16.003 -50.0041466 51.9856339 -s_srs '+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs' -t_srs '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0' -srcnodata "-3.39999999999999996e+38"














library(ncdf4)
library(raster)
library(lubridate)
rastFile = nc_open("/home/sat_ops/goesR/sstClimatology/GOES16_SST_8day.nc")
data_dim <- rastFile$dim

goesRast = raster(ncvar_get(rastFile,'sst') - 273.15)
jdayGoes = format(as.POSIXct(rastFile$dim$time$vals,origin = "1970-01-01",tz = "GMT"),format="%j")
extent(goesRast)=c(min(ncvar_get(rastFile,'latitude')), max(ncvar_get(rastFile, 'latitude')),min(ncvar_get(rastFile,'longitude')), max(ncvar_get(rastFile, 'longitude')))
crs(goesRast) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
goesRast= flip(t(goesRast),2)

clim.nc = nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/aqua_clim_rolling8.nc")
rec_time  <- max(clim.nc$dim$time$vals)
rec_len   <- clim.nc$dim$time$len
sst.c <- ncvar_get(clim.nc, "sst",start = c(1,1,as.numeric(jdayGoes)), count = c(clim.nc$dim$lon$len,clim.nc$dim$lat$len,1))
sst.c = flip(t(raster(sst.c)),2)
proj4string(sst.c) = CRS("+init=epsg:3857")
extent(sst.c)=c(min(ncvar_get(clim.nc,'lon')), max(ncvar_get(clim.nc, 'lon')),min(ncvar_get(clim.nc,'lat')), max(ncvar_get(clim.nc, 'lat')))

rep = projectRaster(sst.c, crs="+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
extent(rep) = extent(sst.c)
rep=resample(rep, goesRast, crs = "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
sst.c[sst.c < 0] = 0
goesRast[goesRast < 0] = 0

sst.anomnc = goesRast - sst.c









