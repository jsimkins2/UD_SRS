library(ncdf4)
library(raster)
library(lubridate)
rastFile = nc_open("/home/sat_ops/goesR/sstClimatology/GOES16_SST_8day.nc")
data_dim <- rastFile$dim

goesRast = raster(ncvar_get(rastFile,'sst') - 273.15)
jdayGoes = format(as.POSIXct(rastFile$dim$time$vals,origin = "1970-01-01",tz = "GMT"),format="%j")
extent(goesRast)=c(min(ncvar_get(rastFile,'longitude')), max(ncvar_get(rastFile, 'longitude')),min(ncvar_get(rastFile,'latitude')), max(ncvar_get(rastFile, 'latitude')))
crs(goesRast) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
goesRast= flip(t(goesRast),2)


clim.nc = nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/aqua_clim_rolling8.nc")
rec_time  <- max(clim.nc$dim$time$vals)
rec_len   <- clim.nc$dim$time$len
sst.c <- ncvar_get(clim.nc, "sst",start = c(1,1,as.numeric(jdayGoes)), count = c(clim.nc$dim$lon$len,clim.nc$dim$lat$len,1))
sst.c = flip(t(raster(sst.c)),2)
crs(sst.c) = "+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs"
extent(sst.c)=c(min(ncvar_get(clim.nc,'lon')), max(ncvar_get(clim.nc, 'lon')),min(ncvar_get(clim.nc,'lat')), max(ncvar_get(clim.nc, 'lat')))

sst.c = projectRaster(sst.c, crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
extent(sst.c) = extent(goesRast)
sst.c=resample(sst.c, goesRast, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sst.c[sst.c < 0] = 0
goesRast[goesRast < 0] = 0

sst.anomnc = goesRast - sst.c
#sst.anomnc = (sst.anomnc * 9/5) + 32
var.list <- list()
var.list[[1]] <- ncdf4::ncvar_def(name="sst", units="Celsius", missval=-999, longname = "Sea Surface Temperature Anomaly", dim=data_dim)
loc.file <- paste0("/home/sat_ops/goesR/sstClimatology/SSTanomalyGoesAqua8dayCelsius.nc")
#loc.file <- "Downloads/atest.nc"
#writing all we need to the output file
loc <- ncdf4::nc_create(filename=loc.file, vars=var.list, force_v4 = T)
ncdf4::ncvar_put(nc=loc, "sst", vals=as.matrix(sst.anomnc))
ncdf4::ncatt_put(nc=loc, 0, "Conventions", "CF=1.0")
ncdf4::ncatt_put(nc=loc, 0,"creator_name", "James Simkins")
ncdf4::ncatt_put(nc=loc, 0, "creator_email", "simkins@udel.edu")
ncdf4::ncatt_put(nc=loc, 0, "institution", "University of Delaware Ocean Exploration, Remote Sensing and Biogeography Group (ORB)")
ncdf4::ncatt_put(nc=loc, 0, "url", "http://orb.ceoe.udel.edu/")
ncdf4::ncatt_put(nc=loc, 0, "source", "satellite observation NASA MODIS-Aqua instrument")
ncdf4::ncatt_put(nc=loc, 0, "groundstation", "University of Delaware, Newark, Center for Remote Sensing")
ncdf4::ncatt_put(nc=loc, 0, "software", "0.0")
ncdf4::ncatt_put(nc=loc, 0, "inputMET1", "0.0")
ncdf4::ncatt_put(nc=loc, 0, "inputOZONE1", "0.0")
ncdf4::ncatt_put(nc=loc, 0, "inputCalibrationFile", "0.0")
ncdf4::ncatt_put(nc=loc, 0, "product_list", "sst")
ncdf4::ncatt_put(nc=loc, 0, "summary", "GOES-16 SST Anomaly based on MODIS Aqua 8 day Aggregate; Regridded to 
                 EPSG 4326 - lat/lon projection. Processed at the Univeristy of Delaware")
ncdf4::nc_close(loc)
