library(ncdf4)
library(raster)
library(maptools)

d <- nc_open("/home/sat_ops/goesR/sstClimatology/GOES16_SST_8day.nc")
jdayGoes = format(as.POSIXct(d$dim$time$vals,origin = "1970-01-01",tz = "GMT"),format="%j")
sst <- ncvar_get(d, "sst")
lon <- ncvar_get(d, "longitude")
lat <- ncvar_get(d, "latitude")
nc_close(d)

xy <- cbind(rep(lon, length(lat)), rep(lat, each=length(lon)))
xyv <- na.omit(cbind(xy, as.vector(sst)))
r <- raster(extent(range(lon), range(lat)), res=1/30)
r <- rasterize(xyv[, 1:2], r, xyv[,3], fun=mean) 
r[r < 273.15] = NA
r = r - 273.15

clim.nc = nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/aqua_clim_rolling8.nc")
rec_time  <- max(clim.nc$dim$time$vals)
rec_len   <- clim.nc$dim$time$len
sst.c <- ncvar_get(clim.nc, "sst",start = c(1,1,as.numeric(jdayGoes)), count = c(clim.nc$dim$lon$len,clim.nc$dim$lat$len,1))
lon <- ncvar_get(clim.nc, "lon")
lat <- ncvar_get(clim.nc, "lat")
xy <- cbind(rep(lon, length(lat)), rep(lat, each=length(lon)))
xyv <- na.omit(cbind(xy, as.vector(sst.c)))
r2 <- raster(extent(range(lon), range(lat)), res=1/30)
r2 <- rasterize(xyv[, 1:2], r, xyv[,3], fun=mean) 
r2[r2 < 0] = NA

sstanom = r-r2
lon3 = unique(coordinates(sstanom)[,1])
lat3 = unique(coordinates(sstanom)[,2])
# create and write the netCDF file -- ncdf4 version
# define dimensions
londim <- ncdim_def("lon","degrees_east",lon3) 
latdim <- ncdim_def("lat","degrees_north",lat3) 
timedim <- ncdim_def("time",d$dim$time$units,d$dim$time$vals)

# define variables
fillvalue <- -9999
dlname <- "Sea Surface Temperature Anomaly"
tmp_def <- ncvar_def("tmp","deg_C",list(londim,latdim,timedim),fillvalue,dlname,prec="single")

var.list <- list()
var.list[[1]] <- ncdf4::ncvar_def(name="sst", units="Celsius", missval=-999, longname = "Sea Surface Temperature Anomaly", dim=list(londim,latdim,timedim))
loc.file <- paste0("/home/sat_ops/goesR/sstClimatology/SSTanomalyGoesAqua8dayCelsius.nc")
#writing all we need to the output file
loc <- ncdf4::nc_create(filename=loc.file, vars=var.list)
ncdf4::ncvar_put(nc=loc, "sst", vals=t(as.matrix(sstanom)))
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
ncdf4::ncatt_put(nc=loc, 0, "summary", "GOES-16 SST Anomaly based on MODIS Aqua 8 day Climatology; Regridded to 
                 EPSG 4326 - lat/lon projection. Processed at the Univeristy of Delaware")
ncdf4::nc_close(loc)
