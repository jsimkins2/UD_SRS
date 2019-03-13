# this script is dedicated to repositioning GPT Aqua data to the APS Aqua data lat/lons
# I tried many different avenues but this is going to be the most accurate way to do it

# James Simkins 2019
library(ncdf4)
library(raster)
library(lubridate)
# Set up file output locations
tem_folder = "Downloads/"
finalFolder = paste0("/data/Aqua/1_day/", year(Sys.Date()), "/")
aps_file = "Downloads/aqua.2019001.0101.235959.D.L3.modis.NAT.v09.1000m.nc4"
gpt_file = "Downloads/aqua.2019044.0213.235959.D.L3.modis.NAT.v09.1000m.nc"
# load in the files - once this is set to go we need to change the gpt file to the output of
# gpt_modis_tc_1day.R - the aps file can be the most recent file from our Thredds Aggregation
aps_nc = nc_open(aps_file)
gpt_nc = nc_open(gpt_file)

var_names = names(aps_nc$var)

# save each raster output as it's own netcdf file.
# we'll read these values in and write a new netcdf file with the correct data later
t1 = Sys.time()
for (v in var_names){
  aps_rast = t(raster(ncvar_get(aps_nc,'sst')))
  gpt_rast = t(raster(ncvar_get(gpt_nc,'sst')))
  
  extent(aps_rast)=c(min(ncvar_get(aps_nc,'lon')), max(ncvar_get(aps_nc, 'lon')),min(ncvar_get(aps_nc,'lat')), max(ncvar_get(aps_nc, 'lat')))
  extent(gpt_rast) = c(min(ncvar_get(gpt_nc,'lon')), max(ncvar_get(gpt_nc, 'lon')),min(ncvar_get(gpt_nc,'lat')), max(ncvar_get(gpt_nc, 'lat')))
  
  crs(aps_rast) = "+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs"
  crs(gpt_rast) = "+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs"
  
  x_rast=resample(gpt_rast, aps_rast, crs = "+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs")
  
  writeRaster(x_rast, filename=paste0(tem_folder, v, ".nc"), format="CDF", overwrite=TRUE,
              varname=v, varunit=aps_nc$var[[v]]$units, longname=aps_nc$var[[v]]$longname,
              xname="longitude", yname="latitude")
}
nc_close(aps_nc)
nc_close(gpt_nc)
# runtime is 16 minutes on macbook
Sys.time() - t1
# now, we need to read all of the raster.nc files back in and replace the values within the APS file
aps_file = "Downloads/aqua.2019001.0101.235959.D.L3.modis.NAT.v09.1000m.nc4"
aps_nc = nc_open(aps_file)


# Extract the current file time and write to a dimension in posixct time. Note: this only works for the current storage syntax from basin
dateStr = strsplit(strsplit(gpt_file, 'aqua.')[[1]][2], '.D')[[1]][1]
ncposix <- as.POSIXct(strptime(paste0(substr(dateStr, 1,4),substr(dateStr, 9,12), substr(dateStr, 14,19)), "%Y%m%d%H%M%S", tz="GMT"))
new_time <- ncdf4::ncdim_def(name='now_time', units="seconds since 1970-01-01T00:00:00Z", 
                             vals=as.numeric(ncposix), create_dimvar=TRUE, unlim=TRUE)

aps_nc$dim$time$val  <- as.numeric(ncposix)
data_dim <- aps_nc$dim

var.list = list()
for (v in var_names){
  var.list[[v]] <- ncdf4::ncvar_def(name=as.character(v), 
                                    units=as.character(aps_nc$var[[v]]$units), 
                                    dim=data_dim, missval=NA)
}
nc_close(aps_nc)

#creating the nc4 output file
#loc.file <- paste0(finalFolder,gpt_file)
gpt_file = "Downloads/aqua.2019044.0213.235959.D.L3.modis.NAT.v09.1000m4.nc4"
loc.file <- paste0(gpt_file)

#writing all we need to the output file
loc <- ncdf4::nc_create(filename=loc.file, vars=var.list, force_v4 = TRUE)
for(j in seq_along(var_names)){
  flipped = t(ncvar_get(nc_open(paste0(tem_folder, var_names[[j]], ".nc"))))
  flipped = apply(flipped, 2, rev)
  flippedvar= t(flipped)
  flippedvar[flippedvar == 0] = NA # not sure if this should be like this but the gpt nc files report NAs as 0s (I think)
  ncdf4::ncvar_put(nc=loc, varid=as.character(var_names[[j]]), vals=flippedvar)
  #ncdf4::ncvar_put(nc=loc, varid=as.character(var_names[[j]]), vals=ncvar_get(nc_open(paste0(tem_folder, var_names[[j]], ".nc")), var_names[[j]]))
}
ncdf4::nc_close(loc)



x = nc_open(gpt_file)
r = raster(ncvar_get(x, 'sst'))
extent(r)=c(min(ncvar_get(x,'lon')), max(ncvar_get(x, 'lon')),min(ncvar_get(x,'lat')), max(ncvar_get(x, 'lat')))



gpt_sst= ncvar_get(x, 'sst')
aps_sst = ncvar_get(aps_nc, 'sst')

