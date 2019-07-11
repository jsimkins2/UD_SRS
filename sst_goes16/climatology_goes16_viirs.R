# this script is dedicated to repositioning GPT Aqua data to the APS Aqua data lat/lons
# I tried many different avenues but this is going to be the most accurate way to do it

# James Simkins 2019
library(ncdf4)
library(raster)
library(lubridate)

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