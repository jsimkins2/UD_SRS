# this script is dedicated to repositioning GPT Aqua data to the APS Aqua data lat/lons
# James Simkins
library(ncdf4)
library(raster)

aps_nc = nc_open("Downloads/aqua.2019001.0101.235959.D.L3.modis.NAT.v09.1000m.nc4")
gpt_nc = nc_open("Downloads/aqua.2019044.0213.235959.D.L3.modis.NAT.v09.1000m.nc")#, write = TRUE)



aps_rast = t(raster(ncvar_get(aps_nc,'sst')))
gpt_rast = t(raster(ncvar_get(gpt_nc,'sst')))
extent(aps_rast)=c(min(ncvar_get(aps_nc,'lon')), max(ncvar_get(aps_nc, 'lon')),min(ncvar_get(aps_nc,'lat')), max(ncvar_get(aps_nc, 'lat')))
extent(gpt_rast) = c(min(ncvar_get(gpt_nc,'lon')), max(ncvar_get(gpt_nc, 'lon')),min(ncvar_get(gpt_nc,'lat')), max(ncvar_get(gpt_nc, 'lat')))

crs(aps_rast) = "+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs"
crs(gpt_rast) = "+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs"

x_rast=resample(gpt_rast, aps_rast, crs = "+ellps=WGS84 +proj=merc +lon_0=0.0 +lat_ts=0.0 +units=m +no_defs")
writeRaster(x_rast, filename='Downloads/x_rast.nc', format="CDF", overwrite=TRUE,
            varname="sst", varunit="Kelvin", longname="sst",
            xname="longitude", yname="latitude")