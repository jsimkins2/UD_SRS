# this script downloads the 24h files from the NCEP Stage IV precip database -
# https://nomads.ncep.noaa.gov/pub/data/nccf/com/pcpanl/prod/
# this is the real-time download and parse for THREDDS script

library(lubridate)
library(ncdf4)
datadir = "/home/james/ncep_stageIV/24h/"
Sys.setenv(TZ="UTC")
utcnow = now() - days(1)
utc14 = now() - days(13)
stageUrls = list()
system(paste0("/usr/bin/rm ",paste0(datadir,'notNetcdf/*')))
system(paste0("/usr/bin/rm ",paste0(datadir,'unprojected/*')))
for (d in seq(utc14, utcnow, by="day")){
  print(d)
  dateTime = as.POSIXct.numeric(d,tz='UTC',origin='1970-01-01 00:00:00')
  temFile = paste0("st4_conus.", year(dateTime),sprintf("%02i", month(dateTime)), sprintf("%02i", day(dateTime)),'12.24h.grb2')

  system(paste0("/usr/bin/wget -P /home/james/ncep_stageIV/24h/zipped/ https://nomads.ncep.noaa.gov/pub/data/nccf/com/pcpanl/prod/pcpanl.",year(dateTime),sprintf("%02i", month(dateTime)), sprintf("%02i", day(dateTime)),"/",temFile))
  system(paste0("/usr/bin/mv ",datadir,'zipped/',temFile, " ", datadir,'notNetcdf/'))
  system(paste0("/usr/local/bin/gdal_translate ",datadir,"notNetcdf/",temFile, " ",datadir,'unprojected/', temFile, ".nc -of netcdf"))
  # add these because gdal can't overwrite
  ncTemFile = paste0("ST4.", substr(temFile,11, 20), ".24h")
  if (file.exists(paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",ncTemFile, "r.nc"))){
    system(paste0("/usr/bin/rm ",paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",ncTemFile, "r.nc")))
  }
  system(paste0("/usr/local/bin/gdalwarp ",datadir,'unprojected/', temFile,".nc /data/ncep_stageIV/quality/",year(dateTime),"/",ncTemFile, "r.nc -of netcdf -t_srs epsg:4326"))
  
  f=paste0(ncTemFile,".nc")
  
  timeVal = as.POSIXct(ymd_h(substr(f,5,14)))
  loc = nc_open(paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",ncTemFile, "r.nc"),write = TRUE)
  prec = ncvar_get(nc=loc, varid="Band1")
  timedim = ncdim_def("time", "seconds since 1970-01-01",as.integer(timeVal))
  londim = loc$dim$lon
  latdim = loc$dim$lat
  
  prec_def = ncvar_def(name="Precipitation_Flux", units="mm/day", missval=loc$var$Band1$missval, longname = "Quantitative Precipitation Estimates", dim=list(londim,latdim,timedim))
  ncdf4::ncvar_add(nc=loc, v=prec_def)
  nc_close(loc)
  
  loc = nc_open(paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",ncTemFile, "r.nc"),write = TRUE)
  ncvar_put(nc=loc,varid="Precipitation_Flux", vals=prec)
  nc_close(loc)
  
  system(paste0("/usr/bin/ncks -x -v Band1 ",paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",ncTemFile, "r.nc")," ",paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",ncTemFile, "r.nc"), " -O"))
}
