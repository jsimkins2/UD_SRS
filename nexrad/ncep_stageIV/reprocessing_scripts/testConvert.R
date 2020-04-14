library(lubridate)
library(ncdf4)
datadir = "/home/james/test/"
Sys.setenv(TZ="UTC")
utcnow = now()
utc14 = now() - days(2)
stageUrls = list()
fnames = list.files("/home/james/test/")
for (temFile in fnames){
    system(paste0("/usr/local/bin/gdal_translate ",datadir,temFile, " ",datadir, temFile, ".nc -of netcdf"))
    system(paste0("/usr/local/bin/gdalwarp ",datadir,temFile,".nc ", datadir,temFile,"r.nc -of netcdf -t_srs epsg:4326"))
    
    f=paste0(temFile,"r.nc")
    timeVal = as.POSIXct(ymd_h(substr(f,5,14)))
    loc = nc_open(paste0(datadir,temFile, "r.nc"),write = TRUE)
    prec = ncvar_get(nc=loc, varid="Band1")
    timedim = ncdim_def("time", "seconds since 1970-01-01",as.numeric(timeVal))
    londim = loc$dim$lon
    latdim = loc$dim$lat
    
    prec_def = ncvar_def(name="Precipitation_Flux", units="mm/hr", missval=loc$var$Band1$missval, longname = "Quantitative Precipitation Estimates", dim=list(londim,latdim,timedim))
    ncdf4::ncvar_add(nc=loc, v=prec_def)
    nc_close(loc)
    
    loc = nc_open(paste0(datadir,temFile, "r.nc"),write = TRUE)
    ncvar_put(nc=loc,varid="Precipitation_Flux", vals=prec)
    nc_close(loc)
    
    system(paste0("/usr/bin/ncks -x -v Band1 ",datadir,temFile, "r.nc ",datadir,temFile, "r.nc -O"))
    
  }




