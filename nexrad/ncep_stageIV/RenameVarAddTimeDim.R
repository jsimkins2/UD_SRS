#ST4.2019052818.01hr.nc

library(lubridate)
library(ncdf4)
yearList = seq(2019,2019,by=1)
datadir = "/data/ncep_stageIV/"
for (y in yearList){
  fileNames = list.files(paste0(datadir, y, "/"))
  for (f in fileNames){
    print(f)
    timeVal = as.POSIXct(ymd_h(substr(f,5,14)))
    loc = nc_open(paste0(datadir,y,"/",f),write = TRUE)
    prec = ncvar_get(nc=loc, varid="Band1")
    timedim = ncdim_def("time", "seconds since 1970-01-01",as.numeric(timeVal))
    londim = loc$dim$lon
    latdim = loc$dim$lat

    prec_def = ncvar_def(name="Precipitation_Flux", units="mm/hr", missval=loc$var$Band1$missval, longname = "Quantitative Precipitation Estimates", dim=list(londim,latdim,timedim))
    ncdf4::ncvar_add(nc=loc, v=prec_def)
    nc_close(loc)

    loc = nc_open(paste0(datadir,y,"/",f),write = TRUE)
    ncvar_put(nc=loc,varid="Precipitation_Flux", vals=prec)
    nc_close(loc)

    system(paste0("/usr/bin/ncks -x -v Band1 ",paste0(datadir,y,"/",f)," ",paste0(datadir,y,"/",f), " -O"))
  }
}
