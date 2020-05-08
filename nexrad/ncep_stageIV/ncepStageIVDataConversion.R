library(lubridate)
library(ncdf4)
datadir = "/home/james/ncep_stageIV/"
Sys.setenv(TZ="UTC")
utcnow = now()
utc14 = now() - days(2)
stageUrls = list()
system(paste0("/usr/bin/rm ",paste0(datadir,'notNetcdf/*')))
system(paste0("/usr/bin/rm ",paste0(datadir,'unprojected/*')))
for (d in seq(utc14, utcnow, by="hour")){
  dateTime = as.POSIXct.numeric(d,tz='UTC',origin='1970-01-01 00:00:00')
  temFile = paste0("ST4.", year(dateTime),sprintf("%02i", month(dateTime)), sprintf("%02i", day(dateTime)), sprintf("%02i", hour(dateTime)), '.01h.gz')

  system(paste0("/usr/bin/wget -P /home/james/ncep_stageIV/zipped/ https://nomads.ncep.noaa.gov/pub/data/nccf/com/pcpanl/prod/pcpanl.",year(dateTime),sprintf("%02i", month(dateTime)), sprintf("%02i", day(dateTime)),"/",temFile))
  system(paste0("/usr/bin/gunzip ",datadir,'zipped/', temFile))
  temFile=substr(temFile, 0, nchar(temFile)-3)
  system(paste0("/usr/bin/mv ",datadir,'zipped/',temFile, " ", datadir,'notNetcdf/'))
  system(paste0("/usr/local/bin/gdal_translate ",datadir,"notNetcdf/",temFile, " ",datadir,'unprojected/', temFile, ".nc -of netcdf"))
  if (file.exists(paste0("/data/ncep_stageIV/indiv/",year(dateTime),"/",temFile, "r.nc"))){
    system(paste0("/usr/bin/rm ",paste0("/data/ncep_stageIV/indiv/",year(dateTime),"/",temFile, "r.nc")))
  }
  system(paste0("/usr/local/bin/gdalwarp ",datadir,'unprojected/', temFile,".nc /data/ncep_stageIV/indiv/",year(dateTime),"/",temFile, "r.nc -of netcdf -t_srs epsg:4326"))
  
  f=paste0(temFile,".nc")
  timeVal = as.POSIXct(ymd_h(substr(f,5,14)))
  loc = nc_open(paste0("/data/ncep_stageIV/indiv/",year(dateTime),"/",temFile, "r.nc"),write = TRUE)
  prec = ncvar_get(nc=loc, varid="Band1")
  timedim = ncdim_def("time", "seconds since 1970-01-01",as.numeric(timeVal))
  londim = loc$dim$lon
  latdim = loc$dim$lat
  
  prec_def = ncvar_def(name="Precipitation_Flux", units="mm/hr", missval=loc$var$Band1$missval, longname = "Quantitative Precipitation Estimates", dim=list(londim,latdim,timedim))
  ncdf4::ncvar_add(nc=loc, v=prec_def)
  nc_close(loc)
  
  loc = nc_open(paste0("/data/ncep_stageIV/indiv/",year(dateTime),"/",temFile, "r.nc"),write = TRUE)
  ncvar_put(nc=loc,varid="Precipitation_Flux", vals=prec)
  nc_close(loc)
  
  system(paste0("/usr/bin/ncks -x -v Band1 ",paste0("/data/ncep_stageIV/indiv/",year(dateTime),"/",temFile, "r.nc")," ",paste0("/data/ncep_stageIV/indiv/",year(dateTime),"/",temFile, "r.nc"), " -O"))
}


# this part of the code grabs the last hourly file of each day (23rd hour), grabs the previous 23 files and concatenates them to 
# create a daily file
library(stringr)
utcnow = now()
dateTime = as.POSIXct.numeric(d,tz='UTC',origin='1970-01-01 00:00:00')
yearList = list(year(dateTime))
for (yr in yearList){
  newFiles = list.files(paste0('/data/ncep_stageIV/indiv/',yr, '/'))
  dayFiles = list()
  for (n in newFiles){
    if (substr(n, 13,14) == "23"){
      dayFiles = append(dayFiles, n)
    }
  }
  
  for (day in dayFiles){
    if (!file.exists(paste0('/data/ncep_stageIV/daily/',yr,'/', substr(day,1,12), 'dailycomp.01hr.nc'))){
      print(day)
      y=substr(day,5,8)
      m=substr(day, 9,10)
      n=substr(day, 11,12)
      fnames=list()
      for (h in seq(0,23)){
        fnames = paste(fnames,(paste0('/data/ncep_stageIV/indiv/', yr, '/ST4.',yr,str_pad(m, 2, pad = "0"),str_pad(n, 2, pad = "0"),
                                      str_pad(h, 2, pad = "0"), '.01hr.nc' )))
      }
      system(paste0("/usr/bin/ncks -O --mk_rec_dmn time ",'/data/ncep_stageIV/indiv/', yr, '/ST4.',yr,str_pad(m, 2, pad = "0"),str_pad(n, 2, pad = "0"),
                    str_pad(0, 2, pad = "0"), '.01hr.nc ','/data/ncep_stageIV/indiv/', yr, '/ST4.',yr,str_pad(m, 2, pad = "0"),str_pad(n, 2, pad = "0"),
                    str_pad(0, 2, pad = "0"), '.01hr.nc'))
      system(paste0("/usr/bin/ncrcat -h", fnames, " ",'/data/ncep_stageIV/', 'daily/', yr, '/', '/ST4.',yr,str_pad(m, 2, pad = "0"),str_pad(n, 2, pad = "0"),'dailycomp.01hr.nc'))
    }
  }
}
