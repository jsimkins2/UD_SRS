# Document detailing the downloading/processing/reprojection/aggregation of NCEP Stage IV data

#1) Download the archived data here - https://data.eol.ucar.edu/dataset/21.093
#   (real time data here - https://nomads.ncep.noaa.gov/pub/data/nccf/com/pcpanl/prod/)

#2) Untar all file folders. Go into each folder, remove .txt files, remove *06h* files, remove *24h*, and uncompress all. 
#   Remove duplicate files (20040101 file in a 2003 folder)

#3) 



# Now the script begins. 

library(lubridate)
library(ncdf4)
yearList = seq(2004,2007,by=1)
yearList = append(yearList, 2019)
datadir = "/home/james/ncep_stageIV/backfill_ncep/"

for (y in yearList){
  fileNames = list.files(paste0(datadir, y, "/"))
  for (fi in fileNames){
    temFile = substr(fi, 1,18)
    if (file.exists(paste0(datadir,y,"/",temFile,"r.nc")) == FALSE){

      print(temFile)
      system(paste0("/usr/local/bin/gdal_translate ",datadir,y,"/",temFile, " ",datadir,y,"/", temFile, ".nc -of netcdf"))
      system(paste0("/usr/local/bin/gdalwarp ",datadir,y,"/", temFile,".nc" ," ", datadir,y,"/",temFile, "r.nc -of netcdf -t_srs epsg:4326"))
      
      f=paste0(temFile,"r.nc")
      timeVal = as.POSIXct(ymd_h(substr(f,5,14)))
      loc = nc_open(paste0(datadir,y,"/",temFile, "r.nc"),write = TRUE)
      prec = ncvar_get(nc=loc, varid="Band1")
      timedim = ncdim_def("time", "seconds since 1970-01-01",as.numeric(timeVal))
      londim = loc$dim$lon
      latdim = loc$dim$lat
      
      prec_def = ncvar_def(name="Precipitation_Flux", units="mm/hr", missval=loc$var$Band1$missval, longname = "Quantitative Precipitation Estimates", dim=list(londim,latdim,timedim))
      ncdf4::ncvar_add(nc=loc, v=prec_def)
      nc_close(loc)
      
      loc = nc_open(paste0(datadir,y,"/",temFile, "r.nc"),write = TRUE)
      ncvar_put(nc=loc,varid="Precipitation_Flux", vals=prec)
      nc_close(loc)
      
      system(paste0("/usr/bin/ncks -x -v Band1 ",paste0(datadir,y,"/",temFile, "r.nc")," ",paste0(datadir,y,"/",temFile, "r.nc"), " -O"))
      system(paste0("/usr/bin/rm ",datadir,y,"/",temFile))
      system(paste0("/usr/bin/rm ", datadir,y, "/",temFile,".nc")) 
    }
  }
}


library(stringr)
yearList = seq(2004,2007,by=1)
yearList = append(yearList, 2019)

yearList = list(2019)
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
    

