# Document detailing the downloading/processing/reprojection/aggregation of NCEP Stage IV data

 #1) Download the archived data here - https://data.eol.ucar.edu/dataset/21.093
 #   (real time data here - https://nomads.ncep.noaa.gov/pub/data/nccf/com/pcpanl/prod/)
 
 #2) Untar all file folders. Go into each folder, remove .txt files, remove *06h* files, remove *24h*, and uncompress all. 
 #   Remove duplicate files (20040101 file in a 2003 folder)
 
 #3) 



# Now the script begins. 

library(lubridate)
library(ncdf4)
yearList = seq(2003,2007,by=1)
yearList = append(yearList, 2019)
datadir = "~/Documents/ncep/"

for (y in yearList){
  fileNames = list.files(paste0(datadir, y, "/"))
  for (fi in fileNames){
    temFile = fi
    system(paste0("/usr/local/bin/gdal_translate ",datadir,temFile, " ",datadir, temFile, ".nc -of netcdf"))
    system(paste0("/usr/local/bin/gdalwarp ",datadir, temFile,".nc" , datadir,"/",temFile, ".nc -of netcdf -t_srs epsg:4326"))
    
    f=paste0(temFile,".nc")
    timeVal = as.POSIXct(ymd_h(substr(f,5,14)))
    loc = nc_open(paste0("/data/ncep_stageIV/",year(dateTime),"/",temFile, ".nc"),write = TRUE)
    prec = ncvar_get(nc=loc, varid="Band1")
    timedim = ncdim_def("time", "seconds since 1970-01-01",as.numeric(timeVal))
    londim = loc$dim$lon
    latdim = loc$dim$lat
    
    prec_def = ncvar_def(name="Precipitation_Flux", units="mm/hr", missval=loc$var$Band1$missval, longname = "Quantitative Precipitation Estimates", dim=list(londim,latdim,timedim))
    ncdf4::ncvar_add(nc=loc, v=prec_def)
    nc_close(loc)
    
    loc = nc_open(paste0("/data/ncep_stageIV/",year(dateTime),"/",temFile, ".nc"),write = TRUE)
    ncvar_put(nc=loc,varid="Precipitation_Flux", vals=prec)
    nc_close(loc)
    
    system(paste0("/usr/bin/ncks -x -v Band1 ",paste0("/data/ncep_stageIV/",year(dateTime),"/",temFile, ".nc")," ",paste0("/data/ncep_stageIV/",year(dateTime),"/",temFile, ".nc"), " -O"))
  }
}