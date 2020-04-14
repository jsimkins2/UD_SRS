# this script processes previously downloaded files from the NCEP Stage IV precip database -
# https://data.eol.ucar.edu/cgi-bin/codiac/fgr_form/id=21.093
# this is the backfill script for the 24h files


library(lubridate)
library(ncdf4)

# Unzip all of the downloaded files
#datadir = "/Users/james/Downloads/ncep24Z/"
#dirList = list(200201, 200202, 200301, 200302, 200401, 200402, 200501, 200502, 200601, 200602, 200701,
#               200702, 200801, 200802, 200901, 200902, 201001, 201002, 201101, 201102, 201201, 201202,
#               201301, 201302, 201401, 201402, 201501, 201502, 201601, 201602, 201701, 201702, 201703,
#               201801, 201802, 201803, 201901, 201902, 201903, 202001)

#for (d in dirList){
#  system(paste0("/bin/mv /Users/james/Downloads/",d,"/*24h* /Users/james/Downloads/ncep24/"))
#}


# reproject and add proper dimensions
datadir = "/home/james/ncep_stageIV/backfill_ncep/"
fileList = list.files(paste0(datadir, "notNetcdf"))
for (f in fileList){
  print(f)
  temFile = f 
  dateTime = ymd_hms(paste0(substr(f, 5,12), "000000"))
  # convert to netcdf and reproject
  system(paste0("/usr/local/bin/gdal_translate ",datadir,"notNetcdf/",temFile, " ",datadir,'unprojected/', temFile, ".nc -of netcdf"))
  system(paste0("/usr/local/bin/gdalwarp ",datadir,'unprojected/', temFile,".nc /data/ncep_stageIV/quality/",year(dateTime),"/",temFile, "r.nc -of netcdf -t_srs epsg:4326"))
  
  # open up the file 
  f=paste0(temFile,".nc")
  timeVal = as.POSIXct(ymd_hms(paste0(substr(f, 5,12), "000000")))
  loc = nc_open(paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",temFile, "r.nc"),write = TRUE)
  prec = ncvar_get(nc=loc, varid="Band1")
  timedim = ncdim_def("time", "seconds since 1970-01-01",as.integer(timeVal))
  londim = loc$dim$lon
  latdim = loc$dim$lat
  
  prec_def = ncvar_def(name="Precipitation_Flux", units="mm/hr", missval=loc$var$Band1$missval, longname = "Quantitative Precipitation Estimates", dim=list(londim,latdim,timedim))
  ncdf4::ncvar_add(nc=loc, v=prec_def)
  nc_close(loc)
  
  loc = nc_open(paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",temFile, "r.nc"),write = TRUE)
  ncvar_put(nc=loc,varid="Precipitation_Flux", vals=prec)
  nc_close(loc)
  
  system(paste0("/usr/bin/ncks -x -v Band1 ",paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",temFile, "r.nc")," ",paste0("/data/ncep_stageIV/quality/",year(dateTime),"/",temFile, "r.nc"), " -O"))
}

