# Adding the forecast time varialbe to the nowtime and also saving 180 files in 1 huge matrix 

# testing stuff

in.path = "Downloads/"
in.file = "aqua.2017248.0905.235959.D.L3.modis.NAT.v09.1000m.nc4.360.combined"

aqua_forecast_EOF <- function(inPath, inFile, outPath){
  
  # extract the syntax of the file so we can change the dates and pull out 180 files
  doy <- substr(in.file,10,12)
  
  # get the year and month and date stuff for naming convention 
  ncyear     <- substr(in.file,6,9)
  curentDate <- strptime(paste(ncyear, doy), format="%Y %j", tz = "UTC")

  # Initialize the loop where we are going to go through 180 files and use them to predict the next 3 days using EOF
  for (i in seq(0,179,1)){
    dateTem <- strptime(paste(ncyear, doy - i), format="%Y %j", tz = "UTC")
    tem     <- nc_open(paste0(in.path,"aqua.", ncyear, doy - i, ".", format(dateTem, "%m"), format(dateTem, "%d"),
                              ".235959.D.L3.modis.NAT.v09.1000m.nc4.360.combined"))
    assign(x=paste0("sst",lubridate::yday(dateTem)), value=ncdf4::ncvar_get(tem,"sst"))
  }
  
  # do the EOF stuff now
  
  
  
  
}
