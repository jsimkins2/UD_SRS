# This script serves to add a time variable to the aqua file from http://basin.ceoe.udel.edu/data/dineof/real_time/
# The goal of this is to save this to a Thredds server which can't actually be done at the moment
# We need to make this data indexable by Thredds
# Ideally we do EPOC days or seconds and create the new time dimension
# Update 09/11/2017 - we have added the nowtime and the forecast time dimensions, but we need forecast time values now

# ------------- Begin Code

aqua_add_dims <- function(inPath, inFile, outPath, verbose = FALSE, ...){
  # Open the netcdf file
  aqua.nc <- ncdf4::nc_open(paste0(inPath,inFile))
  
  # list the info for the nc file
  nc.info <- list(var.name = c("sst", "gapfilled_sst", "a_443_qaa", "gapfilled_a_443_qaa"),
                  CF.name = c("sea_surface_temperature", "gapfilled_sea_surface_temperature", "a_443_qaa", "gapfilled_a_443_qaa"),
                  longname = c("Sea Surface Temperature", "Gapfilled Sea Surface Temperature", "
                               Total Absorption at 443 nm, QAA Agorithm", "Gapfilled Total Absorption at 443 nm, QAA Algorithm"), 
                  units = c("Degree_C", "Degree_C", "m^-1", "m^-1"))
  
  # Extract the current file time and write to a dimension in posixct time. Note: this only works for the current storage syntax from basin
  ncyear  <- substr(inFile,6,9)
  ncmon   <- substr(inFile,14,15)
  ncday   <- substr(inFile,16,17)
  nchr    <- substr(inFile,19,20)
  ncmin   <- substr(inFile,21,22)
  ncsec   <- substr(inFile,23,24)
  ncposix <- as.POSIXct(strptime(paste0(ncyear,"-",ncmon,"-",ncday," ", nchr, ":", ncmin, ":",ncsec), "%Y-%m-%d %H:%M:%S", tz="GMT"))
  
  # Create the now_time dimension
  now_time <- ncdim_def(name='now_time', units="seconds since 1970-01-01T00:00:00Z", vals=as.numeric(ncposix), create_dimvar=TRUE, unlim=TRUE)
  
  # Add forecast time dimension for 3 days beyond the now_time, had to do these individually otherwise we had an error for
  # trying to write too many values
  fcst_val1 <- ncposix + 1*(86400)
  fcst_val2 <- ncposix + 2*(86400)
  fcst_val3 <- ncposix + 3*(86400)
  forecast_day1 <- ncdim_def(name='forecast_day1', units="seconds since 1970-01-01T00:00:00Z", vals=as.numeric(fcst_val1), create_dimvar=TRUE, unlim=TRUE)
  forecast_day2 <- ncdim_def(name='forecast_day2', units="seconds since 1970-01-01T00:00:00Z", vals=as.numeric(fcst_val2), create_dimvar=TRUE, unlim=TRUE)
  forecast_day3 <- ncdim_def(name='forecast_day3', units="seconds since 1970-01-01T00:00:00Z", vals=as.numeric(fcst_val3), create_dimvar=TRUE, unlim=TRUE)
  
  # define forecast_day by lat/lon to place in one string
  # Write each dimension to the dim list
  dim <- list(aqua.nc$dim$lon, aqua.nc$dim$lat, aqua.nc$dim$time, now_time, forecast_day1, forecast_day2, forecast_day3)
  
  # Initialize the lists where we are going to store stuff
  dat.list <- list()
  var.list <- list()
  
  # Extract the data from the input file
  for (j in seq_along(nc.info$var.name)){
    dat.list[[j]] <- ncdf4::ncvar_get(aqua.nc, as.character(nc.info$var.name[[j]]))
    var.list[[j]] <- ncdf4::ncvar_def(name=as.character(nc.info$CF.name[[j]]), 
                                     units=as.character(nc.info$units[[j]]), 
                                     dim=dim, missval=-9999, verbose=verbose)
  }
  
  #close the nc file
  ncdf4::nc_close(aqua.nc)
  
  #creating the nc4 output file
  loc.file <- paste0(outPath,"aqua.gapfilled.",ncyear,ncmon,ncday,".nc4")
  
  #writing all we need to the output file
  loc <- ncdf4::nc_create(filename=loc.file, vars=var.list, verbose=verbose)
  for(j in seq_along(nc.info$CF.name)){
    ncdf4::ncvar_put(nc=loc, varid=as.character(nc.info$CF.name[j]), vals=dat.list[[j]])
  }
  ncdf4::nc_close(loc)
}

  
