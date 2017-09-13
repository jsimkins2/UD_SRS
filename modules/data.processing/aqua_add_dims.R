# This script adds now_time and forecast_times to Aqua dineof NC files
# This also separates the forecasted variables from the observations

aqua_add_dims <- function(inPath, inFile, outPath, verbose = FALSE, ...){
  
  # Open the netcdf file
  aqua.nc <- ncdf4::nc_open(paste0(inPath,inFile))
  
  # list the info for the nc file
  nc.info <- list(var.name = c("sst", "gapfilled_sst", "a_443_qaa", "gapfilled_a_443_qaa", "forecast_sst", "forecast_a_443_qaa"),
                  CF.name = c("sea_surface_temperature", "gapfilled_sea_surface_temperature", "a_443_qaa", "gapfilled_a_443_qaa",
                              "forecast_sst", "forecast_a_443_qaa"),
                  longname = c("Sea Surface Temperature", "Gapfilled Sea Surface Temperature", "
                               Total Absorption at 443 nm, QAA Agorithm", "Gapfilled Total Absorption at 443 nm, QAA Algorithm",
                               "3 Day Forecast of SST", "3 Day Forecast of a_443_qaa"), 
                  units = c("Degree_C", "Degree_C", "m^-1", "m^-1", "Degree_C", "m^-1"))

  # Initialize the lists where we are going to store stuff
  dat.list <- list()
  
  # Extract the data from the input file
  for (j in seq_along(nc.info$var.name)){
    if (j < 5){
      dat.list[[j]] <- ncdf4::ncvar_get(aqua.nc, as.character(nc.info$var.name[[j]]))
    } else {
      dat.list[[j]] <- NA
    }
  }
  
  # Moving the forecasted values to a new variable
  names(dat.list) = nc.info$var.name
  dat.list[["forecast_sst"]] = dat.list[["gapfilled_sst"]][,,181:183]
  dat.list[["forecast_a_443_qaa"]] = dat.list[["gapfilled_a_443_qaa"]][,,181:183]
  
  # Resizing the forecasted arrays
  dat.list[["forecast_sst"]] = array(dat.list[["forecast_sst"]], dim = c(102,164,180))
  dat.list[["forecast_a_443_qaa"]] = array(dat.list[["gapfilled_a_443_qaa"]], dim = c(102,164,180))
  
  # Removing the forecasted days from the observations
  for (j in seq_along(dat.list)){
    if (j < 5){
      dat.list[[j]] = dat.list[[j]][,,1:180]
    }
  }
  
  #--------DIMENSIONS-------------
  
  # Extract the current file time and write to a dimension in posixct time. Note: this only works for the current storage syntax from basin
  ncyear  <- substr(inFile,6,9)
  ncmon   <- substr(inFile,14,15)
  ncday   <- substr(inFile,16,17)
  nchr    <- substr(inFile,19,20)
  ncmin   <- substr(inFile,21,22)
  ncsec   <- substr(inFile,23,24)
  ncposix <- as.POSIXct(strptime(paste0(ncyear,"-",ncmon,"-",ncday," ", nchr, ":", ncmin, ":",ncsec), "%Y-%m-%d %H:%M:%S", tz="GMT"))
  
  # Create the now_time dimension
  now_time <- ncdf4::ncdim_def(name='now_time', units="seconds since 1970-01-01T00:00:00Z", vals=as.numeric(ncposix), create_dimvar=TRUE, unlim=TRUE)
  
  # Add forecast time dimension for 3 days beyond the now_time, had to do these individually otherwise we had an error for
  # trying to write too many values
  fcst_val1 <- ncposix + 1*(86400)
  fcst_val2 <- ncposix + 2*(86400)
  fcst_val3 <- ncposix + 3*(86400)
  forecast_day1 <- ncdf4::ncdim_def(name='forecast_day1', units="seconds since 1970-01-01T00:00:00Z", vals=as.numeric(fcst_val1), 
                                    create_dimvar=TRUE, unlim=TRUE)
  forecast_day2 <- ncdf4::ncdim_def(name='forecast_day2', units="seconds since 1970-01-01T00:00:00Z", vals=as.numeric(fcst_val2), 
                                    create_dimvar=TRUE, unlim=TRUE)
  forecast_day3 <- ncdf4::ncdim_def(name='forecast_day3', units="seconds since 1970-01-01T00:00:00Z", vals=as.numeric(fcst_val3), 
                                    create_dimvar=TRUE, unlim=TRUE)
  
  # define forecast_day by lat/lon to place in one string
  # Write each dimension to the dim list
  aqua.nc$dim$time$vals  <- aqua.nc$dim$time$vals[4:183]
  aqua.nc$dim$time$len <- 180
  dim <- list(forecast_day3, forecast_day2, forecast_day1, now_time, aqua.nc$dim$time, aqua.nc$dim$lon, aqua.nc$dim$lat)

  var.list <- list()
  
  for (j in seq_along(nc.info$var.name)){
      var.list[[j]] <- ncdf4::ncvar_def(name=as.character(nc.info$CF.name[[j]]), 
                                        units=as.character(nc.info$units[[j]]), 
                                        dim=dim, missval=NA, verbose=verbose)
  }
  ncdf4::nc_close(aqua.nc)
  
  #creating the nc4 output file
  loc.file <- paste0(outPath,"/aqua.gapfilled.",ncyear,ncmon,ncday,".nc4")
  
  #writing all we need to the output file
  loc <- ncdf4::nc_create(filename=loc.file, vars=var.list, verbose=verbose)
  for(j in seq_along(nc.info$CF.name)){
    ncdf4::ncvar_put(nc=loc, varid=as.character(nc.info$CF.name[j]), vals=dat.list[[j]])
  }
  ncdf4::nc_close(loc)
}

  
