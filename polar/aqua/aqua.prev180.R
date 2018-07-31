# This script talks to basin and finds the most recent time we have data for and creates a netCDF3 file with that time and the previous 180 days
outPath = ""
verbose = FALSE

aqua.file <- ncdf4::nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/Aqua1DayAggregate.nc")
rec_time  <- max(aqua.file$dim$time$vals)
rec_len   <- aqua.file$dim$time$len - 1

# convert epoch seconds to epoch days
epoch_days <- floor(unclass(rec_time)/86400)
epoch_days180 <- epoch_days - 180

# get array value mins for lat/lons
lat.min <- min(which(aqua.file$dim$lat$vals > 38.45)) - 1
lon.min <- max(which(aqua.file$dim$lon$vals < -75.7)) - 1

# get array value maxs for lat/lons
lat.max <- max(which(aqua.file$dim$lat$vals < 39.7)) + 1
lon.max <- min(which(aqua.file$dim$lon$vals > -74.7)) - 1

# Preparing the values for ncdf dimensions
lat.vals  <- seq(aqua.file$dim$lat$vals[lat.min], aqua.file$dim$lat$vals[lat.max], len=(lat.max-lat.min))
lon.vals  <- seq(aqua.file$dim$lon$vals[lon.min], aqua.file$dim$lon$vals[lon.max], len=(lon.max-lon.min))
time.vals <- seq(rec_time - 180*86400, rec_time, 86400)
  
# Dimensions
lat <- ncdf4::ncdim_def(name='latitude', units='degree_north', vals=lat.vals, create_dimvar=TRUE)
lon <- ncdf4::ncdim_def(name='longitude', units='degree_east', vals=lon.vals, create_dimvar=TRUE)
time <- ncdf4::ncdim_def(name='time', units="seconds since 1970-01-01T00:00:00Z", vals=time.vals, create_dimvar=TRUE, unlim=TRUE)
dim<-list(lat,lon,time)

# declare variables of interest
vars.info <- data.frame(var.name = c("sst", "a_443_qaa", "Rrs_555", "chl_oc3"),
                        CF.name  = c("Sea Surface Temperature", "Total Absorption at 443 nm, QAA Agorithm", 
                                     "Remote Sensing Reflectance at 555 nm","Chlorophyll Concentration OC3"),
                        units    = c("Degree_C", "m^-1", "sr^-1", "mg m^-3"))

# Initialize lists for data
dat.list <- list()
var.list <- list()

# Extract the data from Aqua1DayAggregate
for (j in seq_along(vars.info$var.name)){
  dat.list[[j]] <- ncdf4::ncvar_get(aqua.file, as.character(vars.info$var.name[[j]]), start = c(lon.min,lat.min, rec_len - 180), count = c(lon.max-lon.min,lat.max-lat.min,181))
  var.list[[j]] <- ncdf4::ncvar_def(name=as.character(vars.info$CF.name[j]), units=as.character(vars.info$units[j]), dim=dim, missval=-999, verbose=verbose)
}
ncdf4::nc_close(aqua.file)


# jobs done here, let's head out
loc.file <- paste0(outPath,"/aqua.gapfilled.nc3")

#writing all we need to the output file
loc <- ncdf4::nc_create(filename=loc.file, vars=var.list)
for(j in seq_along(vars.info$CF.name)){
  ncdf4::ncvar_put(nc=loc, varid=as.character(vars.info$CF.name[j]), vals=dat.list[[j]])
}
ncdf4::nc_close(loc)
