# this script is designed to process a ncep stage IV file for it's shape / size and resample a ref et value to this size
library(raster)
library(ncdf4)

####################################################################
############ Housekeeping #################
####################################################################
# load in dtime
args = commandArgs(trailingOnly=TRUE)
datetimeET = args[1]
epochDatetime = as.numeric(as.POSIXct(x = paste0(datetimeET, " 23:59:00"), format = "%Y%m%d %H:%M:%S", origin = "1970-01-01", tz="UTC"))
tunits = "seconds since 1970-01-01"

# set the paths
temPath = "/home/sat_ops/deos/temp/"
outPathNC = "/data/DEOS/agriculture/"

# Set variable & units dataframe
var.df = data.frame(names = c('meanTemp', 'maxTemp', 'minTemp', 'HDD', 'CDD', 'meanWS', 'dailyprecip', 'meanDP','energyDens',
                      'refET', 'GDD', 'meanRH', 'maxRH', 'minRH', 'meanST', 'maxST', 'minST', 'meanVWC', 'maxVWC', 'minVWC',
                    'meanSolar', 'meanWD', 'dailyGust', 'dailyMinWC', 'maxHI'),
                    units = c('Kelvin', 'Kelvin', 'Kelvin', ' ', ' ', 'm.s-1', 'mm', 'Kelvin', 'J.m-2',
                                 'mm.day-1', ' ', '%', '%', '%', 'Kelvin', 'Kelvin', 'Kelvin', ' ', ' ', ' ', 'J.m-2',
                                 'Rad', 'm.s-1', 'Kelvin', 'Kelvin'),
                    longname = c("Mean Air Temperature", "Max Air Temperature", "Min Air Temperature", "Heating Degree Day",
                                 "Cooling Degree Day", "Mean Wind Speed", "Daily Aggregated Precipitation Total", "Mean Dew Point",
                                 "Energy Density", "Reference Evapotranspiration", "Growing Degree Day", "Mean Relative Humidity", 
                                 "Max Relative Humidity", "Min Relative Humidity", "Mean Soil Temperature", "Max Soil Temperature",
                                 "Min Soil Temperature", "Mean Volumetric Water Content", "Max Volumetric Water Content",
                                 "Min Volumetric Water Content", "Mean Solar Energy", "Mean Wind Direction", "Daily Wind Gust",
                                 "Daily Min Wind Chill", "Daily Max Heat Index"))

####################################################################
############ Raster Operations #################
####################################################################

for (pyVar in var.df$names){
  # open the most recent ET file 
  if (file.exists(paste0(temPath, pyVar, "_temp.tif"))){
    ETrast = raster(paste0(temPath, pyVar, "_temp.tif"))
  } else { #if the file doesn't exist, create a blank raster to populate the variable
    ETrast = raster(paste0(temPath, "ST4.2020030311.01hr.nc"))
    ETrast = crop(ETrast, extent(-76.5, -74.2, 38.1, 40.8))
    ETrast[ETrast > -1000] = NA
  }

  
  # Open the default Precip File that we want to replicate the grid of 
  precRast = raster(paste0(temPath, "ST4.2020030311.01hr.nc"))
  precRast = crop(precRast, extent(-76.5, -74.2, 38.1, 40.8))
  
  # resample the ETrast to the precRast grid - bilinear by default, works the same way as mean
  resET = resample(ETrast, precRast)
  resET = flip(resET, "y")
  
  # output the raster as netCDF
  writeRaster(resET, paste0(temPath, pyVar, "_temp.nc"),format="CDF", overwrite=TRUE,
              varname=pyVar, varunit=as.character(var.df$units[which(var.df$names == pyVar)]), 
              longname=as.character(var.df$longname[which(var.df$names == pyVar)]),
              xname="longitude", yname="latitude")
}

####################################################################
############ Make the netCDF file CF compliant #################
####################################################################
# name our nc file 
ncfname = paste0(outPathNC, "/", substr(datetimeET, 1,4), "/DEOS_agri_", substr(datetimeET, 1,4),substr(datetimeET, 5,6),substr(datetimeET, 7,8), ".nc")

# create nc file with definitions for each variable
def.list = list()
for (v in seq_along(var.df$names)){
  tempNC = nc_open(paste0(temPath, as.character(var.df$names[v]), "_temp.nc"))
  tempET = ncvar_get(tempNC, as.character(var.df$names[v]))
  dims = tempNC$dim
  londim <- ncdim_def("longitude","degrees_east",as.double(dims$longitude$vals)) 
  latdim <- ncdim_def("latitude","degrees_north",as.double(dims$latitude$vals))
  timedim <- ncdim_def("time",tunits,as.double(epochDatetime))
  fillvalue <- -9999
  dlname <- as.character(var.df$longname[which(var.df$names == as.character(var.df$names[v]))])
  def.list[[v]] <- ncvar_def(as.character(var.df$names[v]),as.character(var.df$units[v]),
                         list(londim,latdim,timedim),fillvalue,dlname,prec="single")
  nc_close(tempNC)
}

# create the netcdf file we're writing to
ncout <- ncdf4::nc_create(filename = ncfname, vars = def.list)

for (v in seq_along(var.df$names)){
  # reopen the nc files we just created
  tempNC = nc_open(paste0(temPath, as.character(var.df$names[v]), "_temp.nc"))
  tempET = ncvar_get(tempNC, as.character(var.df$names[v]))
  # put variables
  ncvar_put(ncout,def.list[[v]],tempET)
}

# add global attributes
ncatt_put(ncout,0,"title",paste0("DEOS Agrictulural Variables Data Product",
                                 " - Spatially Interpolated from DEOS stations using Inverse Distance Weighted Algorithm"))
ncatt_put(ncout,0,"institution","Center for Environmental Monitoring and Analysis - University of Delaware")
ncatt_put(ncout,0,"source","Delaware Environmental Observing System (DEOS)")
ncatt_put(ncout,0,"author","James Simkins")
ncatt_put(ncout,0,attname="proj4", attval=proj4string(resET))
ncatt_put(ncout,0,attname="crs", attval="EPSG:4326")
ncatt_put(ncout,0,"Conventions","CF-1.4")
# close the file, writing data to disk
nc_close(ncout)

