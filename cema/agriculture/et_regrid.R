# this script is designed to process a ncep stage IV file for it's shape / size and resample a ref et value to this size
library(raster)
library(ncdf4)

####################################################################
############ Housekeeping #################
####################################################################

# load in the most recent datetime.txt - this will be used to set the path for the NC outfile
datetimeET = read.table(paste0(temPath, "datetime.txt"),stringsAsFactor = FALSE)[[1]]
epochDatetime = as.numeric(as.POSIXct(x = paste0(datetimeET, " 23:59:00"), format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01", tz="UTC"))
tunits = "seconds since 1970-01-01"

# set the paths
temPath = "/Users/james/Downloads/"
outPathNC = "/Users/james/Downloads/"

####################################################################
############ Raster Operations #################
####################################################################

# open the most recent ET file 
ETrast = raster(paste0(temPath, "ETtemp.tif"))

# Open the default Precip File that we want to replicate the grid of 
precRast = raster(paste0(temppath, "ST4.2020030311.01hr.nc"))
precRast = crop(precRast, extent(-76.5, -74.2, 38.1, 40.8))

# resample the ETrast to the precRast grid - bilinear by default, works the same way as mean
resET = resample(ETrast, precRast)
resET = flip(resET, "y")

# output the raster as netCDF
writeRaster(resET, paste0(temPath, "temp.nc"),format="CDF", overwrite=TRUE,
            varname="refET", varunit="mm day-1", longname="Reference Evapotranspiration",
            xname="longitude", yname="latitude")

####################################################################
############ Make the netCDF file CF compliant #################
####################################################################
# name our nc file 
ncfname = paste0(outPathNC, "DEOS_refET_", substr(datetimeET, 1,4),substr(datetimeET, 6,7),substr(datetimeET, 9,10), ".nc")

# reopen the nc file we just created
tempNC = nc_open(paste0(temPath, "temp.nc"))
tempET = ncvar_get(tempNC, "refET")
dims = tempNC$dim

# define dimensions
londim <- ncdim_def("longitude","degrees_east",as.double(dims$longitude$vals)) 
latdim <- ncdim_def("latitude","degrees_north",as.double(dims$latitude$vals))
timedim <- ncdim_def("time",tunits,as.double(epochDatetime))

fillvalue <- 1e32
dlname <- "Reference Evapotranspiration"
et_def <- ncvar_def("refET","mm day-1",list(londim,latdim,timedim),fillvalue,dlname,prec="single")

# create netCDF file and put arrays
ncout <- nc_create(ncfname,list(et_def),force_v4=TRUE)

# put variables
ncvar_put(ncout,et_def,tempET)

# add global attributes
ncatt_put(ncout,0,"title","DEOS Reference Evapotranspiration - Spatiall Interpolated")
ncatt_put(ncout,0,"institution","Center for Environmental Monitoring and Analysis - University of Delaware")
ncatt_put(ncout,0,"source","Delaware Environmental Observing System (DEOS)")
ncatt_put(ncout,0,"author","James Simkins")
ncatt_put(ncout,0,attname="proj4", attval=proj4string(resET))
ncatt_put(ncout,0,attname="crs", attval="EPSG:4326")
ncatt_put(ncout,0,"Conventions","CF-1.4")
# close the file, writing data to disk
nc_close(ncout)

