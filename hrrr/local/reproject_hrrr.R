library(rgdal)
library(raster)
library(ncdf4)
verbose=FALSE

datadir = '/Users/james/Downloads/hrrr_temp/'
args = commandArgs(trailingOnly=TRUE)
fname = args[1]
hrrrfile = nc_open(paste0(datadir,fname))
x = ncvar_get(hrrrfile, 'x')
y = ncvar_get(hrrrfile,'y')

var <- data.frame(name <- c("temperature_surface", "temperature_925", "temperature_850", "height_500", "height_1000"),
                  longname <- c("air_temperature_at_surface","air_temperature_at_925mb","air_temperature_at_850mb",
                               "isobaric_height_500mb","isobaric_height_1000mb"),
                  units <- c('Kelvin','Kelvin',"Kelvin","gpm","gpm"),stringsAsFactors = FALSE)

var.list <- list()
dat.list <- list()
for(j in 1:length(var$name)){
  varn = ncvar_get(hrrrfile, var$name[j])
  origraster = raster(varn)
  origraster = t(origraster)
  origraster = flip(origraster, direction = 'y')
  extent(origraster)=extent(min(x), max(x), min(y), max(y))
  crs(origraster) = "+ellps=sphere +a=6371229.0 +b=6371229.0 +proj=lcc +lon_0=265.0 +lat_0=25.0 +x_0=0.0 +y_0=0.0 +lat_1=25.0 +no_defs"
  # we use epsg:4326 or ccrs.Geodetic
  newproj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  projected_raster = projectRaster(origraster, crs = newproj)
  writeRaster(projected_raster, filename=paste0(datadir, var$name[j], ".nc"), format="CDF", overwrite=TRUE,
              varname=var$name[j], varunit=var$units[j], longname=var$longname[j],
              xname="longitude", yname="latitude")

  # now that the raster is written to a netcdf, open up the netcdf and pull that information
  ncrast = nc_open(paste0(datadir,var$name[j], ".nc" ))
  temdata = ncvar_get(ncrast, var$name[j])
  if (j==1){
    dim = ncrast$dim
  }
  nc_close(ncrast)
  
  var.list[[j]] <- ncdf4::ncvar_def(name=as.character(var$name[j]), units=as.character(var$units[j]),longname = as.character(var$longname[j]),
                                    dim=dim, missval=-9999.0, verbose=verbose)
  dat.list[[j]] <- temdata
}
nc_close(hrrrfile)

## put data in new file
loc <- nc_create(filename=paste0(datadir, 'reprojected', fname), vars=var.list, verbose=verbose)
for(j in seq_along(var$name)){
  ncdf4::ncvar_put(nc=loc, varid=as.character(var$name[j]), vals=dat.list[[j]])
}
ncatt_put(nc=loc,0,attname="proj4", attval=newproj)
ncatt_put(nc=loc,0,attname="crs", attval="Geodetic")
ncatt_put(nc=loc, 0, "Conventions", "CF=1.0")
nc_close(loc)



