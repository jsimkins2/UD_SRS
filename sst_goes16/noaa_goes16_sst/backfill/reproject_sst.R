library(rgdal)
library(raster)
library(ncdf4)

options(warn=-1)
args = commandArgs(trailingOnly=TRUE)
fname = args[1]
yr = lubridate::year(Sys.Date())
datadir = "/home/sat_ops/goesR/data/noaa_sst/backfill/"

goesfile <- nc_open(paste0(datadir,"raw/", yr, "/",fname))
time.val = ncvar_get(goesfile, "time")

sst = ncvar_get(goesfile, "SST")
sst_raster = raster(sst)
sst_raster = t(sst_raster)
extent(sst_raster)=extent(-5433893, 5433893, -5433893, 5433893)
crs(sst_raster) = "+units=m +lon_0=-75.0 +h=35786023.0 +sweep=x +proj=geos"

# we use epsg:4326 or ccrs.Geodetic
newproj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projected_sst = projectRaster(sst_raster, crs = newproj)
writeRaster(projected_sst, filename=paste0(datadir, "rast_sst/", yr, "/SST", fname), format="CDF", overwrite=TRUE,
          varname="sst", varunit="Kelvin", longname="sea surface temperature",
          xname="longitude", yname="latitude")

dqf = ncvar_get(goesfile, "DQF")
dqf_raster = raster(dqf)
dqf_raster = t(dqf_raster)
extent(dqf_raster)=extent(-5433893, 5433893, -5433893, 5433893)
crs(dqf_raster) = "+units=m +lon_0=-75.0 +h=35786023.0 +sweep=x +proj=geos"
print("writing dqf file")
projected_dqf = projectRaster(dqf_raster, crs = newproj)
writeRaster(projected_dqf, filename=paste0(datadir, "rast_dqf/", yr, "/DQF", fname), format="CDF", overwrite=TRUE,
          varname="dqf", varunit="Kelvin", longname="data quality flags",
          xname="longitude", yname="latitude")

nc_close(goesfile)

options(warn=0)
