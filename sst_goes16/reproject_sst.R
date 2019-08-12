library(rgdal)
library(raster)
library(ncdf4)

options(warn=-1)

endyr = lubridate::year(Sys.Date())
yearseq = seq(2018,endyr)
datadir = "/home/sat_ops/goesR/data/sst/"
for (y in yearseq){
  fnamelist = list.files(paste0(datadir,"raw/", y, "/"))
  for (i in fnamelist){
    if (!file.exists(paste0(datadir, "rast_sst/", y,"/SST", i))){
      print(i)
      goesfile <- nc_open(paste0(datadir,"raw/", y, "/",i))
      tryCatch({
      nc_open(goesfile)
      
      time.val = ncvar_get(goesfile, "time")
      
      sst = ncvar_get(goesfile, "SST")
      sst_raster = raster(sst)
      sst_raster = t(sst_raster)
      extent(sst_raster)=extent(-5433893, 5433893, -5433893, 5433893)
      crs(sst_raster) = "+units=m +lon_0=-75.0 +h=35786023.0 +sweep=x +proj=geos"
      
      # we use epsg:4326 or ccrs.Geodetic
      newproj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      projected_sst = projectRaster(sst_raster, crs = newproj)
      writeRaster(projected_sst, filename=paste0(datadir, "rast_sst/", y, "/SST", i), format="CDF", overwrite=TRUE,
                  varname="sst", varunit="Kelvin", longname="sea surface temperature",
                  xname="longitude", yname="latitude")
      
      dqf = ncvar_get(goesfile, "DQF")
      dqf_raster = raster(dqf)
      dqf_raster = t(dqf_raster)
      extent(dqf_raster)=extent(-5433893, 5433893, -5433893, 5433893)
      crs(dqf_raster) = "+units=m +lon_0=-75.0 +h=35786023.0 +sweep=x +proj=geos"
      print("writing dqf file")
      projected_dqf = projectRaster(dqf_raster, crs = newproj)
      writeRaster(projected_dqf, filename=paste0(datadir, "rast_dqf/", y, "/DQF", i), format="CDF", overwrite=TRUE,
                  varname="dqf", varunit="Kelvin", longname="data quality flags",
                  xname="longitude", yname="latitude")
      
      b15 = ncvar_get(goesfile, "Band15")
      b15_raster = raster(b15)
      b15_raster = t(b15_raster)
      extent(b15_raster)=extent(-5433893, 5433893, -5433893, 5433893)
      crs(b15_raster) = "+units=m +lon_0=-75.0 +h=35786023.0 +sweep=x +proj=geos"
      print("writing b15 file")
      projected_b15 = projectRaster(b15_raster, crs = newproj)
      writeRaster(projected_b15, filename=paste0(datadir, "rast_b15/", y, "/B15", i), format="CDF", overwrite=TRUE,
                  varname="b15", varunit="Kelvin", longname="Band 15 brightness_temperature",
                  xname="longitude", yname="latitude")
      
      nc_close(goesfile)
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
}

options(warn=0)
