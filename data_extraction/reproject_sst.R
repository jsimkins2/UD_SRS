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
    if (!file.exists(paste0(datadir, "reprojected/", y,"/", i))){
      print(i)
      goesfile <- nc_open(paste0(datadir,"raw/", y, "/",i))
      time.val = ncvar_get(goesfile, "time")
      
      sst = ncvar_get(goesfile, "SST")
      sst_raster = raster(sst)
      sst_raster = t(sst_raster)
      extent(sst_raster)=extent(-5433893, 5433893, -5433893, 5433893)
      crs(sst_raster) = "+units=m +lon_0=-75.0 +h=35786023.0 +sweep=x +proj=geos"
      
      # we use epsg:4326 or ccrs.Geodetic
      newproj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      projected_sst = projectRaster(sst_raster, crs = newproj)
      writeRaster(projected_sst, filename=paste0(datadir, "sst.nc"), format="CDF", overwrite=TRUE,
                  varname="sst", varunit="Kelvin", longname="sea surface temperature",
                  xname="longitude", yname="latitude")
      
      dqf = ncvar_get(goesfile, "DQF")
      dqf_raster = raster(dqf)
      dqf_raster = t(dqf_raster)
      extent(dqf_raster)=extent(-5433893, 5433893, -5433893, 5433893)
      crs(dqf_raster) = "+units=m +lon_0=-75.0 +h=35786023.0 +sweep=x +proj=geos"
      print("writing dqf file")
      projected_dqf = projectRaster(dqf_raster, crs = newproj)
      writeRaster(projected_dqf, filename=paste0(datadir, "dqf.nc"), format="CDF", overwrite=TRUE,
                  varname="dqf", varunit="Kelvin", longname="data quality flags",
                  xname="longitude", yname="latitude")
      
      nc_close(goesfile)
    }
  }
}

options(warn=0)
