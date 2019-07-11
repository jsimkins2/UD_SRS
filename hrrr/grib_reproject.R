# example of wget list on sat data
# /home/aps/rtDownload/script.sh

library(raster)
r = raster("Downloads/st4_pr.2019071000.01h")
newproj = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projected_raster = projectRaster(r, crs = newproj)
nc_outfile = "Downloads/ST44.nc"
writeRaster(projected_raster, filename = nc_outfile,
            format="CDF", overwrite=TRUE,varname="precip", longname="NCEP Stage IV Precip",
            xname="lon", yname="lat")


https://nomads.ncep.noaa.gov/pub/data/nccf/com/pcpanl/prod/pcpanl.20190710/
  
dfile = download.file(get("https://nomads.ncep.noaa.gov/pub/data/nccf/com/pcpanl/prod/pcpanl.20190710/ST4.2019071013.01h.gz")
                          , destfile="Downloads/testfile", method = "wget")

get("https://nomads.ncep.noaa.gov/pub/data/nccf/com/pcpanl/prod/pcpanl.20190710/")
