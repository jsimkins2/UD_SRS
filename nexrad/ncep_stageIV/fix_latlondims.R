# this script is intended to fix the lat/lon discrepencies between newer and older files.
# I ran gdal on two different machines here and the precision was different, leading
# ERDDAP to freak out

# this is also a template for fixing dimension issues in R
library(ncdf4)

yearsSeq = seq(2008,2019)
homedir = '/data/ncep_stageIV/'
baseFile = nc_open('/data/ncep_stageIV/2019/ST4.2019091913.01hr.nc')
baseLat = baseFile$dim$lat$vals
baseLon = baseFile$dim$lon$vals
for (y in yearsSeq){
  fnames = list.files(paste0(homedir, y, '/'))
  for (f in fnames){
    temFile = nc_open(paste0(homedir, y, '/', f), write=T)
    temLat=ncvar_get(temFile,'lat')
    temLon = ncvar_get(temFile, 'lon')
    ncvar_put(nc = temFile, varid='lat', vals=baseLat, start=c(1), count=c(length(baseLat)))
    ncvar_put(nc = temFile, varid='lon', vals=baseLon, start=c(1), count=c(length(baseLon)))
    nc_close(temFile)
  }
}
