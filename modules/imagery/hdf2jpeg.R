#this function assumes that the input VIIRS files and output jpeg will be in the same place
library(gdalUtils)
library(raster)
library(jpeg)
library(tiff)
path=as.character(commandArgs(1))
ofile=as.character(commandArgs(2))

sds = get_subdatasets(paste0(path,"/",ofile,".hdf"))
gdal_translate(sds[3], dst_dataset = paste0(path,"/",ofile,".tiff"))
img <- readTIFF(paste0(path,"/",ofile,".tiff"), native=TRUE)
writeJPEG(img, target = paste0(path,"/",ofile,".jpeg"), quality = 1)
  

