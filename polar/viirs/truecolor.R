library(ncdf4)
library(raster)
goesnc = nc_open("http://thredds.demac.udel.edu/thredds/dodsC/GOESR_FD.nc?band01[1][1:962][1:1011],band02[1][1:962][1:1011],band03[1][1:962][1:1011]")

band01 = ncvar_get(goesnc, "band01")
band02 = ncvar_get(goesnc, "band02")
band03 = ncvar_get(goesnc, "band03")
lat = goesnc$dim$hires_lat$vals
lon = goesnc$dim$hires_lon$vals
x = raster(band03)
plot(x)
as.POSIXct(1515673841)

as.POSIXct(1515673841, origin="1970-01-01")

r <- raster(ncol=10, nrow=10)
values(r) <- sample(0:255, ncell(r), replace=TRUE)
ctab <- sample(rainbow(256))
colortable(r) <- ctab
plot(r)
head(colortable(r)) 


x = c(as.character(band01[1,1]), band02[1,1], band03[1,1])
sapply(strsplit(x, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

rgb(band01[1,1], band02[1,1], band03[1,1], maxColorValue = 255)

# make a colored matrix out of these like above and then rasterize that bitch and then plot her up!