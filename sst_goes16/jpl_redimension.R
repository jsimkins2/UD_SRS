library(lubridate)
library(ncdf4)
yearList = seq(2017,2020,by=1)
datadir = "/data/GOES/GOES-R/jpl_goes/"
for (y in yearList){
  fileNames = list.files(paste0(datadir, y, "/"))
  for (f in fileNames){
    print(f)

    loc = nc_open(paste0(datadir,y,"/",f))
    timeVal = loc$dim$time$vals
    timedim = ncdim_def("time", "seconds since 1970-01-01",as.numeric(timeVal))
    londim = loc$dim$lon
    latdim = loc$dim$lat
    
    dim <- list(latdim, londim, timedim)
    var = names(loc$var)
    
    dat.list = list()
    var.list = list()


    for (j in seq_along(var)) {
      dat.list[[j]] <- ncdf4::ncvar_get(loc, var[j])
      var.list[[j]] <- ncdf4::ncvar_def(name = loc$var[var[j]][[1]]$name, 
                                        units = loc$var[var[j]][[1]]$units, 
                                        dim = dim,
                                        missval = loc$var[var[j]][[1]]$missval)
      }

    ## put data in new file
    newnc.file = paste0(datadir,y,"/",substr(f, 1,26), "s.nc")
    newnc <- ncdf4::nc_create(filename = newnc.file, vars = var.list)
    
    for (j in seq_along(var))  {
      ncdf4::ncvar_put(nc = newnc, varid = as.character(var[j]), vals = dat.list[[j]])
    }
    ncdf4::nc_close(loc)
    
    