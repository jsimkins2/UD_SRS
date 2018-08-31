library(ncdf4)
filenames = read.table("/home/james/files.csv", header = FALSE, stringsAsFactors = FALSE, sep = ",")

fnames = as.data.frame(filenames)
for (i in fnames$V1){
  aqua.nc <- nc_open(paste0("/data/Aqua/8_day/2018/",i))
  
  jday = substr(i, 10, 12)              
  data_dim <- aqua.nc$dim
  sst <- ncvar_get(aqua.nc, "sst")
  chl <- ncvar_get(aqua.nc, "chl_oc3")
  aquatime <- aqua.nc$dim$time
  aqualat <- aqua.nc$dim$lat
  aqualon <- aqua.nc$dim$lon
  lon <- ncvar_get(aqua.nc, "lon")
  lat <- ncvar_get(aqua.nc, "lat")
  nc_close(aqua.nc)
  
  sst[sst < 0] = 0
  
  # read in the specific 8 day segment of the 8 day climatology
  clim.nc = nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/aqua_clim_rolling8.nc")
  rec_time  <- max(clim.nc$dim$time$vals)
  rec_len   <- clim.nc$dim$time$len
  
  sst.c <- ncvar_get(clim.nc, "sst",start = c(1,1,as.numeric(jday)), count = c(clim.nc$dim$lon$len,clim.nc$dim$lat$len,1))
  chl.c <- ncvar_get(clim.nc, "chl_oc3",start = c(1,1,as.numeric(jday)), count = c(clim.nc$dim$lon$len,clim.nc$dim$lat$len,1))
  lon.c <- ncvar_get(clim.nc, "lon")
  lat.c <- ncvar_get(clim.nc, "lat")
  nc_close(clim.nc)
  
  sst.c[sst.c < 0] = 0
  
  sst.anomnc = sst - sst.c
  chl.anomnc = chl - chl.c
  var.list <- list()
  var.list[[1]] <- ncdf4::ncvar_def(name="sst", units="Celsius", missval=-999, longname = "Sea Surface Temperature Anomaly", dim=data_dim)
  var.list[[2]] <- ncdf4::ncvar_def(name="chl_oc3", units="mg m^-3", missval=-999, longname = "Chlorophyll Concentration, OC3 Algorithm", dim=data_dim)
  loc.file <- paste0("/data/Aqua/anomaly/SSTanomaly_", i)
  #loc.file <- "Downloads/atest.nc"
  #writing all we need to the output file
  loc <- ncdf4::nc_create(filename=loc.file, vars=var.list, force_v4 = T)
  ncdf4::ncvar_put(nc=loc, "sst", vals=sst.anomnc)
  ncdf4::ncvar_put(nc=loc, "chl_oc3", vals=chl.anomnc)
  ncdf4::ncatt_put(nc=loc, 0, "Conventions", "CF=1.0")
  ncdf4::ncatt_put(nc=loc, 0,"creator_name", "James Simkins")
  ncdf4::ncatt_put(nc=loc, 0, "creator_email", "simkins@udel.edu")
  ncdf4::ncatt_put(nc=loc, 0, "institution", "University of Delaware Ocean Exploration, Remote Sensing and Biogeography Group (ORB)")
  ncdf4::ncatt_put(nc=loc, 0, "url", "http://orb.ceoe.udel.edu/")
  ncdf4::ncatt_put(nc=loc, 0, "source", "satellite observation NASA MODIS-Aqua instrument")
  ncdf4::ncatt_put(nc=loc, 0, "groundstation", "University of Delaware, Newark, Center for Remote Sensing")
  ncdf4::ncatt_put(nc=loc, 0, "software", "0.0")
  ncdf4::ncatt_put(nc=loc, 0, "inputMET1", "0.0")
  ncdf4::ncatt_put(nc=loc, 0, "inputOZONE1", "0.0")
  ncdf4::ncatt_put(nc=loc, 0, "inputCalibrationFile", "0.0")
  ncdf4::ncatt_put(nc=loc, 0, "product_list", "sst")
  ncdf4::ncatt_put(nc=loc, 0, "summary", "MODIS Aqua 8 day Aggregate Anomaly; Regridded to 
                   Mercator lon/lat projection. Processed at the Univeristy of Delaware")
  ncdf4::nc_close(loc)
}
