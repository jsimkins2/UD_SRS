# Delaware Bay Cropping, Gapfilling, and Forecasting Script
library(ncdf4)
library(raster)
verbose = FALSE

###### ********* NEED TO SPECIFY PROPER DIRECTORY FOR RIVER POINTS DATA SET ******** #########
rvr.pts = get(load("/home/james/dineof2/riverpoints.Rdata"))

# some initial parameters
product<-c('sst', 'a_443_qaa', 'Rrs_555', 'chl_oc3')
#startfile='Downloads/aqua.2017342.1208.235959.D.L3.modis.NAT.v09.1000m.nc4'
end=180
forecast=3
work_dir='/data/Aqua/1_day/'


# Read Arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0){
  stop("NEED FILENAME TO PROCESS")
}

# Turn arguments into objects
for(i in 1:length(args)){
  cmd <- unlist(strsplit(args[i], '='))
  if(cmd[1] == "file"){
    print(c('FOUND FILE', cmd[2]))
    startfile <- cmd[2]
  }else if(cmd[1] == 'forecast'){
    forecast <- cmd[2]
    print(c("FOUND FORECAST: ", forecast))
  }else if(cmd[1] == 'days'){
    end <- cmd[2]
    print(c("FOUND DAYS: ", end))
  }
}

# give the user some info as to what is about to go down
print(sprintf("PROCESSING WITH: FILE: %s NUMBER OF DAY: %s FORECAST DAYS: %s", startfile, end, forecast))
end<-as.integer(end)
forecast<-as.integer(forecast)

# split up the start file
file_parts <- unlist(strsplit(startfile, '\\.'))
file_parts <- file_parts[2]
file_year  <- as.integer(paste(unlist(strsplit(file_parts, ''))[1:4], collapse=''))

# list all the files
filelist<-list.files(path=sprintf('%s/%s/', work_dir, file_year), pattern = 'aqua.[0-9].*nc4')
if(file_year > 2002)
  filelist<-c(list.files(path=sprintf('%s/%s/', work_dir, file_year-1), pattern = 'aqua.[0-9].*nc4'), filelist)
if(end > 180 && file_year > 2003)
  filelist<-c(list.files(path=sprintf('%s/%s/', work_dir, file_year-2), pattern = 'aqua.[0-9].*nc4'), filelist)

#grab last 180 files
startIx = which(filelist==startfile)
print(c("STARTFILE IX" , startIx))
#filelist<-filelist[startIx:(startIx+end-1)]
print(c("START", (startIx-end), " END:" , startIx))
filelist<-filelist[(startIx-end+1):startIx]
print(c("FOUND", length(filelist), " FILES TO ADD"))
file_year<-as.integer(paste(unlist(strsplit((unlist(strsplit(filelist[1], '\\.')))[2], ''))[1:4], collapse=''))

# open the current file
#f = nc_open(startfile)
f<-nc_open(sprintf('%s/%s/%s', work_dir, file_year, filelist[1]))
dat.list <- list()
for (i in seq_along(product)){
  dat.list[[i]] <- ncvar_get(f, product[i])
}
lat.a<-ncvar_get(f, 'lat')
lon.a<-ncvar_get(f, 'lon')
nc_close(f)

crop_lon.a <- subset(lon.a, lon.a<=-74.7 & lon.a>=-75.7)
#crop_lon <- lon
crop_lat.a <- subset(lat.a, lat.a<=39.7 & lat.a>=38.45)
#crop_lat <- lat
rows_lon.a <- which(lon.a<=-74.7 & lon.a>=-75.7)
rows_lat.a <- which(lat.a<=39.7 & lat.a>=38.45)
lon_dim <- ncdim_def("lon", "degrees_east", crop_lon.a)
lat_dim <- ncdim_def("lat", "degrees_north", crop_lat.a)

time_data = c();
for(file in 1:length(filelist)){
  day <- unlist(strsplit(filelist[file], '\\.'))[2]
  time_data <- c(time_data, as.integer(as.POSIXlt(day, format='%Y%j', tz='UTC')-as.POSIXlt('1970-01-01', tz='UTC')))
}
if(forecast){
  for(file in 1:forecast){
    time_data <- c(time_data, time_data[length(time_data)] + 1)
    #time_data <- c((time_data[1] + 1), time_data)
  }
}
tm_dim <- ncdim_def("time", "days in between files", time_data, unlim = TRUE)

var_defs<-list()
var_defs[[product[1]]] <- ncvar_def(product[1], 'degrees_C', list(lon_dim, lat_dim, tm_dim), -999, longname = 'sea_surface_temperature', prec="single")
var_defs[[product[2]]] <- ncvar_def(product[2], 'm^-1', list(lon_dim, lat_dim, tm_dim), -999, longname = 'absorption_due_to_phytoplankton', prec="single")
var_defs[[product[3]]] <- ncvar_def(product[3], 'sr^-1', list(lon_dim, lat_dim, tm_dim), -999, longname = 'Remote sensing reflectance at 555 nm', prec="single")
var_defs[[product[4]]] <- ncvar_def(product[4], 'mg m^-3', list(lon_dim, lat_dim, tm_dim), -999, longname = 'Chlorophyll Concentration, OC3 Algorithm', prec="single")

newnc <- nc_create(sprintf('%s.%s.combined', startfile, end) , var_defs, force_v4=FALSE)

for(vars in 1:length(product)){
  ncatt_put(newnc, var_defs[[product[vars]]], 'missing_value', -999)
}

for(file in 1:length(filelist)){
  print(file)
  file_year<-as.integer(paste(unlist(strsplit((unlist(strsplit(filelist[file], '\\.')))[2], ''))[1:4], collapse=''))
  filedate<-as.integer(substr(filelist[file], 14, 17))
  f<-nc_open(sprintf('%s/%s/%s', work_dir, file_year, filelist[file]))
  lat.a<-ncvar_get(f, 'lat')
  lon.a<-ncvar_get(f, 'lon')

  print(sprintf('%s/%s/%s', work_dir, file_year, filelist[file]))
  crop_lon.a <- subset(lon.a, lon.a<=-74.7 & lon.a>=-75.7)
  #crop_lon <- lon
  crop_lat.a <- subset(lat.a, lat.a<=39.7 & lat.a>=38.45)
  #crop_lat <- lat
  rows_lon.a <- which(lon.a<=-74.7 & lon.a>=-75.7)
  rows_lat.a <- which(lat.a<=39.7 & lat.a>=38.45)
  lon_dim <- ncdim_def("lon", "degrees_east", crop_lon.a)
  lat_dim <- ncdim_def("lat", "degrees_north", crop_lat.a)

  
  for(vars in 1:length(product)){
    param <- ncvar_get(f, product[vars])
    
    #open corresponding mursst file
    #g<-nc_open('/data/NOAA/sst/2015/20151231090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc')
    #mur.file = paste0("/data/NOAA/sst/", file_year, "/", file_year, filedate, "090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc")
    mur.file = sprintf("/data/NOAA/sst/%s/%s%04d%s", file_year, file_year, filedate, "090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc" )
    
    print(mur.file)
    # we might be missing a day here or there so now we are going to go back a day and see if it'll work then 
    if (file.exists(mur.file) == TRUE){
      print(c("MUR FILE: ", mur.file)) 
      g <- nc_open(mur.file) #get sst product
      mursst<-ncvar_get(g, 'analysed_sst')
      lat.m = ncvar_get(g, "lat")
      lon.m = ncvar_get(g, "lon")
      rows_lon.m <- which(lon.m<=-74.7 & lon.m>=-75.7)
      rows_lat.m <- which(lat.m<=39.7 & lat.m>=38.45)  
    } else {
      print("MUR FILE DOES NOT EXIST!!!!!!!!!")
      mursst = param
      lon.m = lon.a
      lat.m = lat.a
    }
    if(product[vars] == 'sst'){
      # this step gets rid of all the bad data, basically under the impression that Aqua data is wrong in the river
      val_rows_lat = lat.a[rows_lat.a]
      for (lt in seq_along(val_rows_lat)){
        for (ln in rows_lon.a){
          if (val_rows_lat[lt] > 39.35){
            param[ln,lt] = NA
            print(param[ln,lt])
          }
        }
      }
      
      # find random river points
      rand.rvr.pts = sample(1:(length(rvr.pts)/2), 5) #divide by 2 because there 2 columns in the r dataset
      rand.lon = list()
      rand.lat = list()
      for (r in seq_along(rand.rvr.pts)){
        rand.lon[[r]] = as.numeric(rvr.pts[rand.rvr.pts[r],1])
        rand.lat[[r]] = as.numeric(rvr.pts[rand.rvr.pts[r],2])
      }
      # This for loop finds the nearest latitude and longitude of NA values in Mur SST and places them in the Aqua SST 
      # Also note, the matrices are lon x lat
      for (i in seq_along(rand.rvr.pts)){
        x = rand.lon[[i]]
        y = rand.lat[[i]]
        # convert the random river points to their respective value in the aqua index
        rvr.aqua.lon = rows_lon.a[x]
        rvr.aqua.lat = rows_lat.a[y]
        
        if (is.na(param[rvr.aqua.lon,rvr.aqua.lat]) == TRUE){
          # grab the actual lat/lon value from the specified index so we can figure out where the closest mur value is
          tem.lon = lon.a[rvr.aqua.lon]
          tem.lat = lat.a[rvr.aqua.lat]

          # this is where we set the southern extent of kept data
          if (tem.lat > 39.35){
            param
            # find the closest replacement lon/lat and its index in the mur SST data 
            replace.lon = which(abs(lon.m - tem.lon)==min(abs(lon.m - tem.lon)))
            replace.lat = which(abs(lat.m - tem.lat)==min(abs(lat.m - tem.lat)))
            replace.sst = mursst[replace.lon, replace.lat]
            print(replace.sst)
            param[rvr.aqua.lon, rvr.aqua.lat] = replace.sst
          }
        }
      }
      # bind vectors to make matrix then raster
      xy <- cbind(rep(lon.a, length(lat.a)), rep(lat.a, each=length(lon.a)))
      xyv <- cbind(xy, as.vector(param))
      param = matrix(data = xyv[,3], nrow = length(lon.a), ncol = length(lat.a))
    } ####------------------ end the gapfilling
    if(product[vars] == 'a_443_qaa' || product[vars] == 'chl_oc3'){
      #param <- log10(param)
    }
    param<-replace(param, is.na(param), -999)
    ncvar_put(newnc, var_defs[[product[vars]]], param[rows_lon.a, rows_lat.a], start=c(1,1,file), count=c(length(rows_lon.a), length(rows_lat.a), 1))
  }
  nc_close(f)
}
if(forecast){
  print('ADDING FORECAST DAYS')
  param <- array(-999, dim=c(length(rows_lon.a), length(rows_lat.a)))
  for(file in 1:forecast){
    startIx = c(1,1,(file+length(filelist)))
    for(vars in 1:length(product)){
      print(product[vars])
      ncvar_put(newnc, var_defs[[product[vars]]], param, start=startIx, count=c(length(rows_lon.a), length(rows_lat.a), 1))
    }
  }
}
nc_close(newnc)
print("ALL DONE")
