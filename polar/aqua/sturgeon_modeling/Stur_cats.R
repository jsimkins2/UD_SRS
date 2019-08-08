
#### This code is to categorize the GAM predictions for Atlantic Sturgeon occurrence using DINEOF outputs.

require(magick)
require(RColorBrewer)
require(ncdf4)
require(raster)
require(lubridate)
require(matlab)
require(fields)
working_path<-'/home/mshatley/scripts/sturgeon/'
output_dir <- '/data/Aqua/sturgeon/'

last <- function(x) { return( x[length(x)] ) }
breakpoints<- c(0,0.01, 0.05,1)
colors <- c("green4", "gold", "red")
prob <- 0.95

setwd(working_path)
load("v2_ts_dxdoy_temp_443.Rdata") ### base model as of 6/12/2017
load("DINeof_stuff.Rdata")
#setwd("/home/mwbreece/realtime/future")

#f<- nc_open("http://128.175.28.250:8080/thredds/dodsC/Aqua1DayFilled.nc")
f<- nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/gapfilled_1day_aqua.nc")


tm<- ncvar_get(f, "time")
lon <- ncvar_get(f, "lon") ### loading the longitude from the THREDDS server
#lon<-ncvar_get(f, 'longitude')
lat <- ncvar_get(f, "lat") ### loading the latitude from the THREDDS server
#lat<-ncvar_get(f, 'latitude')
fore_time<- ncvar_get(f, "forecast_time")
fore_sst<- ncvar_get(f,"forecast_sst", start = c(1, 1, 1,length(tm)), count = c(length(lon), length(lat),4,1))
fore_443 <- ncvar_get(f, "forecast_a_443_qaa", start = c(1, 1, 1,length(tm)), count = c(length(lon), length(lat),4,1))
data_status <- ncvar_get(f, 'data_status', start = c(length(tm)), count = c(1))
nc_close(f)

cur_tm <- as.Date(last(tm), origin = "1970-01-01")
fu1_tm <- as.Date(last(tm+1), origin = "1970-01-01")
fu2_tm <- as.Date(last(tm+2), origin = "1970-01-01")
fu3_tm <- as.Date(last(tm+3), origin = "1970-01-01")

###gapfilled #################
DofY <-yday(cur_tm)
yd<- matrix(DofY, length(lon), length(lat))
ydr<- flip(raster(t(yd)),2)
sstr<- flip(raster(t(fore_sst[,,1])),2)
a_443_qaar<- flip(raster(t(fore_443[,,1])),2)
env_stack<- stack(sstr, a_443_qaar, Depthr, ydr) #sturgeon
names(env_stack)<- c("sst", "a_443_qaa", "Depth", "yd") #sturgeon
#setwd("/home/mwbreece/GAMM/plot/img")

p <- predict(env_stack, GAM$gam, filename="file4.img", type="response", fun=predict, overwrite=T, na.rm=T)
gf_p<-fliplr(t(as.matrix(p)))  ##gap filled predictions
p[] <- replace(p[], p[]>0.1, 0.1)

#######categories
#setwd("/home/mwbreece/GAMM/BAM")
load("raster_masks.Rdata")
load("lon_lat_nefop.Rdata")
load("nefop_grid.Rdata")

ind1<-which(sub_lon == min(lon))
ind2<- which(sub_lon == max(lon))
ind3<- which(sub_lat == min(lat))
ind4<- which(sub_lat == max(lat))
nefop_grid[ind1:ind2, ind3:ind4]<- fliplr(t(as.matrix(p)))

p<- flip(raster(t(nefop_grid)),2)

p.lt5_up_mask <- raster::mask(p, d.lt5_up_mask)
p.5.10_up_mask <- raster::mask(p, d.5.10_up_mask)
p.10.15_up_mask <- raster::mask(p, d.10.15_up_mask)
p.gt15_up_mask <- raster::mask(p, d.gt15_up_mask)

p.lt5_mid_mask <- raster::mask(p, d.lt5_mid_mask)
p.5.10_mid_mask <- raster::mask(p, d.5.10_mid_mask)
p.10.15_mid_mask <- raster::mask(p, d.10.15_mid_mask)
p.gt15_mid_mask <- raster::mask(p, d.gt15_mid_mask)

p.lt5_low_mask <- raster::mask(p, d.lt5_low_mask)
p.5.10_low_mask <- raster::mask(p, d.5.10_low_mask)
p.10.15_low_mask <- raster::mask(p, d.10.15_low_mask)
p.gt15_low_mask <- raster::mask(p, d.gt15_low_mask)

p.lt5_river_mask <- raster::mask(p, d.lt5_river_mask)
p.5.10_river_mask <- raster::mask(p, d.5.10_river_mask)
p.10.15_river_mask <- raster::mask(p, d.10.15_river_mask)
p.gt15_river_mask <- raster::mask(p, d.gt15_river_mask)

p.lt5_ocean_mask <- raster::mask(p, d.lt5_ocean_mask)
p.5.10_ocean_mask <- raster::mask(p, d.5.10_ocean_mask)
p.10.15_ocean_mask <- raster::mask(p, d.10.15_ocean_mask)
p.gt15_ocean_mask <- raster::mask(p, d.gt15_ocean_mask)


p.cat_lt5_river<- replace(d.lt5_river_mask, d.lt5_river_mask <= 0, quantile(p.lt5_river_mask,  probs = prob))
p.cat_5.10_river<- replace(d.5.10_river_mask, d.5.10_river_mask <= 0, quantile(p.5.10_river_mask,  probs = prob))
p.cat_10.15_river<- replace(d.10.15_river_mask, d.10.15_river_mask <= 0, quantile(p.10.15_river_mask,  probs = prob))
p.cat_gt15_river<- replace(d.gt15_river_mask, d.gt15_river_mask <= 0, quantile(p.gt15_river_mask,  probs = prob))

p.cat_lt5_up<- replace(d.lt5_up_mask, d.lt5_up_mask <= 0, quantile(p.lt5_up_mask,  probs = prob))
p.cat_5.10_up<- replace(d.5.10_up_mask, d.5.10_up_mask <= 0, quantile(p.5.10_up_mask,  probs = prob))
p.cat_10.15_up<- replace(d.10.15_up_mask, d.10.15_up_mask <= 0, quantile(p.10.15_up_mask,  probs = prob))
p.cat_gt15_up<- replace(d.gt15_up_mask, d.gt15_up_mask <= 0, quantile(p.gt15_up_mask,  probs = prob))

p.cat_lt5_mid<- replace(d.lt5_mid_mask, d.lt5_mid_mask <= 0, quantile(p.lt5_mid_mask,  probs = prob))
p.cat_5.10_mid<- replace(d.5.10_mid_mask, d.5.10_mid_mask <= 0, quantile(p.5.10_mid_mask,  probs = prob))
p.cat_10.15_mid<- replace(d.10.15_mid_mask, d.10.15_mid_mask <= 0, quantile(p.10.15_mid_mask,  probs = prob))
p.cat_gt15_mid<- replace(d.gt15_mid_mask, d.gt15_mid_mask <= 0, quantile(p.gt15_mid_mask,  probs = prob))

p.cat_lt5_low<- replace(d.lt5_low_mask, d.lt5_low_mask <= 0, quantile(p.lt5_low_mask,  probs = prob))
p.cat_5.10_low<- replace(d.5.10_low_mask, d.5.10_low_mask <= 0, quantile(p.5.10_low_mask,  probs = prob))
p.cat_10.15_low<- replace(d.10.15_low_mask, d.10.15_low_mask <= 0, quantile(p.10.15_low_mask,  probs = prob))
p.cat_gt15_low<- replace(d.gt15_low_mask, d.gt15_low_mask <= 0, quantile(p.gt15_low_mask,  probs = prob))

p.cat_lt5_ocean<- replace(d.lt5_ocean_mask, d.lt5_ocean_mask <= 0, quantile(p.lt5_ocean_mask,  probs = prob))
p.cat_5.10_ocean<- replace(d.5.10_ocean_mask, d.5.10_ocean_mask <= 0, quantile(p.5.10_ocean_mask,  probs = prob))
p.cat_10.15_ocean<- replace(d.10.15_ocean_mask, d.10.15_ocean_mask <= 0, quantile(p.10.15_ocean_mask,  probs = prob))
p.cat_gt15_ocean<- replace(d.gt15_ocean_mask, d.gt15_ocean_mask <= 0, quantile(p.gt15_ocean_mask,  probs = prob))



test<- merge(p.cat_lt5_river,p.cat_5.10_river,p.cat_10.15_river, p.cat_gt15_river, p.cat_lt5_up, p.cat_5.10_up, p.cat_10.15_up, p.cat_gt15_up, p.cat_lt5_mid,  p.cat_5.10_mid, p.cat_10.15_mid, p.cat_gt15_mid, p.cat_lt5_low,p.cat_5.10_low, p.cat_10.15_low, p.cat_gt15_low, p.cat_lt5_ocean, p.cat_5.10_ocean, p.cat_10.15_ocean, p.cat_gt15_ocean)

gapfilled <- fliplr(t(as.matrix(test)))[ind1:ind2, ind3:ind4]


#######################future 1 day#################################

DofY <-yday(fu1_tm)
yd<- matrix(DofY, length(lon), length(lat))
ydr<- flip(raster(t(yd)),2)
sstr<- flip(raster(t(fore_sst[,,2])),2)
a_443_qaar<- flip(raster(t(fore_443[,,2])),2)
env_stack<- stack(sstr, a_443_qaar, Depthr, ydr) #sturgeon
names(env_stack)<- c("sst", "a_443_qaa", "Depth", "yd") #sturgeon
#setwd("/home/mwbreece/GAMM/plot/img")

p <- predict(env_stack, GAM$gam, filename="file4.img", type="response", fun=predict, overwrite=T, na.rm=T)
fu1_p<-fliplr(t(as.matrix(p))) ##predictions for 1 day in future
p[] <- replace(p[], p[]>0.1, 0.1)


#######categories 1day future
#setwd("/home/mwbreece/GAMM/BAM")
#already loaded don't need to do again...
#load("raster_masks.Rdata")
#load("lon_lat_nefop.Rdata")
#load("nefop_grid.Rdata")

ind1<-which(sub_lon == min(lon))
ind2<- which(sub_lon == max(lon))
ind3<- which(sub_lat == min(lat))
ind4<- which(sub_lat == max(lat))
nefop_grid[ind1:ind2, ind3:ind4]<- fliplr(t(as.matrix(p)))

p<- flip(raster(t(nefop_grid)),2)

p.lt5_up_mask <- raster::mask(p, d.lt5_up_mask)
p.5.10_up_mask <- raster::mask(p, d.5.10_up_mask)
p.10.15_up_mask <- raster::mask(p, d.10.15_up_mask)
p.gt15_up_mask <- raster::mask(p, d.gt15_up_mask)

p.lt5_mid_mask <- raster::mask(p, d.lt5_mid_mask)
p.5.10_mid_mask <- raster::mask(p, d.5.10_mid_mask)
p.10.15_mid_mask <- raster::mask(p, d.10.15_mid_mask)
p.gt15_mid_mask <- raster::mask(p, d.gt15_mid_mask)

p.lt5_low_mask <- raster::mask(p, d.lt5_low_mask)
p.5.10_low_mask <- raster::mask(p, d.5.10_low_mask)
p.10.15_low_mask <- raster::mask(p, d.10.15_low_mask)
p.gt15_low_mask <- raster::mask(p, d.gt15_low_mask)

p.lt5_river_mask <- raster::mask(p, d.lt5_river_mask)
p.5.10_river_mask <- raster::mask(p, d.5.10_river_mask)
p.10.15_river_mask <- raster::mask(p, d.10.15_river_mask)
p.gt15_river_mask <- raster::mask(p, d.gt15_river_mask)

p.lt5_ocean_mask <- raster::mask(p, d.lt5_ocean_mask)
p.5.10_ocean_mask <- raster::mask(p, d.5.10_ocean_mask)
p.10.15_ocean_mask <- raster::mask(p, d.10.15_ocean_mask)
p.gt15_ocean_mask <- raster::mask(p, d.gt15_ocean_mask)


p.cat_lt5_river<- replace(d.lt5_river_mask, d.lt5_river_mask <= 0, quantile(p.lt5_river_mask,  probs = prob))
p.cat_5.10_river<- replace(d.5.10_river_mask, d.5.10_river_mask <= 0, quantile(p.5.10_river_mask,  probs = prob))
p.cat_10.15_river<- replace(d.10.15_river_mask, d.10.15_river_mask <= 0, quantile(p.10.15_river_mask,  probs = prob))
p.cat_gt15_river<- replace(d.gt15_river_mask, d.gt15_river_mask <= 0, quantile(p.gt15_river_mask,  probs = prob))

p.cat_lt5_up<- replace(d.lt5_up_mask, d.lt5_up_mask <= 0, quantile(p.lt5_up_mask,  probs = prob))
p.cat_5.10_up<- replace(d.5.10_up_mask, d.5.10_up_mask <= 0, quantile(p.5.10_up_mask,  probs = prob))
p.cat_10.15_up<- replace(d.10.15_up_mask, d.10.15_up_mask <= 0, quantile(p.10.15_up_mask,  probs = prob))
p.cat_gt15_up<- replace(d.gt15_up_mask, d.gt15_up_mask <= 0, quantile(p.gt15_up_mask,  probs = prob))

p.cat_lt5_mid<- replace(d.lt5_mid_mask, d.lt5_mid_mask <= 0, quantile(p.lt5_mid_mask,  probs = prob))
p.cat_5.10_mid<- replace(d.5.10_mid_mask, d.5.10_mid_mask <= 0, quantile(p.5.10_mid_mask,  probs = prob))
p.cat_10.15_mid<- replace(d.10.15_mid_mask, d.10.15_mid_mask <= 0, quantile(p.10.15_mid_mask,  probs = prob))
p.cat_gt15_mid<- replace(d.gt15_mid_mask, d.gt15_mid_mask <= 0, quantile(p.gt15_mid_mask,  probs = prob))

p.cat_lt5_low<- replace(d.lt5_low_mask, d.lt5_low_mask <= 0, quantile(p.lt5_low_mask,  probs = prob))
p.cat_5.10_low<- replace(d.5.10_low_mask, d.5.10_low_mask <= 0, quantile(p.5.10_low_mask,  probs = prob))
p.cat_10.15_low<- replace(d.10.15_low_mask, d.10.15_low_mask <= 0, quantile(p.10.15_low_mask,  probs = prob))
p.cat_gt15_low<- replace(d.gt15_low_mask, d.gt15_low_mask <= 0, quantile(p.gt15_low_mask,  probs = prob))

p.cat_lt5_ocean<- replace(d.lt5_ocean_mask, d.lt5_ocean_mask <= 0, quantile(p.lt5_ocean_mask,  probs = prob))
p.cat_5.10_ocean<- replace(d.5.10_ocean_mask, d.5.10_ocean_mask <= 0, quantile(p.5.10_ocean_mask,  probs = prob))
p.cat_10.15_ocean<- replace(d.10.15_ocean_mask, d.10.15_ocean_mask <= 0, quantile(p.10.15_ocean_mask,  probs = prob))
p.cat_gt15_ocean<- replace(d.gt15_ocean_mask, d.gt15_ocean_mask <= 0, quantile(p.gt15_ocean_mask,  probs = prob))



test<- merge(p.cat_lt5_river,p.cat_5.10_river,p.cat_10.15_river, p.cat_gt15_river, p.cat_lt5_up, p.cat_5.10_up, p.cat_10.15_up, p.cat_gt15_up, p.cat_lt5_mid,  p.cat_5.10_mid, p.cat_10.15_mid, p.cat_gt15_mid, p.cat_lt5_low,p.cat_5.10_low, p.cat_10.15_low, p.cat_gt15_low, p.cat_lt5_ocean, p.cat_5.10_ocean, p.cat_10.15_ocean, p.cat_gt15_ocean)
fu1_day<-fliplr(t(as.matrix(test)))[ind1:ind2, ind3:ind4]

#####future 2 days

DofY <-yday(fu2_tm)
yd<- matrix(DofY, length(lon), length(lat))
ydr<- flip(raster(t(yd)),2)
sstr<- flip(raster(t(fore_sst[,,3])),2)
a_443_qaar<- flip(raster(t(fore_443[,,3])),2)
env_stack<- stack(sstr, a_443_qaar, Depthr, ydr) #sturgeon
names(env_stack)<- c("sst", "a_443_qaa", "Depth", "yd") #sturgeon
#setwd("/home/mwbreece/GAMM/plot/img")

p <- predict(env_stack, GAM$gam, filename="file4.img", type="response", fun=predict, overwrite=T, na.rm=T)
fu2_p<-fliplr(t(as.matrix(p))) ##predictions for 2 days in future
p[] <- replace(p[], p[]>0.1, 0.1)


#######categories 2day future
#setwd("/home/mwbreece/GAMM/BAM")
#load("raster_masks.Rdata")
#load("lon_lat_nefop.Rdata")
#load("nefop_grid.Rdata")

ind1<-which(sub_lon == min(lon))
ind2<- which(sub_lon == max(lon))
ind3<- which(sub_lat == min(lat))
ind4<- which(sub_lat == max(lat))
nefop_grid[ind1:ind2, ind3:ind4]<- fliplr(t(as.matrix(p)))

p<- flip(raster(t(nefop_grid)),2)

p.lt5_up_mask <- raster::mask(p, d.lt5_up_mask)
p.5.10_up_mask <- raster::mask(p, d.5.10_up_mask)
p.10.15_up_mask <- raster::mask(p, d.10.15_up_mask)
p.gt15_up_mask <- raster::mask(p, d.gt15_up_mask)

p.lt5_mid_mask <- raster::mask(p, d.lt5_mid_mask)
p.5.10_mid_mask <- raster::mask(p, d.5.10_mid_mask)
p.10.15_mid_mask <- raster::mask(p, d.10.15_mid_mask)
p.gt15_mid_mask <- raster::mask(p, d.gt15_mid_mask)

p.lt5_low_mask <- raster::mask(p, d.lt5_low_mask)
p.5.10_low_mask <- raster::mask(p, d.5.10_low_mask)
p.10.15_low_mask <- raster::mask(p, d.10.15_low_mask)
p.gt15_low_mask <- raster::mask(p, d.gt15_low_mask)

p.lt5_river_mask <- raster::mask(p, d.lt5_river_mask)
p.5.10_river_mask <- raster::mask(p, d.5.10_river_mask)
p.10.15_river_mask <- raster::mask(p, d.10.15_river_mask)
p.gt15_river_mask <- raster::mask(p, d.gt15_river_mask)

p.lt5_ocean_mask <- raster::mask(p, d.lt5_ocean_mask)
p.5.10_ocean_mask <- raster::mask(p, d.5.10_ocean_mask)
p.10.15_ocean_mask <- raster::mask(p, d.10.15_ocean_mask)
p.gt15_ocean_mask <- raster::mask(p, d.gt15_ocean_mask)


p.cat_lt5_river<- replace(d.lt5_river_mask, d.lt5_river_mask <= 0, quantile(p.lt5_river_mask,  probs = prob))
p.cat_5.10_river<- replace(d.5.10_river_mask, d.5.10_river_mask <= 0, quantile(p.5.10_river_mask,  probs = prob))
p.cat_10.15_river<- replace(d.10.15_river_mask, d.10.15_river_mask <= 0, quantile(p.10.15_river_mask,  probs = prob))
p.cat_gt15_river<- replace(d.gt15_river_mask, d.gt15_river_mask <= 0, quantile(p.gt15_river_mask,  probs = prob))

p.cat_lt5_up<- replace(d.lt5_up_mask, d.lt5_up_mask <= 0, quantile(p.lt5_up_mask,  probs = prob))
p.cat_5.10_up<- replace(d.5.10_up_mask, d.5.10_up_mask <= 0, quantile(p.5.10_up_mask,  probs = prob))
p.cat_10.15_up<- replace(d.10.15_up_mask, d.10.15_up_mask <= 0, quantile(p.10.15_up_mask,  probs = prob))
p.cat_gt15_up<- replace(d.gt15_up_mask, d.gt15_up_mask <= 0, quantile(p.gt15_up_mask,  probs = prob))

p.cat_lt5_mid<- replace(d.lt5_mid_mask, d.lt5_mid_mask <= 0, quantile(p.lt5_mid_mask,  probs = prob))
p.cat_5.10_mid<- replace(d.5.10_mid_mask, d.5.10_mid_mask <= 0, quantile(p.5.10_mid_mask,  probs = prob))
p.cat_10.15_mid<- replace(d.10.15_mid_mask, d.10.15_mid_mask <= 0, quantile(p.10.15_mid_mask,  probs = prob))
p.cat_gt15_mid<- replace(d.gt15_mid_mask, d.gt15_mid_mask <= 0, quantile(p.gt15_mid_mask,  probs = prob))

p.cat_lt5_low<- replace(d.lt5_low_mask, d.lt5_low_mask <= 0, quantile(p.lt5_low_mask,  probs = prob))
p.cat_5.10_low<- replace(d.5.10_low_mask, d.5.10_low_mask <= 0, quantile(p.5.10_low_mask,  probs = prob))
p.cat_10.15_low<- replace(d.10.15_low_mask, d.10.15_low_mask <= 0, quantile(p.10.15_low_mask,  probs = prob))
p.cat_gt15_low<- replace(d.gt15_low_mask, d.gt15_low_mask <= 0, quantile(p.gt15_low_mask,  probs = prob))

p.cat_lt5_ocean<- replace(d.lt5_ocean_mask, d.lt5_ocean_mask <= 0, quantile(p.lt5_ocean_mask,  probs = prob))
p.cat_5.10_ocean<- replace(d.5.10_ocean_mask, d.5.10_ocean_mask <= 0, quantile(p.5.10_ocean_mask,  probs = prob))
p.cat_10.15_ocean<- replace(d.10.15_ocean_mask, d.10.15_ocean_mask <= 0, quantile(p.10.15_ocean_mask,  probs = prob))
p.cat_gt15_ocean<- replace(d.gt15_ocean_mask, d.gt15_ocean_mask <= 0, quantile(p.gt15_ocean_mask,  probs = prob))



test<- merge(p.cat_lt5_river,p.cat_5.10_river,p.cat_10.15_river, p.cat_gt15_river, p.cat_lt5_up, p.cat_5.10_up, p.cat_10.15_up, p.cat_gt15_up, p.cat_lt5_mid,  p.cat_5.10_mid, p.cat_10.15_mid, p.cat_gt15_mid, p.cat_lt5_low,p.cat_5.10_low, p.cat_10.15_low, p.cat_gt15_low, p.cat_lt5_ocean, p.cat_5.10_ocean, p.cat_10.15_ocean, p.cat_gt15_ocean)

fu2_day<- fliplr(t(as.matrix(test)))[ind1:ind2, ind3:ind4]
################## future 3 day #################

DofY <-yday(fu3_tm)
yd<- matrix(DofY, length(lon), length(lat))
ydr<- flip(raster(t(yd)),2)
sstr<- flip(raster(t(fore_sst[,,4])),2)
a_443_qaar<- flip(raster(t(fore_443[,,4])),2)
env_stack<- stack(sstr, a_443_qaar, Depthr, ydr) #sturgeon
names(env_stack)<- c("sst", "a_443_qaa", "Depth", "yd") #sturgeon
#setwd("/home/mwbreece/GAMM/plot/img")

p <- predict(env_stack, GAM$gam, filename="file4.img", type="response", fun=predict, overwrite=T, na.rm=T)
fu3_p<-fliplr(t(as.matrix(p)))  ##predictions for 3 days in future

p[] <- replace(p[], p[]>0.1, 0.1)

#######categories 3day future
#setwd("/home/mwbreece/GAMM/BAM")
#load("raster_masks.Rdata")
#load("lon_lat_nefop.Rdata")
#load("nefop_grid.Rdata")

ind1<-which(sub_lon == min(lon))
ind2<- which(sub_lon == max(lon))
ind3<- which(sub_lat == min(lat))
ind4<- which(sub_lat == max(lat))
nefop_grid[ind1:ind2, ind3:ind4]<- fliplr(t(as.matrix(p)))

p<- flip(raster(t(nefop_grid)),2)

p.lt5_up_mask <- raster::mask(p, d.lt5_up_mask)
p.5.10_up_mask <- raster::mask(p, d.5.10_up_mask)
p.10.15_up_mask <- raster::mask(p, d.10.15_up_mask)
p.gt15_up_mask <- raster::mask(p, d.gt15_up_mask)

p.lt5_mid_mask <- raster::mask(p, d.lt5_mid_mask)
p.5.10_mid_mask <- raster::mask(p, d.5.10_mid_mask)
p.10.15_mid_mask <- raster::mask(p, d.10.15_mid_mask)
p.gt15_mid_mask <- raster::mask(p, d.gt15_mid_mask)

p.lt5_low_mask <- raster::mask(p, d.lt5_low_mask)
p.5.10_low_mask <- raster::mask(p, d.5.10_low_mask)
p.10.15_low_mask <- raster::mask(p, d.10.15_low_mask)
p.gt15_low_mask <- raster::mask(p, d.gt15_low_mask)

p.lt5_river_mask <- raster::mask(p, d.lt5_river_mask)
p.5.10_river_mask <- raster::mask(p, d.5.10_river_mask)
p.10.15_river_mask <- raster::mask(p, d.10.15_river_mask)
p.gt15_river_mask <- raster::mask(p, d.gt15_river_mask)

p.lt5_ocean_mask <- raster::mask(p, d.lt5_ocean_mask)
p.5.10_ocean_mask <- raster::mask(p, d.5.10_ocean_mask)
p.10.15_ocean_mask <- raster::mask(p, d.10.15_ocean_mask)
p.gt15_ocean_mask <- raster::mask(p, d.gt15_ocean_mask)


p.cat_lt5_river<- replace(d.lt5_river_mask, d.lt5_river_mask <= 0, quantile(p.lt5_river_mask,  probs = prob))
p.cat_5.10_river<- replace(d.5.10_river_mask, d.5.10_river_mask <= 0, quantile(p.5.10_river_mask,  probs = prob))
p.cat_10.15_river<- replace(d.10.15_river_mask, d.10.15_river_mask <= 0, quantile(p.10.15_river_mask,  probs = prob))
p.cat_gt15_river<- replace(d.gt15_river_mask, d.gt15_river_mask <= 0, quantile(p.gt15_river_mask,  probs = prob))

p.cat_lt5_up<- replace(d.lt5_up_mask, d.lt5_up_mask <= 0, quantile(p.lt5_up_mask,  probs = prob))
p.cat_5.10_up<- replace(d.5.10_up_mask, d.5.10_up_mask <= 0, quantile(p.5.10_up_mask,  probs = prob))
p.cat_10.15_up<- replace(d.10.15_up_mask, d.10.15_up_mask <= 0, quantile(p.10.15_up_mask,  probs = prob))
p.cat_gt15_up<- replace(d.gt15_up_mask, d.gt15_up_mask <= 0, quantile(p.gt15_up_mask,  probs = prob))

p.cat_lt5_mid<- replace(d.lt5_mid_mask, d.lt5_mid_mask <= 0, quantile(p.lt5_mid_mask,  probs = prob))
p.cat_5.10_mid<- replace(d.5.10_mid_mask, d.5.10_mid_mask <= 0, quantile(p.5.10_mid_mask,  probs = prob))
p.cat_10.15_mid<- replace(d.10.15_mid_mask, d.10.15_mid_mask <= 0, quantile(p.10.15_mid_mask,  probs = prob))
p.cat_gt15_mid<- replace(d.gt15_mid_mask, d.gt15_mid_mask <= 0, quantile(p.gt15_mid_mask,  probs = prob))

p.cat_lt5_low<- replace(d.lt5_low_mask, d.lt5_low_mask <= 0, quantile(p.lt5_low_mask,  probs = prob))
p.cat_5.10_low<- replace(d.5.10_low_mask, d.5.10_low_mask <= 0, quantile(p.5.10_low_mask,  probs = prob))
p.cat_10.15_low<- replace(d.10.15_low_mask, d.10.15_low_mask <= 0, quantile(p.10.15_low_mask,  probs = prob))
p.cat_gt15_low<- replace(d.gt15_low_mask, d.gt15_low_mask <= 0, quantile(p.gt15_low_mask,  probs = prob))

p.cat_lt5_ocean<- replace(d.lt5_ocean_mask, d.lt5_ocean_mask <= 0, quantile(p.lt5_ocean_mask,  probs = prob))
p.cat_5.10_ocean<- replace(d.5.10_ocean_mask, d.5.10_ocean_mask <= 0, quantile(p.5.10_ocean_mask,  probs = prob))
p.cat_10.15_ocean<- replace(d.10.15_ocean_mask, d.10.15_ocean_mask <= 0, quantile(p.10.15_ocean_mask,  probs = prob))
p.cat_gt15_ocean<- replace(d.gt15_ocean_mask, d.gt15_ocean_mask <= 0, quantile(p.gt15_ocean_mask,  probs = prob))



test<- merge(p.cat_lt5_river,p.cat_5.10_river,p.cat_10.15_river, p.cat_gt15_river, p.cat_lt5_up, p.cat_5.10_up, p.cat_10.15_up, p.cat_gt15_up, p.cat_lt5_mid,  p.cat_5.10_mid, p.cat_10.15_mid, p.cat_gt15_mid, p.cat_lt5_low,p.cat_5.10_low, p.cat_10.15_low, p.cat_gt15_low, p.cat_lt5_ocean, p.cat_5.10_ocean, p.cat_10.15_ocean, p.cat_gt15_ocean)

fu3_day<- fliplr(t(as.matrix(test)))[ind1:ind2, ind3:ind4]

#layer 1 is gapfilled
#layer 2 is fu1_day
#layer 3 is fu2_day
#layer 4 is fu3_day

#day is last(tm)
#longitude is sub_lon
#latitude is  sub_lat
########





## fixed to here so that lat and lon are the same as the original netcdf 




tm_dim <- ncdim_def("time", "Days since 1970-01-01 00:00:00", c(as.integer(last(tm))), unlim = TRUE)
mdltm_dim <- ncdim_def("forecast_time", "Days since 1970-01-01 00:00:00", c(last(tm), last(tm)+1, last(tm)+2, last(tm)+3))
#cat_lon_dim <- ncdim_def("categories_lon", "degrees_east", lon)
#cat_lat_dim <- ncdim_def("categories_lat", "degrees_north", lat)
lon_dim <- ncdim_def("lon", "degrees_east", lon)
lat_dim <- ncdim_def("lat", "degrees_north", lat)
var_defs = list()
var_defs[['sturgeon_predictions']] <- ncvar_def('Sturgeon_Predictions', 'Probability of Sturgeon Presence', list(lon_dim, lat_dim, mdltm_dim, tm_dim), -999, longname = 'Sturgeon Model Predictions', prec="single", compression=5)
var_defs[['sturgeon_categories']] <- ncvar_def('Sturgeon_Categories', 'Sturgeon Categories', list(lon_dim, lat_dim, mdltm_dim, tm_dim), -999, longname = 'Sturgeon Categories', prec="single", compression=5)
var_defs[['data_status']] <- ncvar_def(name='data_status', units='NA', list(tm_dim), -999, longname = 'Data Status Flag', prec="integer", compression=5)

print("start nc file")
newnc <- nc_create(sprintf('%s/model_predictions.%s.nc4',output_dir, cur_tm), var_defs, force_v4=TRUE)

ncvar_put(newnc, var_defs[['sturgeon_predictions']], replace(gf_p, is.na(gf_p), -999), start=c(1,1,1, 1), count=c(lon_dim$len, lat_dim$len, 1, 1))
ncvar_put(newnc, var_defs[['sturgeon_predictions']], replace(fu1_p, is.na(fu1_p), -999), start=c(1,1,2, 1), count=c(lon_dim$len, lat_dim$len, 1, 1))
ncvar_put(newnc, var_defs[['sturgeon_predictions']], replace(fu2_p, is.na(fu2_p), -999), start=c(1,1,3, 1), count=c(lon_dim$len, lat_dim$len, 1, 1))
ncvar_put(newnc, var_defs[['sturgeon_predictions']], replace(fu3_p, is.na(fu3_p), -999), start=c(1,1,4, 1), count=c(lon_dim$len, lat_dim$len, 1, 1))

ncvar_put(newnc, var_defs[['sturgeon_categories']], replace(gapfilled, is.na(gapfilled), -999), start=c(1,1,1, 1), count=c(lon_dim$len, lat_dim$len, 1, 1))
ncvar_put(newnc, var_defs[['sturgeon_categories']], replace(fu1_day, is.na(fu1_day), -999), start=c(1,1,2, 1), count=c(lon_dim$len,  lat_dim$len, 1, 1))
ncvar_put(newnc, var_defs[['sturgeon_categories']], replace(fu2_day, is.na(fu2_day), -999), start=c(1,1,3, 1), count=c(lon_dim$len,  lat_dim$len, 1, 1))
ncvar_put(newnc, var_defs[['sturgeon_categories']], replace(fu3_day, is.na(fu3_day), -999), start=c(1,1,4, 1), count=c(lon_dim$len,  lat_dim$len, 1, 1))

ncatt_put(newnc, 'data_status', 'Forecast_Status_Flag', '1 indicates the data may be degraded due to a lack of satellite coverage in the 3 days leading up to the forecast. Caution should be used')
ncvar_put(newnc, var_defs[['data_status']], data_status)

nc_close(newnc)


