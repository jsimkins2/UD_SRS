library(ncdf4)
library(raster)
library(rasterVis)
library(lubridate)
library(maps)
library(mapdata)
library(maptools)
library(grid)
library(png)
library(magick)

copyright.draw <- function(label, image, x, y, size, ...) {
  lab <- textGrob(label = label,
                  x = unit(x, "npc"), y = unit(y, "npc"),
                  just = c("left", "centre"), gp = gpar(...))
  logo <- rasterGrob(image = image,
                     x = unit(x, "npc") + unit(1, "grobwidth", lab), y = unit(y, "npc"),
                     width = unit(size, "cm"), height = unit(size, "cm"),
                     just = c("left", "centre"), gp = gpar(...))
  grid.draw(lab)
  grid.draw(logo)
}
im1 = readPNG("/home/sat_ops/goes_r/nexrad/cema38.png")
im2 = readPNG("/home/sat_ops/goes_r/nexrad/udel38.png")


args     <- commandArgs(TRUE)
lastline <- args[1]

# this loop goes through the available files in the Aqua data folder and compares what hasn't been done 
# in the logfile.txt. We go through the last 15 files and the next 15 files because sometimes we miss a 
# day because processing errors take place
for (j in seq(-15,15,1)){
  llyr       <- as.numeric(substr(lastline, 7, 10))
  lljday     <- as.numeric(substr(lastline, 11, 13)) + j
  if (lljday < 1){
    lljday <- ifelse(leap_year(llyr), 366 + lljday, 365 + lljday)
    llyr   <- llyr - 1
  }
  # Figure out Day Of Interest
  doi <- as.Date(as.numeric(lljday) - 1, origin = paste0(llyr, "-01-01"))
  
  # We have to do this or else it'll go back through all the data
  if (doi < as.Date('2018-01-01')){
    llyr   <- 2018
    lljday <- 1
    doi    <- as.Date(as.numeric(lljday), origin = paste0(llyr, "-01-01"))
  }
  days.in.year = ifelse(leap_year(llyr), 366, 365)
  
  ejday = as.numeric(lljday)
  eyr = year(doi)
  emn = month(as.Date(ejday, origin = paste0(eyr, "-01-01")))
  edom = day(as.Date(ejday, origin = paste0(eyr, "-01-01"))) - 1

  jday = ifelse(ejday + 7 > days.in.year, ejday + 7 - days.in.year, ejday + 7)
  yr = ifelse(jday < 7, year(doi) + 1, year(doi))
  mn = month(as.Date(jday, origin = paste0(yr, "-01-01")))
  dom = day(as.Date(jday, origin = paste0(yr, "-01-01"))) - 1

  
  # add 0s where needed to match the filename string
  ejday = ifelse(ejday > 99, ejday, ifelse(ejday > 9, paste0("0", ejday), paste0("00", ejday)))
  edom = ifelse(edom > 9, edom, paste0("0", edom))
  emn = ifelse(emn > 9, emn, paste0("0", emn))
  jday = ifelse(jday < 10, paste0("00",jday), ifelse(jday < 100, paste0("0",jday), jday))
  dom = ifelse(dom > 9, dom, paste0("0",dom)) #dom = day of month
  mn = ifelse(mn > 9, mn, paste0("0",mn))
  
  
  #fname represents the day of interest
  fname = paste0("/aqua.",eyr, ejday, ".", emn, edom,
                 "-", yr, jday, ".", mn, dom, ".D.L3.modis.NAT.v09.1000m.nc4")
  
  # now here is the full file
  aqua.file = paste0("/data/Aqua/8_day/",eyr,fname)
  if (file.exists(aqua.file) == TRUE){
    # if the file hasn't been processed, place it in logfile so we don't run it again next hour
    if (any(grep(fname, readLines("/home/james/anomalies/logfile.txt"))) == FALSE){
      cat(paste0(fname,"\n"), file = "/home/james/anomalies/logfile.txt", append = TRUE)
      break #since the file exists and we haven't run it yet, we need to exit the loop
    }
  }
  
  # if no file exists, then we better stop the whole operation because it'll crash anyway
  if (j == 15 && file.exists(aqua.file) == FALSE){
    stop("No Files Exists")
  }
}

  
  # read in the aqua file
  aqua.nc <- nc_open(aqua.file)
  #aqua.file = "/Users/leathers/Downloads/aqua.2017304.1031-2017311.1107.D.L3.modis.NAT.v09.1000m.nc4"
  sst <- ncvar_get(aqua.nc, "sst")
  lon <- ncvar_get(aqua.nc, "lon")
  lat <- ncvar_get(aqua.nc, "lat")
  nc_close(aqua.nc)
  
  sst[sst < 0] = 0
  
  # create new vectors
  aqua.xy <- cbind(rep(lon, length(lat)), rep(lat, each=length(lon)))
  aqua.xyv <- cbind(aqua.xy, as.vector(sst))
  
  #rasterize
  aqua.raster <- raster(extent(range(lon), range(lat)), res=0.00980196)
  aqua.raster <- rasterize(aqua.xyv[, 1:2], aqua.raster, aqua.xyv[,3], fun=mean) 
  
  # read in the specific 8 day segment of the 8 day climatology
  clim.nc = nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/aqua_clim_rolling8.nc")
  rec_time  <- max(clim.nc$dim$time$vals)
  rec_len   <- clim.nc$dim$time$len

  #clim.nc = nc_open("/data/Aqua/climatology/8_day/8_day_climatology.nc4")
  #e.day <- ceiling(((as.numeric(jday) + as.numeric(ejday)) / 2)/8)
  #clim.nc = nc_open("/Users/leathers/Downloads/8_day_climatology.nc4")
  #e.day <- 39
  sst.c <- ncvar_get(clim.nc, "sst",start = c(1,1,as.numeric(jday)), count = c(clim.nc$dim$lon$len,clim.nc$dim$lat$len,1))
  lon.c <- ncvar_get(clim.nc, "lon")
  lat.c <- ncvar_get(clim.nc, "lat")
  nc_close(clim.nc)
  
  sst.c[sst.c < 0] = 0
  
  # vectorize the cliamtology data
  clim.xy <- cbind(rep(lon, length(lat)), rep(lat, each=length(lon)))
  clim.xyv <- cbind(clim.xy, as.vector(sst.c))
  
  # rasterize the climatology data 
  clim.raster <- raster(extent(range(lon), range(lat)), res=0.00980196)
  clim.raster <- rasterize(clim.xyv[, 1:2], clim.raster, clim.xyv[,3], fun=mean) 
  
  # compute the difference and plot
  diff = aqua.raster - clim.raster
  
  # write out diff raster to netcdf 
  outname = substr(fname, 2, 90)
  print(fname)
  writeRaster(diff, paste0("/home/james/anomalies/data/SSTanomaly_", outname),  overwrite=TRUE, format="CDF", varname="SST", varunit="Celsius", 
              longname="Sea Surface Temperature", xname="lon", yname="lat")
  ###### PLOTTING ###### 
  
  #Color palette creation
  color0=c(0,0,100)
  color1=c(5,5,130)
  color2=c(20,20,165)
  color3=c(30,42,195)
  color4=c(33,70,225)
  color5=c(37,110,249)
  color6=c(48,153,255)
  color7=c(75,200,255)
  color8=c(255,255,255) #white
  color9=c(255,210,30)
  color10=c(255,160,10)
  color11=c(250,105,4)
  color12=c(240,53,1)
  color13=c(210,16,0)
  color14=c(165,3,0)
  color15=c(135,0,0)
  color16=c(110,0,0)
  
  clr.list = list()
  for (i in 0:16){
    tem = eval(parse(text = paste0("color",i)))
    clr.list[i+1] = rgb(tem[1]/255,tem[2]/255,tem[3]/255)
  }
  
  # Getting all of the maps together
  ext <- as.vector(extent(diff))
  
  usa = map("worldHires","usa", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
  usa <- map2SpatialPolygons(usa, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  can = map("worldHires", "Canada", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(can$names, ":"), function(x) x[1])
  can <- map2SpatialPolygons(can, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  mex = map("worldHires", "Mexico", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(mex$names, ":"), function(x) x[1])
  mex <- map2SpatialPolygons(mex, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  cub = map("worldHires", "cuba", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(cub$names, ":"), function(x) x[1])
  cub <- map2SpatialPolygons(cub, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  pr  = map("worldHires", "puerto rico", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(pr$names, ":"), function(x) x[1])
  pr <- map2SpatialPolygons(pr, IDs=IDs,
                            proj4string=CRS(projection(diff)))
  
  jam = map("worldHires", "jamaica", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(jam$names, ":"), function(x) x[1])
  jam <- map2SpatialPolygons(jam, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  bah = map("worldHires", "bahamas", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(bah$names, ":"), function(x) x[1])
  bah <- map2SpatialPolygons(bah, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  hai = map("worldHires", "haiti", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(hai$names, ":"), function(x) x[1])
  hai <- map2SpatialPolygons(hai, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  cay = map("worldHires", "cayman islands", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(cay$names, ":"), function(x) x[1])
  cay <- map2SpatialPolygons(cay, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  dr  = map("worldHires", "Dominican Republic", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(dr$names, ":"), function(x) x[1])
  dr <- map2SpatialPolygons(dr, IDs=IDs,
                            proj4string=CRS(projection(diff)))
  
  vi  = map("worldHires", "virgin islands", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(vi$names, ":"), function(x) x[1])
  vi <- map2SpatialPolygons(vi, IDs=IDs,
                            proj4string=CRS(projection(diff)))
  
  tc  = map("worldHires", "turks and caicos", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(tc$names, ":"), function(x) x[1])
  tc <- map2SpatialPolygons(tc, IDs=IDs,
                            proj4string=CRS(projection(diff)))
  
  ang = map("worldHires", "anguilla", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(ang$names, ":"), function(x) x[1])
  ang <- map2SpatialPolygons(ang, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  bel = map("worldHires", "belize", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(bel$names, ":"), function(x) x[1])
  bel <- map2SpatialPolygons(bel, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  hon = map("worldHires", "honduras", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(hon$names, ":"), function(x) x[1])
  hon <- map2SpatialPolygons(hon, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  gua = map("worldHires", "Guatemala", fill=TRUE,xlim=ext[1:2], ylim=ext[3:4],plot=FALSE)
  IDs <- sapply(strsplit(gua$names, ":"), function(x) x[1])
  gua <- map2SpatialPolygons(gua, IDs=IDs,
                             proj4string=CRS(projection(diff)))
  
  
  diverge0 <- function(p) {
    # p: a trellis object resulting from rasterVis::levelplot
    # ramp: the name of an RColorBrewer palette (as character), a character 
    #       vector of colour names to interpolate, or a colorRampPalette.
    require(RColorBrewer)
    require(rasterVis)
    ramp <- colorRampPalette(clr.list)
    rng <- range(p$legend[[1]]$args$key$at)
    s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
    i <- findInterval(rng[which.min(abs(rng))], s)
    zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
    p$legend[[1]]$args$key$at <- s[zlim]
    p$par.settings$regions$col <- ramp(1000)[zlim[-length(zlim)]]
    p
  }
  
  # Begin the PNG
  # max pixels = range lat / .00980196 * range lon / .00980196
  png(filename = "/home/james/anomalies/AquaAtlantic_SST_anomaly.png",width = 8, height = 6, units = "in", res = 200)
    lp = levelplot(diff, margin = FALSE, maxpixels = 18075029, xlab = "Longitude", ylab = "Latitude", 
                   main=paste0("MODIS-Aqua 8-Day SST Anomalies ", emn, "/", edom, "/", eyr, "-", mn, "/", dom, "/", yr), cex=2) +
    #layer(sp.polygons(usa, fill='gray16')) + layer(sp.polygons(can, fill='gray16')) + layer(sp.polygons(mex, fill='gray16')) + 
    #layer(sp.polygons(cub, fill='gray16')) + layer(sp.polygons(pr, fill='gray16')) + layer(sp.polygons(jam, fill='gray16')) + 
    #layer(sp.polygons(bel, fill='gray16')) + layer(sp.polygons(bah, fill='gray16')) + layer(sp.polygons(hai, fill='gray16')) + 
    #layer(sp.polygons(cay, fill='gray16')) + layer(sp.polygons(dr, fill='gray16')) + layer(sp.polygons(vi, fill='gray16')) +
    #layer(sp.polygons(tc, fill='gray16'))  + layer(sp.polygons(ang, fill='gray16')) + layer(sp.polygons(gua, fill = 'gray16')) +
    layer(sp.polygons(hon, fill='gray16'))
    diverge0(lp)
    copyright.draw("", im1, .93, .96, 1, fontsize = 8)
    copyright.draw("", im2, .02, .96, 1, fontsize = 8) 
    trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
    grid::grid.text('DegC', y=unit(-.015, "npc"), 
                    x=unit(.26, "npc"), gp=gpar(fontsize=8))
    trellis.unfocus()
  garbage <- dev.off()

  
  
  
  ##### --------------- AQUA PASS AND CLIMO PLOTS ------------------ #####
  library(ggplot2)
  
  usa <- map_data("usa")
  canada <- map_data("worldHires", "canada")
  mexico <- map_data("worldHires", "mexico")# we already did this, but we can do it again
  cuba <- map_data("worldHires", "cuba", col = 'gray16')
  nicaragua <- map_data("worldHires", "nicaragua", col = 'gray16')
  pr <- map_data("worldHires", "puerto rico", col = 'gray16')
  jam <- map_data("worldHires", "jamaica", col = 'gray16')
  cr <- map_data("worldHires", "Costa Rica", col = 'gray16')
  bah <- map_data("worldHires", "bahamas", col = 'gray16')
  hai <- map_data("worldHires", "haiti", col = 'gray16')
  cay <- map_data("worldHires", "cayman islands", col = 'gray16')
  dr <- map_data("worldHires", "Dominican Republic", col = 'gray16')
  vi <- map_data("worldHires", "virgin islands", col = 'gray16')
  tc <- map_data("worldHires", "turks and caicos", col = 'gray16')
  ang <- map_data("worldHires", "anguilla", col = 'gray16')
  bel <- map_data("worldHires", "belize", col = 'gray16')
  els <- map_data("worldHires", "el salvador", col = 'gray16')
  hon <- map_data("worldHires", "honduras", col = 'gray16')
  pan <- map_data("worldHires", "panama", col = 'gray16')
  guat <- map_data("worldHires", "Guatemala", col = 'gray16')
  
  color0=c(0,0,100)
  color1=c(5,5,130)
  color2=c(20,20,165)
  color3=c(30,42,195)
  color4=c(33,70,225)
  color5=c(37,110,249)
  color6=c(48,153,255)
  color7=c(75,200,255)
  color8=c(140,235,255)
  color9=c(200,250,255)
  color10=c(255,250,170)
  color11=c(255,237,80)
  color12=c(255,210,30)
  color13=c(255,160,10)
  color14=c(250,105,4)
  color15=c(240,53,1)
  color16=c(210,16,0)
  color17=c(165,3,0)
  color18=c(135,0,0)
  color19=c(110,0,0)
  
  clr.list = list()
  for (i in 0:19){
    tem = eval(parse(text = paste0("color",i)))
    clr.list[i+1] = rgb(tem[1]/255,tem[2]/255,tem[3]/255)
  }
  clr.ramp = colorRampPalette(clr.list)
  
  # Aqua Pass
  map.p = rasterToPoints(aqua.raster)
  df <- data.frame(map.p)
  #Make appropriate column headings
  colnames(df) <- c("Longitude", "Latitude", "MAP")

  png(filename = "/home/james/anomalies/AquaAtlantic_SST.png",width = 8, height = 6, units = "in", res = 200)
  gg = ggplot(data=df, aes(y=Latitude, x=Longitude)) +
    geom_raster(aes(fill=MAP)) +
    ggtitle(paste0("MODIS-Aqua 8-Day SST ", emn, "/", edom, "/", eyr, "-", mn, "/", dom, "/", yr)) +
    #geom_polygon(data = land, aes(x=long, y = lat, group = group), 
    #fill = "gray16") +
    
    theme_bw() +
    coord_equal() +
    scale_fill_gradientn("SST (DegC)", limits=c(0,40), na.value = "white", colours = clr.ramp(20)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16, angle=90),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          legend.key = element_blank(),
          legend.title = element_text(size=10),
          legend.key.height = unit(1, "in"),
          legend.key.width = unit(.5, "in")
    ) + geom_polygon(data = usa, aes(x=long, y = lat, group = group), 
                     fill = "gray16", 
                     color="black") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") + 
    geom_polygon(data = mexico, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = cuba, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = nicaragua, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = pr, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = jam, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = cr, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = bah, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = hai, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = cay, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = dr, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = vi, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = tc, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = ang, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = bel, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = els, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = hon, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = pan, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = guat, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    coord_fixed(xlim = c(-100, -55),  ylim = c(18, 48), ratio = 1.2)
  
  garbage <- dev.off()
  
  cema <- image_read("/home/sat_ops/goes_r/nexrad/cema38.png")
  udel <- image_read("/home/sat_ops/goes_r/nexrad/udel38.png")
  # add logos to ggplot
  plot <- image_read("/home/james/anomalies/AquaAtlantic_SST.png")
  final_plot  <- image_composite(plot, cema, offset = "+1540+5")
  final_plot  <- image_composite(final_plot, udel, offset = "+5+5")
  image_write(final_plot, "/home/james/anomalies/AquaAtlantic_SST.png")
  
  
  # CLIM
  map.p = rasterToPoints(clim.raster)
  df <- data.frame(map.p)
  #Make appropriate column headings
  colnames(df) <- c("Longitude", "Latitude", "MAP")
  
  png(filename = "/home/james/anomalies/AquaAtlantic_SST_climatology.png",width = 8, height = 6, units = "in", res = 200)
  ggplot(data=df, aes(y=Latitude, x=Longitude)) +
    geom_raster(aes(fill=MAP)) +
    ggtitle(paste0("MODIS-Aqua 8-Day SST Climatology ", emn, "/", edom, "/", eyr, "-", mn, "/", dom, "/", yr)) +
    #geom_polygon(data = land, aes(x=long, y = lat, group = group), 
    #fill = "gray16") +
    
    theme_bw() +
    coord_equal() +
    scale_fill_gradientn("SST (DegC)", limits=c(0,40), na.value = "white", colours = clr.ramp(20)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16, angle=90),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          legend.key = element_blank(),
          legend.title = element_text(size=10),
          legend.key.height = unit(1, "in"),
          legend.key.width = unit(.5, "in")
    ) + geom_polygon(data = usa, aes(x=long, y = lat, group = group), 
                     fill = "gray16", 
                     color="black") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") + 
    geom_polygon(data = mexico, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = cuba, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = nicaragua, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = pr, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = jam, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = cr, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = bah, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = hai, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = cay, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = dr, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = vi, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = tc, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = ang, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = bel, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = els, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = hon, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = pan, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    geom_polygon(data = guat, aes(x=long, y = lat, group = group), 
                 fill = "gray16", color="gray16") +
    coord_fixed(xlim = c(-100, -55),  ylim = c(18, 48), ratio = 1.2)
  
  garbage <- dev.off()
  
  # add logos to ggplot
  plot <- image_read("/home/james/anomalies/AquaAtlantic_SST_climatology.png")
  final_plot  <- image_composite(plot, cema, offset = "+1540+5")
  final_plot  <- image_composite(final_plot, udel, offset = "+5+5")
  image_write(final_plot, "/home/james/anomalies/AquaAtlantic_SST_climatology.png")
  
  diff.hist = hist(diff, plot = FALSE)
  
  png(filename = "/home/james/anomalies/AquaAtlantic_SST_histogram.png",width = 8, height = 6, units = "in", res = 200)
  plot(diff.hist, main = paste0(emn, "/", edom, "/", eyr, "-", mn, "/", dom, "/", yr, " MODIS-Aqua Anomaly Histogram"), xlab = "Sea Surface Temperature (DegC)")
  copyright.draw("", im1, .93, .96, 1, fontsize = 8)
  copyright.draw("", im2, .02, .96, 1, fontsize = 8) 
  garbage <- dev.off()
  
  cat(1)
