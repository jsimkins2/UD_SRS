# viirs comparison
require("mapdata")
require("ncdf4")
require("raster")
require("RColorBrewer")
require(colorRamps)

sst = raster("Documents/Delaware/processing/viirs/test/l2.ocbg.mosaic.nc", varname = "sst")
newext <- c(-90, -60, 28, 50) 
sst.c <- crop(sst, newext)

ocbg = raster("Downloads/V2013092172400.L2_SNPP_SST.nc", varname = "geophysical_data/sst")

newext <- c(-90, -60, 28, 50) 
ocbg.c <- crop(ocbg, newext)

sst.diff = sst.c - ocbg.c

clrs = rev(grDevices::rainbow(200, start = 0, end = .95, s = 1, v = 1))
clrs[1] = "#000000"
r.range <- c(minValue(sst.c), maxValue(sst.c))

png(filename = "~/Documents/Delaware/processing/viirs/test/ocbg.sst.greatherthan1.png",width = 8, height = 8, units = "in", res = 200)
plot(sst.c, col=clrs, legend=FALSE, axes=FALSE)
plot(sst.c, legend.only=TRUE, col=clrs,
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=seq(round(r.range[1],1), round(r.range[2],1), length = 20),
                    labels=round(seq(round(r.range[1],1), round(r.range[2],1), length = 20),0), 
                    cex.axis=0.6),
     legend.args=list(text='SST (Celsius)', side=4, font=2, line=2.5, cex=0.8))
map('worldHires', c('USA', 'Canada'), xlim=c(newext[[1]], newext[[2]]), ylim=c(newext[[3]], newext[[4]]), 
    col = "gray47", fill=T, add = T)
map('rivers', add = T, lwd = .4)
dev.off()