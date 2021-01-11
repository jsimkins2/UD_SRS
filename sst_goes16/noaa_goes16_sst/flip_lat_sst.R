# flip the latitude for ERDDAP viewing 
library(ncdf4)
args = commandArgs(trailingOnly=TRUE)
fname = args[1]
yr = lubridate::year(Sys.Date())

sstfile = nc_open(fname, write=TRUE)
lat=ncvar_get(sstfile, "latitude")
lat=rev(lat)

sstvar = ncvar_get(sstfile, "SST")
sstflipped = t(sstvar)
sstflipped = apply(sstflipped, 2, rev)
sst = t(sstflipped)


dqfvar = ncvar_get(sstfile, "DQF")
dqfflipped = t(dqfvar)
dqfflipped = apply(dqfflipped, 2, rev)
dqf = t(dqfflipped)


ncvar_put(sstfile, varid = 'DQF', vals = dqf)
ncvar_put(sstfile, varid = 'SST', vals = sst)
ncvar_put(sstfile, varid = 'latitude', vals = lat)
nc_close(sstfile)
print(paste0("Done flipping ", fname))
