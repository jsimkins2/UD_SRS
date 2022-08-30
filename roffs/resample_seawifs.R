# running on basin
library(raster)
library(ncdf4)
library(stringr)

# generate file list
sealist <- Sys.glob("/data/seawifs/daily_9km/*.nc4")
# set the paths
temPath = "/home/james/temp/"
outPathNC = "/data/seawifs/daily_50km/"

# Set variable & units dataframe
var.df = data.frame(names = c('chl_ocx','chlor_a','a_412_giop','a_443_giop','a_510_giop','a_555_giop',
                              'adg_443_giop','adg_s_giop','aph_443_giop','aph_unc_443_giop','bb_443_giop',
                              'bb_510_giop','bb_670_giop','kd_490','par','poc','rrs_412','rrs_443','rrs_490',
                              'rrs_510','rrs_555','rrs_670','Zeu_lee'), 
                    units=c('mg m^-3','mg m^-3','m^-1','m^-1','m^-1','m^-1','m^-1','m^-1 nm^-1','m^-1','m^-1',
                             'm^-1','m^-1','m^-1','m^-1','einstein m^-2 day^-1','mg m^-3','sr^-1','sr^-1','sr^-1',
                             'sr^-1','sr^-1','sr^-1','m'), 
                    longname=c('chl_ocx','chlor_a','a_412_giop','a_443_giop','a_510_giop','a_555_giop',
                               'adg_443_giop','adg_s_giop','aph_443_giop','aph_unc_443_giop','bb_443_giop',
                               'bb_510_giop','bb_670_giop','kd_490','par','poc','rrs_412','rrs_443','rrs_490',
                               'rrs_510','rrs_555','rrs_670','Zeu_lee'))

####################################################################
############ Raster Operations #################
####################################################################
for (f in sealist){
  # name our nc file 
  ncfname = paste0(outPathNC, str_split(f, "/")[[1]][5])
  if (file.exists(ncfname) == FALSE){
    seaTem = nc_open(f)
    tval = seaTem$dim$time$vals
    tunits = seaTem$dim$time$units
    nc_close(seaTem)
    print(tval)
    for (pyVar in var.df$names){
    # open the most recent ET file 
      seaRast = raster(f, varname=pyVar, band=1)
      extent(seaRast) = c(-180, 180, -89.96, 89.96)
      # 9km / 50km = 0.18
      # 2160 lats * 0.18 = ~ 389, 4320 lons * 0.18 = ~ 778
      newRast = raster(nrow = 389, ncol = 778)
      extent(newRast) = c(-180, 180, -89.96, 89.96)
      resRast = resample(x=seaRast, y=newRast, method='bilinear') # can be set to nearest neighbor using 'ngb' method
      
      # output the raster as netCDF
      writeRaster(resRast, paste0(temPath, pyVar, "_temp.nc"),format="CDF", overwrite=TRUE,
                  varname=pyVar,
                  xname="longitude", yname="latitude")
    }
    # endT is 118 seconds * 4509 files = 532062 is 147 hours
    ####################################################################
    ############ Make the netCDF file CF compliant #################
    ####################################################################
    # create nc file with definitions for each variable
    def.list = list()
    for (v in seq_along(var.df$names)){
      tempNC = nc_open(paste0(temPath, as.character(var.df$names[v]), "_temp.nc"))
      tempET = ncvar_get(tempNC, as.character(var.df$names[v]))
      dims = tempNC$dim
      londim <- ncdim_def("longitude","degrees_east",as.double(dims$longitude$vals)) 
      latdim <- ncdim_def("latitude","degrees_north",as.double(dims$latitude$vals))
      timedim <- ncdim_def("time",tunits,tval)
      fillvalue <- -9999
      dlname <- as.character(var.df$longname[which(var.df$names == as.character(var.df$names[v]))])
      def.list[[v]] <- ncvar_def(as.character(var.df$names[v]),as.character(var.df$units[v]),
                                 list(londim,latdim,timedim),fillvalue,dlname,prec="single")
      nc_close(tempNC)
    }
    
    # create the netcdf file we're writing to
    ncout <- ncdf4::nc_create(filename = ncfname, vars = def.list)
    
    for (v in seq_along(var.df$names)){
      # reopen the nc files we just created
      tempNC = nc_open(paste0(temPath, as.character(var.df$names[v]), "_temp.nc"))
      tempET = ncvar_get(tempNC, as.character(var.df$names[v]))
      # put variables
      ncvar_put(ncout,def.list[[v]],tempET)
    }
    
    # add global attributes
    ncatt_put(ncout,0,"title",paste0("Seawifs"," - coarsened to spatial resolution of 50km using bilinear interpolation"))
    ncatt_put(ncout,0,"institution","Center for Environmental Monitoring and Analysis - University of Delaware")
    ncatt_put(ncout,0,"source","Seawifs 9km - repository found on Basin at /data/seawifs/daily_9km/")
    ncatt_put(ncout,0,"author","James Simkins")
    ncatt_put(ncout,0,attname="proj4", attval=proj4string(resRast))
    ncatt_put(ncout,0,attname="crs", attval="EPSG:4326")
    ncatt_put(ncout,0,"Conventions","CF-1.4")
    # close the file, writing data to disk
    nc_close(ncout)
  }
}
