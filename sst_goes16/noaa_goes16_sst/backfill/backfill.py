import xarray as xr
from netCDF4 import Dataset
import metpy
from datetime import datetime
import cartopy.feature as cfeature
from siphon.catalog import TDSCatalog
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib import patheffects
import metpy
from metpy.plots import colortables
import xarray as xr
from xarray.backends import NetCDF4DataStore
from netCDF4 import Dataset, num2date
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.io.img_tiles as cimgt
from dateutil import tz
import time
from time import mktime
import os.path
import os
import numpy as np
from os import listdir
from os.path import isfile, join
import calendar

# ignore this, we aren't using this data we are just using the projection found here
nowdate = datetime.utcnow()
cat = TDSCatalog('https://thredds-test.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/SeaSurfaceTemperature/FullDisk/' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')
dataset_name = sorted(cat.datasets.keys())[-1]
dataset = cat.datasets[dataset_name]

ds = dataset.remote_access(service='OPENDAP')
orig_d = ds

def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

def JulianDate_to_MMDDYYY(y,jd):
    month = 1
    day = 0
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month = month + 1
    return month,jd,y

#datadir = "/data/GOES/GOES-R/backfill/"
datadir = "/data/GOES/GOES-R/backfill/2018/"
sst_data = sorted([f for f in listdir(datadir) if isfile(join(datadir, f))])


for fname in sst_data:
    ftime = fname.split('_')[3][1:-3]
    gdatetime=datetime.strptime(ftime, '%Y%j%H%M')
    dataset = fname
    # rename in case they added an M6 instead of M3
    if str(dataset).split('_')[1] != 'ABI-L2-SSTF-M3':
        dataset_name = str(dataset).split('_')[0] + '_ABI-L2-SSTF-M3_G16' + str(dataset).split('G16')[1]
    else:
        dataset_name = str(dataset)
    if os.path.isfile("/data/GOES/GOES-R/sst/" + str(gdatetime.year) + "/" + str(dataset_name)) == False:
        if os.path.isfile("/data/GOES/GOES-R/suspect/"+ str(dataset_name)) == False:
            print("at the top" + dataset_name)
            ds = xr.open_dataset(datadir + fname)
            dat = ds.metpy.parse_cf('SST')
            proj = dat.metpy.cartopy_crs
            dat_dqf = ds.metpy.parse_cf('DQF')
        
            dat = dat.where(dat > -1)
            dat.values[np.isnan(dat.values)] = -999
        
            f = Dataset("/home/sat_ops/goesR/data/noaa_sst/backfill/raw/" + str(nowdate.year) + "/" + str(fname),'w', format='NETCDF4') #'w' stands for write
            # dimensions
            f.createDimension('x', dat['x'].size)
            f.createDimension('y', dat['y'].size)
            f.createDimension('time', 1)
        
            # variables
        
            x = f.createVariable('x', 'f4', 'x')
            x.standard_name = dat['x'].standard_name
            x.units = dat['x'].units
        
            y = f.createVariable('y', 'f4', 'y')
            y.standard_name = dat['y'].standard_name
            y.units = dat['y'].units
        
            time = f.createVariable('time', 'f8', 'time')
            time.standard_name = dat['t'].standard_name
            time.long_name = dat['t'].long_name
        
            proj = f.createVariable('goes_imager_projection', 'f4')
            proj.grid_mapping_name = orig_d.variables['goes_imager_projection'].grid_mapping_name
            proj.perspective_point_height = orig_d.variables['goes_imager_projection'].perspective_point_height
            proj.semi_major_axis = orig_d.variables['goes_imager_projection'].semi_major_axis
            proj.semi_minor_axis = orig_d.variables['goes_imager_projection'].semi_minor_axis
            proj.inverse_flattening = orig_d.variables['goes_imager_projection'].inverse_flattening
            proj.latitude_of_projection_origin = orig_d.variables['goes_imager_projection'].latitude_of_projection_origin
            proj.longitude_of_projection_origin = orig_d.variables['goes_imager_projection'].longitude_of_projection_origin
            proj.sweep_angle_axis = orig_d.variables['goes_imager_projection'].sweep_angle_axis
        
            sst = f.createVariable('SST', 'f4', ('time', 'y', 'x'))
            sst.long_name = orig_d.variables['SST'].long_name
            sst.standard_name = orig_d.variables['SST'].standard_name
            sst.units = orig_d.variables['SST'].units
            sst.resolution = orig_d.variables['SST'].resolution
            sst.missing_value = orig_d.variables['SST']._FillValue
            sst.grid_mapping = 'goes_imager_projection'
        
            dqf = f.createVariable('DQF', 'f4', ('time', 'y', 'x'))
            dqf.long_name = orig_d.variables['DQF'].long_name
            dqf.standard_name = orig_d.variables['DQF'].standard_name
            dqf.units = orig_d.variables['DQF'].units
            dqf.flag_values = orig_d.variables['DQF'].flag_values
            dqf.flag_meanings = orig_d.variables['DQF'].flag_meanings
            dqf.grid_mapping = 'goes_imager_projection'
        
            # data
            x[:] = dat['x'].values
            y[:] = dat['y'].values
            sst[:]= dat.values
            dqf[:]=dat_dqf
            time[:] = dat['t'].values
            f.close()
            
            print("completed initial saving of raw file")
            reproject = "Rscript /home/sat_ops/goesR/data/noaa_sst/backfill/reproject_sst.R " + str(fname)
            os.system(reproject)
            
            print("ALL DONE REPROJECTING")
            nowdate = datetime.utcnow()
            
        
            d = Dataset("/home/sat_ops/goesR/data/noaa_sst/backfill/raw/" + str(nowdate.year) + "/" + str(dataset_name))
            ds = NetCDF4DataStore(d)
            ds = xr.open_dataset(ds)
            dat = ds.metpy.parse_cf('SST')
            proj = dat.metpy.cartopy_crs
            print("/home/sat_ops/goesR/data/noaa_sst/backfill/rast_sst/" + str(nowdate.year) + "/SST" + str(dataset))
            dat_dqf  = xr.open_dataset("/home/sat_ops/goesR/data/noaa_sst/backfill/rast_dqf/" + str(nowdate.year) + "/DQF" + str(dataset))
            new = xr.open_dataset("/home/sat_ops/goesR/data/noaa_sst/backfill/rast_sst/" + str(nowdate.year) + "/SST" + str(dataset))
    
            new = new.where((new['latitude']>16) & (new['latitude']<52) & (new['longitude']>-100) & (new['longitude']<-50), drop=True)
            dat_dqf = dat_dqf.where((dat_dqf['latitude']>16) & (dat_dqf['latitude']<52) & (dat_dqf['longitude']>-100) & (dat_dqf['longitude']<-50), drop=True)
            # now write it all to netcdf!
            print("WRITING FINAL FILE TO NETCDF")
            f = Dataset("/data/GOES/GOES-R/sst/" + str(gdatetime.year) + "/" + str(dataset_name),'w', format='NETCDF4') #'w' stands for write
            # dimensions
            f.createDimension('longitude', new['longitude'].size)
            f.createDimension('latitude', new['latitude'].size)
            f.createDimension('time', 1)
    
            # variables
            longitude = f.createVariable('longitude', 'f4', 'longitude')
            longitude.standard_name = new['longitude'].long_name
            longitude.units = new['longitude'].units
    
            latitude = f.createVariable('latitude', 'f4', 'latitude')
            latitude.standard_name = new['latitude'].long_name
            latitude.units = new['latitude'].units
    
            time = f.createVariable('time', 'f8', 'time')
            time.standard_name = "time"
            time.long_name = "EPOCH Time"
            time.units = "seconds since 1970-01-01T00:00:00Z"
    
            proj = f.createVariable('projection', 'f4')
            proj.proj4_string = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
            proj.epsg_code = 4326
    
            sst = f.createVariable('SST', 'f4', ('time', 'latitude', 'longitude'))
            sst.long_name = d.variables['SST'].long_name
            sst.standard_name = d.variables['SST'].standard_name
            sst.units = d.variables['SST'].units
            sst.resolution = d.variables['SST'].resolution
            sst.missing_value = -999.0
    
            dqf = f.createVariable('DQF', 'f4', ('time', 'latitude', 'longitude'))
            dqf.long_name = d.variables['DQF'].long_name
            dqf.standard_name = d.variables['DQF'].standard_name
            dqf.units = d.variables['DQF'].units
            dqf.flag_values = d.variables['DQF'].flag_values
            dqf.flag_meanings = d.variables['DQF'].flag_meanings
    
            latitude[:] = new['latitude'].values
            longitude[:] = new['longitude'].values
            sst[:]= new['sst'].values
            dqf[:]=dat_dqf['dqf'].values
            time[:] = (dat['time'].values.astype('uint64') / 1e9).astype('uint32')
    
            # metadata
            f.creator_name = "James Simkins"
            f.creator_email = "simkins@udel.edu"
            f.institution = "University of Delaware Ocean Exploration, Remote Sensing and Biogeography Group (ORB)"
            f.url = "http://orb.ceoe.udel.edu/"
            f.source = "satellite observation NASA MODIS-Aqua instrument"
            f.groundstation = "University of Delaware, Newark, Center for Remote Sensing"
            f.software = 0.0
            f.inputMET1 = 0.0
            f.inputOZONE1 = 0.0
            f.inputCalibrationFile = 0.0
            f.product_list = "SST, DQF"
            f.summary = "GOES16 SST product, reprojected to EPSG:4326."
    
            f.close()
            # add netcdf file compression.  -L 5 means level 5 compression, -O means overwrite
            fname2 = "/data/GOES/GOES-R/sst/" + str(gdatetime.year) + "/" + str(dataset_name)
            temp_str = "ncks " + fname2 + " " + fname2 + " -L 5 -O"
            os.system(temp_str)
            
            print("flipping lat")
            flip_lat = "Rscript /home/sat_ops/goesR/data/noaa_sst/backfill/flip_lat_sst.R " + fname2
            os.system(flip_lat)
            os_dat1  = "rm /home/sat_ops/goesR/data/noaa_sst/backfill/rast_dqf/" + str(nowdate.year) + "/DQF" + str(dataset)
            os.system(os_dat1)
            os_dat2 = "rm /home/sat_ops/goesR/data/noaa_sst/backfill/rast_sst/" + str(nowdate.year) + "/SST" + str(dataset)
            print(os_dat2)
            os.system(os_dat2)
            os_dat3 = "rm /home/sat_ops/goesR/data/noaa_sst/backfill/raw/" + str(nowdate.year) + "/" + str(fname)
            os.system(os_dat3)





