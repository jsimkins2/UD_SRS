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
from os import listdir
from os.path import isfile, join

nowdate = datetime.utcnow()

datadir = "/home/sat_ops/goesR/data/sst/raw/" + str(nowdate.year) + "/"
dataset = [f for f in listdir(datadir) if isfile(join(datadir, f))][-1]

if os.path.isfile("/home/sat_ops/goesR/data/sst/reprojected/" + str(nowdate.year) + "/" + str(dataset)) == False:
    print('saving reprojected file')
    d = Dataset(datadir + dataset)
    ds = NetCDF4DataStore(d)
    ds = xr.open_dataset(ds)
    dat = ds.metpy.parse_cf('SST')
    proj = dat.metpy.cartopy_crs
    dat_dqf  = xr.open_dataset("/home/sat_ops/goesR/data/sst/dqf.nc")
    new = xr.open_dataset("/home/sat_ops/goesR/data/sst/sst.nc")
    
    
    
    f = Dataset("/home/sat_ops/goesR/data/sst/reprojected/" + str(nowdate.year) + "/" + str(dataset),'w', format='NETCDF4') #'w' stands for write
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
    
    # data 
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
else:
    print('all caught up!')
