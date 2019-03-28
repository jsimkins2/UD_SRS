import xarray as xr
from netCDF4 import Dataset
import numpy as np
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
import os
import sys

nowdate = datetime.utcnow()

datadir = "/home/sat_ops/goesR/data/sst/rast_sst/" + str(nowdate.year) + "/"
filenames = sorted([f for f in listdir(datadir) if isfile(join(datadir, f))]) #

for dataset in filenames:
    dataset = dataset[3:]
    # rename in case they added an M6 instead of M3
    if str(dataset).split('_')[1] != 'ABI-L2-SSTF-M3':
        dataset_name = str(dataset).split('_')[0] + '_ABI-L2-SSTF-M3_G16' + str(dataset).split('G16')[1]
    else:
        dataset_name = str(dataset)
    if os.path.isfile("/data/GOES/GOES-R/sst/" + str(nowdate.year) + "/" + str(dataset_name)) == False:
        if os.path.isfile("/data/GOES/GOES-R/suspect/"+ str(dataset_name)) == False:
            print(str(dataset))
            d = Dataset("/home/sat_ops/goesR/data/sst/raw/" + str(nowdate.year) + "/" + dataset)
            ds = NetCDF4DataStore(d)
            ds = xr.open_dataset(ds)
            dat = ds.metpy.parse_cf('SST')
            proj = dat.metpy.cartopy_crs
            dat_dqf  = xr.open_dataset("/home/sat_ops/goesR/data/sst/rast_dqf/" + str(nowdate.year) + "/DQF" + str(dataset))
            new = xr.open_dataset("/home/sat_ops/goesR/data/sst/rast_sst/" + str(nowdate.year) + "/SST" + str(dataset))
            
            dat15 = xr.open_dataset("/home/sat_ops/goesR/data/sst/rast_b15/" + str(nowdate.year) + "/B15" + str(dataset))
            
            new = new.where((new['latitude']>16) & (new['latitude']<52) & (new['longitude']>-100) & (new['longitude']<-50), drop=True)
            dat_dqf = dat_dqf.where((dat_dqf['latitude']>16) & (dat_dqf['latitude']<52) & (dat_dqf['longitude']>-100) & (dat_dqf['longitude']<-50), drop=True)
            dat15 = dat15.where((dat15['latitude']>16) & (dat15['latitude']<52) & (dat15['longitude']>-100) & (dat15['longitude']<-50), drop=True)
    
            # now write it all to netcdf!
            f = Dataset("/data/GOES/GOES-R/sst/" + str(nowdate.year) + "/" + dataset_name,'w', format='NETCDF4') #'w' stands for write
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
            
        
            band15 = f.createVariable('Band15', 'f4', ('time', 'latitude', 'longitude'))
            band15.long_name = "GOES-16 Band 15 Brightness Temperature"
            band15.standard_name = "brightness_temperature"
            band15.units = "kelvin"
            band15.valid_min = 0
            band15.valid_max = 4095
            
            # data 
            latitude[:] = new['latitude'].values
            longitude[:] = new['longitude'].values
            sst[:]= new['sst'].values
            dqf[:]=dat_dqf['dqf'].values
            band15[:]=dat15['b15'].values
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
            fname = "/data/GOES/GOES-R/sst/" + str(nowdate.year) + "/" + str(dataset_name)
            temp_str = "ncks " + fname + " " + fname + " -L 5 -O"
            os.system(temp_str)
            
            flip_lat = "Rscript /home/sat_ops/goesR/data/sst/scripts/flip_lat_sst.R " + fname
            os.system(flip_lat)
    
        else:
            print('all caught up!')


# quality control 
ds = xr.open_dataset("/data/GOES/GOES-R/sst/" + str(nowdate.year) + "/" + dataset_name)
sst = ds.metpy.parse_cf("SST")
sst = sst.sel(longitude=slice(-88, -50), latitude=slice(16,18))
sst = sst.where(sst.values > 270, np.nan)
sst = sst.fillna(-999)
if np.percentile(sst.values, 60) == -999:
    print(sst.time.values)
    move_files = "mv " + "/data/GOES/GOES-R/sst/" + str(nowdate.year) + "/" + dataset_name + " /data/GOES/GOES-R/suspect/"
    os.system(move_files)
    print("moving bad files")
if np.percentile(sst.values, 60) != -999:
    dqf = ds.metpy.parse_cf("DQF")
    dqf = dqf.sel(longitude=slice(-65, -50), latitude=slice(37,45))
    dqf = dqf.where(dqf.values > 1, np.nan)
    dqf = dqf.fillna(-999)
    if np.percentile(dqf.values, 60) == -999:
        move_files = "mv " + "/data/GOES/GOES-R/sst/" + str(nowdate.year) + "/" + dataset_name + " /data/GOES/GOES-R/suspect/"
        os.system(move_files)


# make a current daily compsite of the most recent file
outpath = "/data/GOES/GOES-R/daily_composite/"
today = pd.datetime.today()
goes_nc = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
goes_nc = goes_nc.sel(time=datetime.strftime(today.date(), '%Y-%m-%d'))
goes_nc = goes_nc.drop('Band15')
for t in range(len(goes_nc.time.values)):
    x = goes_nc['SST'][t]
    goes_nc['SST'][t] = x.where(goes_nc['DQF'][t] == 0)
    x = goes_nc['DQF'][t]
    goes_nc['DQF'][t] = x.where(goes_nc['DQF'][t] == 3)

goes_nc = goes_nc.drop(['DQF'])
goes_nc['sst'] = goes_nc['SST']
goes_nc = goes_nc.drop(['SST'])
# resample to a daily composite
goes_nc = goes_nc.resample(time='1D').mean('time')
goes_nc.to_netcdf(path=outpath + '/' + str(datelist[d].year) + '/GOES16_SST_dailycomposite_' + str(datelist[d].year) + str("{0:0=3d}".format(datelist[d].dayofyear)) + '.nc', format='NETCDF3_CLASSIC')



