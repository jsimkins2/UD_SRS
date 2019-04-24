import xarray as xr
from netCDF4 import Dataset
import metpy
from datetime import datetime, timedelta
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
import pandas as pd

def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

nowdate = datetime.utcnow()
nowdate = nowdate.replace(second=0, microsecond=0)
sstFile = dataset_names[1]

# Use both servers to double check
mainURL = 'https://thredds.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/SeaSurfaceTemperature/FullDisk/'
testURL = 'https://thredds-test.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/SeaSurfaceTemperature/FullDisk/'
urlstr = [mainURL, testURL]

b15data = [f for f in listdir("/home/sat_ops/goesR/data/fulldisk/") if isfile(join("/home/sat_ops/goesR/data/fulldisk/", f))]
ldatetime = []
for t in b15data:
    t = t.split('_s')[1][:-4]
    ldatetime.append(datetime.strptime(t, '%Y%j%H%M%S'))

# Begin loop
for url in urlstr:
    for d in range(0,1):
        temday = nowdate - timedelta(days=d)
        # UCAR Thredds Server call
        cat = TDSCatalog(url + str(temday.year) + str("%02d"%temday.month) + str("%02d"%temday.day) + '/catalog.xml')
        dataset_names = sorted(cat.datasets.keys())

        # check to see if we already have the dataset name saved, note that we change the name to M3 so ERDDAP can sort the data
        for sstFile in dataset_names:
            # rename in case they added an M6 instead of M3
            if str(sstFile).split('_')[1] != 'ABI-L2-SSTF-M3':
                fname = str(sstFile).split('_')[0] + '_ABI-L2-SSTF-M3_G16' + str(sstFile).split('G16')[1]
            else:
                fname = str(sstFile)
            #make sure the file doesn't already exist
            if os.path.isfile("/home/sat_ops/goesR/data/sst/raw/" + str(nowdate.year) + "/" + str(fname)) == False:
                if os.path.isfile("/data/GOES/GOES-R/sst/" + str(nowdate.year) + "/" + str(fname)) == False:
                    if os.path.isfile("/data/GOES/GOES-R/suspect/" + str(fname)) == False:
                        print('downloading ' + fname + ' to raw folder')

                        dataset = cat.datasets[sstFile]
                        ds = dataset.remote_access(service='OPENDAP')
                        d = ds
                        ds = NetCDF4DataStore(ds)
                        ds = xr.open_dataset(ds)
                        dat = ds.metpy.parse_cf('SST')

                        proj = dat.metpy.cartopy_crs
                        dat_dqf = ds.metpy.parse_cf('DQF')

                        dat = dat.where(dat > -1)
                        dat.values[np.isnan(dat.values)] = -999

                        # Now grab band 15
                        gdatetime=datetime.strptime(sstFile.split('_')[3][1:-1], '%Y%j%H%M%S')
                        b15index= ldatetime.index(nearest(ldatetime, gdatetime))

                        ds = Dataset("/home/sat_ops/goesR/data/fulldisk/" + b15data[b15index])
                        ds = NetCDF4DataStore(ds)
                        ds = xr.open_dataset(ds)
                        d2 = ds
                        dat15 = d2.metpy.parse_cf('CMI_C15')

                        f = Dataset("/home/sat_ops/goesR/data/sst/raw/" + str(nowdate.year) + "/" + str(fname),'w', format='NETCDF4') #'w' stands for write
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
                        proj.grid_mapping_name = d.variables['goes_imager_projection'].grid_mapping_name
                        proj.perspective_point_height = d.variables['goes_imager_projection'].perspective_point_height
                        proj.semi_major_axis = d.variables['goes_imager_projection'].semi_major_axis
                        proj.semi_minor_axis = d.variables['goes_imager_projection'].semi_minor_axis
                        proj.inverse_flattening = d.variables['goes_imager_projection'].inverse_flattening
                        proj.latitude_of_projection_origin = d.variables['goes_imager_projection'].latitude_of_projection_origin
                        proj.longitude_of_projection_origin = d.variables['goes_imager_projection'].longitude_of_projection_origin
                        proj.sweep_angle_axis = d.variables['goes_imager_projection'].sweep_angle_axis

                        sst = f.createVariable('SST', 'f4', ('time', 'y', 'x'))
                        sst.long_name = d.variables['SST'].long_name
                        sst.standard_name = d.variables['SST'].standard_name
                        sst.units = d.variables['SST'].units
                        sst.resolution = d.variables['SST'].resolution
                        sst.missing_value = d.variables['SST']._FillValue
                        sst.grid_mapping = d.variables['SST'].grid_mapping

                        dqf = f.createVariable('DQF', 'f4', ('time', 'y', 'x'))
                        dqf.long_name = d.variables['DQF'].long_name
                        dqf.standard_name = d.variables['DQF'].standard_name
                        dqf.units = d.variables['DQF'].units
                        dqf.flag_values = d.variables['DQF'].flag_values
                        dqf.flag_meanings = d.variables['DQF'].flag_meanings
                        dqf.grid_mapping = d.variables['DQF'].grid_mapping

                        band15 = f.createVariable('Band15', 'f4', ('time', 'y', 'x'))
                        band15.long_name = dat15.long_name
                        band15.standard_name = dat15.standard_name
                        band15.units = dat15.units

                        # data
                        x[:] = dat['x'].values
                        y[:] = dat['y'].values
                        sst[:]= dat.values
                        dqf[:]=dat_dqf
                        band15[:]=dat15.values
                        time[:] = dat['t'].values
                        f.close()
                    else:
                        print('File Already Exists in Suspect Folder')
                else:
                    print('File Already Exists in Aggregation')
            else:
                print('File Already Exists in Raw Folder')
