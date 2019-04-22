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

nowdate = datetime.utcnow()

def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

def JulianDate_to_MMDDYYY(y,jd):
    month = 1
    day = 0
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month = month + 1
    return month,jd,y

datadir = "/home/sat_ops/goesR/data/sst/temp2/"
sst_data = sorted([f for f in listdir(datadir) if isfile(join(datadir, f))])
datadir = "/home/sat_ops/goesR/data/sst/temp/"
b15data = sorted([f for f in listdir(datadir) if isfile(join(datadir, f))]) #

t = "OR_ABI-L2-MCMIPF-M3_G16_s20191081020223.nc"
ldatetime = []
for t in b15data:
    jday = t.split('_s')[1][:-3][4:7]
    year = t.split('_s')[1][:-3][0:4]
    mdy = JulianDate_to_MMDDYYY(int(year),int(jday))
    hms = t.split('_s')[1][:-3][7:13]
    t = str(mdy[2]) + jday + hms
    ldatetime.append(datetime.strptime(t, '%Y%j%H%M%S'))

ldatetime = sorted(ldatetime)
for fname in sst_data:
    ds = xr.open_dataset(fname)
    dat = ds.metpy.parse_cf('SST')
    proj = dat.metpy.cartopy_crs
    dat_dqf = ds.metpy.parse_cf('DQF')

    dat = dat.where(dat > -1)
    dat.values[np.isnan(dat.values)] = -999
    # Now grab band 15
    ftime = fname.split('_')[3][1:-3]

    gdatetime=datetime.strptime(ftime, '%Y%j%H%M')
    b15index= ldatetime.index(nearest(ldatetime, gdatetime))

    ds = Dataset("/home/sat_ops/goesR/data/sst/temp/" + b15data[b15index])
    print(str(b15data[b15index]) + " matched with " + str(fname))
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
