import xarray as xr
import metpy
from datetime import datetime
from xarray.backends import NetCDF4DataStore
from netCDF4 import Dataset, num2date
from siphon.catalog import TDSCatalog
from pyproj import Proj
import cartopy.crs as ccrs
import os.path
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import os.path
import os
import sys

datadir = '/home/sat_ops/goesR/radar/prectype/hrrr_temp/'
nowdate = datetime.utcnow()
cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/latest.xml')
dataset_name = sorted(cat.datasets.keys())[-1]
dataset = cat.datasets[dataset_name]


if os.path.isfile(datadir + dataset_name[:-5] + 'nc') == False:
    os.system("rm " + datadir + "*")
    print('downloading ' + dataset_name)
    ds = dataset.remote_access(service='OPENDAP')
    ds = NetCDF4DataStore(ds)
    ds = xr.open_dataset(ds)
    tempiso = ds.metpy.parse_cf('Temperature_isobaric')
    hproj = tempiso.metpy.cartopy_crs
    hproj = Proj(hproj.proj4_init)
    
    
    # Write the netcdf file
    
    f = Dataset(datadir + dataset_name[:-5] + 'nc','w', format='NETCDF4') #'w' stands for write
    # dimensions
    f.createDimension('x', tempiso['x'].size)
    f.createDimension('y', tempiso['y'].size)
    f.createDimension('time', 1)
    
    # variables
    
    x = f.createVariable('x', 'f4', 'x')
    x.standard_name = tempiso['x'].standard_name
    x.units = tempiso['x'].units
    
    y = f.createVariable('y', 'f4', 'y')
    y.standard_name = tempiso['y'].standard_name
    y.units = tempiso['y'].units
    
    time = f.createVariable('time', 'f8', 'time')
    time.standard_name = tempiso['time'][1].standard_name
    time.long_name = tempiso['time'][1].long_name
    
    proj = f.createVariable('projection', 'f4')
    proj.crs = 'lambert_conformal_conic'
    proj.srs = hproj.srs
    
    tsurf = f.createVariable('temperature_surface', 'f4', ('time', 'y', 'x'))
    tsurf.long_name = 'Temperature @ Surface'
    tsurf.standard_name = 'temperature_surface'
    tsurf.units = 'Kelvin'
    
    t925 = f.createVariable('temperature_925', 'f4', ('time', 'y', 'x'))
    t925.long_name = 'Temperature @ 925mb'
    t925.standard_name = 'temperature_925'
    t925.units = 'Kelvin'
    
    t850 = f.createVariable('temperature_850', 'f4', ('time', 'y', 'x'))
    t850.long_name = 'Temperature @ 850mb'
    t850.standard_name = 'temperature_850'
    t850.units = 'Kelvin'
    
    ht500 = f.createVariable('height_500', 'f4', ('time', 'y', 'x'))
    ht500.long_name = 'Geopotential_height_isobaric @ 500mb'
    ht500.standard_name = 'height_500'
    ht500.units = 'gpm'
    
    ht1000 = f.createVariable('height_1000', 'f4', ('time', 'y', 'x'))
    ht1000.long_name = 'Geopotential_height_isobaric @ 1000mb'
    ht1000.standard_name = 'height_1000'
    ht1000.units = 'gpm'
    
    # fill the data
    x[:] = tempiso['x'].values
    y[:] = tempiso['y'].values
    time[:] = tempiso['time'][1].values
    tsurf[:]= ds.metpy.parse_cf('Temperature_height_above_ground')[1][0].values
    t925[:]=tempiso[1][3].values
    t850[:]=tempiso[1][2].values
    ht500[:]=ds.metpy.parse_cf("Geopotential_height_isobaric")[1][0].values
    ht1000[:]=ds.metpy.parse_cf("Geopotential_height_isobaric")[1][3].values
    f.close()
    
    # now run the reprojecting script 
    print('beginning the R reprojection')
    os.system("Rscript /home/sat_ops/goesR/github/UD_SRS/hrrr/reproject_hrrr.R " + dataset_name[:-5] + 'nc')

if os.path.isfile(datadir + 'gridsurf.npy') == False:
    print('regridding the hrrr')
    os.system("python3 /home/sat_ops/goesR/github/UD_SRS/hrrr/regrid_hrrr.py")
    