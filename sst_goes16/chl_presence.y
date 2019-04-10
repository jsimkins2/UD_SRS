from datetime import datetime, timedelta
import cartopy.feature as cfeature
from siphon.catalog import TDSCatalog
import matplotlib.pyplot as plt
from matplotlib import patheffects
import metpy
from metpy.plots import colortables
import xarray as xr
from xarray.backends import NetCDF4DataStore
from netCDF4 import Dataset, num2date
from dateutil import tz
import time
from time import mktime
import os.path
import os
from os import listdir
from os.path import isfile, join
import cartopy.crs as ccrs
import pandas as pd

# compare with viirs blue band which is 486_qaa
goesdata1 = xr.open_dataset('Downloads/ABI-L1b-RadC_2019_091_16_OR_ABI-L1b-RadC-M3C01_G16_s20190911602156_e20190911604529_c20190911604573.nc')
band1 = goesdata1.metpy.parse_cf("Rad")
proj = band1.metpy.cartopy_crs
newproj = ccrs.Mercator()
ts = pd.to_datetime(band1.t.values)

# compare with viirs green band which is 862_qaa
goesdata3 = xr.open_dataset('Downloads/ABI-L1b-RadC_2019_091_16_OR_ABI-L1b-RadC-M3C03_G16_s20190911602156_e20190911604529_c20190911604571.nc')
band3 = goesdata3.metpy.parse_cf("Rad")
#band3 = band3.where(band1 < 120)


goesdata2 = xr.open_dataset('Downloads/ABI-L1b-RadC_2019_091_16_OR_ABI-L1b-RadC-M3C02_G16_s20190911602156_e20190911604529_c20190911604566.nc')
band2 = goesdata3.metpy.parse_cf("Rad")
#band2 = band2.where(band1 < 120)

#viirs = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/Aqua1DayAggregate.nc')
viirs = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc')
viirs = viirs.sel(time=slice('2019-04-01', '2019-04-01'))
chl = viirs['chl_oc3']
tschl = pd.to_datetime(chl.time.values)


fig = plt.figure(figsize=[12,12], dpi=70)
ax1 = fig.add_subplot(2,2,1, projection=newproj)
im1 = ax1.pcolormesh(band1['x'], band1['y'], band1, transform=proj, cmap='jet', vmin=70, vmax = 100)
ax1.set_extent((-74, -76, 38, 40))
ax1.set_title('GOES16 Blue Band ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=10)
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax1.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax1.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
plt.colorbar(im1, fraction=0.046, pad=0.04, label = "L1b Radiances W m-2 sr-1 um-1")

ax2 = fig.add_subplot(2,2,2, projection=newproj)
im2 = ax2.pcolormesh(band1['x'], band1['y'], band3, transform=proj, cmap='jet', vmin=0, vmax = 20)
ax2.set_extent((-74, -76, 38, 40))
ax2.set_title('GOES16 Veggie Band ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=10)
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax2.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax2.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
plt.colorbar(im2, fraction=0.046, pad=0.04, label = "L1b Radiances W m-2 sr-1 um-1")

ax = fig.add_subplot(2,2,3, projection=newproj)
im = ax.pcolormesh(band1['x'], band1['y'], band2, transform=proj, cmap='jet', vmin=0, vmax = 20)
ax.set_extent((-74, -76, 38, 40))
ax.set_title('GOES16 Red Band ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=10)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
plt.colorbar(im, fraction=0.046, pad=0.04, label = "L1b Radiances W m-2 sr-1 um-1")
#plt.savefig("Downloads/redband.png", dpi=100)


ax = fig.add_subplot(2,2,4, projection=newproj)
im = ax.pcolormesh(chl.lon, chl.lat, chl[0], transform=ccrs.PlateCarree(), cmap='jet', vmin=0, vmax = 15)
ax.set_extent((-74, -76, 38, 40))
ax.set_title('VIIRS CHL_OC3 ' + '2019-04-01' + ' UTC', size=10)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
plt.colorbar(im, fraction=0.046, pad=0.04, label = "chl_oc3 mg m^-3")
plt.savefig("Downloads/chl_presence_vmaxlow.png", dpi=100)





tschl.tolist().strftime('%Y-%m-%d')


