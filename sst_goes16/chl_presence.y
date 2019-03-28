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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.io.img_tiles as cimgt
from dateutil import tz
import time
from time import mktime
import os.path
import os
from os import listdir
from os.path import isfile, join
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

goesdata1 = xr.open_dataset('Downloads/ABI-L1b-RadC_2019_086_17_OR_ABI-L1b-RadC-M3C01_G16_s20190861702162_e20190861704535_c20190861704580.nc')
band1 = goesdata1.metpy.parse_cf("Rad")

goesdata3 = xr.open_dataset('Downloads/ABI-L1b-RadC_2019_086_17_OR_ABI-L1b-RadC-M3C03_G16_s20190861702162_e20190861704535_c20190861704582.nc')
band3 = goesdata3.metpy.parse_cf("Rad")
proj = band1.metpy.cartopy_crs
newproj = ccrs.Mercator()
ts = pd.to_datetime(band3.t.values)


band3 = band3.where(band1 < 120)
fig = plt.figure(figsize=[12,12], dpi=100)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(band1['x'], band1['y'], band3, transform=proj, cmap='jet', vmin=0, vmax = 50)
ax.set_extent((-74, -76, 38, 40))

ax.set_title('GOES16 Veggie Band ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=24)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '10m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))

plt.colorbar(im, fraction=0.046, pad=0.04, label = "L1b Radiances W m-2 sr-1 um-1")
plt.savefig("Downloads/chl_goes16_veggie_band.png", dpi=100)


goesdata2 = xr.open_dataset('Downloads/ABI-L1b-RadC_2019_086_17_OR_ABI-L1b-RadC-M3C02_G16_s20190861702162_e20190861704535_c20190861704576.nc')
band2 = goesdata3.metpy.parse_cf("Rad")
band2 = band2.where(band1 < 120)

fig = plt.figure(figsize=[12,12], dpi=100)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(band1['x'], band1['y'], band2, transform=proj, cmap='jet', vmin=0, vmax = 90)
ax.set_extent((-74, -76, 38, 40))

ax.set_title('GOES16 Red Band ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=24)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '10m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))

plt.colorbar(im, fraction=0.046, pad=0.04, label = "L1b Radiances W m-2 sr-1 um-1")
plt.savefig("Downloads/chl_goes16_red_band.png", dpi=100)





viirs = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc')
viirs = viirs.sel(time=slice('2019-03-27', '2019-03-27'))
chl = viirs['chl_oc3']
proj = ccrs.PlateCarree()

fig = plt.figure(figsize=[12,12], dpi=100)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(chl.lon, chl.lat, chl[0], transform=proj, cmap='jet', vmin=0, vmax = 5)
ax.set_extent((-74, -76, 38, 40))

ax.set_title('VIIRS CHL_OC3 ' + str(chl.time.values[0]) + ' UTC', size=24)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '10m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))

plt.colorbar(im, fraction=0.046, pad=0.04, label = "chl_oc3")
plt.savefig("Downloads/chl_viirs.png", dpi=100)






