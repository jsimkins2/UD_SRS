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

goesdata1 = xr.open_dataset('Downloads/c01.nc')
band1 = goesdata1.metpy.parse_cf("Rad")
proj = band1.metpy.cartopy_crs
newproj = ccrs.Mercator()
ts = pd.to_datetime(band1.t.values)


goesdata3 = xr.open_dataset('Downloads/c03.nc')
band3 = goesdata3.metpy.parse_cf("Rad")
#band3 = band3.where(band1 < 120)


goesdata2 = xr.open_dataset('Downloads/c02.nc')
band2 = goesdata3.metpy.parse_cf("Rad")
#band2 = band2.where(band1 < 120)

viirs = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc')
viirs = viirs.sel(time=slice('2019-04-03', '2019-04-03'))
chl = viirs['chl_oc3']
proj = ccrs.PlateCarree()


fig = plt.figure(figsize=[7,7], dpi=70)
ax = fig.add_subplot(1,4,1, projection=newproj)
im = ax.pcolormesh(band1['x'], band1['y'], band1, transform=proj, cmap='jet', vmin=50, vmax = 150)
ax.set_extent((-74, -76, 38, 40))
ax.set_title('GOES16 Blue Band ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=8)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))

ax = fig.add_subplot(1,4,2, projection=newproj)
im = ax.pcolormesh(band1['x'], band1['y'], band3, transform=proj, cmap='jet', vmin=0, vmax = 50)
ax.set_extent((-74, -76, 38, 40))
ax.set_title('GOES16 Veggie Band ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=8)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))


ax = fig.add_subplot(1,4,3, projection=newproj)
im = ax.pcolormesh(band1['x'], band1['y'], band2, transform=proj, cmap='jet', vmin=0, vmax = 90)
ax.set_extent((-74, -76, 38, 40))
ax.set_title('GOES16 Red Band ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=24)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
#plt.colorbar(im, fraction=0.046, pad=0.04, label = "L1b Radiances W m-2 sr-1 um-1")


ax = fig.add_subplot(1,4,4, projection=newproj)
im = ax.pcolormesh(chl.lon, chl.lat, chl[0], transform=proj, cmap='jet', vmin=0, vmax = 5)
ax.set_extent((-74, -76, 38, 40))
ax.set_title('VIIRS CHL_OC3 ' + str(chl.time.values[0]) + ' UTC', size=24)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
#plt.colorbar(im, fraction=0.046, pad=0.04, label = "chl_oc3")
#plt.savefig("Downloads/chl_presence.png", dpi=100)






