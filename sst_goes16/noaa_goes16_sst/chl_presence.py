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

# compare with viirs blue band which is 486_qaa becuase goes blue band is wavelength of .47 microns
goesdata1 = xr.open_dataset('/Users/james/Downloads/ABI-L1b-RadC_2019_085_16_OR_ABI-L1b-RadC-M3C03_G16_s20190851617164_e20190851619537_c20190851619582.nc')
band1 = goesdata1.metpy.parse_cf("Rad")
proj = band1.metpy.cartopy_crs
newproj = ccrs.Mercator()
ts = pd.to_datetime(band1.t.values)

# open the Rrs_448
#aqua = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/Aqua1DayAggregate.nc')
viirs = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc')
#viirs = xr.open_dataset("/Users/james/Downloads/Viirs.V2019091.mosaic.nc4")
#viirs = viirs.sel(time=slice('2019-04-01', '2019-04-01'))
viirs = viirs.sel(time=slice(ts.strftime('%Y-%m-%d'),ts.strftime('%Y-%m-%d')))
r486 = viirs['Rrs_486']
timev = pd.to_datetime(r486.time.values)
#pandas format to strftime
#str(timev.format(formatter=lambda x: x.strftime('%Y-%m-%d')))


# compare with viirs green band which is 862_qaa because goes green band wavelength of .86 microns
goesdata3 = xr.open_dataset('/Users/james/Downloads/ABI-L1b-RadC_2019_085_16_OR_ABI-L1b-RadC-M3C03_G16_s20190851617164_e20190851619537_c20190851619582.nc')
band3 = goesdata3.metpy.parse_cf("Rad")
r862 = viirs['Rrs_862']



fig = plt.figure(figsize=[12,12], dpi=70)
ax1 = fig.add_subplot(2,2,1, projection=newproj)
im1 = ax1.pcolormesh(band1['x'], band1['y'], band1, transform=proj, cmap='jet', vmin=0, vmax = 60)
ax1.set_extent((-72, -77.7, 37, 41.5))
ax1.set_title('GOES16 Blue Band (470nm) ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=10)
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax1.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax1.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
plt.colorbar(im1, fraction=0.046, pad=0.04, label = "L1b Radiances W m-2 sr-1 um-1")

ax2 = fig.add_subplot(2,2,2, projection=newproj)
im = ax2.pcolormesh(r486.lon, r486.lat, r486[0], transform=ccrs.PlateCarree(), cmap='jet', vmin=0, vmax = 0.008)
ax2.set_extent((-72, -77.7, 37, 41.5))
ax2.set_title('VIIRS Rrs 486 ' + ts.strftime('%Y-%m-%d') + ' UTC', size=10)
ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax2.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax2.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
plt.colorbar(im, fraction=0.046, pad=0.04, label = "sr^-1")

ax = fig.add_subplot(2,2,3, projection=newproj)
im = ax.pcolormesh(band3['x'], band3['y'], band3, transform=proj, cmap='jet', vmin=1, vmax = 10)
ax.set_extent((-72, -77.7, 37, 41.5))
ax.set_title('GOES16 Veggie Band (860nm) ' + ts.strftime('%Y-%m-%d %H:%M') + ' UTC', size=10)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
plt.colorbar(im, fraction=0.046, pad=0.04, label = "L1b Radiances W m-2 sr-1 um-1")



ax = fig.add_subplot(2,2,4, projection=newproj)
im = ax.pcolormesh(r862.lon, r862.lat, r862[0], transform=ccrs.PlateCarree(), cmap='jet', vmin=0, vmax = 0.0006)
ax.set_extent((-72, -77.7, 37, 41.5))
ax.set_title('VIIRS Rrs_862 ' + ts.strftime('%Y-%m-%d')  + ' UTC', size=10)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='white',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
plt.colorbar(im, fraction=0.046, pad=0.04, label = "sr^-1")
plt.savefig("/Users/james/Downloads/goes16_viirs_" + ts.strftime('%m%d') + "_l1b.png", dpi=100, bbox_inches='tight')
