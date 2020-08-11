# prepare goes16 sst data for dineof
from datetime import datetime, timedelta
import cartopy.feature as cfeature
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
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import matplotlib.font_manager as font_manager
from scipy.ndimage.filters import gaussian_filter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl
import random
mpl.rcParams["contour.negative_linestyle"] = 'solid'



# load our area, in this case area 2
area = [36.5, 40, -75.25, -72.0]
extArea = [area[2], area[3], area[0], area[1]]

# load SST
today = datetime.utcnow() - timedelta(days=7)
goes_thredds = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/GOESR_SST_DAILY.nc")
goes_thredds = goes_thredds.sel(latitude=slice(area[1],area[0]), longitude=slice(area[2], area[3]), time=datetime.strftime(today, '%Y-%m-%d'))
sst = goes_thredds.metpy.parse_cf('sst')[-1]
sst = sst.where(sst.values > 0)
#sst = sst - 273.15
sst = sst*(9/5) + 32
sstproj = sst.metpy.cartopy_crs

# load chlorophyll
viirs = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc')
viirs = viirs.sel(lat=slice(area[0],area[1]), lon=slice(area[2], area[3]), time=datetime.strftime(today, '%Y-%m-%d'))
r486 = viirs.metpy.parse_cf('Rrs_486')[0]
r486 = r486.where(r486.values > 0) #blue
r862 = viirs.metpy.parse_cf('Rrs_862')[0] #green
r862 = r862.where(r862.values > 0) 
viirsproj = ccrs.PlateCarree()

ocolor = np.zeros((r486.shape[0], r486.shape[1]))
for i in range(ocolor.shape[0]):
    for j in range(ocolor.shape[1]):
        if np.isnan(r486[i,j]) == True or np.isnan(r862[i,j]) == True:
            ocolor[i,j] == np.nan
        else:
            if r486[i,j] > r862[i,j] + np.nanstd(r486):
                ocolor[i,j] = 0 #b
            if r862[i,j] > r486[i,j] + np.nanstd(r862):
                ocolor[i,j] = 1 #g
            else:
                ocolor[i,j] = 2 #bg

# load bathymetry data
bath = xr.open_dataset("Documents/Delaware/bath.nc")
bath=bath.sel(lat=slice(area[0],area[1]), lon=slice(area[2], area[3]))
bath = bath.metpy.parse_cf('Band1')
bathproj = ccrs.PlateCarree()


# Plot it all up
import matplotlib.ticker as mticker
fig = plt.figure(figsize=[17,14], dpi=120)
ax1 = fig.add_subplot(1,1,1, projection=ccrs.Mercator())
ax1.set_extent((extArea),crs=ccrs.PlateCarree())
im1 = ax1.pcolormesh(sst['longitude'], sst['latitude'], sst, transform=sstproj, cmap='rainbow')
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=1, linestyle='--')
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xpadding = -40
gl.ypadding = -50
bbox = dict(boxstyle="round", ec="white", fc="white", alpha=0.5, pad=0.01)
gl.xlabel_style = {'size': 11, 'color': 'black', 'weight': 'bold', 'bbox':bbox}
gl.ylabel_style = {'size':11, 'color': 'black', 'weight': 'bold', 'bbox':bbox}

sstcontours = ax1.contour(sst['longitude'], sst['latitude'], gaussian_filter(sst,0.75),10, transform=sstproj, colors='black',linewidths=2)
hfont = {'fontname':'monospace'}#'style':'oblique'}
degsym = '\u00b0'
for (i,j), z in np.ndenumerate(sst.values):
    divi = random.choice(range(20,30))
    if i % divi == 0 and i > 14 and i < (sst.values.shape[0] - 5):
        if j % divi == 0 and j > 14 and j < (sst.values.shape[1] - 5):
            if np.isnan(sst.values[i,j]) == False:
                ax1.text(sst['longitude'].values[j],sst['latitude'].values[i],str('{:2.0f}'.format(z) + degsym), ha='center', va='center', size=16, 
                transform=sstproj,weight='bold',color='black',**hfont)


bathcontours = ax1.contour(bath['lon'], bath['lat'],bath,20, transform=bathproj, colors='black',linestyle='solid',linewidths=0.4)
for (i,j), z in np.ndenumerate(bath.values):
    divi = random.choice(range(20,40))
    if i % divi == 0 and i > 14 and i < (bath.values.shape[0] - 25):
        if j % divi == 0 and j > 14 and j < (bath.values.shape[1] - 25):
            if np.isnan(bath.values[i,j]) == False:
                ax1.text(bath['lon'].values[j],bath['lat'].values[i],(z*-1), ha='center', va='center', size=8, 
                transform=bathproj,weight='bold',color='black',rotation=-10,**hfont)
ax1.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '10m',
                                edgecolor='black', facecolor='gray',linewidth=0.5))
for (i,j), z in np.ndenumerate(r486):
    divi = random.choice(range(60,70))
    if i % divi == 0 and i > 14 and i < (ocolor.shape[0] - 20):
        if j % divi == 0 and j > 14 and j < (ocolor.shape[1] - 20):
            if np.isnan(ocolor[i,j]) == False:
                if r486[i,j] > (r862[i,j] + (r862[i,j] *.15)):
                    val = 'BLUE'
                else:
                    if r862[i,j] > (r486[i,j] + (r486[i,j] *.15)):
                        val = 'GREEN'
                    else:
                        val = 'BLUE \nGREEN'
                ax1.text(r862['lon'].values[j],r862['lat'].values[i],val, ha='center', va='center', size=10, 
                transform=viirsproj,weight='bold',color='black',rotation=30,**hfont)



