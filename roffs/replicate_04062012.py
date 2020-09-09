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
#area = [36.5, 40, -75.25, -72.0]
area = [27, 29, -91.5, -88]
extArea = [area[2], area[3], area[0], area[1]]

# load SST
today = datetime.utcnow() - timedelta(days=7)
#goes_thredds = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/GOESR_SST_DAILY.nc")
goes_thredds = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc")
goes_thredds = goes_thredds.sel(lat=slice(area[0],area[1]), lon=slice(area[2], area[3]), time="04-06-2012")
sst = goes_thredds.metpy.parse_cf('sst')[-1]
sst = sst.where(sst.values > 0)
#sst = sst - 273.15
sst = sst*(9/5) + 32
sstproj = sst.metpy.cartopy_crs

# load chlorophyll
viirs = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/viirs_1day_aggregate.nc')
viirs = viirs.sel(lat=slice(area[0],area[1]), lon=slice(area[2], area[3]), time="04-06-2012")
chl_oc3 = viirs.metpy.parse_cf('chl_oc3')[0]
viirsproj = ccrs.PlateCarree()

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
im1 = ax1.pcolormesh(sst['lon'], sst['lat'], sst, transform=sstproj, cmap='rainbow')
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=1, linestyle='--')
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xpadding = -40
gl.ypadding = -50
bbox = dict(boxstyle="round", ec="white", fc="white", alpha=0.5, pad=0.01)
gl.xlabel_style = {'size': 11, 'color': 'black', 'weight': 'bold', 'bbox':bbox}
gl.ylabel_style = {'size':11, 'color': 'black', 'weight': 'bold', 'bbox':bbox}

sstcontours = ax1.contour(sst['lon'], sst['lat'], gaussian_filter(sst,0.5),10, transform=sstproj, colors='black',linewidths=1)
hfont = {'fontname':'monospace'}#'style':'oblique'}
degsym = '\u00b0'
for i in range(5,sst.shape[0] - 5, int(np.round(sst.shape[0]/10))):
    for j in range(5,sst.shape[1] - 5, int(np.round(sst.shape[1]/10))):
        if np.isnan(sst.values[i,j]) == False:
            ax1.text(sst['lon'].values[j],sst['lat'].values[i],str('{:2.0f}'.format(z) + degsym), ha='center', va='center', size=16, 
            transform=sstproj,weight='bold',color='black',**hfont)

bathcontours = ax1.contour(bath['lon'], bath['lat'],bath,20, transform=bathproj, colors='black',linestyle='solid',linewidths=0.5)
for (i,j), z in np.ndenumerate(bath.values):
    divi = random.choice(range(20,40))
    if i % divi == 0 and i > 14 and i < (bath.values.shape[0] - 25):
        if j % divi == 0 and j > 14 and j < (bath.values.shape[1] - 25):
            if np.isnan(bath.values[i,j]) == False:
                ax1.text(bath['lon'].values[j],bath['lat'].values[i],(z*-1), ha='center', va='center', size=8, 
                transform=bathproj,weight='bold',color='black',rotation=-10,**hfont)
ax1.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '10m',
                                edgecolor='black', facecolor='gray',linewidth=0.5))
for i in range(8,chl_oc3.shape[0] - 8, int(np.round(chl_oc3.shape[0]/8))):
    for j in range(8,chl_oc3.shape[1] - 8, int(np.round(chl_oc3.shape[1]/8))):
        if np.isnan(chl_oc3[i,j]) == False:
            if chl_oc3[i,j].values > 0 and chl_oc3[i,j].values < 0.12:
                val = 'BLUE'
            elif chl_oc3[i,j].values > 0.12 and chl_oc3[i,j].values < 0.15:
                val = 'BLENDED \nBLUE'
            elif chl_oc3[i,j].values > 0.15 and chl_oc3[i,j].values < 0.2:
                val = 'BLUE \nGREEN'
            else:
                val = 'GREEN'
            ax1.text(chl_oc3['lon'].values[j],chl_oc3['lat'].values[i],val, ha='center', va='center', size=16, 
            transform=viirsproj,weight='bold',color='black',rotation=00,**hfont)








for (i,j), z in np.ndenumerate(chl_oc3):
    divi = random.choice(range(20,40))
    if i % divi == 0 and i > 14 and i < (chl.shape[0] - 20):
        if j % divi == 0 and j > 14 and j < (chl.shape[1] - 20):
            if np.isnan(chl_oc3[i,j]) == False:
                if chl_oc3[i,j].values > 0 and chl_oc3[i,j].values < 0.12:
                    val = 'BLUE'
                elif chl_oc3[i,j].values > 0.12 and chl_oc3[i,j].values < 0.15:
                    val = 'BLENDED \nBLUE'
                elif chl_oc3[i,j].values > 0.15 and chl_oc3[i,j].values < 0.2:
                    val = 'BLUE \nGREEN'
                else:
                    val = 'GREEN'
                ax1.text(chl_oc3['lon'].values[j],chl_oc3['lat'].values[i],val, ha='center', va='center', size=16, 
                transform=viirsproj,weight='bold',color='black',rotation=00,**hfont)





