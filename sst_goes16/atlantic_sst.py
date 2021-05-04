# this script makes an image of SST over the NW Atlantic
from datetime import datetime, timedelta
import cartopy.feature as cfeature
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
import pandas as pd
import cartopy
import matplotlib.ticker as mticker
import cmocean
# open the main dataset
main_sst = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/GOESR_SST_DAILY.nc")
main_sst = main_sst.reindex(latitude=list(reversed(main_sst.latitude)))

dat = main_sst.metpy.parse_cf('sst')[-1]
dat = dat * (9/5) + 32
proj = dat.metpy.cartopy_crs


timestamp = pd.to_datetime(dat.time.values)
#from_zone = tz.gettz('UTC')
#to_zone = tz.gettz('America/New_York')
#utc = timestamp.replace(tzinfo=from_zone)
#local = utc.astimezone(to_zone)

local = str(timestamp.month).zfill(2) + '-' + str(timestamp.day).zfill(2) + '-' + str(timestamp.year)

#######################################################################
#######################################################################
####################### CONUS Plotting ################################
#######################################################################
#######################################################################
# Define Plotting Locations


toprecx = 0.125
toprecy = 0.752
toptext = toprecy + 0.008
toptextleft = toprecx + 0.005
toptextright = toprecx + 0.69
bottomtextleft = toprecx + 0.005
bottomtextheight = 0.15
imgdir = "/var/www/html/imagery/ocean/"


fig = plt.figure(figsize=[16, 16], dpi=100)
ax = fig.add_subplot(1,1,1, projection=proj)
ax.set_extent((dat['longitude'].min(), dat['longitude'].max(), dat['latitude'].min(), dat['latitude'].max()), crs=proj)

im = ax.pcolormesh(dat['longitude'], dat['latitude'], dat, cmap=cmocean.cm.thermal, vmin=32, vmax=90)
cbaxes = inset_axes(ax, width="100%", height="3%", loc='lower center', borderpad=0) 
cb1 = fig.colorbar(im, orientation='horizontal', cax=cbaxes, ticks=[32, 40, 50, 60, 70, 80, 90])
cb1.ax.set_xticklabels(['','40' +  u"\u00b0" +  'F',  '50' +  u"\u00b0" +  'F',  '60' +  u"\u00b0" +  'F',
                        '70'+  u"\u00b0" +  'F',  '80'+  u"\u00b0" +  'F', ''],
                        fontweight='bold')
cb1.outline.set_visible(False) # Remove the colorbar outline
cb1.ax.tick_params(width = 0) # Remove the colorbar ticks 
cb1.ax.xaxis.set_tick_params(pad=-16) # Put the colobar labels inside the colorbar
gl = ax.gridlines(crs=proj, linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
gl.xlabels_bottom = False
gl.ylabels_left = False
ax.set_title("")
ax.add_feature(cartopy.feature.LAND,facecolor='lightslategray')
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=1))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.774,0.022,
                              fill=True, color='darkslateblue', alpha=1, zorder=1000,
                              transform=fig.transFigure, figure=fig)])

title = 'GOES16 Sea Surface Temperature Daily Composite [Fahrenheit]'
timestr = local

fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000, fontweight='bold')
fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=10, zorder=2000, fontweight='bold')

fig.savefig(imgdir + "G16_atlantic" + '.png', dpi=100, bbox_inches='tight')

#######################################################################
#######################################################################
####################### Mid Atlantic Plotting ################################
#######################################################################
#######################################################################
# Define Plotting Locations


toprecx = 0.125
toprecy = 0.828
toptext = toprecy + 0.008
toptextleft = toprecx + 0.005
toptextright = toprecx + 0.69
bottomtextleft = toprecx + 0.005
bottomtextheight = 0.15
imgdir = "/var/www/html/imagery/ocean/"


fig = plt.figure(figsize=[12, 12], dpi=100)
ax = fig.add_subplot(1,1,1, projection=proj)
ax.set_extent((-80, -68, 33, 44), crs=proj)

im = ax.pcolormesh(dat['longitude'], dat['latitude'], dat, cmap=cmocean.cm.thermal, vmin=32, vmax=90)
cbaxes = inset_axes(ax, width="100%", height="3%", loc='lower center', borderpad=0) 
cb1 = fig.colorbar(im, orientation='horizontal', cax=cbaxes, ticks=[32, 40, 50, 60, 70, 80, 90])
cb1.ax.set_xticklabels(['','40' +  u"\u00b0" +  'F',  '50' +  u"\u00b0" +  'F',  '60' +  u"\u00b0" +  'F',
                        '70'+  u"\u00b0" +  'F',  '80'+  u"\u00b0" +  'F', ''],
                        fontweight='bold')
cb1.outline.set_visible(False) # Remove the colorbar outline
cb1.ax.tick_params(width = 0) # Remove the colorbar ticks 
cb1.ax.xaxis.set_tick_params(pad=-16) # Put the colobar labels inside the colorbar
gl = ax.gridlines(crs=proj, linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
gl.xlabels_bottom = False
gl.ylabels_left = False
ax.set_title("")
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=1))
#ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m',
#                               edgecolor='black', facecolor='gray',linewidth=1))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='lightslategray',linewidth=0.5))
                                
fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.774,0.022,
                              fill=True, color='darkslateblue', alpha=1, zorder=1000,
                              transform=fig.transFigure, figure=fig)])

title = 'GOES16 Sea Surface Temperature Daily Composite [Fahrenheit]'
timestr = local

fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000, fontweight='bold')
fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=10, zorder=2000, fontweight='bold')

fig.savefig(imgdir + "G16_delaware" + '.png', dpi=100, bbox_inches='tight')






