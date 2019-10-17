import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects, ticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as image
from mpl_toolkits.axes_grid1 import make_axes_locatable

import xarray as xr
from netCDF4 import Dataset, num2date
import xarray as xr
from xarray.backends import NetCDF4DataStore
import numpy as np
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import metpy
import pandas as pd

workdir = "/home/james/"
datadir = "/data/ncep_stageIV/indiv/2019"
file_names = sorted([f for f in listdir(datadir) if isfile(join(datadir, f))])
file_names = file_names[-24:]


ens_list = []
for num in range(len(file_names)):
     ens = 'ens%d' % num
     ens_list.append(xr.open_mfdataset(os.path.join(datadir, file_names[num])))
     
ds = xr.concat(ens_list, dim='time')

dat = ds.metpy.parse_cf("Precipitation_Flux")
proj = dat.metpy.cartopy_crs
newproj = ccrs.Mercator()
newtimestamp = dat.time.values[-1]
dat = dat.sum(dim='time')

#######################################################################
#######################################################################
####################### CONUS Plotting ################################
#######################################################################
#######################################################################
# Define Plotting Locations
fs_x = 16
fs_y = 12
dpi = 100
toptext = 0.794
toptextleft = 0.13
toptextright = 0.65
bottomtextleft = 0.13
bottomtextheight = 0.187
toprecx = 0.125
toprecy = 0.786
bottomrecx = 0.125
bottomrecy = 0.178
symbol = u'$\u26A1$'

fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(dat['lon'], dat['lat'], dat.values,
     cmap=metpy.plots.colortables.get_colortable('precipitation'),transform=proj, vmin=0, vmax=25)

ax.set_extent((-65, -128, 21, 47))
ax.set_title("")
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))

divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="3%", pad=-0.6, axes_class=plt.Axes)
fig.add_axes(ax_cb)
cbar = plt.colorbar(im, cax=ax_cb)
cbar.set_label('millimeters precipitation', rotation=90)

# top rectangle
fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7,0.025,
                              fill=True, alpha=1, facecolor='darkslateblue', zorder=3,
                              transform=fig.transFigure, figure=fig)])
fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.7,0.025,
                      fill=True, facecolor='darkslateblue', zorder=3, alpha=1,
                      transform=fig.transFigure, figure=fig)])
title = 'NCEP Stage IV Precipitation Totals - 24hr Aggregation - Powered by CEMA'
timestr = str(pd.to_datetime(newtimestamp)) + ' UTC'
clabeltext = "Last 24hr Precipitation Totals (mm)"
fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=14, zorder=2000)
fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=14, zorder=2000)
fig.text(bottomtextleft,bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=14, zorder=2000)
im1 = image.imread("/home/james/ncep_stageIV/zfolder/udelcemagoes38.png")
plt.figimage(im1, 20, 58, zorder=1)
ax.outline_patch.set_visible(False)
ax.background_patch.set_visible(False)
fig.savefig(workdir + "/ncep_stageIV/ncepStageIV_conus.png", dpi=dpi, bbox_inches='tight', transparent=False)
plt.close()

#######################################################################
#######################################################################
####################### Mid Atlantic ##################################
#######################################################################
#######################################################################
fs_x = 8
fs_y = 8
dpi = 100
toptext = 0.861
textleft = 0.137
toptextright = 0.69
bottomtextleft = 0.139
bottomtextheight = 0.115
toprecx = 0.1352
toprecy = 0.854
bottomrecx = 0.1352
bottomrecy = 0.11

fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(dat['lon'], dat['lat'], dat.values,
     cmap=metpy.plots.colortables.get_colortable('precipitation'),transform=proj, vmin=0, vmax=25)

ax.set_extent((-69, -81, 34.5, 44))
ax.set_title("")
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))

divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="3%", pad=-0.6, axes_class=plt.Axes)
fig.add_axes(ax_cb)
cbar = plt.colorbar(im, cax=ax_cb)
cbar.set_label('millimeters precipitation', rotation=90)

# top rectangle
fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7538,0.025,
                              fill=True, alpha=1, facecolor='darkslateblue', zorder=3,
                              transform=fig.transFigure, figure=fig)])
fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.7538,0.025,
                      fill=True, facecolor='darkslateblue', zorder=3, alpha=1,
                      transform=fig.transFigure, figure=fig)])
title = 'NCEP Stage IV Precipitation Totals - 24hr Aggregation - Powered by CEMA'
timestr = str(pd.to_datetime(newtimestamp)) + ' UTC'
clabeltext = "Last 24hr Precipitation Totals (mm)"

fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
fig.text(textleft, toptext,title,horizontalalignment='left', color = 'white', size=9, zorder=2000)

fig.text(bottomtextleft, bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
im1 = image.imread("/home/james/ncep_stageIV/zfolder/udelcemagoes24.png")
plt.figimage(im1, 15, 34,zorder=1)

ax.outline_patch.set_visible(False)
ax.background_patch.set_visible(False)

fig.savefig(workdir + "/ncep_stageIV/ncepStageIV_midatlantic.png", dpi=dpi, bbox_inches='tight', transparent=False)
plt.close()