import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects, ticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as image
from mpl_toolkits.axes_grid1 import make_axes_locatable
from os.path import isfile, join
import os
from os import listdir
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
from metpy.units import units
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from metpy.plots import colortables
from matplotlib.colors import BoundaryNorm

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
dat.values = (dat.values * units.mm).to(units.inch)
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
toptextright = 0.73
bottomtextleft = 0.13
bottomtextheight = 0.209
toprecx = 0.125

toprecy = 0.789
bottomrecx = 0.1295
bottomrecy = 0.2016
symbol = u'$\u26A1$'
cmap = metpy.plots.colortables.get_colortable('precipitation')
levels = [0,0.01,0.1,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,6,8,10,15]
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(dat['lon'], dat['lat'], dat.values,cmap=cmap, norm=norm,transform=proj, vmin=0, vmax=15)

ax.set_extent((-65, -128, 21, 47))
ax.set_title("")
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))


cbaxes = inset_axes(ax, width="100%", height="3%", loc=3) 
cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks = [0,0.01,0.1,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,6,8,10,15])
cb.ax.set_xticklabels(['0','0.01','0.1','0.3','0.4','0.5','0.75','1','1.5','2','3','4','5','6','8','10','15'])

#cbar = plt.colorbar(im, location = 'bottom',size="2%", pad='-1.95%',ticks = [5,10,15,20,25,30,35,40,45,50,55,60])
#cbar.set_label('millimeters precipitation')
# top rectangle
fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7752,0.0232,
                              fill=True, alpha=1, facecolor='darkslateblue', zorder=3,
                              transform=fig.transFigure, figure=fig)])
fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.7752,0.0232,
                      fill=True, facecolor='darkslateblue', zorder=3, alpha=1,
                      transform=fig.transFigure, figure=fig)])
title = 'NCEP Stage IV Precipitation Totals - 24hr Aggregation'
timestr = str(pd.to_datetime(newtimestamp)) + ' UTC'
clabeltext = "Last 24hr Precipitation Totals (inches)"
fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=14, zorder=2000)
fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=14, zorder=2000)
fig.text(bottomtextleft,bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=14, zorder=2000)
im1 = image.imread("/home/james/ncep_stageIV/zfolder/udelcemagoes38.png")
plt.figimage(im1, 20, 90, zorder=1)
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
toptext = 0.8645
textleft = 0.148
toptextright = 0.67
bottomtextleft = 0.148
bottomtextheight = 0.147
toprecx = 0.1445
toprecy = 0.859
bottomrecx = 0.1445
bottomrecy = 0.1417

fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(dat['lon'], dat['lat'], dat.values,cmap=cmap, norm=norm,transform=proj, vmin=0, vmax=15)

ax.set_extent((-69, -81, 34.5, 44))
ax.set_title("")
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))

cbaxes = inset_axes(ax, width="100%", height="3%", loc=3) 
cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal',ticks = [0,0.01,0.1,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,6,8,10,15])
cb.ax.set_xticklabels(['0','0.01','0.1','0.3','0.4','0.5','0.75','1','1.5','2','3','4','5','6','8','10','15'])


# top rectangle
fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7545,0.0232,
                              fill=True, alpha=1, facecolor='darkslateblue', zorder=3,
                              transform=fig.transFigure, figure=fig)])
fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.7545,0.0232,
                      fill=True, facecolor='darkslateblue', zorder=3, alpha=1,
                      transform=fig.transFigure, figure=fig)])
title = 'NCEP Stage IV Precipitation Totals - 24hr Aggregation'
timestr = str(pd.to_datetime(newtimestamp)) + ' UTC'
clabeltext = "Last 24hr Precipitation Totals (inches)"

fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
fig.text(textleft, toptext,title,horizontalalignment='left', color = 'white', size=10, zorder=2000)

fig.text(bottomtextleft, bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
im1 = image.imread("/home/james/ncep_stageIV/zfolder/udelcemagoes24.png")
plt.figimage(im1, 20, 75,zorder=1)

ax.outline_patch.set_visible(False)
ax.background_patch.set_visible(False)

fig.savefig(workdir + "/ncep_stageIV/ncepStageIV_midatlantic.png", dpi=dpi, bbox_inches='tight', transparent=False)
plt.close()


