import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects, ticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as image

import xarray as xr
from netCDF4 import Dataset, num2date
import xarray as xr
from xarray.backends import NetCDF4DataStore
import numpy as np
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import metpy

workdir = "/home/sat_ops/"
datadir = "/data/ncep_stageIV/indiv/2019"
file_names = sorted([f for f in listdir(datadir) if isfile(join(datadir, f))])
file_names = file_names[-24:]

ncep_nc = xr.open_dataset("http://thredds.demac.udel.edu/thredds/dodsC/ncep_stage_iv.nc")
# grab the last 24 hours of sst dataset
ncep_nc = ncep_nc.isel(time=range(-24, 0))
dat = ncep_nc.metpy.parse_cf("Precipitation_Flux")
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
toptextright = 0.76
bottomtextleft = 0.13
bottomtextheight = 0.183
toprecx = 0.125
toprecy = 0.786
bottomrecx = 0.125
bottomrecy = 0.178
symbol = u'$\u26A1$'

fig = plt.figure(figsize=[fs_x, fs_y], dpi=dpi)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(dat['lon'], dat['lat'], dat.values, 
     cmap=metpy.plots.colortables.get_colortable('precipitation'),transform=proj)

ax.set_extent((-65, -128, 21, 47))
ax.set_title("")
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=0.5))
                                

# top rectangle
fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.7745,0.025,
                              fill=True, alpha=1, facecolor='darkslateblue', zorder=3,
                              transform=fig.transFigure, figure=fig)])
fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.7745,0.025,
                      fill=True, facecolor='darkslateblue', zorder=3, alpha=1,
                      transform=fig.transFigure, figure=fig)])
title = 'NCEP Stage IV Precipitation Flux - 24hr Aggregation - Powered by Cema'
timestr = newtimestamp

fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=14, zorder=2000)
fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=14, zorder=2000)

fig.text(bottomtextleft,bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=14, zorder=2000)
im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes38.png")
plt.figimage(im1, 23, 58, zorder=1)
ax.outline_patch.set_visible(False)
ax.background_patch.set_visible(False)
output_file = workdir + "ltng_conus/" + ABI_datetime[n] + ".png"
fig.savefig(output_file, dpi=dpi, bbox_inches='tight', transparent=False)
fig.savefig(workdir + "img_conus/" + ABI_datetime[n] + ".png", dpi=dpi, bbox_inches='tight', transparent=False)
plt.close()