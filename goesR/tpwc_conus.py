from datetime import datetime
import cartopy.feature as cfeature
from siphon.catalog import TDSCatalog
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
# defining sizing for plotting stuff
nowdate = datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/TotalPrecipitableWater/CONUS/' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')
dataset_name = sorted(cat.datasets.keys())[-1]
dataset = cat.datasets[dataset_name]
ds = dataset.remote_access(service='OPENDAP')
ds = NetCDF4DataStore(ds)
ds = xr.open_dataset(ds)

variable = 'TPW' 
vmin = 0
vmax = 60
cmap = "jet"

#dqf = ds.metpy.parse_cf('DQF')
dat = ds.metpy.parse_cf(variable)
proj = dat.metpy.cartopy_crs

dqf = ds.metpy.parse_cf('DQF_Overall')

dat = dat.where(dqf == 0)
dat = dat.where(dat.variable > vmin)
dat = dat.where(dat.variable < vmax)
#dat = dat - 273.15
# Plot in Mercator
import cartopy.crs as ccrs
newproj = ccrs.Mercator()

ymd = ds.time_coverage_end.split("T")[0]
hms = ds.time_coverage_end.split("T")[1][:-3]
timestamp = ymd + " " + hms
timestamp = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
from_zone = tz.gettz('UTC')
to_zone = tz.gettz('America/New_York')
utc = timestamp.replace(tzinfo=from_zone)
local = utc.astimezone(to_zone)

lt = time.localtime()
dst = lt.tm_isdst
lt = time.localtime()
dst = lt.tm_isdst

if dst == 0:
    et = "EDT"
else:
    et = "EST"




toptext = 0.779
toptextleft = 0.13
toptextright = 0.8
bottomtextleft = 0.13
bottomtextheight = 0.22
toprecx = 0.125
toprecy = 0.776
bottomrecx = 0.125
bottomrecy = 0.214

#137776

fig = plt.figure(figsize=[16, 12], dpi=100)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(dat['x'], dat['y'], dat, cmap='jet', transform=proj, vmin=vmin, vmax=vmax)
ax.set_extent((dat['x'].min() + 900000, dat['x'].max() - 50000, dat['y'].min() + 500000, dat['y'].max() - 50000), crs=proj)

cbaxes = inset_axes(ax, width="100%", height="3.1%", loc='lower center', borderpad=0) 
cb1 = fig.colorbar(im, orientation='horizontal', cax=cbaxes, ticks=[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
cb1.ax.set_xticklabels(['5', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55'])
cb1.outline.set_visible(False) # Remove the colorbar outline
cb1.ax.tick_params(width = 0) # Remove the colorbar ticks 
cb1.ax.xaxis.set_tick_params(pad=-15.5) # Put the colobar labels inside the colorbar
ax.set_title("")
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, edgecolor='gray')
ax.add_feature(cfeature.STATES, linestyle=':', edgecolor='gray')

# top rectangle
fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.774,0.02,
                              fill=True, color='darkslateblue', alpha=1, zorder=1000,
                              transform=fig.transFigure, figure=fig)])
# bottom rectangle
fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.774,0.02,
                              fill=True, color='darkslateblue', alpha=1, zorder=1000,
                              transform=fig.transFigure, figure=fig)])
request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/dark_nolabels/{z}/{x}/{y}.png")

clabeltext='Total Precipitable Water [mm]'
title = 'GOES16 Total Precipitable Water'
timestr = local.strftime('%Y-%m-%d %H:%M ') + et

ax.add_image(request, 7, zorder=0, interpolation='none')
fig.text(toptextleft, bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=10, zorder=2000)

fig.savefig("/home/sat_ops/goesR/products/output/tpwc.png", dpi=100, bbox_inches='tight')




