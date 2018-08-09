from datetime import datetime
import cartopy.feature as cfeature
from siphon.catalog import TDSCatalog
import matplotlib.pyplot as plt
from matplotlib import patheffects
import metpy
from metpy.plots import colortables
import xarray as xr
from xarray.backends import NetCDF4DataStore
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# defining sizing for plotting stuff
toptext = 0.857
toptextleft = 0.13
toptextright = 0.7
bottomtextleft = 0.13
bottomtextheight = 0.155
toprecx = 0.125
toprecy = 0.78
bottomrecx = 0.125
bottomrecy = 0.245

%matplotlib inline

nowdate = datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/SeaSurfaceTemperature/FullDisk/' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')
dataset_name = sorted(cat.datasets.keys())[-1]
dataset = cat.datasets[dataset_name]
ds = dataset.remote_access(service='OPENDAP')
ds = NetCDF4DataStore(ds)
ds = xr.open_dataset(ds)



dat = ds.metpy.parse_cf('SST')
proj = dat.metpy.cartopy_crs

dat = dat.where(dat.variable > 274)
dat = dat.where(dat.variable < 310)
dat = dat - 273.15
# Plot in Mercator
import cartopy.crs as ccrs
newproj = ccrs.Mercator()

fig = plt.figure(figsize=[16, 12], dpi=100)
ax = fig.add_subplot(1,1,1, projection=newproj)
im = ax.pcolormesh(dat['x'], dat['y'], dat, cmap='jet', transform=proj, vmin=-2, vmax=38)
ax.set_extent((dat['x'].min() + 3000000, dat['x'].max()- 1000000, dat['y'].min()+ 5500000, dat['y'].max()- 1000000), crs=proj)
cbaxes = inset_axes(ax, width="100%", height="4%", loc='lower center', borderpad=0) 
cb1 = fig.colorbar(im, orientation='horizontal', cax=cbaxes, ticks=[0, 5, 10, 15, 20, 25, 30, 35])
cb1.ax.set_xticklabels(['0', '5', '10', '15', '20', '25', '30', '35'])
cb1.outline.set_visible(False) # Remove the colorbar outline
cb1.ax.tick_params(width = 0) # Remove the colorbar ticks 
cb1.ax.xaxis.set_tick_params(pad=-15.5) # Put the colobar labels inside the colorbar
ax.set_title("")
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, edgecolor='gray')
ax.add_feature(cfeature.STATES, linestyle=':', edgecolor='gray')

# top rectangle
fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.775,0.03,
                              fill=True, color='darkslateblue', alpha=1, zorder=1000,
                              transform=fig.transFigure, figure=fig)])
# bottom rectangle
fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.775,0.03,
                              fill=True, color='darkslateblue', alpha=1, zorder=1000,
                              transform=fig.transFigure, figure=fig)])
                              
fig.text(toptextleft, bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=10, zorder=2000)

fig.savefig("Downloads/sst.png", dpi=100, bbox_inches='tight')







