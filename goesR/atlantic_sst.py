
from datetime import datetime, timedelta
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
import os.path
import os
from os import listdir
from os.path import isfile, join

# first, go ahead and run the download/reproject scripts
import os
import sys
os.system("sh /home/sat_ops/goesR/data/sst/scripts/sh_reproject_sst.sh")


# now get everything ready for plotting
nowdate = datetime.utcnow()
datadir = '/home/sat_ops/goesR/data/sst/raw/' + str(nowdate.year) + '/'
imgdir = '/home/sat_ops/goesR/products/img_sst/'
file_names = [f for f in listdir(datadir) if isfile(join(datadir, f))][-24:]

for dataset in file_names:
    dataset = dataset[:-3]
    sstimg = imgdir + dataset + '.png'  # GOES16 East
    if os.path.isfile(sstimg) == False:
        print(dataset)
        ds = Dataset(datadir + dataset + '.nc')
        ds = NetCDF4DataStore(ds)
        ds = xr.open_dataset(ds)
        
        
        dqf = ds.metpy.parse_cf('DQF')[0]
        dat = ds.metpy.parse_cf('SST')[0]
        proj = dat.metpy.cartopy_crs

        import cartopy.crs as ccrs
        newproj = ccrs.Mercator()
        
        #dat = dat.where(dqf == 0)
        dat = dat.where(dat.variable > 271)
        dat = dat.where(dat.variable < 313)
        dat = dat - 273.15
        # Plot in Mercator

        yr = dataset[55:59]
        jday = dataset[59:62]
        yrjday = yr + jday
        mon = datetime.strptime(yrjday, '%Y%j').month
        day = datetime.strptime(yrjday, '%Y%j').day
    
        ymd = yr + str(mon) + str(day)
        hms = dataset[62:68]
        timestamp = ymd + " " + hms
        timestamp = datetime.strptime(timestamp, "%Y%m%d %H%M%S")
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
        
        
        
        
        toptext = 0.867
        toptextleft = 0.13
        toptextright = 0.76
        bottomtextleft = 0.13
        bottomtextheight = 0.15
        toprecx = 0.125
        toprecy = 0.86
        bottomrecx = 0.125
        bottomrecy = 0.146
        
        fig = plt.figure(figsize=[12, 12], dpi=100)
        ax = fig.add_subplot(1,1,1, projection=newproj)
        im = ax.pcolormesh(dat['x'], dat['y'], dat, cmap='jet', transform=proj, vmin=-2, vmax=40)
        ax.set_extent((dat['x'].min() + 4000000, dat['x'].max()- 3200000, dat['y'].min()+ 5500000, dat['y'].max()- 650000), crs=proj)
        
        cbaxes = inset_axes(ax, width="100%", height="3%", loc='lower center', borderpad=0) 
        cb1 = fig.colorbar(im, orientation='horizontal', cax=cbaxes, ticks=[0, 5, 10, 15, 20, 25, 30, 35])
        cb1.ax.set_xticklabels(['0', '5', '10', '15', '20', '25', '30', '35'])
        cb1.outline.set_visible(False) # Remove the colorbar outline
        cb1.ax.tick_params(width = 0) # Remove the colorbar ticks 
        cb1.ax.xaxis.set_tick_params(pad=-16) # Put the colobar labels inside the colorbar
        ax.set_title("")
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                        edgecolor='black', facecolor='none',linewidth=1))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                        edgecolor='black', facecolor='none',linewidth=1))
        ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                        edgecolor='black', facecolor='none',linewidth=1))
        # top rectangle
        fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.774,0.022,
                                      fill=True, color='darkslateblue', alpha=1, zorder=1000,
                                      transform=fig.transFigure, figure=fig)])
        # bottom rectangle
        fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.774,0.022,
                                      fill=True, color='darkslateblue', alpha=1, zorder=1000,
                                      transform=fig.transFigure, figure=fig)])
        
        clabeltext='Sea Surface Temperature [Celsius]'
        title = 'GOES16 Raw Sea Surface Temperature'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        
        fig.text(toptextleft, bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
        fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
        fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=10, zorder=2000)
        
        fig.savefig(imgdir + dataset + '.png', dpi=100, bbox_inches='tight')
        plt.close()
import imageio
import numpy as np

img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-24:]

imglen = len(img_names)
images = []
dur_vals = []
for i in range(0,imglen - 1):
    dur_vals.append(.09)
dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goesR/products/output/' + 'raw_sst_goes16.gif', images, duration=dur_vals)






























