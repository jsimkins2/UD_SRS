import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from siphon.catalog import TDSCatalog
import urllib
from netCDF4 import Dataset, num2date
from matplotlib import ticker
import matplotlib as mpl
import pyart
from dateutil import tz
import time
from time import mktime
import os.path
import numpy as np
import matplotlib.image as image
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import calendar
from pyproj import Proj     
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, ListedColormap # Linear interpolation for color maps
import imageio
############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/lightning/"
datadir = "/home/sat_ops/goesR/data/glm/"
site = 'KDOX'
################ Grab the Lat/Lon of the site we want ####################
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]
    
# Define projections
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3100, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)

DH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-4.5,llcrnrlon=lon0-5.5,
    urcrnrlat=lat0+5,urcrnrlon=lon0+5.5,resolution='h') 

GLM_files = sorted([f for f in listdir(datadir) if isfile(join(datadir, f))])[-76:]
GLM_names = []

for i in range(0,len(GLM_files)):
    fname = str(GLM_files[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    GLM_names.append(fname)

fnamelist = sorted(GLM_names, key=int)
ABI_datetime = []
for i in fnamelist:    
    if os.path.isfile(workdir + "conus/" + str(i) + ".png") == False:
        ABI_datetime.append(i)

cmapl = LinearSegmentedColormap.from_list('this', ['darkred','darkred','red', 'red','red', 'orangered', 'orangered' ,'orangered', 'darkorange', 'darkorange', 'darkorange', 'yellow', 'yellow'], N = len(fnamelist)/2)

col_list = []
col_u = []
for i in range(cmapl.N):
    rgb = cmapl(i)[:3] # will return rgba, we take only first 3 so we get rgb
    col_u.append(matplotlib.colors.rgb2hex(rgb))

for i in range(len(col_u)):
    col_list.append(col_u[i].encode("utf-8"))

col_list = [val for val in col_list for _ in (0, 1)]
ltng_lat = {}
ltng_lon = {}
ltng_date = {}


if len(ABI_datetime) > 0:
    for n in range(0, len(fnamelist)):
        C_file = datadir + 'OR_GLM-L2-LCFA_G16_s' + str(fnamelist[n]) + '.nc'  # GOES16 East
        C = Dataset(C_file, 'r')
        ltng_lat[n] = C.variables['group_lat'][:]
        ltng_lon[n] = C.variables['group_lon'][:]
        add_seconds = C.time_coverage_start
        add_seconds = add_seconds.encode('ascii', 'ignore')
        ltng_date[n] = datetime.strptime(add_seconds, '%Y-%m-%dT%H:%M:%S.0Z')
    
    for lt in range(0, 15):
        if lt==0:
            conus_flash_count = 0
        custom = lt + len(ltng_lon) - 15
        ylt = len(ltng_lat[custom])
        conus_flash_count = conus_flash_count + ylt
    
    for lt in range(0, 15):
        subset = []
        llat=DH.latmin
        hlat=DH.latmax
        llon=DH.lonmin
        rlon=DH.lonmax
        if lt==0:
            midatl_flash_count = 0
        custom = lt + len(ltng_lon) - 15
        for v in range(0,len(ltng_lat[custom])):
            latval=ltng_lat[custom][v]
            lonval=ltng_lon[custom][v]
            if latval > llat and latval < hlat and lonval > llon and lonval < rlon:
                subset.append(latval)
        midatl_flash_count = midatl_flash_count + len(subset)
        
    abi_time = ltng_date[len(fnamelist) - 1]
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz('America/New_York')
    utc = abi_time.replace(tzinfo=from_zone)
    local = utc.astimezone(to_zone)
    
    lt = time.localtime()
    dst = lt.tm_isdst
    lt = time.localtime()
    dst = lt.tm_isdst
    if dst == 0:
        et = "EDT"
    else:
        et = "EST"
    

    
    
    rec_height = 120000
    rec_width = mH.xmax
    # Lightning plotting time 
    plt.figure(figsize=[16, 12], dpi=100)
    # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
    mH.drawcoastlines(color = 'cyan')
    mH.drawcountries(color = 'cyan')
    mH.drawmapboundary(fill_color='darkblue')
    mH.fillcontinents(color = 'black')
    mH.drawstates(color = 'cyan')
    
    symbol = u'$\u26A1$'
    #group_x, group_y = mH(group_lon, group_lat)
    for g in range(0, len(ltng_lat)):
        group_x, group_y = mH(ltng_lon[g], ltng_lat[g])
        mH.scatter(group_x, group_y, s=6, marker=symbol, c=col_list[g], zorder=3, edgecolor=col_list[g], lw=0)
    
    title = 'NOAA GOES-16 Lightning Progression'
    timestr = local.strftime('%Y-%m-%d %H:%M ') + et
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0, mH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
    plt.text(4440000, 3200000,timestr,horizontalalignment='left', color = 'black', size=14)
    plt.text(9000, 3200000,title,horizontalalignment='left', color = 'black', size=14)
    clabeltext = 'Latest 5-min Flash Count=' + str(conus_flash_count)
    currentAxis.add_patch(Rectangle((0, 0), 1000000000, rec_height * 0.8, alpha=1, zorder=3, facecolor='darkslateblue'))
    plt.text(9000, 15000,clabeltext,horizontalalignment='left', color = 'white', size=14)
    
    im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes38.png")
    plt.figimage(im1, 15, 38, zorder=1)
    output_file = workdir + "conus/" + fnamelist[-1] + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    
    ########################################################################
    ################# NOW PLOT MIDATLANTIC DOMAIN ########################
    ########################################################################
    imgdir = "/home/sat_ops/goesR/lightning/midatl/"
    rec_height = 40000
    rec_width = DH.xmax
    
    plt.figure(figsize=[8, 8], dpi=100)
    DH.drawcoastlines(color = 'cyan')
    DH.drawcountries(color = 'cyan')
    DH.drawmapboundary(fill_color='darkblue')
    DH.fillcontinents(color = 'black')
    DH.drawstates(color = 'cyan')
    DH.drawcounties(color = 'cyan')
    symbol = u'$\u26A1$'
    for g in range(0, len(ltng_lat)):
        group_x, group_y = DH(ltng_lon[g], ltng_lat[g])
        DH.scatter(group_x, group_y, s=6, marker=symbol, c=col_list[g], zorder=3, edgecolor=col_list[g], lw=0)

    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0, DH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
    plt.text(674000, 1025000,timestr,horizontalalignment='left', color = 'black', size=10)
    plt.text(7000, 1025000,title,horizontalalignment='left', color = 'black', size=10)
    clabeltext = 'Latest 5-min Flash Count=' + str(midatl_flash_count)
    currentAxis.add_patch(Rectangle((0, 0), DH.xmax, rec_height * 0.8, alpha=1, zorder=3, facecolor='darkslateblue'))
    plt.text(5000, 8000,clabeltext,horizontalalignment='left', color = 'white', size=10)
    
    im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes24.png")

    plt.figimage(im1, 15, 32,zorder=1)
    output_file = workdir + "midatl/" + fnamelist[-1] + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    

######################## Lightning Progression Gifs ########################
######################## ######################## ######################## 
imgdir = '/home/sat_ops/goesR/lightning/conus/'
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-36:]

imglen = len(img_names)
images = []
dur_vals = []
for i in xrange(1,imglen):
    if i != imglen:
        dur_vals.append(.07)
dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave(workdir + 'lightning_progression_conus.gif', images, duration=dur_vals)

# now for the midatlantic
imgdir = "/home/sat_ops/goesR/lightning/midatl/"
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-36:]

imglen = len(img_names)
images = []
dur_vals = []
for i in xrange(1,imglen):
    if i != imglen:
        dur_vals.append(.07)
dur_vals.append(2)

for i in img_names:
    input_file=imgdir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave(workdir + 'lightning_progression_midatlantic.gif', images, duration=dur_vals)
