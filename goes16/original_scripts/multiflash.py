# script designed for basin.ceoe.udel.edu
# James Simkins
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
from pyproj import Proj     
from scipy import ndimage
from scipy import stats
import os.path
import matplotlib.image as image
from dateutil import tz
import time
from time import mktime
from collections import OrderedDict
from matplotlib.colors import LinearSegmentedColormap, ListedColormap # Linear interpolation for color maps
import pyart
import boto

# extract the radar location we want
site = 'KDOX'
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]

DH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-4,llcrnrlon=lon0-5,
    urcrnrlat=lat0+4.5,urcrnrlon=lon0+5,resolution='h')   

with open("/home/sat_ops/goes_r/lightning/logfile.txt") as f:
    file_names = f.readlines()

file_names = [x.strip() for x in file_names] 
fnamelist=[]

for i in xrange(1,len(file_names)):
    fname = str(file_names[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    fnamelist.append(fname)
    
fnamelist = sorted(fnamelist, key=int)

# this isn't being used right now but hopefully we can figure this out later on
ABI_datetime = []
for i in fnamelist:    
    if os.path.isfile("/home/sat_ops/goes_r/lightning/multi_con/" + str(i) + ".png") == False:
        ABI_datetime.append(i)

# later on we want to plot all files leading to the file 
ABI_datetime.append(9999)
# Make a new map object for the HRRR model domain map projection
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)

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
for elem in fnamelist:
    ltng_lat[elem] = []
    ltng_lon[elem] = []
    ltng_date[elem] = []

if len(ABI_datetime) > 0:
    for n in xrange(0, len(fnamelist)):
        # C is for Conus File OR_ABI-L2-CMIPC-M3C02_G16_s20180601912.nc
        # RED BAND
        C_file = '/home/sat_ops/goes_r/lightning/polished_data/OR_GLM-L2-LCFA_G16_s' + str(fnamelist[n]) + '.nc'  # GOES16 East
        if os.path.isfile(C_file) == False:
            import smtplib
            dfile=open("/home/sat_ops/goes_r/cloud_prod/noaa_format/udelsatinfo.txt", "r")
            dw = dfile.readline()
            server = smtplib.SMTP('smtp.gmail.com', 587)
            server.starttls()
            server.login("goessatelliteudel@gmail.com", dw)
            msg = "TC CONUS IS BREAKING"
            server.sendmail("goessatelliteudel@gmail.com", "simkins@udel.edu", msg)
            server.quit()
        C = Dataset(C_file, 'r')
        # Load the RGB arrays and apply a gamma correction (square root)
        #event_lat = C.variables['event_lat'][:]
        #event_lon = C.variables['event_lon'][:]
        ltng_lat[n] = C.variables['group_lat'][:]
        ltng_lon[n] = C.variables['group_lon'][:]
        
        #flash_lat = C.variables['flash_lat'][:]
        #flash_lon = C.variables['flash_lon'][:]
        
        add_seconds = C.time_coverage_start
        add_seconds = add_seconds.encode('ascii', 'ignore')
        ltng_date[n] = datetime.strptime(add_seconds, '%Y-%m-%dT%H:%M:%S.0Z')
        # end the Abi_datetime plotting for that file
    
    # Now Plot
    plt.figure(figsize=[16, 12], dpi=200)
    mH.drawcoastlines(color = 'cyan')
    mH.drawcountries(color = 'cyan')
    mH.drawmapboundary(fill_color='midnightblue')
    mH.fillcontinents(color = 'black')
    mH.drawstates(color = 'cyan')
    #map.draw
    symbol = u'$\u26A1$'
    #symbol = u'$\u2193$'
    # Plot flashes as small red dots
    #group_x, group_y = mH(group_lon, group_lat)
    for n in xrange(0, len(fnamelist)):
        group_x, group_y = mH(ltng_lon[n], ltng_lat[n])
        mH.scatter(group_x, group_y, s=6, marker=symbol, c=col_list[n], zorder=3, edgecolor=col_list[n], lw=0)
    
    
    # METHOD 1: Hardcode zones:
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
    
    plt.title('NOAA GOES-16 Lightning Mapper\n%s' % local.strftime('%B %d, %Y %H:%M ') + et)
    # add logo
    im1 = image.imread("/home/sat_ops/goes_r/nexrad/cema38.png")
    im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
    #im[:, :, -1] = 0.5
    plt.figimage(im1, 1210, 745, zorder=1)
    plt.figimage(im2, 15, 745, zorder=1)
    # save file
    output_file = '/home/sat_ops/goes_r/lightning/multi_con/' + fnamelist[n] + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()


    # Now Plot Delaware
    ax = plt.figure(figsize=[8, 8], dpi=100)
    DH.drawcoastlines(linewidth=0.7, color = 'cyan')
    DH.drawcountries(linewidth=0.7, color = 'cyan')
    DH.drawmapboundary(fill_color='midnightblue')
    DH.fillcontinents(color = 'black')
    DH.drawstates(linewidth=0.7, color = 'cyan')
    #map.draw
    symbol = u'$\u26A1$'
    #symbol = u'$\u2193$'
    # Plot flashes as small red dots
    #group_x, group_y = mH(group_lon, group_lat)
    for n in xrange(0, len(fnamelist)):
        group_x, group_y = DH(ltng_lon[n], ltng_lat[n])
        DH.scatter(group_x, group_y, s=6, marker=symbol, c=col_list[n], zorder=3, edgecolor=col_list[n], lw=0)
    
    
    # METHOD 1: Hardcode zones:
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
    
    plt.title('NOAA GOES-16 Lightning Mapper\n%s' % local.strftime('%B %d, %Y %H:%M ') + et)
    # add logo
    im1 = image.imread("/home/sat_ops/goes_r/nexrad/cema38.png")
    im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
    #im[:, :, -1] = 0.5
    plt.figimage(im1, 560, 633, zorder=1)
    plt.figimage(im2, 13, 633, zorder=1)
    # save file
    output_file = '/home/sat_ops/goes_r/lightning/multi_mid/' + fnamelist[n] + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
