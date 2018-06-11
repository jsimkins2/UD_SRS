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
# Make a new map object for the HRRR model domain map projection
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)
            
            
for n in xrange(0, len(fnamelist)):
    print n
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
    group_lat = C.variables['group_lat'][:]
    group_lon = C.variables['group_lon'][:]
    
    #flash_lat = C.variables['flash_lat'][:]
    #flash_lon = C.variables['flash_lon'][:]
    
    add_seconds = C.time_coverage_start
    add_seconds = add_seconds.encode('ascii', 'ignore')
    DATE = datetime.strptime(add_seconds, '%Y-%m-%dT%H:%M:%S.0Z')

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
    group_x, group_y = mH(group_lon, group_lat)
    mH.scatter(group_x, group_y, s=24, marker=symbol, c="yellow", zorder=3, edgecolor="yellow", lw=.2)
    
    
    # METHOD 1: Hardcode zones:
    abi_time = DATE
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
    output_file = '/home/sat_ops/goes_r/lightning/tem_images/' + str(n) + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    C.close()
    
import imageio
import numpy as np
images = []
dur_vals = []
for i in xrange(1,len(fnamelist) - 1):
    dur_vals.append(.1)
    
dur_vals.append(2)
#print dur_vals
new_seq = range(1,len(fnamelist))
from collections import OrderedDict
#new_seq = sorted(new_seq, key=int, reverse=True)

for i in new_seq:
    input_file='/home/sat_ops/goes_r/lightning/tem_images/' + str(i) + '.png'
    images.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/lightning/lightning_conus.gif', images, duration=dur_vals)

