# script designed for basin.ceoe.udel.edu
# James Simkins
# Parts of this script used code develop by Brian Blaylock - https://github.com/blaylockbk/pyBKB_v2  - Thanks Brian!
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
import matplotlib.image as image
from dateutil import tz
import time
from time import mktime

#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)


# initial block mean allows us to downsample the Red Band which is a higher resolution than the others
def block_mean(ar, fact):
    assert isinstance(fact, int), type(fact)
    sx, sy = ar.shape
    X, Y = np.ogrid[0:sx, 0:sy]
    regions = sy/fact * (X/fact) + Y/fact
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
    res.shape = (sx/fact, sy/fact)
    return res


import pyart
import boto

# extract the radar location we want
site = 'KDOX'
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]


# create a logfile with most recent 36 files (3 hours)
with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C01_logfile.txt") as f:
    C01_names = f.readlines()

C01_names = [x.strip() for x in C01_names] 

with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C02_logfile.txt") as f:
    C02_names = f.readlines()

C02_names = [x.strip() for x in C02_names] 

with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C03_logfile.txt") as f:
    C03_names = f.readlines()

C03_names = [x.strip() for x in C03_names] 

with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C13_logfile.txt") as f:
    C13_names = f.readlines()

C13_names = [x.strip() for x in C13_names] 

with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/ltng_logfile.txt") as f:
    GLM_files = f.readlines()

GLM_files = [x.strip() for x in GLM_files] 
sname1 = []
sname2 = []
sname3 = []
sname4 = []
GLM_names = []
# we're going to set C02 as the default here because it's the most important band
for i in xrange(1,len(C01_names)):
    fname = str(C01_names[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    sname1.append(fname)

for i in xrange(1,len(C02_names)):
    fname = str(C02_names[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    sname2.append(fname)

for i in xrange(1,len(C03_names)):
    fname = str(C03_names[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    sname3.append(fname)

for i in xrange(1,len(C13_names)):
    fname = str(C13_names[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    sname4.append(fname)

for i in xrange(1,len(GLM_files)):
    fname = str(GLM_files[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    GLM_names.append(fname)
    
# we only want to do timestamps where we have all 3 bands
from collections import OrderedDict
match = set(sname1) & set(sname2) & set(sname3) & set(sname4)
match = sorted(match, key=int)

ABI_datetime = []
for i in match:    
    if os.path.isfile("/home/sat_ops/goes_r/lightning/add_clouds/im_tc_con/" + str(i) + ".png") == False:
        ABI_datetime.append(i)

# convert our lightning file list and our ABI file list to integers so we can match the lightning files with the ABI files
GLM_int = []
ABI_int = []
for st in xrange(0, len(ABI_datetime)):
    ABI_int.append(int(ABI_datetime[st] + '000'))

for st in xrange(0, len(GLM_names)):
    GLM_int.append(int(GLM_names[st]))


def between(l1,low,high):
    l2 = []
    for i in l1:
        if(i >= low and i < high):
            l2.append(i)
    if len(l2) > 15:
        l2=l2[-15:]
    return l2

def below(l1, val):
    l2 = []
    for i in l1:
       if(i <= val):
            l2.append(i)
    l2=l2[-15:]
    return l2

def above(l1, val):
    l2 = []
    for i in l1:
       if(i <= val):
            l2.append(i)
    l2=l2[-15:]
    return l2
    
# set up a dictionary for the lightning files that match up with each ABI file
abi_glm = {}

# ok, the way this works is we need to grab all GLMs within 5 minutes of the ABI file
# the thought here is place the GLM files in a dictionary corresponding with the ABI file (0,1,2 etc.)
# if there is only 1 ABI_int, then all we need to do is find the GLM ints for that one file, hence 
# we only use the below function
# however, if there is more than one file, reset abi_glm so we have a blank dict, then find the last 5 minutes
# of GLM for each ABI datetime file. 
# Then, we need to do the same for the last one since there won't be a between for that file because there isn't
# a higher value to use the between function. Therefore, we use the below function
print GLM_int
if ABI_int > 0:
    # if there is only 1 ABI_datetime file then just fill the dictionary for that time slot
    abi_glm[len(ABI_int) -1] = below(GLM_int, ABI_int[len(ABI_int) - 1])
    print len(ABI_int)-1
    if ABI_int > 1:
        abi_glm = {}
        for a in xrange(1, len(ABI_int)):
            if a == 1:
                abi_glm[0] = below(GLM_int, ABI_int[0])
            if a > 1:
                abi_glm[a - 1] = between(GLM_int,ABI_int[a-2] , ABI_int[a-1])
        #print len(ABI_int)-1
        #abi_glm[len(ABI_int) -1] = below(GLM_int, ABI_int[len(ABI_int) - 1])

# Make a new map object for the HRRR model domain map projection
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)

DH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-4,llcrnrlon=lon0-5,
    urcrnrlat=lat0+4.5,urcrnrlon=lon0+5,resolution='h')   

# begin the loop that makes the images
if len(abi_glm) > 0:
    for n in xrange(0, len(abi_glm)):
        print ABI_int[n]
        # Below code has been patched together with code created by Brian Blaylock - https://github.com/blaylockbk/pyBKB_v2
        # C is for Conus File OR_ABI-L2-CMIPC-M3C02_G16_s20180601912.nc
        # RED BAND
        C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M3C02_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
        if os.path.isfile(C_file) == False:
            C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M4C02_G16_s' + str(ABI_datetime[n]) + '.nc'

        C = Dataset(C_file, 'r')
        # Load the RGB arrays and apply a gamma correction (square root)
        R = C.variables['CMI'][:].data # Band 2 is red (0.64 um)
        R[R < 0] = np.nan
        R = np.sqrt(block_mean(R, 2)) # Red band is twice the res of B and G, and 4 times that of IR bands
        R = block_mean(R, 2)
        C.close()
        
        # GREEN BAND
        C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M3C03_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
        if os.path.isfile(C_file) == False:
            C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M4C03_G16_s' + str(ABI_datetime[n]) + '.nc'
            if os.path.isfile(C_file) == False:
                print "No file exists in data folder"
        C = Dataset(C_file, 'r')
        # Load the RGB arrays and apply a gamma correction (square root)
        G = C.variables['CMI'][:].data
        G[G < 0] = np.nan
        # Band 3 is "green" (0.865 um)
        G = np.sqrt(block_mean(G, 2))
        C.close()
        
        # BLUE BAND 
        C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M3C01_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
        if os.path.isfile(C_file) == False:
            C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M4C01_G16_s' + str(ABI_datetime[n]) + '.nc'
            if os.path.isfile(C_file) == False:
                print "No file exists in data folder"
        C = Dataset(C_file, 'r')
        # Load the RGB arrays and apply a gamma correction (square root)
        B = C.variables['CMI'][:].data # Band 1 is blue (0.47 um)
        B[B < 0] = np.nan
        B = np.sqrt(block_mean(B, 2))
        C.close()
        
        # LW BAND
        C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M3C13_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
        if os.path.isfile(C_file) == False:
            C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M4C13_G16_s' + str(ABI_datetime[n]) + '.nc'
            if os.path.isfile(C_file) == False:
                print "No file exists in data folder"
        
        C = Dataset(C_file, 'r')
        # Load the RGB arrays and apply a gamma correction (square root)
        b13 = C.variables['CMI'][:] # Band 1 is blue (0.47 um)
        
        
        # "True Green" is some linear interpolation between the three channels
        G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G
        # The final RGB array :)
        RGB = np.dstack([R, G_true, B])
        add_seconds = C.variables['t'][0]
        DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=add_seconds)
        # Satellite height
        sat_h = C.variables['goes_imager_projection'].perspective_point_height
        # Satellite longitude
        sat_lon = C.variables['goes_imager_projection'].longitude_of_projection_origin
        # Satellite sweep
        sat_sweep = C.variables['goes_imager_projection'].sweep_angle_axis
        # The projection x and y coordinates equals
        # the scanning angle (in radians) multiplied by the satellite height (http://proj4.org/projections/geos.html)
        X = C.variables['x'][:] * sat_h
        Y = C.variables['y'][:] * sat_h
        
        ########## LIGHTNING ##############
        # Now that we have all the satellite ABIs read in, we need to call the lightning files 
        # create dictionaries for the lat/lons
        ltng_lat = {}
        ltng_lon = {}
        print abi_glm
        for lt in xrange(0, len(abi_glm[n])):
            print lt
        for lt in xrange(0, len(abi_glm[n])):
            ltfile = abi_glm[n][lt]
            L_file = '/home/sat_ops/goes_r/lightning/polished_data/OR_GLM-L2-LCFA_G16_s' + str(ltfile) + '.nc'  # GOES16 East
            L = Dataset(L_file, 'r')
            ltng_lat[lt] = L.variables['group_lat'][:]
            ltng_lon[lt] = L.variables['group_lon'][:]

        
        
        # The geostationary projection is perhaps the easiest way, and we don't need to use the Proj object.
        # Essentially, we are stretching the image across a map with the same projection and dimensions.
        
        # map object with pyproj
        p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
        # Convert map points to latitude and longitude with the magic provided by Pyproj
        XX, YY = np.meshgrid(X, Y)
        lons, lats = p(XX, YY, inverse=True)
        xH, yH = mH(lons, lats)
        
        # Create a color tuple for pcolormesh
        rgb = RGB[:,:-1,:] # Using one less column is very imporant, else your image will be scrambled! (This is the stange nature of pcolormesh)
        rgb = np.minimum(rgb, 1) # Force the maximum possible RGB value to be 1 (the lowest should be 0).
        colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]), 3) # flatten array, becuase that's what pcolormesh wants.
        colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster?? according to stackoverflow.
        
        # adding this additional line here to see if it helps
        # for some reason this helps
        colorTuple[colorTuple < 0] = 0
        colorTuple[colorTuple > 1] = 1
        # change the alphas to show the b13 below the tc image
        # colorTuple[:,3][colorTuple[:,2] < 0.3] = 0.3
        colorTuple[:,3][colorTuple[:,2] < 0.2] = 0.0
        colorTuple[:,3][colorTuple[:,2] < 0.1] = 0.0

        # Now we can plot the GOES data on the HRRR map domain and projection
        # first, plot the b13 stuff
        plt.figure(figsize=[16, 12], dpi=200)
        m = mH.pcolormesh(xH, yH, b13, cmap='Greys', vmax=280, vmin=180)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
        newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
        mH.drawcoastlines(color = 'cyan')
        mH.drawcountries(color = 'cyan')
        mH.drawstates(color = 'cyan')
        #map.draw
        symbol = u'$\u26A1$'
        #symbol = u'$\u2193$'
        # Plot flashes as small red dots
        #group_x, group_y = mH(group_lon, group_lat)
        for g in xrange(0, len(ltng_lat)):
            group_x, group_y = mH(ltng_lon[g], ltng_lat[g])
            mH.scatter(group_x, group_y, s=2, marker=symbol, c='yellow', zorder=3, edgecolor='yellow', lw=0)
        
        
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
        output_file = '/home/sat_ops/goes_r/lightning/add_clouds/im_tc_con/' + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
    
    
        # Now Plot Delaware

        plt.figure(figsize=[8, 8], dpi=100)
        D = DH.pcolormesh(xH, yH, b13, cmap='Greys', vmax=280, vmin=180)
        xH, yH = DH(lons, lats)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        newmap = DH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
        newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
        DH.drawcoastlines(linewidth=0.7, color = 'cyan')
        DH.drawcountries(linewidth=0.7, color = 'cyan')
        DH.drawstates(linewidth=0.7, color = 'cyan')
        #map.draw
        symbol = u'$\u26A1$'
        #symbol = u'$\u2193$'
        # Plot flashes as small red dots
        #group_x, group_y = mH(group_lon, group_lat)
        for g in xrange(0, len(abi_glm[n])):
            group_x, group_y = DH(ltng_lon[g], ltng_lat[g])
            DH.scatter(group_x, group_y, s=8, marker=symbol, c='yellow', zorder=3, edgecolor='yellow', lw=0)
        
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
        plt.figimage(im1, 560, 633, zorder=1)
        plt.figimage(im2, 13, 633, zorder=1)
        # save file
        output_file = "/home/sat_ops/goes_r/lightning/add_clouds/im_tc_mid/" + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
else:
    print "Up to Date"
    
    
