site = 'KDOX'
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
import pyart
import boto
import os
import tempfile
import matplotlib as mpl
from dateutil import tz
import time
from time import mktime
import matplotlib.image as image

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



# open the logfiles for each ABI of the most recent 6 hours of data
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

sname1 = []
sname2 = []
sname3 = []
sname4 = []
# parse the name strings so we can organize them better
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

# we only want to do timestamps where we have all 3 bands
from collections import OrderedDict
match = set(sname1) & set(sname2) & set(sname3) & set(sname4)
match = sorted(match, key=int)

# if the file already exists, don't copy it on the to do list
ABI_datetime = []
for i in match:    
    if os.path.isfile("/home/sat_ops/goes_r/cloud_prod/noaa_format/image_midatlantic/" + str(i) + ".png") == False:
        ABI_datetime.append(i)

print ABI_datetime

#get the radar location (this is used to set up the basemap and plotting grid)

loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]

# before we go in the loop for the images, let's do some stuff now before we go into the loop
lt = time.localtime()
dst = lt.tm_isdst
lt = time.localtime()
dst = lt.tm_isdst
if dst == 0:
    et = "EDT"
else:
    et = "EST"

# Make a new map object for the HRRR model domain map projection
mH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-3,llcrnrlon=lon0-4,
    urcrnrlat=lat0+3.5,urcrnrlon=lon0+4,resolution='h')

# begin the loop that makes the images
for n in xrange(0, len(ABI_datetime)):
    
    print n
    # C is for Conus File OR_ABI-L2-CMIPC-M3C02_G16_s20180601912.nc
    # C is for Conus File OR_ABI-L2-CMIPC-M3C02_G16_s20180601912.nc
    C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M3C02_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
    if os.path.isfile(C_file) == False:
        C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M4C02_G16_s' + str(ABI_datetime[n]) + '.nc'
        if os.path.isfile(C_file) == False:
                import smtplib
                dfile=open("/home/sat_ops/goes_r/cloud_prod/noaa_format/udelsatinfo.txt", "r")
                dw = dfile.readline()
                server = smtplib.SMTP('smtp.gmail.com', 587)
                server.starttls()
                server.login("goessatelliteudel@gmail.com", dw)
                msg = "MidAtlantic IS BREAKING"
                server.sendmail("goessatelliteudel@gmail.com", "simkins@udel.edu", msg)
                server.quit()

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
    
    # Now we can plot the GOES data on the HRRR map domain and projection
    # change the alphas to show the b13 below the tc image
    colorTuple[:,3][colorTuple[:,2] < 0.2] = 0.2
    colorTuple[:,3][colorTuple[:,2] < 0.15] = 0.1
    colorTuple[:,3][colorTuple[:,2] < 0.1] = 0.0
    
    # Now we can plot the GOES data on the HRRR map domain and projection
    plt.figure(figsize=[8, 8])
    m = mH.pcolormesh(xH, yH, b13, cmap='Greys', vmax=280, vmin=180)
    
    # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
    newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
    newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
    mH.drawstates()
    mH.drawcoastlines(linewidth=0.7,color='k')
    mH.drawcounties(linewidth=0.1,color='k')

    # convert from UTC to local time
    abi_time = DATE
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz('America/New_York')
    utc = abi_time.replace(tzinfo=from_zone)
    local = utc.astimezone(to_zone)
    # add logo
    im1 = image.imread("/home/sat_ops/goes_r/nexrad/cema38.png")
    im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
    #im[:, :, -1] = 0.5
    plt.figimage(im1, 680, 760, zorder=1)
    plt.figimage(im2, 13, 760, zorder=1)
    # save file
    plt.title('NOAA GOES-16\n%s' % local.strftime('%B %d, %Y %H:%M ') + et)
    output_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/image_midatlantic/' + str(ABI_datetime[n]) + ".png"
    plt.savefig(output_file, dpi=120, bbox_inches='tight')
    plt.close()
    C.close()
