### James Simkins

## View NEXRAD Level II Data

## Download from Amazon Web Service
## https://s3.amazonaws.com/noaa-nexrad-level2/index.html

# Download URL Example: https://noaa-nexrad-level2.s3.amazonaws.com/2016/08/05/KCBX/KCBX20160805_205859_V06

# working in pyg27 environment, where AWIPS isn't fucking anything up
import os
import numpy as np
from numpy import ma
from datetime import datetime
import matplotlib.pyplot as plt
import sys
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from pyproj import Geod #this is used to convert range/azimuth to lat/lon
from mpl_toolkits.basemap import Basemap

from metpy.io import Level2File
from metpy.plots import add_metpy_logo, ctables

SAVEDIR = '/Users/leathers/Documents/Delaware/goesR/pyBKB/'
import pyart

# have to download this data with aws via these commands
# aws s3 ls --recursive s3://noaa-nexrad-level2/2018/02/12/KDOX - finds what files are in the folder
# aws s3 cp s3://noaa-nexrad-level2/2018/02/12/KDOX/KDOX20180212_171032_V06 /. - copies the file we want to cd

filename = "KDOX20180212_171032_V06"
f = Level2File("Documents/Delaware/goesR/nexraddata/KDOX20180212_171032_V06")

rLAT = f.sweeps[0][0][1].lat
rLON = f.sweeps[0][0][1].lon
lats = f.sweeps[:][:][1].lat
rDT = f.dt # this is in local time



fig = plt.figure(1,figsize=[8,8])
ax  = fig.add_subplot(111)
top_right_lat = rLAT+1.25
top_right_lon = rLON+1.25
bot_left_lat = rLAT-1
bot_left_lon = rLON-1

mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=bot_left_lon, lat_2=bot_left_lat, \
            lat_0=top_right_lon, lon_0=top_right_lat)
            
            
%matplotlib inline
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
from pyproj import Proj     
from scipy import ndimage
from scipy import stats

# initial block mean allows us to downsample the Red Band which is a higher resolution than the others
def block_mean(ar, fact):
    assert isinstance(fact, int), type(fact)
    sx, sy = ar.shape
    X, Y = np.ogrid[0:sx, 0:sy]
    regions = sy/fact * (X/fact) + Y/fact
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
    res.shape = (sx/fact, sy/fact)
    return res





# C is for Conus File
C_file = 'Documents/Delaware/goesR/satdata/OR_ABI-L2-CMIPC-M3C02_G16_s20180431702207_e20180431704580_c20180431705095.nc'  # GOES16 East
C = Dataset(C_file, 'r')
# Load the RGB arrays and apply a gamma correction (square root)
R = C.variables['CMI'][:].data # Band 2 is red (0.64 um)
R = block_mean(R, 2)

C_file = 'Documents/Delaware/goesR/satdata/OR_ABI-L2-CMIPC-M3C03_G16_s20180431702207_e20180431704580_c20180431705056.nc'  # GOES16 East
C = Dataset(C_file, 'r')
# Load the RGB arrays and apply a gamma correction (square root)
G = np.sqrt(C.variables['CMI'][:].data) # Band 3 is "green" (0.865 um)

C_file = 'Documents/Delaware/goesR/satdata/OR_ABI-L2-CMIPC-M3C01_G16_s20180431702207_e20180431704580_c20180431705058.nc'  # GOES16 East
C = Dataset(C_file, 'r')
# Load the RGB arrays and apply a gamma correction (square root)
B = np.sqrt(C.variables['CMI'][:].data) # Band 1 is blue (0.47 um)

# "True Green" is some linear interpolation between the three channels
G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G

# The final RGB array :)
RGB = np.dstack([R, G_true, B])

add_seconds = C.variables['t'][0]
DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=add_seconds)


plt.figure(figsize=[10, 8])
plt.imshow(RGB)
plt.title(DATE)



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


# The geostationary projection is perhaps the easiest way, and we don't need to use the Proj object.
# Essentially, we are stretching the image across a map with the same projection and dimensions.
plt.figure(figsize=[10, 8])
m = Basemap(projection='geos', lon_0=sat_lon,
            resolution='i', area_thresh=1000,
            llcrnrx=X.min(),llcrnry=Y.min(),
            urcrnrx=X.max(),urcrnry=Y.max())
m.imshow(np.flipud(RGB)) # Remember, "images" are upside down, so flip up/down
m.drawcoastlines()
m.drawcountries()
m.drawstates()
plt.title(DATE)


# map object with pyproj
p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
# Convert map points to latitude and longitude with the magic provided by Pyproj
XX, YY = np.meshgrid(X, Y)
lons, lats = p(XX, YY, inverse=True)

# Make a new map object for the HRRR model domain map projection
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=bot_left_lon, lat_2=bot_left_lat, \
            lat_0=top_right_lon, lon_0=top_right_lat)

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
plt.figure(figsize=[20, 16])

# The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.

mH.drawstates()
mH.drawcountries()
mH.drawcoastlines()

plt.title('GOES-16 True Color\n%s' % DATE.strftime('%B %d, %Y %H:%M UTC'))





