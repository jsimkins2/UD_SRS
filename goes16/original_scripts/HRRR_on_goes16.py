#pyKBK github guy 

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
C_file = 'Downloads/OR_ABI-L2-CMIPC-M3C02_G16_s20180381702190_e20180381704563_c20180381705075.nc'  # GOES16 East
C = Dataset(C_file, 'r')
# Load the RGB arrays and apply a gamma correction (square root)
R = C.variables['CMI'][:].data # Band 2 is red (0.64 um)
R = block_mean(R, 2)

C_file = 'Downloads/OR_ABI-L2-CMIPC-M3C03_G16_s20180381702190_e20180381704563_c20180381705033.nc'  # GOES16 East
C = Dataset(C_file, 'r')
# Load the RGB arrays and apply a gamma correction (square root)
G = np.sqrt(C.variables['CMI'][:].data) # Band 3 is "green" (0.865 um)

C_file = 'Downloads/OR_ABI-L2-CMIPC-M3C01_G16_s20180381702190_e20180381704563_c20180381705038.nc'  # GOES16 East
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
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)

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






##### HRRR stuff now 
import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
import numpy as np
grib='Downloads/hrrr.t17z.wrfsfcf00.grib2' # Set the file name of your input GRIB file
grbs=pygrib.open(grib)

Hvalues, Hlat, Hlon = grbs[1].data()
validDATE = grbs[1].validDate
anlysDATE = grbs[1].analDate
msg = str(grbs[1])

# Load some pre-downloaded HRRR 500 mb heights data that I have in a numpy dictionary.
grib2='Downloads/hrrr.t17z.wrfprsf00.grib2' # Set the file name of your input GRIB file
gr =pygrib.open(grib2)
H500 = grbs[11] #get record number 38
H1000 = grbs[38]


H500values = H500.values
H500lat,H500lon = H500.latlons()
H1000values = H1000.values
H1000values = H500values - H1000values

# Mask points with no reflectivity
dBZ = Hvalues
dBZ = np.ma.array(dBZ)
dBZ[dBZ == -10] = np.ma.masked

plt.figure(figsize=[20,16])
# Plot the GOES image
newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
newmap.set_array(None) # without this, the linewidth is set to zero, but the RGB color is ignored

# Plot the HRRR reflectivity
mH.pcolormesh(Hlon, Hlat, dBZ, latlon=True, cmap='gist_ncar', vmax=80, vmin=0)
cb = plt.colorbar(orientation='horizontal', shrink=.7, pad=.01)
cb.set_label('Simulated Radar Reflectivity (dBZ)')

# Plot the HRRR 500 mb height
cs = mH.contour(H500lon, H500lat, H1000values, latlon=True, colors='k', levels=range(5020, 5860, 60))
cs2 = mH.contour(H500lon, H500lat, H1000values, latlon = True, colors='magenta', linestyles="dotted",linewidths=4,levels=[5400])
plt.clabel(cs, fmt = '%1.0f')
plt.clabel(cs2, fmt = '%1.0f')


# Plot other map elements
mH.drawstates()
mH.drawcountries()
mH.drawcoastlines()

# Title
date_fmt = '%B %d, %Y %H:%M'
plt.title('GOES-16 True Color, HRRR Reflectivity, HRRR 1000-500 mb Thickness\nGOES date: %s\n HRRR date: %s' % (DATE.strftime(date_fmt), validDATE.strftime(date_fmt)))
plt.savefig("Documents/deos_projects/thickness_goes16_HRRR_0207.png")




