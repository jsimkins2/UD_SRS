# script designed for basin.ceoe.udel.edu
# James Simkins

#from IPython import get_ipython
#get_ipython().run_line_magic('matplotlib', 'inline')
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



#ls |tail -7 > logfile.txt
with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C01_logfile.txt") as f:
    C01_names = f.readlines()
C01_names = [x.strip() for x in C01_names] 

with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C02_logfile.txt") as f:
    C02_names = f.readlines()
C02_names = [x.strip() for x in C02_names] 
with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C03_logfile.txt") as f:
    C03_names = f.readlines()
C03_names = [x.strip() for x in C03_names] 

for n in range(7):
    
    # C is for Conus File
    C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/' + C02_names[n]  # GOES16 East
    C = Dataset(C_file, 'r')
    # Load the RGB arrays and apply a gamma correction (square root)
    R = C.variables['CMI'][:].data # Band 2 is red (0.64 um)
    R = np.sqrt(block_mean(R, 2))
    
    C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/' + C03_names[n]  # GOES16 East
    C = Dataset(C_file, 'r')
    # Load the RGB arrays and apply a gamma correction (square root)
    G = np.sqrt(C.variables['CMI'][:].data) # Band 3 is "green" (0.865 um)
    
    C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/' + C01_names[n]  # GOES16 East
    C = Dataset(C_file, 'r')
    # Load the RGB arrays and apply a gamma correction (square root)
    B = np.sqrt(C.variables['CMI'][:].data) # Band 1 is blue (0.47 um)
    
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
    
    
    # The geostationary projection is perhaps the easiest way, and we don't need to use the Proj object.
    # Essentially, we are stretching the image across a map with the same projection and dimensions.
    
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
    plt.figure(figsize=[10, 8])
    
    # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
    newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
    newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
    
    mH.drawstates()
    mH.drawcountries()
    mH.drawcoastlines()
    
    plt.title('GOES-16 True Color\n%s' % DATE.strftime('%B %d, %Y %H:%M UTC'))
    output_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/conus/g16tc' + str(n) + '.png'
    plt.savefig(output_file, dpi=100)


import imageio
images = []
seq = [1,2,3,4,5,6,7]

for i in range(7):
    input_file='/home/sat_ops/goes_r/cloud_prod/noaa_format/conus/g16tc' + str(i) + '.png'
    images.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/cloud_prod/noaa_format/conus_goes16.gif', images, duration=.4)

'''


import calendar
import time
import datetime

min_list = [02,07,12,17,22,27,32,37,42,47,52,57]
band = [1,2,3]
b=1

now_time = datetime.datetime.now()
tt = now_time.timetuple()
jday = tt.tm_yday
# have to add a 0 in front of the single digit months
mon=tt.tm_mon
hr=tt.tm_hour
minute=tt.tm_min
minute = min(min_list, key=lambda x:abs(x-int(minute)))

if len(str(tt.tm_mon)) == 1 :
    mon="0" + str(tt.tm_mon)
if len(str(tt.tm_hour)) == 1 :
    hr="0" + str(tt.tm_hour)
if len(str(tt.tm_min)) == 1 :
    minute="0" + str(tt.tm_min)
if len(str(jday)) == 2:
    jday="0" + str(jday)
    if len(str(jday)) == 1:
        jday="0" 


fname="OR_ABI-L2-CMIPC-M3C0" + str(b) + "_G16_s" + str(tt.tm_year) + str(jday) + str(hr) + str(minute) +".nc"

# from here we can have a logfile that finds the closest value to this maf which should actually work if we read a logfile



'''




