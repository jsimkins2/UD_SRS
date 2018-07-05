# This script creates true color imagery for CONUS and Mid Atlantic & Lightning/TC over CONUS and Mid Atlantic
# 4 different images altogether 
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
import matplotlib.image as image
from dateutil import tz
import time
from time import mktime
from collections import OrderedDict
import pyart
import boto
#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)


############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/truecolor/"
datadir = "/home/sat_ops/goesR/data/mcmipc/"
site = 'KDOX'

############## Define Functions that will be used later #########################
# initial block mean allows us to downsample the Red Band which is a higher resolution than the others
def block_mean(ar, fact):
    assert isinstance(fact, int), type(fact)
    sx, sy = ar.shape
    X, Y = np.ogrid[0:sx, 0:sy]
    regions = sy/fact * (X/fact) + Y/fact
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
    res.shape = (sx/fact, sy/fact)
    return res

def contrast_correction(color, contrast):
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.minimum(COLOR, 1)
    COLOR = np.maximum(COLOR, 0)
    return COLOR

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

################ Grab the Lat/Lon of the site we want ####################
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]


# create a logfile with most recent files created in the shell script
with open(workdir + "logfile.txt") as f:
    file_names = f.readlines()
file_names = [x.strip() for x in file_names] 

fnamelist = []
for i in xrange(1,len(file_names)):
    fname = str(file_names[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    fnamelist.append(fname)

with open("/home/sat_ops/goesR/truecolor/ltng_logfile.txt") as f:
    GLM_files = f.readlines()

GLM_files = [x.strip() for x in GLM_files] 
GLM_names = []

for i in xrange(1,len(GLM_files)):
    fname = str(GLM_files[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    GLM_names.append(fname)
    
# sort through the files and make sure that they don't exist before we process them
fnamelist = sorted(fnamelist, key=int)
ABI_datetime = []
for i in fnamelist:    
    if os.path.isfile(workdir + "image_conus/" + str(i) + ".png") == False:
        ABI_datetime.append(i)

# Sort the lightning files and grab the corresponding lightning files with each ABI scan 
GLM_int = []
ABI_int = []
abi_glm = {}

for st in xrange(0, len(ABI_datetime)):
    ABI_int.append(int(ABI_datetime[st]))

for st in xrange(0, len(GLM_names)):
    GLM_int.append(int(GLM_names[st]))

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

# Make a new map object for the HRRR model domain map projection
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)

DH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-4,llcrnrlon=lon0-5,
    urcrnrlat=lat0+4.5,urcrnrlon=lon0+5,resolution='h') 
    

# begin the loop that makes the images
if len(ABI_datetime) > 0:
    for n in xrange(0, len(ABI_datetime)):
        print ABI_datetime[n]
        # C is for Conus File OR_ABI-L2-CMIPC-M3C02_G16_s20180601912.nc
        # RED BAND
        C_file = datadir + 'OR_ABI-L2-MCMIPC-M3_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
        if os.path.isfile(C_file) == False:
            C_file = '/home/sat_ops/goes_r/night_scan/polished_data/OR_ABI-L2-MCMIPC-M4_G16_s' + str(ABI_datetime[n]) + '.nc'
            if os.path.isfile(C_file) == False:
                import smtplib
                dfile=open("/home/sat_ops/goes_r/noaa_format/udelsatinfo.txt", "r")
                dw = dfile.readline()
                server = smtplib.SMTP('smtp.gmail.com', 587)
                server.starttls()
                server.login("goessatelliteudel@gmail.com", dw)
                msg = "TC CONUS IS BREAKING"
                server.sendmail("goessatelliteudel@gmail.com", "simkins@udel.edu", msg)
                server.quit()
        
        # below is code from Brian Blaylock of Univ Utah, thanks Brian!
        Cnight = Dataset(C_file, 'r')
        
        # Load the RGB arrays
        R = Cnight.variables['CMI_C02'][:].data
        G = Cnight.variables['CMI_C03'][:].data
        B = Cnight.variables['CMI_C01'][:].data
        
        # Turn empty values into nans
        R[R==-1] = np.nan
        G[G==-1] = np.nan
        B[B==-1] = np.nan
        
        # Apply range limits for each channel becuase RGB values must be between 0 and 1
        R = np.maximum(R, 0)
        R = np.minimum(R, 1)
        G = np.maximum(G, 0)
        G = np.minimum(G, 1)
        B = np.maximum(B, 0)
        B = np.minimum(B, 1)
        
        # Apply the gamma correction
        gamma = 0.4
        R = np.power(R, gamma)
        G = np.power(G, gamma)
        B = np.power(B, gamma)
        
        # Calculate the "True" Green
        G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G
        G_true = np.maximum(G_true, 0)
        G_true = np.minimum(G_true, 1)
        
        # Grab the IR / Apply range limits for clean IR channel/Normalize the channel between a range/invert colors/lessen the brightness
        cleanIR = Cnight.variables['CMI_C13'][:].data
        cleanIR[cleanIR==-1] = np.nan
        cleanIR = np.maximum(cleanIR, 90)
        cleanIR = np.minimum(cleanIR, 313)
        cleanIR = (cleanIR-90)/(313-90)
        cleanIR = 1 - cleanIR
        cleanIR = cleanIR/1.5
        
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
            ltng_lat[lt] = L.variables['flash_lat'][:]
            ltng_lon[lt] = L.variables['flash_lon'][:]


        # Modify the RGB color contrast
        contrast = 125
        RGB_contrast = contrast_correction(np.dstack([R, G_true, B]), contrast)
        RGB_contrast_IR = np.dstack([np.maximum(RGB_contrast[:,:,0], cleanIR), np.maximum(RGB_contrast[:,:,1], cleanIR), np.maximum(RGB_contrast[:,:,2], cleanIR)])
        
        # Satellite Date, height, lon, sweep
        add_seconds = Cnight.variables['t'][0]
        DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=add_seconds)
        sat_h = Cnight.variables['goes_imager_projection'].perspective_point_height
        sat_lon = Cnight.variables['goes_imager_projection'].longitude_of_projection_origin
        sat_sweep = Cnight.variables['goes_imager_projection'].sweep_angle_axis
        
        # Configure EST/EDT depending on time of year
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
        
        # The projection x and y coordinates equals
        # the scanning angle (in radians) multiplied by the satellite height (http://proj4.org/projections/geos.html)
        X = Cnight.variables['x'][:] * sat_h
        Y = Cnight.variables['y'][:] * sat_h
        
        # map object with pyproj
        p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
        # Convert map points to latitude and longitude with the magic provided by Pyproj
        XX, YY = np.meshgrid(X, Y)
        lons, lats = p(XX, YY, inverse=True)
        lats[np.isnan(R)] = np.nan
        lons[np.isnan(R)] = np.nan
        xH, yH = mH(lons, lats)
        
        # Create a color tuple for pcolormesh
        rgb = RGB_contrast_IR[:,:-1,:] # Using one less column is very imporant, else your image will be scrambled! (This is the stange nature of pcolormesh)
        rgb = np.minimum(rgb, 1)
        colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]),3) # flatten array, becuase that's what pcolormesh wants.
        colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster?? according to stackoverflow.
        
        # Now we can plot the GOES data on the HRRR map domain and projection
        plt.figure(figsize=[16, 12], dpi=100)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
        newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
        mH.drawstates()
        mH.drawcountries()
        mH.drawcoastlines(linewidth=0.7,color='k')
        
        # add logo
        plt.title('NOAA GOES-16\n%s' % local.strftime('%B %d, %Y %H:%M ') + et)
        im1 = image.imread("/home/sat_ops/goes_r/nexrad/cema38.png")
        im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
        plt.figimage(im1, 1210, 745, zorder=1)
        plt.figimage(im2, 15, 745, zorder=1)
        
        # save file
        output_file = workdir + "img_conus/" + str(ABI_datetime[n]) + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        
        symbol = u'$\u26A1$'
        #group_x, group_y = mH(group_lon, group_lat)
        for g in xrange(0, len(abi_glm[n])):
            group_x, group_y = mH(ltng_lon[g], ltng_lat[g])
            mH.scatter(group_x, group_y, s=14, marker=symbol, c='yellow', zorder=3, edgecolor='yellow', lw=0)
        
        output_file = workdir + "ltng_conus/" + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        
        ########################################################################
        ################# NOW PLOT MIDATLANTIC DOMAIN ########################
        ########################################################################
        plt.figure(figsize=[8, 8], dpi=100)
        xH, yH = DH(lons, lats)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        newmap = DH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
        newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
        DH.drawcoastlines(linewidth=0.7, color = 'k')
        DH.drawcountries(linewidth=0.7, color = 'k')
        DH.drawstates(linewidth=0.7, color = 'k')
        
        plt.title('NOAA GOES-16 Lightning Mapper\n%s' % local.strftime('%B %d, %Y %H:%M ') + et)
        # add logo
        im1 = image.imread("/home/sat_ops/goes_r/nexrad/cema38.png")
        im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
        #im[:, :, -1] = 0.5
        plt.figimage(im1, 535, 633, zorder=1)
        plt.figimage(im2, 13, 633, zorder=1)
        # save file
        output_file = workdir + "img_mid/" + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')

        symbol = u'$\u26A1$'
        #symbol = u'$\u2193$'
        # Plot flashes as small red dots
        #group_x, group_y = mH(group_lon, group_lat)
        for g in xrange(0, len(abi_glm[n])):
            group_x, group_y = DH(ltng_lon[g], ltng_lat[g])
            DH.scatter(group_x, group_y, s=14, marker=symbol, c='yellow', zorder=3, edgecolor='yellow', lw=0)
        
        output_file = workdir + "ltng_mid/" + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        Cnight.close()
    
else:
    print "Up to Date"
    
