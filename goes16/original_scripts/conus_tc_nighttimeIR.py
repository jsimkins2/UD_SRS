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



# create a logfile with most recent 36 files (3 hours)
with open("/home/sat_ops/goes_r/night_scan/logfile.txt") as f:
    file_names = f.readlines()
file_names = [x.strip() for x in file_names] 

fnamelist = []
# we're going to set C02 as the default here because it's the most important band
for i in xrange(1,len(file_names)):
    fname = str(file_names[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    fnamelist.append(fname)


# we only want to do timestamps where we have all 3 bands
from collections import OrderedDict
fnamelist = sorted(fnamelist, key=int)

ABI_datetime = []
for i in fnamelist:    
    if os.path.isfile("/home/sat_ops/goes_r/night_scan/image_conus/" + str(i) + ".png") == False:
        ABI_datetime.append(i)

# Make a new map object for the HRRR model domain map projection
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)
                    
# begin the loop that makes the images
if len(ABI_datetime) > 0:
    for n in xrange(0, len(ABI_datetime)):
        print n
        # C is for Conus File OR_ABI-L2-CMIPC-M3C02_G16_s20180601912.nc
        # RED BAND
        C_file = '/home/sat_ops/goes_r/night_scan/polished_data/OR_ABI-L2-MCMIPC-M3_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
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
        
        cleanIR = Cnight.variables['CMI_C13'][:].data
        cleanIR[cleanIR==-1] = np.nan
        
        # Apply range limits for clean IR channel
        cleanIR = np.maximum(cleanIR, 90)
        cleanIR = np.minimum(cleanIR, 313)
        
        # Normalize the channel between a range
        cleanIR = (cleanIR-90)/(313-90)
        
        # Invert colors
        cleanIR = 1 - cleanIR
        
        # Lessen the brightness of the coldest clouds so they don't appear so bright near the day/night line
        cleanIR = cleanIR/1.5
        #cleanIR = np.flipud(cleanIR)
        # numpy d stack it all for the final dataset
        RGB_IR = np.dstack([np.maximum(R, cleanIR), np.maximum(G_true, cleanIR), np.maximum(B, cleanIR)])
        # RGB_IR = np.flipud(RGB_IR)
        def contrast_correction(color, contrast):
            """
            Modify the contrast of an R, G, or B color channel
            See: #www.dfstudios.co.uk/articles/programming/image-programming-algorithms/image-processing-algorithms-part-5-contrast-adjustment/
        
            Input:
                C - contrast level
            """
            F = (259*(contrast + 255))/(255.*259-contrast)
            COLOR = F*(color-.5)+.5
            COLOR = np.minimum(COLOR, 1)
            COLOR = np.maximum(COLOR, 0)
            return COLOR
        
        # Modify the RGB color contrast
        contrast = 125
        RGB_contrast = contrast_correction(np.dstack([R, G_true, B]), contrast)
        RGB_contrast_IR = np.dstack([np.maximum(RGB_contrast[:,:,0], cleanIR), np.maximum(RGB_contrast[:,:,1], cleanIR), np.maximum(RGB_contrast[:,:,2], cleanIR)])
            
        
        add_seconds = Cnight.variables['t'][0]
        DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=add_seconds)

        
        # Satellite height
        sat_h = Cnight.variables['goes_imager_projection'].perspective_point_height
        
        # Satellite longitude
        sat_lon = Cnight.variables['goes_imager_projection'].longitude_of_projection_origin
        
        # Satellite sweep
        sat_sweep = Cnight.variables['goes_imager_projection'].sweep_angle_axis
        
        # The projection x and y coordinates equals
        # the scanning angle (in radians) multiplied by the satellite height (http://proj4.org/projections/geos.html)
        X = Cnight.variables['x'][:] * sat_h
        Y = Cnight.variables['y'][:] * sat_h
        
        
        # The geostationary projection is perhaps the easiest way, and we don't need to use the Proj object.
        # Essentially, we are stretching the image across a map with the same projection and dimensions.
        
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
        #rgb = np.flipud(rgb) # Force the maximum possible RGB value to be 1 (the lowest should be 0).
        colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]),3) # flatten array, becuase that's what pcolormesh wants.
        colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster?? according to stackoverflow.
        
        # Now we can plot the GOES data on the HRRR map domain and projection
        plt.figure(figsize=[16, 12], dpi=100)
        #newmap2 = mH.pcolormesh(xH, yH, cleanIR, cmap='Greys_r')
        #newmap2.set_array(None)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
        newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
        
        #mH.imshow(np.flipud(RGB_IR))
        mH.drawstates()
        mH.drawcountries()
        mH.drawcoastlines(linewidth=0.7,color='k')
        
        
    
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
        
        plt.title('NOAA GOES-16\n%s' % local.strftime('%B %d, %Y %H:%M ') + et)
        # add logo
        im1 = image.imread("/home/sat_ops/goes_r/nexrad/cema38.png")
        im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
        #im[:, :, -1] = 0.5
        plt.figimage(im1, 1210, 745, zorder=1)
        plt.figimage(im2, 15, 745, zorder=1)
        # save file
        output_file = '/home/sat_ops/goes_r/night_scan/image_conus/' + str(ABI_datetime[n]) + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        Cnight.close()
    
    else:
        print "Up to Date"
    
