# This script creates true color imagery for CONUS and Mid Atlantic & Lightning/TC over CONUS and Mid Atlantic
# 4 different images altogether 
# James Simkins
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
import pyart
#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)


############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/truecolor/"
datadir = "/home/sat_ops/goesR/data/mcmipc/"
ltngdir = "/home/sat_ops/goesR/data/glm/"
imgdir = "/home/sat_ops/goesR/truecolor/img_conus/"
site = 'KDOX'

############## Define Functions that will be used later #########################
# initial block mean allows us to downsample the Red Band which is a higher resolution than the others
def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

def contrast_correction(color, contrast):
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.minimum(COLOR, 1)
    COLOR = np.maximum(COLOR, 0)
    return COLOR

def JulianDate_to_MMDDYYY(y,jd):
    month = 1
    day = 0
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month = month + 1
    return month,jd,y

################ Grab the Lat/Lon of the site we want ####################
# Note, this is for the regional map
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]


# create a logfile with most recent files created in the shell script
file_names = [f for f in listdir(datadir) if isfile(join(datadir, f))]

fnamelist = []
for i in range(0,len(file_names)):
    fname = str(file_names[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    fnamelist.append(fname)

fnamelist = sorted(fnamelist, key=int)[-3:]
GLM_files = [f for f in listdir(ltngdir) if isfile(join(ltngdir, f))]
GLM_names = []

for i in range(0,len(GLM_files)):
    fname = str(GLM_files[i].split("_", 4)[3])
    fname = str(fname.split(".", 2)[0])
    fname = fname[1:]
    GLM_names.append(fname)

lnamelist = sorted(GLM_names, key=int)[-501:]
# sort through the files and make sure that they don't exist before we process them

ABI_datetime = []
for i in fnamelist:    
    if os.path.isfile(imgdir + str(i) + ".png") == False:
        ABI_datetime.append(i)


ldatetime = []
for t in lnamelist:
    jday = t[4:7]
    year = t[0:4]
    mdy = JulianDate_to_MMDDYYY(int(year),int(jday))
    hms = t[7:13]
    t = str(mdy[2]) + jday + hms
    ldatetime.append(datetime.strptime(t, '%Y%j%H%M%S'))

# Make a new map object for the HRRR model domain map projection
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3100, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)

DH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-4.5,llcrnrlon=lon0-5.5,
    urcrnrlat=lat0+5,urcrnrlon=lon0+5.5,resolution='h') 
    

# begin the loop that makes the images
if len(ABI_datetime) > 0:
    for n in range(0, len(ABI_datetime)):
        t = ABI_datetime[n]
        gdatetime=datetime.strptime(t, '%Y%j%H%M%S')
        ltng_index = ldatetime.index(nearest(ldatetime, gdatetime))
        ltng_files = lnamelist[ltng_index - 15: ltng_index]
        
        C_file = datadir + 'OR_ABI-L2-MCMIPC-M3_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
        if os.path.isfile(C_file) == False:
            C_file = datadir + 'OR_ABI-L2-MCMIPC-M4_G16_s' + str(ABI_datetime[n]) + '.nc'
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
        
        for lt in range(0, len(ltng_files)):
            ltfile = ltng_files[lt]
            L_file = ltngdir + 'OR_GLM-L2-LCFA_G16_s' + str(ltfile) + '.nc'  # GOES16 East
            L = Dataset(L_file, 'r')
            ltng_lat[lt] = L.variables['flash_lat'][:]
            ltng_lon[lt] = L.variables['flash_lon'][:]
        
        for lt in range(0, len(ltng_lat)):
            if lt==0:
                conus_flash_count = 0
            ylt = len(ltng_lat[lt])
            conus_flash_count = conus_flash_count + ylt

        for lt in range(0, len(ltng_lat)):
            subset = []
            llat=DH.latmin
            hlat=DH.latmax
            llon=DH.lonmin
            rlon=DH.lonmax
            if lt==0:
                midatl_flash_count = 0
            for v in range(0,len(ltng_lat[lt])):
                latval=ltng_lat[lt][v]
                lonval=ltng_lon[lt][v]
                if latval > llat and latval < hlat and lonval > llon and lonval < rlon:
                    subset.append(latval)
            midatl_flash_count = midatl_flash_count + len(subset)

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
        
        rec_height = 120000
        rec_width = mH.xmax
        # Now we can plot the GOES data on the HRRR map domain and projection
        plt.figure(figsize=[16, 12], dpi=100)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
        newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
        mH.drawstates(linewidth=0.7,color='k')
        mH.drawcountries(linewidth=0.7,color='k')
        mH.drawcoastlines(linewidth=0.7,color='k')

        title = 'NOAA GOES-16 True Color'
        timestr = local.strftime('%B %d, %Y %H:%M ') + et
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        currentAxis = plt.gca()
        currentAxis.add_patch(Rectangle((0, mH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
        plt.text(4440000, 3200000,timestr,horizontalalignment='left', color = 'black', size=14)
        plt.text(9000, 3200000,title,horizontalalignment='left', color = 'black', size=14)
        
        # add logo
        im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes38.png")
        plt.figimage(im1, 15, 15, zorder=1)
        
        # save file
        output_file = workdir + 'img_conus/' + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        
        # Lightning plotting time 
        plt.figure(figsize=[16, 12], dpi=100)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
        newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
        mH.drawstates(linewidth=0.7,color='k')
        mH.drawcountries(linewidth=0.7,color='k')
        mH.drawcoastlines(linewidth=0.7,color='k')
        
        symbol = u'$\u26A1$'
        #group_x, group_y = mH(group_lon, group_lat)
        for g in range(0, len(ltng_files)):
            group_x, group_y = mH(ltng_lon[g], ltng_lat[g])
            mH.scatter(group_x, group_y, s=18, marker=symbol, c='red', zorder=3, edgecolor='red', lw=0)
        
        currentAxis = plt.gca()
        currentAxis.add_patch(Rectangle((0, mH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
        plt.text(4440000, 3200000,timestr,horizontalalignment='left', color = 'black', size=14)
        plt.text(9000, 3200000,title + ' & GLM Lightning Flashes',horizontalalignment='left', color = 'black', size=14)
        clabeltext = 'Flash Count=' + str(conus_flash_count)
        currentAxis.add_patch(Rectangle((0, 0), 1000000000, rec_height * 0.8, alpha=1, zorder=3, facecolor='darkslateblue'))
        plt.text(9000, 15000,clabeltext,horizontalalignment='left', color = 'white', size=14)
        
        plt.figimage(im1, 15, 38, zorder=1)
        output_file = workdir + "ltng_conus/" + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        
        ########################################################################
        ################# NOW PLOT MIDATLANTIC DOMAIN ########################
        ########################################################################
        imgdir = "/home/sat_ops/goesR/truecolor/img_mid/"
        rec_height = 40000
        rec_width = DH.xmax
        
        plt.figure(figsize=[8, 8], dpi=100)
        xH, yH = DH(lons, lats)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        newmap = DH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
        newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
        DH.drawcoastlines(linewidth=0.7, color = 'k')
        DH.drawcountries(linewidth=0.7, color = 'k')
        DH.drawstates(linewidth=0.7, color = 'k')
        
        title = 'NOAA GOES-16 True Color'
        timestr = local.strftime('%Y-%m-%d %H:%M ') + et
        currentAxis = plt.gca()
        currentAxis.add_patch(Rectangle((0, DH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
        plt.text(674000, 1025000,timestr,horizontalalignment='left', color = 'black', size=10)
        plt.text(7000, 1025000,title,horizontalalignment='left', color = 'black', size=10)


        # add logo
        im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes24.png")
        plt.figimage(im1, 15, 12,zorder=1)
        
        # save file
        output_file = workdir + "img_mid/" + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        
        plt.figure(figsize=[8, 8], dpi=100)
        xH, yH = DH(lons, lats)
        # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
        newmap = DH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
        newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.
        DH.drawcoastlines(linewidth=0.7, color = 'k')
        DH.drawcountries(linewidth=0.7, color = 'k')
        DH.drawstates(linewidth=0.7, color = 'k')
        
        symbol = u'$\u26A1$'
        #symbol = u'$\u2193$'
        # Plot flashes as small red dots
        #group_x, group_y = mH(group_lon, group_lat)
        for g in range(0, len(ltng_files)):
            group_x, group_y = DH(ltng_lon[g], ltng_lat[g])
            DH.scatter(group_x, group_y, s=18, marker=symbol, c='red', zorder=3, edgecolor='red', lw=0)

        currentAxis = plt.gca()
        currentAxis.add_patch(Rectangle((0, DH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
        plt.text(674000, 1025000,timestr,horizontalalignment='left', color = 'black', size=10)
        plt.text(7000, 1025000,title + ' & GLM Lightning Flashes',horizontalalignment='left', color = 'black', size=10)
        clabeltext = 'Flash Count=' + str(midatl_flash_count)
        currentAxis.add_patch(Rectangle((0, 0), DH.xmax, rec_height * 0.8, alpha=1, zorder=3, facecolor='darkslateblue'))
        plt.text(5000, 8000,clabeltext,horizontalalignment='left', color = 'white', size=10)
        
        plt.figimage(im1, 15, 32,zorder=1)
        output_file = workdir + "ltng_mid/" + ABI_datetime[n] + ".png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        Cnight.close()
    
else:
    print "Up to Date"
    

######################## TRUE COLOR GIFS ########################
######################## ######################## ######################## 
import imageio
import numpy as np
imgdir = '/home/sat_ops/goesR/truecolor/img_conus/'
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-72:]

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
imageio.mimsave(workdir + 'truecolor_conus.gif', images, duration=dur_vals)

# now for the midatlantic
imgdir = '/home/sat_ops/goesR/truecolor/img_mid/'
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-72:]

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
imageio.mimsave(workdir + 'truecolor_midatlantic.gif', images, duration=dur_vals)

######################## TRUE COLOR WITH LIGHTNING GIFS ########################
######################## ######################## ######################## 
imgdir = '/home/sat_ops/goesR/truecolor/ltng_conus/'
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-72:]

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
imageio.mimsave(workdir + 'lightning_truecolor_conus.gif', images, duration=dur_vals)

# now for the midatlantic
imgdir = "/home/sat_ops/goesR/truecolor/ltng_mid/"
img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
img_names = sorted(img_list)[-72:]

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
imageio.mimsave(workdir + 'lightning_truecolor_midatlantic.gif', images, duration=dur_vals)
