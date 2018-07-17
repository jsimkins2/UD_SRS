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
from cpt_convert import loadCPT # Import the CPT convert function
from matplotlib.colors import LinearSegmentedColormap, ListedColormap # Linear interpolation for color maps
from collections import OrderedDict
from matplotlib.patches import Rectangle
import pyart
import boto

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

# 16 bands

nums = [01,02,03,04,05,06,07,010,011,10,11,12,13,14,15,16]
abi_no = []
for n in nums:
    abi_no.append("%02d"%n)

channel_list = [None] * 16
channel_list[0] = r"1 - Blue Band 0.47 $\mu$m"
channel_list[1]= r'2 - Red Band 0.64 $\mu$m'
channel_list[2]= r'3 - Veggie Band 0.86 $\mu$m'
channel_list[3]= r'4 - Cirrus Band 1.37 $\mu$m'
channel_list[4]= r'5 - Snow/Ice Band 1.6 $\mu$m'
channel_list[5]= r'6 - Cloud Particle Size Band 2.2 $\mu$m'
channel_list[6]= r'7 - Shortwave Window Band 3.9 $\mu$m'
channel_list[7]= r'8 - Upper-Level Tropo. WV Band 6.2 $\mu$m'
channel_list[8]= r'9 - Mid-Level Tropo. WV Band 6.9 $\mu$m'
channel_list[9]= r'10 - Low-Level WV Band 7.3 $\mu$m'
channel_list[10]= r'11 - Cloud-Top Phase Band 8.4 $\mu$m'
channel_list[11]= r'12 - Ozone Band 9.6 $\mu$m'
channel_list[12]= r'13 - Clean IR Longwave Band 10.3 $\mu$m'
channel_list[13]= r'14 - IR Longwave Band 11.2 $\mu$m'
channel_list[14]= r'15 - Dirty Longwave Band 12.3 $\mu$m'
channel_list[15]= r'16 - CO2 Longwave IR 13.3 $\mu$m'

#get the radar location (this is used to set up the basemap and plotting grid)
site = 'KDOX'
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]


# Make a CONUS map
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)
# make a MidAtlantic Map
DH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-4,llcrnrlon=lon0-5,
    urcrnrlat=lat0+4.5,urcrnrlon=lon0+5,resolution='h')   
    
for Band in abi_no:
    # create a logfile with most recent 24 files (2 hours): the logfile.txt is created in a shell script abi_ind.sh
    with open("/home/sat_ops/goes_r/ind_bands/logfiles/C" + Band + "_logfile.txt") as f:
        C_names = f.readlines()
    
    C_names = [x.strip() for x in C_names] 
    sname = []
    # we only want the date strings
    for i in xrange(1,len(C_names)):
        fname = str(C_names[i].split("_", 4)[3])
        fname = str(fname.split(".", 2)[0])
        fname = fname[1:]
        sname.append(fname)
    
    # have to sort the match here
    match = set(sname)
    match = sorted(match, key=int)
    
    ABI_datetime = []
    for i in match:    
        if os.path.isfile("/home/sat_ops/goes_r/ind_bands/conus/band" + Band + "/" + str(i) + ".png") == False:
            ABI_datetime.append(i)
    
    Bandint = int(Band)
    if Bandint <= 6:
    # Converts a CPT file to be used in Python
        cpt_convert = 'Greys_r'
        v_min = 0
        v_max = 100
        cblabel = 'Albedo [%]'
        state_col = 'cyan'
    elif Bandint == 7:
        # Converts a CPT file to be used in Python
        colorscheme = 'SVGAIR2_TEMP.cpt'
        v_min = -50.15
        v_max = 120
        cblabel = "Brightness Temperature [DegC]"
        cpt = loadCPT('/home/sat_ops/goes_r/ind_bands/colortables/' + colorscheme)
        cpt_convert = LinearSegmentedColormap('cpt', cpt)
        state_col = 'cyan'
    elif Bandint > 7 and Bandint < 11:
        # Converts a CPT file to be used in Python
        colorscheme = 'SVGAWVX_TEMP.cpt'
        v_min = -112.15
        v_max = 56.85
        cblabel = "Brightness Temperature [DegC]"
        #cpt_convert = LinearSegmentedColormap.from_list('this', ['darkgreen', 'green', 'lightgreen', 'white', 'blue', 'yellow', 'red', 'k'])
        # NOAAs color scheme for water vapor
        cpt_convert = LinearSegmentedColormap.from_list('this', ['#a71613', '#c8954f', '#3f3f3f', '#d2d2d2','#d2d2d2', '#08657f', '#ebdd0e', '#de810a', '#dd030a', '#262626', '#ffffff'], N=256)
        state_col = 'black'
    elif Bandint > 10:
        colorscheme = 'IR4AVHRR6.cpt'
        v_min = -103
        v_max = 84
        cblabel = "Brightness Temperature [DegC]"
        cpt = loadCPT('/home/sat_ops/goes_r/ind_bands/colortables/' + colorscheme)
        cpt_convert = LinearSegmentedColormap('cpt', cpt) 
        state_col = 'black'
    # begin the loop that makes the images
    if len(ABI_datetime) > 0:
        for n in xrange(0, len(ABI_datetime)):
            print n
            # C is for Conus File OR_ABI-L2-CMIPC-M3C02_G16_s20180601912.nc
            # RED BAND
            C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M3C' + Band + '_G16_s' + str(ABI_datetime[n]) + '.nc'  # GOES16 East
            if os.path.isfile(C_file) == False:
                C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M4C'+ Band + '_G16_s' + str(ABI_datetime[n]) + '.nc'

            C = Dataset(C_file, 'r')
            data = C.variables['CMI'][:].data
            
            if Bandint <= 6:
                data[data==-1] = np.nan
                data = data*100 #(so we are in % albedo)
                data = np.maximum(data,0)
                data = np.minimum(data,100)
            else:
                data = data - 273.15
            
            # Datetime organization
            add_seconds = C.variables['t'][0]
            DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=add_seconds)
            # Generate Local Time
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
            
            
            # Satellite height, lon, sweep 
            sat_h = C.variables['goes_imager_projection'].perspective_point_height
            sat_lon = C.variables['goes_imager_projection'].longitude_of_projection_origin
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

            # Now we can plot the GOES data on CONUS LCC projection
            ax = plt.figure(figsize=[16, 12], dpi=100)
            mH.drawstates(linewidth=0.7,color=state_col)
            mH.drawcountries(linewidth=0.7,color=state_col)
            mH.drawcoastlines(linewidth=0.7,color=state_col)
            mH.pcolormesh(xH, yH, data, cmap=cpt_convert,vmin = v_min, vmax=v_max)
        
            # plot colorbar
            if Bandint <= 6:
                # Insert the colorbar at the bottom
                cb = mH.colorbar(location='bottom', size = '2%', pad = '-4%',  ticks=[20, 40, 60, 80]) 
                cb.ax.set_xticklabels(['20', '40', '60', '80'])
            else:
                # Insert the colorbar at the bottom
                cb = mH.colorbar(location='bottom', size = '2%', pad = '-4%')

            
            cb.outline.set_visible(False) # Remove the colorbar outline
            cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
            cb.ax.xaxis.set_tick_params(pad=-14.5) # Put the colobar labels inside the colorbar
            cb.ax.tick_params(axis='x', colors='black', labelsize=8) # Change the color and size of the colorbar labels
            
            # plot black rectangle, the extents don't work so had to hardcode one in
            extent = [mH.llcrnrlon, mH.llcrnrlat, mH.urcrnrlon, mH.urcrnrlat]
            lon_difference = (extent[2] - extent[0]) # Max Lon - Min Lon
            lat_difference = (extent[3] - extent[1]) # Max lat - Min lat
            currentAxis = plt.gca()
            currentAxis.add_patch(Rectangle((0, 0), 1000000000, 120000, alpha=1, zorder=3, facecolor='black'))
            
            # Add the image description inside the black rectangle 
            Title = " NOAA GOES-16 ABI-CMI:  " + channel_list[Bandint - 1] + "      (" + cblabel + ")"
            # Insert the institution name
            Institution = "University of Delaware CEMA"
            atime = local.strftime('%H:%M ') + et + local.strftime('         %m-%d-%Y')
            plt.text(1000, 10000,Title,horizontalalignment='left', color = 'white', size=10)
            plt.text(4500000,10000,atime, horizontalalignment='left', color = 'yellow', size=10)
            
            # add UDEL, CEMA, and GOES logos
            im1 = image.imread("/home/sat_ops/goes_r/nexrad/combined_logo.png")
            plt.figimage(im1, 20, 45, zorder=1)
            
            # save the conus file 
            output_file = "/home/sat_ops/goes_r/ind_bands/conus/band" + Band + "/" + str(ABI_datetime[n]) + ".png"
            print output_file
            plt.savefig(output_file, dpi=100, bbox_inches='tight')
            plt.close()
            
            ################### Now Plot MidAtlantic Domain #######################
            xH, yH = DH(lons, lats)
            ax = plt.figure(figsize=[8, 8], dpi=100)
            DH.drawstates(linewidth=0.7,color=state_col)
            DH.drawcountries(linewidth=0.7,color=state_col)
            DH.drawcoastlines(linewidth=0.7,color=state_col)
            DH.pcolormesh(xH, yH, data, cmap=cpt_convert,vmin = v_min, vmax=v_max)
            
            # plot colorbar
            if Bandint <= 6:
                # Insert the colorbar at the bottom
                cb = mH.colorbar(location='bottom', size = '2%', pad = '-6%',  ticks=[20, 40, 60, 80]) 
                cb.ax.set_xticklabels(['20', '40', '60', '80'])
            else:
                # Insert the colorbar at the bottom
                cb = mH.colorbar(location='bottom', size = '2%', pad = '-6%')
            
            
            cb.outline.set_visible(False) # Remove the colorbar outline
            cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
            cb.ax.xaxis.set_tick_params(pad=-12.2) # Put the colobar labels inside the colorbar
            cb.ax.tick_params(axis='x', colors='black', labelsize=8) # Change the color and size of the colorbar labels
            
            # plot black rectangle, the extents don't work so had to hardcode one in
            extent = [mH.llcrnrlon, mH.llcrnrlat, mH.urcrnrlon, mH.urcrnrlat]
            lon_difference = (extent[2] - extent[0]) # Max Lon - Min Lon
            lat_difference = (extent[3] - extent[1]) # Max lat - Min lat
            currentAxis = plt.gca()
            currentAxis.add_patch(Rectangle((0, 0), 1000000, 55000, alpha=1, zorder=3, facecolor='black'))
            
            # Add the image description inside the black rectangle 
            Title = "NOAA GOES-16 ABI-CMI:  " + channel_list[Bandint - 1]
            # Insert the institution name
            Institution = "University of Delaware CEMA"
            atime = local.strftime('%H:%M ') + et + local.strftime('   %m-%d-%Y')
            clabeltext =  "("+ cblabel + ")"
            plt.text(1000, 20000,Title,horizontalalignment='left', color = 'white', size=10)
            plt.text(1000,1000,atime, horizontalalignment='left', color = 'yellow', size=10)
            plt.text(265000,3500,clabeltext, horizontalalignment='left', color = 'gray', size=9)
                        
            # add UDEL, CEMA, and GOES logos
            im1 = image.imread("/home/sat_ops/goes_r/nexrad/combined_logo_small.png")
            plt.figimage(im1, 490, 8, zorder=1)
            
            # save the conus file
            output_file = "/home/sat_ops/goes_r/ind_bands/midAtl/band" + Band + "/" + str(ABI_datetime[n]) + ".png"
            print output_file
            plt.savefig(output_file, dpi=100, bbox_inches='tight')
            plt.close()
            C.close()
        
        else:
            print "Up to Date"


