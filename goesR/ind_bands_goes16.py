# script designed for basin.ceoe.udel.edu
# James Simkins
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap, ListedColormap # Linear interpolation for color maps
import urllib
from netCDF4 import Dataset, num2date
from matplotlib import ticker
import matplotlib as mpl
import pyart
from dateutil import tz
import time
from time import mktime
import os.path
from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.image as image
from datetime import datetime, timedelta
from pyproj import Proj
from matplotlib.patches import Rectangle
import warnings
from metpy.plots.ctables import registry
from cpt_convert import loadCPT # Import the CPT convert function
import imageio
import numpy as np
warnings.simplefilter("ignore", category=DeprecationWarning)

############# Initial Set Up ##################
############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/indbands/"
datadir = "/home/sat_ops/goesR/data/mcmipc/"
site = 'KDOX'


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

# Define projections
mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3100, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)

DH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
    llcrnrlat=lat0-4.5,llcrnrlon=lon0-5.5,
    urcrnrlat=lat0+5,urcrnrlon=lon0+5.5,resolution='h') 
    

mcmipc_list = [f for f in listdir(datadir) if isfile(join(datadir, f))]
fnamelist = sorted(mcmipc_list)[-3:]

dataset_list = []
for r in fnamelist:
    if os.path.isfile(workdir + 'conus/imgband01/' + r.split('.')[0] + ".png") == False:
        dataset_list.append(r)
        
if len(dataset_list) > 0:       
    for n in dataset_list:
        Cnight = Dataset(datadir + n, 'r')
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
        xH, yH = mH(lons, lats)
        xD, yD = DH(lons, lats)
        for Band in abi_no:
            # Declare initial information for the band we want
            Bandint = int(Band)
            if Bandint <= 6:
            # Converts a CPT file to be used in Python
                cpt_convert = 'Greys_r'
                v_min = 0
                v_max = 100
                clabeltext = 'Albedo [%]'
                state_col = 'cyan'
            elif Bandint == 7:
                # Converts a CPT file to be used in Python
                colorscheme = 'SVGAIR2_TEMP.cpt'
                v_min = -50.15
                v_max = 120
                clabeltext = "Brightness Temperature [DegC]"
                cpt = loadCPT('/home/sat_ops/goesR/indbands/colortables/' + colorscheme)
                cpt_convert = LinearSegmentedColormap('cpt', cpt)
                state_col = 'cyan'
            elif Bandint > 7 and Bandint < 11:
                # Converts a CPT file to be used in Python
                v_min = -112.15
                v_max = 56.85
                clabeltext = "Brightness Temperature [DegC]"
                # NOAAs color scheme for water vapor
                wv_cmap = registry.get_colortable('WVCIMSS')
                state_col = 'black'
            elif Bandint > 10:
                colorscheme = 'IR4AVHRR6.cpt'
                v_min = -103
                v_max = 84
                clabeltext = "Brightness Temperature [DegC]"
                cpt = loadCPT('/home/sat_ops/goesR/indbands/colortables/' + colorscheme)
                cpt_convert = LinearSegmentedColormap('cpt', cpt) 
                state_col = 'black'
            
            # Load the Band Data
            data = Cnight.variables['CMI_C'+ Band][:].data
            if Bandint <= 6:
                data[data==-1] = np.nan
                data = data*100 #(so we are in % albedo)
                data = np.maximum(data,0)
                data = np.minimum(data,100)
            else:
                data = data - 273.15

            ########################################################################
            ################# CONUS DOMAIN FIRST ########################
            ########################################################################
            rec_height = 120000
            rec_width = mH.xmax
            
            plt.figure(figsize=[16, 12], dpi=100)
            mH.drawstates(linewidth=0.7,color=state_col)
            mH.drawcountries(linewidth=0.7,color=state_col)
            mH.drawcoastlines(linewidth=0.7,color=state_col)
            if Bandint > 7 and Bandint < 11:
                mH.pcolormesh(xH, yH, data, cmap=wv_cmap,vmin = v_min, vmax=v_max)
            else:
                mH.pcolormesh(xH, yH, data, cmap=cpt_convert,vmin = v_min, vmax=v_max)
    
            if Bandint <= 6:
                # Insert the colorbar at the bottom
                cb = mH.colorbar(location='bottom', size = '2%', pad = '-2.15%',  ticks=[20, 40, 60, 80]) 
                cb.ax.set_xticklabels(['20', '40', '60', '80'], color = 'orangered')
            else:
                # Insert the colorbar at the bottom
                cb = mH.colorbar(location='bottom', size = '2%', pad = '-1.95%')
            
            cb.outline.set_visible(False) # Remove the colorbar outline
            cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
            cb.ax.xaxis.set_tick_params(pad=-13.75) # Put the colobar labels inside the colorbar
            cb.ax.tick_params(axis='x', colors='black', labelsize=10) # Change the color and size of the colorbar labels
            
            title = "Band " + channel_list[Bandint - 1]
            timestr = local.strftime('%B %d, %Y %H:%M ') + et
            timestr = local.strftime('%Y-%m-%d %H:%M ') + et
            currentAxis = plt.gca()
            currentAxis.add_patch(Rectangle((0, 0), 1000000000, rec_height * 1.3, alpha=1, zorder=3, facecolor='darkslateblue'))
            plt.text(9000, 90000,clabeltext,horizontalalignment='left', color = 'white', size=11)
            # btw, 5400000 comes from the width of mH basemap
            currentAxis.add_patch(Rectangle((0, mH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
            plt.text(4440000, 3200000,timestr,horizontalalignment='left', color = 'black', size=14)
            plt.text(9000, 3200000,title,horizontalalignment='left', color = 'black', size=14)
    
            # add logo
            im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes38.png")
            plt.figimage(im1, 15, 50, zorder=1)
            
            # save file
            output_file = workdir + 'conus/imgband' + Band + '/' + n.split('.')[0]
            plt.savefig(output_file, dpi=100, bbox_inches='tight')
            web_file = '/var/www/html/imagery/static_goes/conus/band' + Band + '/' + n.split('.')[0]
            plt.savefig(web_file, dpi=100, bbox_inches='tight')
            plt.close()
            
            ########################################################################
            ################# NOW PLOT MIDATLANTIC DOMAIN ########################
            ########################################################################
            rec_height = 40000
            rec_width = DH.xmax
            
            plt.figure(figsize=[8, 8], dpi=100)
            DH.drawstates(linewidth=0.7,color=state_col)
            DH.drawcountries(linewidth=0.7,color=state_col)
            DH.drawcoastlines(linewidth=0.7,color=state_col)
            if Bandint > 7 and Bandint < 11:
                DH.pcolormesh(xD, yD, data, cmap=wv_cmap, vmin = v_min, vmax = v_max)
            else:
                DH.pcolormesh(xD, yD, data, cmap=cpt_convert,vmin = v_min, vmax=v_max)
            
            if Bandint <= 6:
                cb = DH.colorbar(location='bottom', size = '2%', pad = '-1.95%',  ticks=[20, 40, 60, 80]) 
                cb.ax.set_xticklabels(['20', '40', '60', '80'], color = 'orangered')
            else:
                cb = DH.colorbar(location='bottom', size = '2%', pad = '-1.95%')
            
            cb.outline.set_visible(False) # Remove the colorbar outline
            cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
            cb.ax.xaxis.set_tick_params(pad=-10.5) # Put the colobar labels inside the colorbar
            cb.ax.tick_params(axis='x', colors='yellow', labelsize=7) # Change the color and size of the colorbar labels
            title = "Band " + channel_list[Bandint - 1]
            timestr = local.strftime('%B %d, %Y %H:%M ') + et
            timestr = local.strftime('%Y-%m-%d %H:%M ') + et
            currentAxis = plt.gca()
            currentAxis.add_patch(Rectangle((0, 0), DH.xmax, rec_height * 1.2, alpha=1, zorder=3, facecolor='darkslateblue'))
            plt.text(5000, 30000,clabeltext,horizontalalignment='left', color = 'white', size=8)
            # btw, 5400000 comes from the width of mH basemap
            currentAxis.add_patch(Rectangle((0, DH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
            plt.text(674000, 1025000,timestr,horizontalalignment='left', color = 'black', size=10)
            plt.text(7000, 1025000,title,horizontalalignment='left', color = 'black', size=10)
            
            # add logo
            im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemagoes24.png")
            plt.figimage(im1, 15, 42,zorder=1)
    
            output_file = workdir + 'midatl/imgband' + Band + '/' + n.split('.')[0]
            plt.savefig(output_file, dpi=100, bbox_inches='tight')
            plt.close()
        
        Cnight.close()


######################## Make gifs for conus ########################
######################## ######################## ######################## 

for Band in abi_no:
    imgdir = workdir + 'conus/imgband' + Band + '/'
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)[-20:]
    
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
    imageio.mimsave(workdir + 'gifs/' + 'Band' + Band + '_conus.gif', images, duration=dur_vals)

######################## ######################## ######################## 
# now for the midatlantic
######################## ######################## ######################## 
for Band in abi_no:

    imgdir = workdir + 'midatl/imgband' + Band + '/'
    
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)[-20:]
    
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
    imageio.mimsave(workdir + 'gifs/' + 'Band' + Band + '_midatlantic.gif', images, duration=dur_vals)
