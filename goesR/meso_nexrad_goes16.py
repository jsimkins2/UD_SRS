site = 'KDIX'
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
import pyart.io
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

# define function that finds nearest time given a list of datetime strings
def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

# create a list from logfiles - ABI bands 1,2,3 or C01,C02,C03
with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C01_logfile.txt") as f:
    C01_names = f.readlines()

# remove extra /n that comes with f.readlines
C01_names = [x.strip() for x in C01_names] 

with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C02_logfile.txt") as f:
    C02_names = f.readlines()

# remove extra /n that comes with f.readlines
C02_names = [x.strip() for x in C02_names] 

with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C03_logfile.txt") as f:
    C03_names = f.readlines()

# remove extra /n that comes with f.readlines
C03_names = [x.strip() for x in C03_names] 
sname1 = []
sname2 = []
sname3 = []

# we're going to set C02 as the default here because it's the 'most important' band
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

# we only want to do timestamps where we have all 3 bands, otherwise we can't make a true color
from collections import OrderedDict
match = set(sname1) & set(sname2) & set(sname3)
match = sorted(match, key=int)

# if the image already exists, we don't want to duplicate it, so let's not add it to the list
ABI_datetime = []
for i in match:    
    if os.path.isfile("/home/sat_ops/goes_r/cloud_prod/noaa_format/meso/nex_image/" + str(i) + ".png") == False:
        ABI_datetime.append(i)

# convert strings to datetime object
goes_date = []
for i in xrange(0, len(ABI_datetime)):
    year = str(ABI_datetime[i])[0:4]
    jday = str(ABI_datetime[i])[4:7]
    hr = str(ABI_datetime[i])[7:9]
    mt = str(ABI_datetime[i])[9:11]
    goes_date.append(datetime(int(year), 1, 1, int(hr), int(mt)) + timedelta(int(jday) -1))

# grab the locally stored nexrad files
with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/meso/data_nexrad.txt") as f:
    nex_names = f.readlines()

# remove extra /n that comes with f.readlines
nex_names = [x.strip() for x in nex_names] 

nex_dates = []
abi_match = []
nex_match = []

# convert to datetime objects
for i in xrange(0,len(nex_names)):
    nex = str(nex_names[i])[5:]
    tem = datetime.strptime(nex, '%Y%m%d_%H%M')
    nex_dates.append(tem)
    
# find closest matching dates
for i in xrange(0, len(nex_dates)):
    a = nearest(goes_date, nex_dates[i])
    adex = goes_date.index(a)
    abi_match.append(adex)
    
abi_match = sorted(list(set(abi_match)))

for i in xrange(0, len(abi_match)):
    print goes_date[abi_match[i]]
    n = nearest(nex_dates, goes_date[abi_match[i]])
    ndex=nex_dates.index(n)
    print nex_dates[ndex]
    nex_match.append(ndex)
    
nex_match = sorted(list(set(nex_match)))

# clear out abi_match and re-make so we don't have duplicates
abi_match = []
for i in xrange(0, len(nex_match)):
    a = nearest(goes_date, nex_dates[nex_match[i]])
    adex = goes_date.index(a)
    abi_match.append(adex)
    
if len(abi_match) > 0:
    # before we go in the loop, we have ot extract some info about KDOX
    
    #get the radar location (this is used to set up the basemap and plotting grid)
    loc = pyart.io.nexrad_common.get_nexrad_location(site)
    lon0 = loc[1] ; lat0 = loc[0]
    
    # need to convert to local time and grab daylight savings time info
    lt = time.localtime()
    dst = lt.tm_isdst
    if dst == 0:
        et = "EDT"
    else:
        et = "EST"
    
    for i in xrange(0, len(abi_match)):
        abi = abi_match[i]
        nex = nex_match[i]
        print goes_date[abi]
        print nex_dates[nex]
        
        diff = nex_dates[nex] - goes_date[abi]
        time_diff = divmod(diff.days * 86400 + diff.seconds, 60)
        
        if abs(time_diff[0]) < 11:
            
            # read in the stored radar data
            radar = pyart.io.read_cfradial('/home/sat_ops/goes_r/cloud_prod/noaa_format/meso/data/' + nex_names[nex])
            display = pyart.graph.RadarMapDisplay(radar)
            x,y = display._get_x_y(0,True,None)
            # C is for Conus File OR_ABI-L2-CMIPC-M3C02_G16_s20180601912.nc
            C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M3C02_G16_s' + str(ABI_datetime[abi]) + '.nc'  # GOES16 East
            C = Dataset(C_file, 'r')
            # Load the RGB arrays and apply a gamma correction (square root)
            R = C.variables['CMI'][:].data # Band 2 is red (0.64 um)
            R = np.sqrt(block_mean(R, 2))
            
            C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M3C03_G16_s' + str(ABI_datetime[abi]) + '.nc' # GOES16 East
            C = Dataset(C_file, 'r')
            # Load the RGB arrays and apply a gamma correction (square root)
            G = np.sqrt(C.variables['CMI'][:].data) # Band 3 is "green" (0.865 um)
            
            C_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR_ABI-L2-CMIPC-M3C01_G16_s' + str(ABI_datetime[abi]) + '.nc' # GOES16 East
            C = Dataset(C_file, 'r')
            # Load the RGB arrays and apply a gamma correction (square root)
            B = np.sqrt(C.variables['CMI'][:].data) # Band 1 is blue (0.47 um)
            
            # "True Green" is some linear interpolation between the three channels
            # note that I've added some multiplying factors here to enhance contrast
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
            
            # Make a new map object for the HRRR model domain map projection
            mH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
                        llcrnrlat=lat0-3,llcrnrlon=lon0-4,
                        urcrnrlat=lat0+3.5,urcrnrlon=lon0+4,resolution='h')
            
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
            plt.figure(figsize=[7, 7])
            
            # The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
            newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
            newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.

            fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(7,7),dpi=200)
            #set up a basemap with a lambert conformal projection centered 
            # on the radar location, extending 1 degree in the meridional direction
            # and 1.5 degrees in the longitudinal in each direction away from the 
            # center point.
            mH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
                        llcrnrlat=lat0-3,llcrnrlon=lon0-4,
                        urcrnrlat=lat0+3.5,urcrnrlon=lon0+4,resolution='h')
                       
            newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
            newmap.set_array(None) # without this, the linewidth is set to zero, but the RGB color is ignored
            
            #get the plotting grid into lat/lon coordinates
            x0,y0 = mH(lon0,lat0)
            glons,glats = mH((x0+x*1000.), (y0+y*1000.),inverse=True)
            #read in the lowest scan angle reflectivity field in the NEXRAD file 
            refl = np.squeeze(radar.get_field(sweep=0,field_name='reflectivity'))
            # Mask points with no reflectivity
            dBZ = refl
            dBZ = np.ma.array(dBZ)
            dBZ[dBZ == -10] = np.ma.masked
            #set up the plotting parameters (NWSReflectivity colormap, contour levels,
            # and colorbar tick labels)
            cmap = 'pyart_NWSRef'
            levs = np.linspace(0,80,41,endpoint=True)
            ticks = np.linspace(0,80,9,endpoint=True)
            label = 'Radar Reflectivity Factor ($\mathsf{dBZ}$)'
            #define the plot axis to the be axis defined above
            ax = axes
            #normalize the colormap based on the levels provided above
            norm = mpl.colors.BoundaryNorm(levs,256)
            cs = mH.pcolormesh(glons,glats,dBZ,norm=norm,cmap=cmap,ax=ax,latlon=True)
            #add geographic boundaries and lat/lon labels
            mH.drawparallels(np.arange(20,70,0.5),labels=[1,0,0,0],fontsize=12,
                            color='k',ax=ax,linewidth=0.001)
            mH.drawmeridians(np.arange(-150,-50,1),labels=[0,0,1,0],fontsize=12,
                           color='k',ax=ax,linewidth=0.001)
            mH.drawcounties(linewidth=0.1,color='k',ax=ax)
            mH.drawstates(linewidth=1,color='k',ax=ax)
            mH.drawcoastlines(linewidth=0.7,color='k',ax=ax)
            #mark the radar location with a black dot
            mH.scatter(lon0,lat0,marker='o',s=20,color='k',ax=ax,latlon=True)
            mH.scatter(-75.7506,39.6780,marker='*',s=3,color='k',ax=ax,latlon=True) # UDEL
            #add the colorbar axes and create the colorbar based on the settings above
            cax = fig.add_axes([0.075,0.075,0.85,0.025])
            cbar = plt.colorbar(cs,ticks=ticks,norm=norm,cax=cax,orientation='horizontal')
            cbar.set_label(label,fontsize=12)
            cbar.ax.tick_params(labelsize=11)
            #add a title to the figure

            abi_time = goes_date[abi]
            from_zone = tz.gettz('UTC')
            to_zone = tz.gettz('America/New_York')
            utc = abi_time.replace(tzinfo=from_zone)
            local = utc.astimezone(to_zone)

            # get the kdox zulu time, and convert it to local time
            kdox_t1 = nex_dates[nex]
            kdox_newtime = kdox_t1.replace(tzinfo=from_zone)
            kdox_local = kdox_newtime.astimezone(to_zone)
    
            #fig.text(0.5,0.95, site + ' (0.5$^{\circ}$) Reflectivity ' + 
                    #'at ' + kdox_local.strftime('%Y-%m-%d at %H:%M ') + et,horizontalalignment='center',fontsize=14)
                    # should be .94 below
            fig.text(0.5,0.94, 'GOES-16 True Color & ' + site + ' Reflectivity \n' +
                 local.strftime('%Y-%m-%d %H:%M ') + et,horizontalalignment='center',fontsize=12)
            # add logo
            im = image.imread("/home/sat_ops/goes_r/nexrad/xsmall_ud_cema.png")
            plt.figimage(im, origin = 'upper', zorder=1)
            # save file
            output_file = '/home/sat_ops/goes_r/cloud_prod/noaa_format/meso/nex_image/' + str(ABI_datetime[abi]) + ".png"
            fig.savefig(output_file, dpi=120, bbox_inches='tight')
            plt.close()
        
    else:
        print "Up to Date"






