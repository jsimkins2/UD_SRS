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

site = 'KDOX'
# first we need to check if KDOX is up and running

# the abi string is not in reverse, so we gotta change this
#abi = int(i)*-1
# kdox is in Dover
#get the radar location (this is used to set up the basemap and plotting grid)
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]
#use boto to connect to the AWS nexrad holdings directory
s3conn = boto.connect_s3()
bucket = s3conn.get_bucket('noaa-nexrad-level2')
#create a datetime object for the current time in UTC and use the
# year, month, and day to drill down into the NEXRAD directory structure.
now = datetime.utcnow()
date = ("{:4d}".format(now.year) + '/' + "{:02d}".format(now.month) + '/' +
        "{:02d}".format(now.day) + '/')
print date
#get the bucket list for the selected date
#Note: this returns a list of all of the radar sites with data for 
# the selected date
ls = bucket.list(prefix=date,delimiter='/')
for key in ls:
    #only pull the data and save the arrays for the site we want
    if site in key.name.split('/')[-2]:
        #set up the path to the NEXRAD files
        path = date + site + '/' + site
        #grab the last file in the file list
        fname = bucket.get_all_keys(prefix=path)[-1]
        print fname
        #get the file 
        s3key = bucket.get_key(fname)
        #save a temporary file to the local host
        localfile = tempfile.NamedTemporaryFile(delete=False)
        #write the contents of the NEXRAD file to the temporary file
        s3key.get_contents_to_filename(localfile.name)
        #use the read_nexrad_archive function from PyART to read in NEXRAD file
        radar = pyart.io.read_nexrad_archive(localfile.name)
        #get the date and time from the radar file for plot enhancement
        ktime = radar.time['units'].split(' ')[-1].split('T')
        print(site + ': ' + ktime[0] + ' at ' + ktime[1] )
        
        checktime = str(ktime[0] + ' ' + ktime[1][:-1])
        checktime = datetime.strptime(checktime, '%Y-%m-%d %H:%M:%S')
        kdifftime = now - checktime
        
        if divmod(kdifftime.days * 86400 + kdifftime.seconds, 60)[0] > 16 == True:
            site = 'KDIX'
            print "KDOX is down, using KDIX now"


# if KDIX is down, well we gotta use something I guess 
s3conn = boto.connect_s3()
bucket = s3conn.get_bucket('noaa-nexrad-level2')
#create a datetime object for the current time in UTC and use the
# year, month, and day to drill down into the NEXRAD directory structure.
now = datetime.utcnow()
date = ("{:4d}".format(now.year) + '/' + "{:02d}".format(now.month) + '/' +
        "{:02d}".format(now.day) + '/')
print date
#get the bucket list for the selected date
#Note: this returns a list of all of the radar sites with data for 
# the selected date
ls = bucket.list(prefix=date,delimiter='/')
for key in ls:
    #only pull the data and save the arrays for the site we want
    if site in key.name.split('/')[-2]:
        #set up the path to the NEXRAD files
        path = date + site + '/' + site
        #grab the last file in the file list
        fname = bucket.get_all_keys(prefix=path)[-1]
        print fname
        #get the file 
        s3key = bucket.get_key(fname)
        #save a temporary file to the local host
        localfile = tempfile.NamedTemporaryFile(delete=False)
        #write the contents of the NEXRAD file to the temporary file
        s3key.get_contents_to_filename(localfile.name)
        #use the read_nexrad_archive function from PyART to read in NEXRAD file
        radar = pyart.io.read_nexrad_archive(localfile.name)
        #get the date and time from the radar file for plot enhancement
        ktime = radar.time['units'].split(' ')[-1].split('T')
        print(site + ': ' + ktime[0] + ' at ' + ktime[1] )
        
        checktime = str(ktime[0] + ' ' + ktime[1][:-1])
        checktime = datetime.strptime(checktime, '%Y-%m-%d %H:%M:%S')
        kdifftime = now - checktime
        
        if divmod(kdifftime.days * 86400 + kdifftime.seconds, 60)[0] > 16 == True:
            site = 'KLWX'
            print "KDIX is down, using KLWX now"




lt = time.localtime()
dst = lt.tm_isdst
if dst == 0:
    et = "EDT"
else:
    et = "EST"

seq = [-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1]
for i in seq:
    # the abi string is not in reverse, so we gotta change this
    #abi = int(i)*-1
    # kdox is in Dover
    #select the radar site
    #get the radar location (this is used to set up the basemap and plotting grid)
    loc = pyart.io.nexrad_common.get_nexrad_location(site)
    lon0 = loc[1] ; lat0 = loc[0]
    #use boto to connect to the AWS nexrad holdings directory
    s3conn = boto.connect_s3()
    bucket = s3conn.get_bucket('noaa-nexrad-level2')
    #create a datetime object for the current time in UTC and use the
    # year, month, and day to drill down into the NEXRAD directory structure.
    now = datetime.utcnow()
    date = ("{:4d}".format(now.year) + '/' + "{:02d}".format(now.month) + '/' +
            "{:02d}".format(now.day) + '/')
    print date
    #get the bucket list for the selected date
    #Note: this returns a list of all of the radar sites with data for 
    # the selected date
    ls = bucket.list(prefix=date,delimiter='/')
    for key in ls:
        #only pull the data and save the arrays for the site we want
        if site in key.name.split('/')[-2]:
            #set up the path to the NEXRAD files
            path = date + site + '/' + site
            #grab the last file in the file list
            fname = bucket.get_all_keys(prefix=path)[i]
            print fname
            #get the file 
            s3key = bucket.get_key(fname)
            #save a temporary file to the local host
            localfile = tempfile.NamedTemporaryFile(delete=False)
            #write the contents of the NEXRAD file to the temporary file
            s3key.get_contents_to_filename(localfile.name)
            #use the read_nexrad_archive function from PyART to read in NEXRAD file
            radar = pyart.io.read_nexrad_archive(localfile.name)
            #get the date and time from the radar file for plot enhancement
            ktime = radar.time['units'].split(' ')[-1].split('T')
            print(site + ': ' + ktime[0] + ' at ' + ktime[1] )
            #set up the plotting grid for the data
            display = pyart.graph.RadarMapDisplay(radar)
            x,y = display._get_x_y(0,True,None)
    fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(7,7),dpi=200)
    #set up a basemap with a lambert conformal projection centered 
    # on the radar location, extending 1 degree in the meridional direction
    # and 1.5 degrees in the longitudinal in each direction away from the 
    # center point.
    mH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
               llcrnrlat=lat0-2,llcrnrlon=lon0-3,
               urcrnrlat=lat0+2.5,urcrnrlon=lon0+3,resolution='h')
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
    mH.drawcounties(linewidth=0.3,color='k',ax=ax)
    mH.drawstates(linewidth=1,color='k',ax=ax)
    mH.drawcoastlines(linewidth=0.7,color='k',ax=ax)
    #mark the radar location with a black dot
    mH.scatter(lon0,lat0,marker='o',s=20,color='k',ax=ax,latlon=True)
    mH.scatter(-75.7506,39.6780,marker='*',s=10,color='k',ax=ax,latlon=True) # UDEL
    #add the colorbar axes and create the colorbar based on the settings above
    cax = fig.add_axes([0.075,0.075,0.85,0.025])
    cbar = plt.colorbar(cs,ticks=ticks,norm=norm,cax=cax,orientation='horizontal')
    cbar.set_label(label,fontsize=12)
    cbar.ax.tick_params(labelsize=11)
    #add a title to the figure
    # get the kdox zulu time, and convert it to local time
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz('America/New_York')
    kd_year = str(ktime[0])[0:4]
    kd_month = str(ktime[0])[5:7]
    kd_day = str(ktime[0])[8:10]
    kd_hr = str(ktime[1])[0:2]
    kd_min= str(ktime[1])[3:5]
    kd_sec = str(ktime[1])[6:8]
    kdox_t1 = datetime(int(kd_year), int(kd_month), int(kd_day), int(kd_hr), int(kd_min), int(kd_sec))
    utc = kdox_t1.replace(tzinfo=from_zone)
    local = utc.astimezone(to_zone)

    kdox_newtime = kdox_t1.replace(tzinfo=from_zone)
    kdox_local = kdox_newtime.astimezone(to_zone)
    kdox_dt = datetime.fromtimestamp(mktime(kdox_local.timetuple()))
    fig.text(0.5,0.95, 'DEOS Radar Reflectivity\n'
            + kdox_dt.strftime('%Y-%m-%d at %H:%M ') + et,horizontalalignment='center',fontsize=16)
    # display the figure
    im1 = image.imread("/home/sat_ops/goes_r/nexrad/cema38.png")
    im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
    #im[:, :, -1] = 0.5
    plt.figimage(im1, 720, 805, zorder=1)
    plt.figimage(im2, 13, 805, zorder=1)
    
    output_file = '/home/sat_ops/goes_r/nexrad/image_kdox/' + str(abs(i)) + ".png"
    fig.savefig(output_file, dpi=120, bbox_inches='tight')
    plt.close()
    # close and delete the temporary file holding the radar data
    localfile.close()
    os.remove(localfile.name)


import imageio
import numpy as np
images = []
dur_vals = []
for i in xrange(1,12):
    dur_vals.append(.1)
    
dur_vals.append(2)
#print dur_vals
new_seq = range(1,13)
from collections import OrderedDict
new_seq = sorted(new_seq, key=int, reverse=True)

for i in new_seq:
    input_file='/home/sat_ops/goes_r/nexrad/image_kdox/' + str(i) + '.png'
    images.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/nexrad/kdox_radar.gif', images, duration=dur_vals)







