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

def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))


# create a logfile with most recent 36 files (3 hours)
with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C01_logfile.txt") as f:
    C01_names = f.readlines()

C01_names = [x.strip() for x in C01_names] 

with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C02_logfile.txt") as f:
    C02_names = f.readlines()

C02_names = [x.strip() for x in C02_names] 

with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/C03_logfile.txt") as f:
    C03_names = f.readlines()

C03_names = [x.strip() for x in C03_names] 

sname1 = []
sname2 = []
sname3 = []
# we're going to set C02 as the default here because it's the most important band
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

# we only want to do timestamps where we have all 3 bands
from collections import OrderedDict
match = set(sname1) & set(sname2) & set(sname3)
match = sorted(match, key=int)

# need to go to the last 10 files or so here, so we need to do seq -10 

ABI_datetime = []
for i in match:    
    if os.path.isfile("/home/sat_ops/goes_r/nexrad/image_nxrd_goes/" + str(i) + ".png") == False:
        ABI_datetime.append(i)

# begin the loop that makes the images
seq = range(1, len(ABI_datetime) + 40)

nxlist = []
nex_goes_match = []
nex_names = []
nex_dates = []
abi_match = []
nex_match = []
select = []
for s in seq:
    nxlist.append(s*-1)

nxlist=sorted(nxlist, key=int)

goes_date = []
if len(ABI_datetime) > 0:
    for i in xrange(0, len(ABI_datetime)):
        year = str(ABI_datetime[i])[0:4]
        jday = str(ABI_datetime[i])[4:7]
        hr = str(ABI_datetime[i])[7:9]
        mt = str(ABI_datetime[i])[9:11]
        goes_date.append(datetime(int(year), 1, 1, int(hr), int(mt)) + timedelta(int(jday) -1))
    for j in xrange(0, len(nxlist)):
        # the abi string is not in reverse, so we gotta change this
        #abi = int(i)*-1
        # kdox is in Dover
        #select the radar site
        site = 'KDOX'
        #get the radar location (this is used to set up the basemap and plotting grid)
        loc = pyart.io.nexrad_common.get_nexrad_location(site)
        lon0 = loc[1] ; lat0 = loc[0]
        #use boto to connect to the AWS nexrad holdings directory
        s3conn = boto.connect_s3()
        bucket = s3conn.get_bucket('noaa-nexrad-level2')
        #create a datetime object for the current time in UTC and use the
        # year, month, and day to drill down into the NEXRAD directory structure.
        date = ("{:4d}".format(goes_date[len(goes_date)-1].year) + '/' + "{:02d}".format(goes_date[len(goes_date)-1].month) + '/' +
                "{:02d}".format(goes_date[len(goes_date)-1].day) + '/')
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
                nex_names.append(bucket.get_all_keys(prefix=path)[nxlist[j]])
    # turn the names into datetime object            
    for i in xrange(0,len(nex_names)):
        nex = str(nex_names[i])[-20:-7]
        tem = datetime.strptime(nex, '%Y%m%d_%H%M')
        nex_dates.append(tem)
    
    for i in xrange(0, len(goes_date)):
        n = nearest(nex_dates, goes_date[i])
        select.append(n)
    
    for i in xrange(0, len(select)):
        a = nearest(goes_date, select[i])
        adex = goes_date.index(a)
        abi_match.append(adex)
        n = nearest(nex_dates, select[i])
        ndex=nex_dates.index(n)
        nex_match.append(nxlist[ndex])

    for i in xrange(0, len(abi_match)):
        abi = abi_match[i]
        nex = nex_match[i]
        print ABI_datetime[abi]
        print nex
        site = 'KDOX'
        #get the radar location (this is used to set up the basemap and plotting grid)
        loc = pyart.io.nexrad_common.get_nexrad_location(site)
        lon0 = loc[1] ; lat0 = loc[0]
        #use boto to connect to the AWS nexrad holdings directory
        s3conn = boto.connect_s3()
        bucket = s3conn.get_bucket('noaa-nexrad-level2')
        #create a datetime object for the current time in UTC and use the
        # year, month, and day to drill down into the NEXRAD directory structure.
        date = ("{:4d}".format(goes_date[abi].year) + '/' + "{:02d}".format(goes_date[abi].month) + '/' +
                "{:02d}".format(goes_date[abi].day) + '/')
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
                fname = bucket.get_all_keys(prefix=path)[nex]
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
                    llcrnrlat=lat0-2,llcrnrlon=lon0-3,
                    urcrnrlat=lat0+2.5,urcrnrlon=lon0+3,resolution='h')
        
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
        
        mH.drawstates()
        mH.drawcountries()
        mH.drawcoastlines()
        
        
        # now plot the nexrad and the goes
        plt.title('GOES-16 True Color\n%s' % DATE.strftime('%B %d, %Y %H:%M UTC'))

        fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(7,7),dpi=200)
        #set up a basemap with a lambert conformal projection centered 
        # on the radar location, extending 1 degree in the meridional direction
        # and 1.5 degrees in the longitudinal in each direction away from the 
        # center point.
        mH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
                    llcrnrlat=lat0-2,llcrnrlon=lon0-3,
                    urcrnrlat=lat0+2.5,urcrnrlon=lon0+3,resolution='h')
                   
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
        mH.drawcounties(linewidth=0.5,color='k',ax=ax)
        mH.drawstates(linewidth=1.5,color='k',ax=ax)
        mH.drawcoastlines(linewidth=1.5,color='k',ax=ax)
        #mark the radar location with a black dot
        mH.scatter(lon0,lat0,marker='o',s=20,color='k',ax=ax,latlon=True)
        mH.scatter(-75.7506,39.6780,marker='*',s=3,color='k',ax=ax,latlon=True) # UDEL
        #add the colorbar axes and create the colorbar based on the settings above
        cax = fig.add_axes([0.075,0.075,0.85,0.025])
        cbar = plt.colorbar(cs,ticks=ticks,norm=norm,cax=cax,orientation='horizontal')
        cbar.set_label(label,fontsize=12)
        cbar.ax.tick_params(labelsize=11)
        #add a title to the figure
        # need to convert to local time and grab daylight savings time info
        from dateutil import tz
        import time
        from time import mktime

        abi_time = goes_date[abi]
        from_zone = tz.gettz('UTC')
        to_zone = tz.gettz('America/New_York')
        utc = abi_time.replace(tzinfo=from_zone)
        local = utc.astimezone(to_zone)
        lt = time.localtime()
        dst = lt.tm_isdst
        if dst == 0:
            et = "EDT"
        else:
            et = "EST"

        # get the kdox zulu time, and convert it to local time
        kd_year = str(ktime[0])[0:4]
        kd_month = str(ktime[0])[5:7]
        kd_day = str(ktime[0])[8:10]
        kd_hr = str(ktime[1])[0:2]
        kd_min= str(ktime[1])[3:5]
        kd_sec = str(ktime[1])[6:8]
        kdox_t1 = datetime(int(kd_year), int(kd_month), int(kd_day), int(kd_hr), int(kd_min), int(kd_sec))
        kdox_newtime = kdox_t1.replace(tzinfo=from_zone)
        kdox_local = kdox_newtime.astimezone(to_zone)
        kdox_dt = datetime.fromtimestamp(mktime(kdox_local.timetuple()))
        fig.text(0.5,0.95, site + ' (0.5$^{\circ}$) Reflectivity ' + 
                'at ' + kdox_dt.strftime('%Y-%m-%d at %H:%M ') + et,horizontalalignment='center',fontsize=14)
                # should be .94 below
        fig.text(0.5,0.92, 'GOES-16 True Color & KDOX Reflectivity %s' % local.strftime('%Y-%m-%d at %H:%M ') + et,horizontalalignment='center',fontsize=10)
        #display the figure
        output_file = '/home/sat_ops/goes_r/nexrad/image_nxrd_goes/' + str(ABI_datetime[abi]) + ".png"
        fig.savefig(output_file, dpi=120, bbox_inches='tight')
        plt.close()
        # close and delete the temporary file holding the radar data
        localfile.close()
        os.remove(localfile.name)
    else:
        print "Up to Date"









