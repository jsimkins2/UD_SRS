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

seq = [-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1]
for i in seq:
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
    now = datetime.utcnow()
    date = ("{:4d}".format(now.year) + '/' + "{:02d}".format(now.month) + '/' +
            "{:02d}".format(now.day) + '/')
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
            #get the file 
            s3key = bucket.get_key(fname)
            #save a temporary file to the local host
            localfile = tempfile.NamedTemporaryFile(delete=False)
            #write the contents of the NEXRAD file to the temporary file
            s3key.get_contents_to_filename(localfile.name)
            #use the read_nexrad_archive function from PyART to read in NEXRAD file
            radar = pyart.io.read_nexrad_archive(localfile.name)
            #get the date and time from the radar file for plot enhancement
            time = radar.time['units'].split(' ')[-1].split('T')
            print(site + ': ' + time[0] + ' at ' + time[1] )
    
            #set up the plotting grid for the data
            display = pyart.graph.RadarMapDisplay(radar)
            x,y = display._get_x_y(0,True,None)

    fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(9,9),dpi=100)
    #set up a basemap with a lambert conformal projection centered 
    # on the radar location, extending 1 degree in the meridional direction
    # and 1.5 degrees in the longitudinal in each direction away from the 
    # center point.
    mH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
               llcrnrlat=lat0-1.4,llcrnrlon=lon0-1.75,
               urcrnrlat=lat0+1.5,urcrnrlon=lon0+1.75,resolution='h')
               

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
    mH.scatter(-75.7506,39.6780,marker='*',s=20,color='k',ax=ax,latlon=True) # UDEL
    #add the colorbar axes and create the colorbar based on the settings above
    cax = fig.add_axes([0.075,0.075,0.85,0.025])
    cbar = plt.colorbar(cs,ticks=ticks,norm=norm,cax=cax,orientation='horizontal')
    cbar.set_label(label,fontsize=12)
    cbar.ax.tick_params(labelsize=11)
    #add a title to the figure
    fig.text(0.5,0.92, site + ' (0.5$^{\circ}$) Reflectivity\n ' + 
            time[0] + ' at ' + time[1],horizontalalignment='center',fontsize=16)
    #display the figure
    output_file = '/home/sat_ops/goes_r/nexrad/image_kdox/' + str(abs(i)) + ".png"
    fig.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()
    
    
    # close and delete the temporary file holding the radar data
    localfile.close()
    os.remove(localfile.name)


import imageio
import numpy as np
images = []
dur_vals = []
for i in xrange(1,12):
    dur_vals.append(.15)
dur_vals.append(2)
#print dur_vals
new_seq = range(1,13)
from collections import OrderedDict
new_seq = sorted(new_seq, key=int, reverse=True)

for i in new_seq:
    input_file='/home/sat_ops/goes_r/nexrad/image_kdox/' + str(i) + '.png'
    images.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/nexrad/kdox_radar.gif', images, duration=dur_vals)







