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
import datetime
############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/radar/"
imgdir = "imgconus/"


# Go to the Unidata Thredds Server for the Current Day
nowdate = datetime.datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/grib/nexrad/composite/unidata/NEXRAD_Unidata_Reflectivity-' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')
# run through the last 5 files for now and we'll see whether or not we've already created them or not
raw_list = []
raw_list = list(cat.catalog_refs)[-5:-1]
raw_list.append(str(cat.catalog_refs[-1]))

dataset_list = []
for r in raw_list:
    if os.path.isfile(workdir + imgdir + str(r.split('.')[0]) + ".png") == False:
        dataset_list.append(r)

for i in range(0,len(dataset_list)):
    cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/grib/nexrad/composite/unidata/NEXRAD_Unidata_Reflectivity-' + \
                      str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/' + dataset_list[i] + '/catalog.xml')
    
    dataset_name = sorted(cat.datasets.keys())[-1]
    print dataset_name
    dataset = cat.datasets[dataset_name]
    nc = dataset.remote_access()
    if i==0:
        refltime = len(nc.variables['reftime'][:]) - len(dataset_list) + i
        refltimelen = len(nc.variables['reftime'][:])
    else:
        refltime = refltime + 1
    
    print refltime
    geoy = np.array(nc.variables['y'][:]) * 1000.
    geox = np.array(nc.variables['x'][:]) * 1000.
    refl = np.array(nc.variables['Base_reflectivity_surface_layer'][refltime,:,:])
    proj_var = nc.variables['LambertConformal_Projection']
    time_var = nc.variables['time']
    timestamp = num2date(time_var[:].squeeze(), time_var.units)
    print timestamp[refltime]
    
    timestamp = timestamp[refltime]
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz('America/New_York')
    utc = timestamp.replace(tzinfo=from_zone)
    local = utc.astimezone(to_zone)
    
    lt = time.localtime()
    dst = lt.tm_isdst
    lt = time.localtime()
    dst = lt.tm_isdst
    
    if dst == 0:
        et = "EDT"
    else:
        et = "EST"
    
    # need to convert the NaNs to 0 and then mask them in order for them to all be hidden in the image
    dBZ = refl
    dBZ[np.isnan(dBZ)] = 0
    dBZ = np.ma.array(dBZ)
    dBZ[dBZ < 10] = np.ma.masked
    
    
    '''
    mH = Basemap(projection='lcc', lon_0=proj_var.longitude_of_central_meridian, lat_0=proj_var.latitude_of_projection_origin, \
                 lat_1=proj_var.standard_parallel, rsphere=proj_var.earth_radius, \
                 width=geox[-1]-geox[0], height=geoy[-1]-geoy[0], resolution='h', area_thresh=1500)
    '''
    mH = Basemap(resolution='i', projection='lcc', area_thresh=1500, \
            width=1800*3000, height=1060*3000, \
            lat_1=38.5, lat_2=38.5, \
            lat_0=38.5, lon_0=-97.5)
    # need to find the central projection coordinates give the central lat lon and add those to the geox, geoy arrays to match the new projection with original
    x0,y0 = mH(260,40)
    geox = geox + x0
    geoy = geoy + y0
    
    # create a meshgrid so that they are the same size and we can convert the projection coordinates to lat lon so basemap can understand them
    XX, YY = np.meshgrid(geox, geoy)
    lons, lats = mH(XX, YY, inverse=True)
    
    # Begin the plotting 
    plt.figure(figsize=[16, 12], dpi=100)
    # Plot the GOES image
    cmap = 'pyart_NWSRef'
    levs = np.linspace(0,80,41,endpoint=True)
    norm = mpl.colors.BoundaryNorm(levs,256)
    # Plot the HRRR reflectivity
    mH.pcolormesh(lons, lats, dBZ, latlon=True,
                  cmap=cmap,
                  vmax=80, vmin=0)
    cb = mH.colorbar(location='left', size = '2%', pad = '-2%', ticks=[5, 15, 25, 35, 45, 55, 65, 75])
    cb.ax.set_yticklabels(['5', '15', '25', '35', '45', '55', '65', '75'])
    
    cb.outline.set_visible(False) # Remove the colorbar outline
    cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb.ax.yaxis.set_tick_params(pad=-18.5) # Put the colobar labels inside the colorbar
    cb.ax.tick_params(axis='y', colors='black', labelsize=10) 
    cb.ax.set_ylabel('Reflectivity [dBZ]')
    #cb.ax.text(1.1, 0.55, 'Reflectivity [dBZ]', rotation=90)
    
    # Plot other map elements
    mH.drawcountries()
    mH.drawcounties(linewidth=0.3, color='lightgray')
    mH.drawcoastlines(linewidth=0.7,color='k')
    mH.drawstates(color='k')
    
    plt.title('NEXRAD Unidata Reflectivity\n%s' % local.strftime('%B %d, %Y %H:%M ') + et)
    im1 = image.imread("/home/sat_ops/goes_r/nexrad/cema38.png")
    im2 = image.imread("/home/sat_ops/goes_r/nexrad/udel38.png")
    plt.figimage(im1, 1210, 745, zorder=1)
    plt.figimage(im2, 15, 745, zorder=1)
    # save file
    output_file = workdir + imgdir + str(dataset_name.split('.')[0]) + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    
