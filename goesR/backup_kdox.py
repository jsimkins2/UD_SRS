import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from siphon.catalog import TDSCatalog, get_latest_access_url
import urllib
from netCDF4 import Dataset, num2date
from matplotlib import ticker
import matplotlib as mpl
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
import pyart
from siphon.radarserver import get_radarserver_datasets, RadarServer
import numpy.ma as ma
############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/radar/"
site = 'KDOX'

###!!!!!!!!!!! USING KDIX BECAUSE KDOX SUCKS !!!!!!!!!!!!!!!!!#
#site = 'KDIX'
############## Import Functions for Steiner Smith Algorithm ############## Thanks to Jure Zbontar for python-ing these 
##### https://github.com/jzbontar/steiner_smith/blob/master/steiner_smith.py #####
import numpy as np
import scipy.signal
import scipy.ndimage.morphology

def ipol_nearest(src, trg, data):
    tree = scipy.spatial.cKDTree(src)
    dists, ix = tree.query(trg, k=1)
    return data[ix]

def compute_spinchange(data, window=(11, 21)):
    spin = np.abs(np.diff(data, axis=1)) > 2
    spin = np.column_stack((np.zeros(data.shape[0]), spin))
    kernel = np.ones(window)
    spin = scipy.signal.fftconvolve(spin, kernel, 'same')
    possible = np.ones_like(spin)
    possible[:,0] = 0
    possible = scipy.signal.fftconvolve(possible, kernel, 'same')
    return spin / possible

def steiner_smith(radar, refl_thresh=5, spin_thresh_a=8, spin_thresh_b=40, spin_thresh_c=15, grad_thresh=20):
    data0 = radar.get_field(0, 'reflectivity')
    x0, y0, _ = radar.get_gate_x_y_z(0)
    data2 = radar.get_field(2, 'reflectivity')
    x2, y2, _ = radar.get_gate_x_y_z(2)
    src = np.column_stack((x2.ravel(), y2.ravel()))
    trg = np.column_stack((x0.ravel(), y0.ravel()))
    data2 = ipol_nearest(src, trg, data2.ravel()).reshape(data0.shape)
    spin_thresh = (spin_thresh_a - (data0.filled(0) - spin_thresh_b) / spin_thresh_c) * 0.01
    zpixel = data0.filled(0) >= refl_thresh
    echotop = data2.filled(0) >= refl_thresh
    echotop = scipy.ndimage.morphology.binary_dilation(echotop)
    spinchange = compute_spinchange(data0.filled(0)) >= spin_thresh
    elevation_diff = np.median(radar.get_elevation(2)) - np.median(radar.get_elevation(0))
    vertgrad = np.abs(data0 - data2).filled(0) > grad_thresh * elevation_diff
    r1 = ~zpixel
    r2 =  zpixel & ~echotop
    r3 =  zpixel &  echotop & ~spinchange
    r4 =  zpixel &  echotop &  spinchange & ~vertgrad
    r5 =  zpixel &  echotop &  spinchange &  vertgrad
    return r1 | r2 | r5


def togrid(polar, x, y, gridsize=1024, lim=460):
    src = np.column_stack((x.ravel(), y.ravel()))
    grid = np.linspace(-lim * 1000, lim * 1000, gridsize)
    mgrid = np.meshgrid(grid, grid[::-1])
    trg = np.column_stack((mgrid[0].ravel(), mgrid[1].ravel()))
    grid = ipol_nearest(src, trg, polar.ravel())
    return grid.reshape(gridsize, gridsize)

def dump_ref_cmap(fname, data, vmin=-32., vmax=94.5, cmap='pyart_NWSRef'):
    cmap = plt.get_cmap(cmap)
    data = (data - vmin) / (vmax - vmin)
    data_bgr = cmap(data)[:,:,2::-1] * 255
    data_bgr[data.mask] = 0


##############################################################
##############################################################
# Grab the Radar
nowdate = datetime.utcnow()

rdr = get_radarserver_datasets('http://thredds.ucar.edu/thredds/')
url = rdr['NEXRAD Level II Radar from IDD'].follow().catalog_url
rs = RadarServer(url)
query = rs.query()
query.stations('KDOX').time(datetime.utcnow())
rs.validate_query(query)
catalog = rs.get_catalog(query)
data_str = str(catalog.datasets[0]).split('_')
data_time = datetime.strptime(data_str[2] + data_str[3].split('.')[0], '%Y%m%d%H%M')
d =  data_time - nowdate

if d.total_seconds() > 3600:
    site='KDIX'
    query.stations('KDOX').time(datetime.utcnow())
    rs.validate_query(query)
    catalog = rs.get_catalog(query)
    data_str = str(catalog.datasets[0]).split('_',)
    data_time = datetime.strptime(data_str[2] + data_str[3].split('.')[0], '%Y%m%d%H%M')
    d =  data_time - nowdate
    if d.total_seconds() > 3600:
        site='KDIX'
        query.stations('KDOX').time(datetime.utcnow())
        rs.validate_query(query)
        catalog = rs.get_catalog(query)

ds = list(catalog.datasets.values())[0]
# grab the lat/lon of site, use basemap
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]
kdoxH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
           llcrnrlat=lat0-2,llcrnrlon=lon0-3,
           urcrnrlat=lat0+2.5,urcrnrlon=lon0+3,resolution='h')

if os.path.isfile(workdir + 'ref' + site + '/' + str(ds) + ".png") == False:
    # grab the radar file
    radar = pyart.io.read_nexrad_cdm(ds.access_urls['OPENDAP'])
    # create timestamp
    timestamp = radar.time['units'].split(' ')[-1].split('T')
    timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
    timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
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
        
    # Handle data 
    display = pyart.graph.RadarMapDisplay(radar)
    x,y = display._get_x_y(0,True,None)
    x0,y0 = kdoxH(lon0,lat0)
    glons,glats = kdoxH((x0+x*1000.), (y0+y*1000.),inverse=True)
    
    # Reflectivity
    ref = radar.get_field(0, 'reflectivity')
    x, y, _ = radar.get_gate_x_y_z(0)
    grid = togrid(ref, x, y, gridsize=256, lim=250)
    clutter = steiner_smith(radar)
    dBZ = np.ma.masked_where(clutter, ref)
    
    # Velocity
    vel = np.squeeze(radar.get_field(sweep=1,field_name='velocity'))
    
    ###################################################
    ###############      Plotting      ################
    ###################################################
    rec_height = 20000
    rec_width = kdoxH.xmax
    cmap = 'pyart_NWSRef'
    plt.figure(figsize=[8, 8], dpi=100)
    kdoxH.drawcounties(linewidth=0.5, color = 'k')
    kdoxH.pcolormesh(glons, glats, dBZ, latlon=True,
          cmap='pyart_NWSRef',
          vmax=80, vmin=0)
    kdoxH.drawcoastlines(linewidth=0.7, color = 'k')
    kdoxH.drawstates(linewidth=1.1, color = 'k')
    cb = kdoxH.colorbar(location='bottom', size = '4%', pad = '-4%', ticks=[5, 15, 25, 35, 45, 55, 65, 75])
    cb.ax.set_xticklabels(['5', '15', '25', '35', '45', '55', '65', '75'])
    cb.outline.set_visible(False) # Remove the colorbar outline
    cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb.ax.xaxis.set_tick_params(pad=-16) # Put the colobar labels inside the colorbar
    cb.ax.tick_params(axis='x', colors='black', labelsize=10) # Change the color and size of the colorbar labels
    kdoxH.scatter(lon0,lat0,marker='o',s=20,color='k',latlon=True)
    kdoxH.scatter(-75.7506,39.6780,marker='*',s=10,color='k',latlon=True) # UDEL
    
    clabeltext='Base Reflectivity [dBZ]'
    title = 'CEMA Radar Reflectivity'
    timestr = local.strftime('%Y-%m-%d %H:%M ') + et
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0, 0), kdoxH.xmax, rec_height * 1.8, alpha=1, zorder=3, facecolor='darkslateblue'))
    plt.text(5000, rec_height * 1.25,clabeltext,horizontalalignment='left', color = 'white', size=10)
    # btw, 5400000 comes from the width of mH basemap
    currentAxis.add_patch(Rectangle((0, kdoxH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
    plt.text(357000, kdoxH.ymax - 15000,timestr,horizontalalignment='left', color = 'black', size=12)
    plt.text(7000, kdoxH.ymax - 15000,title,horizontalalignment='left', color = 'black', size=12)
    
    
    # add logo
    im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemaunidata24.png")
    plt.figimage(im1, 15, 57,zorder=1)
    # save file
    
    output_file = workdir + 'ref' + site + '/' + str(ds) + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    
    
    ###################### Plot Velocity ######################
    
    vel = np.ma.array(vel)
    
    plt.figure(figsize=[8, 8], dpi=100)
    kdoxH.drawcounties(linewidth=0.5, color = 'k')
    kdoxH.pcolormesh(glons, glats, vel, latlon=True,
          cmap='pyart_NWSVel', vmax=20, vmin=-20)
    kdoxH.drawcoastlines(linewidth=0.7, color = 'k')
    kdoxH.drawstates(linewidth=1.1, color = 'k')
    cb = kdoxH.colorbar(location='bottom', size = '4%', pad = '-4%', ticks=[-25, -15, -5, 5, 15, 25])
    cb.ax.set_xticklabels(['-25', '-15', '-5', '5', '15', '25'])
    cb.outline.set_visible(False) # Remove the colorbar outline
    cb.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb.ax.xaxis.set_tick_params(pad=-16) # Put the colobar labels inside the colorbar
    cb.ax.tick_params(axis='x', colors='black', labelsize=10) # Change the color and size of the colorbar labels
    kdoxH.scatter(lon0,lat0,marker='o',s=20,color='k',latlon=True)
    kdoxH.scatter(-75.7506,39.6780,marker='*',s=10,color='k',latlon=True) # UDEL
    
    clabeltext='Radial Velocity [m/s]'
    title = 'CEMA Radial Velocity'
    timestr = local.strftime('%Y-%m-%d %H:%M ') + et
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0, 0), kdoxH.xmax, rec_height * 1.8, alpha=1, zorder=3, facecolor='darkslateblue'))
    plt.text(5000, rec_height * 1.25,clabeltext,horizontalalignment='left', color = 'white', size=10)
    # btw, 5400000 comes from the width of mH basemap
    currentAxis.add_patch(Rectangle((0, kdoxH.ymax - rec_height), rec_width, rec_height , alpha=1, zorder=3, edgecolor='black',facecolor='white'))
    plt.text(357000, kdoxH.ymax - 15000,timestr,horizontalalignment='left', color = 'black', size=12)
    plt.text(7000, kdoxH.ymax - 15000,title,horizontalalignment='left', color = 'black', size=12)
    # add logo
    im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemaunidata24.png")
    plt.figimage(im1, 15, 57,zorder=1)
    # save file
    
    output_file = workdir + 'vel' + site + '/' + str(ds) + ".png"
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    
    ###################################################
    ###############   Create the Gifs  ################
    ###################################################
    
    
    # Now create the gif
    import imageio
    imgdir = workdir + 'vel' + site + '/' 
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)[-15:]
    
    imglen = len(img_names)
    images = []
    dur_vals = []
    for i in range(0,imglen -1):
        if i != imglen:
            dur_vals.append(.07)
    dur_vals.append(2)
    
    for i in img_names:
        input_file=imgdir + str(i)
        images.append(imageio.imread(input_file))
    imageio.mimsave(workdir + 'kdox_velocity.gif', images, duration=dur_vals)
    
    
    ###################################################################
    imgdir = workdir + 'ref' + site + '/' 
    img_list = [f for f in listdir(imgdir) if isfile(join(imgdir, f))]
    img_names = sorted(img_list)[-15:]
    
    imglen = len(img_names)
    images = []
    dur_vals = []
    for i in range(0,imglen -1):
        if i != imglen:
            dur_vals.append(.07)
    dur_vals.append(2)
    
    for i in img_names:
        input_file=imgdir + str(i)
        images.append(imageio.imread(input_file))
    imageio.mimsave(workdir + 'kdox_radar.gif', images, duration=dur_vals)


