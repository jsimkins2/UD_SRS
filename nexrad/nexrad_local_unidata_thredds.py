
# Velocity Script
# By James Simkins

# The velocity is very shakey, so I need to try and fix tha that issue
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
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
from matplotlib.patches import Rectangle
import pyart
from siphon.radarserver import RadarServer, get_radarserver_datasets
import numpy.ma as ma
import netCDF4
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import shapely.geometry as sgeom
import xarray
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from urllib.error import HTTPError
import imageio

# monkey patch from nightmare that is january 10th
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh


def plot_velocity(radar, dataset, imgdir):
    print('Plotting dataset: ', dataset)
    my_gf = pyart.filters.GateFilter(radar)
    my_gf.exclude_below('reflectivity', 12)
    my_ds_gf = pyart.correct.despeckle_field(radar, 'velocity', gatefilter=my_gf)
    
    timestamp = radar.time['units'].split(' ')[-1].split('T')
    timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
    timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S') + timedelta(minutes=4)
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
    
    lats = radar.gate_latitude
    lons = radar.gate_longitude
    min_lon = -78.24357604980469 #lons['data'].min() + 2.5
    min_lat = 36.69588851928711 #lats['data'].min() + 2
    max_lat = 40.95521545410156 #lats['data'].max() - 2
    max_lon = -72.63585662841797 #lons['data'].max() - 2.5
    display = pyart.graph.RadarMapDisplay(radar)
    lat0 = display.loc[0]
    lon0 = display.loc[1]
    
    projection = ccrs.Mercator(
                    central_longitude=lon0,
                    min_latitude=min_lat, max_latitude=max_lat)
                    
    fig = plt.figure(figsize=[8, 8], dpi=100)
    my_ax = plt.subplot(projection = ccrs.Mercator())
    
    display.plot_ppi_map(
        'velocity', 1,
        projection=projection, colorbar_flag=False,
        min_lon=min_lon, max_lon=max_lon, min_lat=min_lat, max_lat=max_lat,
        vmin=-20, vmax=20, cmap=pyart.graph.cm.NWSVel, title_flag=False, lat_lines = None, lon_lines = None, embelish = False, gatefilter=my_ds_gf)
    
    display.ax.set_extent([min_lon, max_lon, min_lat, max_lat])
    display.plot_point(lon0, lat0, color= 'black', markersize=5)
    display.plot_point(-75.7506,39.6780,marker='*',markersize=2,color='white')
    #extent = [min_lon, max_lon, min_lat,max_lat]
    #my_ax.set_extent(extent)
    cbaxes = inset_axes(my_ax, width="100%", height="4%", loc='lower center', borderpad=0) 
    cb1 = plt.colorbar(display.plots[0], orientation='horizontal', cax=cbaxes, ticks=[-15, -10, -5, 0, 5, 10, 15])
    cb1.ax.set_xticklabels(['-15', '-10', '-5', '0', '5', '10', '15'])
    #cb1.set_label('Meshed reflectivity (dBZ)')
    cb1.outline.set_visible(False) # Remove the colorbar outline
    cb1.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb1.ax.xaxis.set_tick_params(pad=-15.5) # Put the colobar labels inside the colorbar
    my_ax.set_title("")
    political_boundaries = cartopy.feature.NaturalEarthFeature(category='cultural',
                                   name='admin_0_boundary_lines_land',
                                   scale='10m', facecolor='none')
    
    states = cartopy.feature.NaturalEarthFeature(category='cultural',
                                   name='admin_1_states_provinces_lines',
                                   scale='10m', facecolor='none')
    
    coast = cartopy.feature.NaturalEarthFeature(category='physical', scale='10m',
                                facecolor='none', name='coastline')
    reader = shpreader.Reader('/home/sat_ops/goesR/zfolder/countyl010g_shp_nt00964/countyl010g.shp')
    counties = list(reader.geometries())
    COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())
    
    my_ax.add_feature(COUNTIES, facecolor='none', edgecolor='darkslategray')
    my_ax.add_feature(political_boundaries, linestyle='-', edgecolor='lightgray', linewidth=1)
    my_ax.add_feature(states, linestyle='-', edgecolor='lightgray',linewidth=1)
    my_ax.add_feature(coast, linestyle='-', edgecolor='lightgray',linewidth=1)
    #request = cimgt.GoogleTiles(url="https://cartodb-basemaps-{s}.global.ssl.fastly.net/dark_all/{z}/{x}/{y}.png")
    request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/dark_nolabels/{z}/{x}/{y}.png")

    clabeltext='Radial Velocity [m/s]'
    title = 'CEMA Base Radial Velocity'
    timestr = local.strftime('%Y-%m-%d %H:%M ') + et
    
    my_ax.add_image(request, 7, zorder=0, interpolation='none')
    # top rectangle
    fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.775,0.03,
                                  fill=True, color='darkslateblue', alpha=1, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
    # bottom rectangle
    fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.775,0.03,
                                  fill=True, color='darkslateblue', alpha=1, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
    fig.text(toptextleft, bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
    fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
    fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=10, zorder=2000)
    gl = display.ax.gridlines(draw_labels=True,
                              linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylabels_left = False

    im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemaunidata24.png")
    plt.figimage(im1, 15, logoy,zorder=3000, alpha=0.8)
    #output_file = workdir + 'vel' + site + '/' + str(dataset) + ".png"
    plt.savefig(imgdir + str(dataset) + '.png', dpi=100, bbox_inches='tight')
    plt.close()

def plot_reflectivity(radar, dataset, imgdir):
    print('Plotting dataset: ', dataset)
    my_gf = pyart.filters.GateFilter(radar)
    my_gf.exclude_above('differential_reflectivity', 4)
    my_gf.exclude_below('reflectivity', 10)
    my_gf.exclude_below('cross_correlation_ratio', .8)
    my_ds_gf = pyart.correct.despeckle_field(radar, 'reflectivity', gatefilter=my_gf)
    timestamp = radar.time['units'].split(' ')[-1].split('T')
    timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
    timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S') + timedelta(minutes=4)
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
    
    lats = radar.gate_latitude
    lons = radar.gate_longitude
    min_lon = -78.24357604980469 #lons['data'].min() + 2.5
    min_lat = 36.69588851928711 #lats['data'].min() + 2
    max_lat = 40.95521545410156 #lats['data'].max() - 2
    max_lon = -72.63585662841797 #lons['data'].max() - 2.5
    display = pyart.graph.RadarMapDisplay(radar)
    lat0 = display.loc[0]
    lon0 = display.loc[1]
    
    projection = ccrs.Mercator(
                    central_longitude=lon0,
                    min_latitude=min_lat, max_latitude=max_lat)
                    
    fig = plt.figure(figsize=[8, 8], dpi=100)
    my_ax = plt.subplot(projection = ccrs.Mercator())
    display.plot_ppi_map(
        'reflectivity', 0,
        projection=projection, colorbar_flag=False,
        min_lon=min_lon, max_lon=max_lon, min_lat=min_lat, max_lat=max_lat,
        vmin=0, vmax=64, cmap=pyart.graph.cm.NWSRef, title_flag=False, lat_lines = None, lon_lines = None, embelish = False, gatefilter=my_ds_gf)
    
    display.ax.set_extent([min_lon, max_lon, min_lat, max_lat])
    display.plot_point(lon0, lat0, color= 'black', markersize=5)
    display.plot_point(-75.7506,39.6780,marker='*',markersize=2,color='white')
    #extent = [min_lon, max_lon, min_lat,max_lat]
    #my_ax.set_extent(extent)
    cbaxes = inset_axes(my_ax, width="100%", height="4%", loc='lower center', borderpad=0) 
    cb1 = plt.colorbar(display.plots[0], orientation='horizontal', cax=cbaxes, ticks=[5, 15, 25, 35, 45, 55, 65])
    cb1.ax.set_xticklabels(['5', '15', '25', '35', '45', '55', '65'])
    #cb1.set_label('Meshed reflectivity (dBZ)')
    cb1.outline.set_visible(False) # Remove the colorbar outline
    cb1.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb1.ax.xaxis.set_tick_params(pad=-15.5) # Put the colobar labels inside the colorbar
    my_ax.set_title("")
    political_boundaries = cartopy.feature.NaturalEarthFeature(category='cultural',
                                   name='admin_0_boundary_lines_land',
                                   scale='10m', facecolor='none')
    
    states = cartopy.feature.NaturalEarthFeature(category='cultural',
                                   name='admin_1_states_provinces_lines',
                                   scale='10m', facecolor='none')
    
    coast = cartopy.feature.NaturalEarthFeature(category='physical', scale='10m',
                                facecolor='none', name='coastline')
    reader = shpreader.Reader('/home/sat_ops/goesR/zfolder/countyl010g_shp_nt00964/countyl010g.shp')
    counties = list(reader.geometries())
    COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())
    
    my_ax.add_feature(COUNTIES, facecolor='none', edgecolor='darkslategray')
    my_ax.add_feature(political_boundaries, linestyle='-', edgecolor='lightgray', linewidth=1)
    my_ax.add_feature(states, linestyle='-', edgecolor='lightgray',linewidth=1)
    my_ax.add_feature(coast, linestyle='-', edgecolor='lightgray',linewidth=1)
    
    #request = cimgt.GoogleTiles(url="https://cartodb-basemaps-{s}.global.ssl.fastly.net/dark_all/{z}/{x}/{y}.png")
    request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/dark_nolabels/{z}/{x}/{y}.png")
    
    
    clabeltext='Reflectivtity [dBZ]'
    title = 'CEMA Base Reflectivity'
    timestr = local.strftime('%Y-%m-%d %H:%M ') + et
    
    my_ax.add_image(request, 7, zorder=0, interpolation='none')
    # top rectangle
    fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.775,0.03,
                                  fill=True, color='darkslateblue', alpha=1, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
    # bottom rectangle
    fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.775,0.03,
                                  fill=True, color='darkslateblue', alpha=1, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
    fig.text(toptextleft, bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
    fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
    fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=10, zorder=2000)
    gl = display.ax.gridlines(draw_labels=True,
                              linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylabels_left = False

    im1 = image.imread("/home/sat_ops/goesR/zfolder/udelcemaunidata24.png")
    plt.figimage(im1, 15, logoy,zorder=3000, alpha=0.8)

    #output_file = workdir + 'vel' + site + '/' + str(dataset) + ".png"
    plt.savefig(imgdir + str(dataset) + '.png', dpi=100, bbox_inches='tight')
    plt.close()

def create_gif(workdir, imgdir, gifname):
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
    imageio.mimsave(workdir + gifname, images, duration=dur_vals)





try:
    site = 'KDOX'
    nowtime = datetime.utcnow().replace(second=0, microsecond=0)
    cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/nexrad/level2/' + site + '/' + str(nowtime.year) + str(nowtime.month).zfill(2) + str(nowtime.day).zfill(2) + '/catalog.xml')
    dataset = list(cat.datasets.values())[1]
    data_ymd = str(dataset).split('_')[2]
    data_hm = str(str(dataset).split('_')[3]).split('.')[0]
    data_datetime = datetime.strptime(data_ymd + data_hm,'%Y%m%d%H%M')
    time_diff = nowtime - data_datetime
    if time_diff.total_seconds() > 1600:
        raise ValueError('Datetimes too far apart, moving to next site')
except:
    pass
    try:
        site = 'KLWX'
        nowtime = datetime.utcnow().replace(second=0, microsecond=0)
        cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/nexrad/level2/' + site + '/' + str(nowtime.year) + str(nowtime.month).zfill(2) + str(nowtime.day).zfill(2) + '/catalog.xml')
        dataset = list(cat.datasets.values())[1]
        data_ymd = str(dataset).split('_')[2]
        data_hm = str(str(dataset).split('_')[3]).split('.')[0]
        data_datetime = datetime.strptime(data_ymd + data_hm,'%Y%m%d%H%M')
        time_diff = nowtime - data_datetime
        if time_diff.total_seconds() > 1600:
            raise ValueError('Datetimes too far apart, moving to next site')
    except:
        pass
        site = 'KDIX'
        nowtime = datetime.utcnow().replace(second=0, microsecond=0)
        cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/nexrad/level2/' + site + '/' + str(nowtime.year) + str(nowtime.month).zfill(2) + str(nowtime.day).zfill(2) + '/catalog.xml')
        dataset = list(cat.datasets.values())[1]
        data_ymd = str(dataset).split('_')[2]
        data_hm = str(str(dataset).split('_')[3]).split('.')[0]
        data_datetime = datetime.strptime(data_ymd + data_hm,'%Y%m%d%H%M')
        time_diff = nowtime - data_datetime
        if time_diff.total_seconds() > 1600:
            raise ValueError('Datetimes too far apart, moving to next site')

workdir = '/home/sat_ops/goesR/radar/'

# defining sizing for plotting stuff

toptext = 0.857
toptextleft = 0.13
toptextright = 0.7
bottomtextleft = 0.13
bottomtextheight = 0.155
toprecx = 0.124
toprecy = 0.849
bottomrecx = 0.124
bottomrecy = 0.145
logoy =  58.8


if os.path.isfile(workdir + 'ref' + site + '/' + str(dataset) + ".png") == False:
    # open the radar data
    radar = pyart.io.read_nexrad_cdm(dataset.access_urls['OPENDAP'])
    
    # plot the reflectivity
    imgdir = workdir + 'ref' + site + '/'
    plot_reflectivity(radar=radar, dataset=dataset, imgdir=imgdir)
    create_gif(workdir=workdir, imgdir=imgdir, gifname="kdox_radar.gif")
    
    # plot the velocity
    imgdir = workdir + 'vel' + site + '/'
    plot_velocity(radar=radar, dataset=dataset, imgdir=imgdir)
    create_gif(workdir=workdir, imgdir=imgdir, gifname="kdox_velocity.gif")

