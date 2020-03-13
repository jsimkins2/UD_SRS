# deos geopandas to be run at home 
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import geopandas as gpd
from scipy.interpolate import griddata
import rioxarray
import xarray as xr
import pyproj
from pyproj import Proj
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import geopandas
from shapely.geometry import box, mapping
import matplotlib.colors as clr
from matplotlib.colors import BoundaryNorm
import matplotlib as mpl
import pandas as pd
import matplotlib.patheffects as path_effects
import time
import matplotlib.image as image
from datetime import datetime, timedelta
# declare paths
shapePaths = "/Users/James/Downloads/mapLayers/"
colorPaths = "/Users/James/Downloads/colorramps/"
## define custom functions
def remove_nan_observations(x, y, z):
    r"""Remove all x, y, and z where z is nan.
    Will not destroy original values.
    Parameters
    ----------
    x: array_like
        x coordinate
    y: array_like
        y coordinate
    z: array_like
        observation value
    Returns
    -------
    x, y, z
        List of coordinate observation pairs without
        nan valued observations.
    """
    x_ = x[~np.isnan(z)]
    y_ = y[~np.isnan(z)]
    z_ = z[~np.isnan(z)]

    return x_, y_, z_
def check_crs(crs):
    """Checks if the crs represents a valid grid, projection or ESPG string.
    Examples
    --------
    >>> p = check_crs('+units=m +init=epsg:26915')
    >>> p.srs
    '+units=m +init=epsg:26915 '
    >>> p = check_crs('wrong')
    >>> p is None
    True
    Returns
    -------
    A valid crs if possible, otherwise None
    """
    if isinstance(crs, pyproj.Proj) or isinstance(crs, Grid):
        out = crs
    elif isinstance(crs, dict) or isinstance(crs, string_types):
        try:
            out = pyproj.Proj(crs)
        except RuntimeError:
            try:
                out = pyproj.Proj(init=crs)
            except RuntimeError:
                out = None
    else:
        out = None
    return out
def proj_to_cartopy(proj):
    """Converts a pyproj.Proj to a cartopy.crs.Projection
    Parameters
    ----------
    proj: pyproj.Proj
        the projection to convert
    Returns
    -------
    a cartopy.crs.Projection object
    """

    import cartopy.crs as ccrs

    proj = check_crs(proj)

    #if proj.is_latlong():
        #return ccrs.PlateCarree()

    srs = proj.srs

    km_proj = {'lon_0': 'central_longitude',
               'lat_0': 'central_latitude',
               'x_0': 'false_easting',
               'y_0': 'false_northing',
               'k': 'scale_factor',
               'zone': 'zone',
               }
    km_globe = {'a': 'semimajor_axis',
                'b': 'semiminor_axis',
                }
    km_std = {'lat_1': 'lat_1',
              'lat_2': 'lat_2',
              }
    kw_proj = dict()
    kw_globe = dict()
    kw_std = dict()
    for s in srs.split('+'):
        s = s.split('=')
        if len(s) != 2:
            continue
        k = s[0].strip()
        v = s[1].strip()
        try:
            v = float(v)
        except:
            pass
        if k == 'proj':
            if v == 'tmerc':
                cl = ccrs.TransverseMercator
            if v == 'lcc':
                cl = ccrs.LambertConformal
            if v == 'merc':
                cl = ccrs.Mercator
            if v == 'utm':
                cl = ccrs.UTM
        if k in km_proj:
            kw_proj[km_proj[k]] = v
        if k in km_globe:
            kw_globe[km_globe[k]] = v
        if k in km_std:
            kw_std[km_std[k]] = v

    globe = None
    if kw_globe:
        globe = ccrs.Globe(**kw_globe)
    if kw_std:
        kw_proj['standard_parallels'] = (kw_std['lat_1'], kw_std['lat_2'])

    # mercatoooor
    if cl.__name__ == 'Mercator':
        kw_proj.pop('false_easting', None)
        kw_proj.pop('false_northing', None)

    return cl(globe=globe, **kw_proj)

#### custom colormap  
import matplotlib as mpl
startcolor = '#8B4513'
midcolor = '#FFFFFF'
endcolor = '#008000'
own_cmap1 = mpl.colors.LinearSegmentedColormap.from_list( 'own2', [startcolor, midcolor, endcolor] )


######################################################################
######################################################################
# Load in the deos boundaries 
deos_boundarys = gpd.read_file(shapePaths + 'deoscounties.shp')
bigdeos = gpd.read_file(shapePaths + 'TRISTATE_OVERVIEW.shp')
inland_bays = gpd.read_file(shapePaths + 'InlandBays.shp')
state_outline = gpd.read_file(shapePaths + 'tristateMultiaddedPACo.shp')

# create cartopy instance of obscure projection
c = Proj('+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
oldproj = proj_to_cartopy(c)

# Ugly lat/lon lists for plotting purposes 
lons = [-75.437285,-75.842194,-75.682511,-75.062685,-75.838585,-75.750292,
 -75.518358,-75.729665,-75.436712,-75.564711,-75.592728,-75.700788,-75.454989,-75.247235,-75.748172,-75.593015,-75.588603,-75.726915,
 -75.439749,-75.319427,-75.727202,-75.581097,-75.703596,-75.740838,-75.367499,-75.607224,-75.455562,-75.731098,-75.554627,-76.04611,
 -75.777507,-75.726113,-75.615819,-75.716946,-75.988069,-75.118033,-75.639711,-75.234744]
lats = [39.088212,39.864455,39.871675,38.546423,39.709929,39.669478,39.73279,39.037735,38.594724,38.921882,38.541496,39.78871,38.636091,
38.679178,39.391078,39.169629,38.720602,39.605879,38.881661,38.628356,39.819994,39.277231,38.652478,39.24165,38.480304,39.802519,
 39.809165,39.72918,39.939513,39.74167,40.06946,40.164743,40.085846,39.942779,39.923355,38.616095,39.060538,38.542241]
 
######################################################################
timeframes = ["daily", "weekly", "monthly"]
for tf in timeframes:
    if tf == "daily":
        nearday = datetime.utcnow() - timedelta(days=2)
        farday = datetime.utcnow() - timedelta(days=2)
        vmin=-10
        vmax=10
    if tf == "weekly":
        nearday = datetime.utcnow() - timedelta(days=2)
        farday = datetime.utcnow() - timedelta(days=9)
        vmin=-30
        vmax=30
    if tf == "monthly":
        nearday = datetime.utcnow() - timedelta(days=2)
        farday = datetime.utcnow() - timedelta(days=32)
        vmin=-100
        vmax=100
    
    nearstr = str(str(nearday.year) + "-" + "{:02d}".format(nearday.month) + "-" + "{:02d}".format(nearday.day))
    farstr = str(str(farday.year) + "-" + "{:02d}".format(farday.month) + "-" + "{:02d}".format(farday.day))
    # load in the ET data
    et_data = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSRefET.nc")
    et_data = et_data.sel(time=slice(nearstr, farstr))
    et_data.refET.values = et_data.refET.values * 0.1 #convert to cm 
    #prec_data = xr.open_dataset("")
        
        
    ## plot ET data 
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())
    ax.set_extent([-76.15, -75.03, 38.44, 40.26], crs=ccrs.PlateCarree())
    for ind in range(0,len(bigdeos)):
            ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                          facecolor='silver', edgecolor='black')
    im=ax.pcolormesh(et_data['longitude'].values,et_data['latitude'].values,et_data.refET.values[0],cmap=own_cmap1,transform=ccrs.PlateCarree(),zorder=2, vmin=vmin, vmax=vmax)
    
    for l in range(0,len(lons)):
        if lons[l] != -76.35 and lons[l] != -74.68 and lons[l] != -75.062685 and lons[l] != -75.118033 and lons[l] != -75.640685 and lons[l] != -75.527755 and lons[l] != -75.118033 and lons[l] != -75.148629 and lons[l] != -75.727202:
                text = plt.text(lons[l],lats[l],str(round(et_data.refET.sel(longitude=lons[l], latitude=lats[l], method="nearest").values[0], 2)), size=9,weight='bold',verticalalignment='center',
                horizontalalignment='center',transform=ccrs.PlateCarree(),zorder=5)
                text.set_path_effects([path_effects.Stroke(linewidth=4, foreground='white'),path_effects.Normal()])
    for ind in range(0,len(deos_boundarys)):
        ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='gray', zorder=3, linewidth=1.5)
    for ind in range(0,len(inland_bays)):
        ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                          facecolor='white', edgecolor='black',zorder=3, linewidth=1.5)
    ax.add_geometries([state_outline['geometry'][74]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=1.5)
    ax.add_geometries([bigdeos['geometry'][121]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=1.5)
    plt.title("(Precipitation - Ref ET) Proof of Concept")
    cbar = plt.colorbar(im)
    cbar.set_label('cm day-1')
    plt.text(-76.02, 38.503, 'Ref. ET',horizontalalignment='left',color='black',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())
    plt.text(-76.1, 38.473, str(farstr + " to "),horizontalalignment='left',color='black',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())
    plt.text(-76.1, 38.45, str(nearstr),horizontalalignment='left',color='black',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())

    im1 = image.imread("Downloads/maplayers/deos_logo.png")
    plt.figimage(im1, 350, 140 ,zorder=30, alpha=1)
    plt.savefig("Downloads/P_et.png")



