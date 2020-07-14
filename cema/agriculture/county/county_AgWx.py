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
import matplotlib.image as image
import time
from datetime import datetime, timedelta

# declare paths
shapePaths = "/home/james/mapLayers/"
colorPaths = "/home/james/colorramps/"
my_dpi = 100
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
def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))
    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256*4)
    return cmap
def linear_rbf(x, y, z, xi, yi):
    dist = distance_matrix(x,y, xi,yi)
    # Mutual pariwise distances between observations
    internal_dist = distance_matrix(x,y, x,y)
    # Now solve for the weights such that mistfit at the observations is minimized
    weights = np.linalg.solve(internal_dist, z)
    # Multiply the weights for each interpolated point by the distances
    zi =  np.dot(dist.T, weights)
    return zi

def distance_matrix(x0, y0, x1, y1):
    obs = np.vstack((x0, y0)).T
    interp = np.vstack((x1, y1)).T
    # Make a distance matrix between pairwise observations
    # Note: from <http://stackoverflow.com/questions/1871536>
    # (Yay for ufuncs!)
    d0 = np.subtract.outer(obs[:,0], interp[:,0])
    d1 = np.subtract.outer(obs[:,1], interp[:,1])
    return np.hypot(d0, d1)


# read in deos special shapefiles
deos_boundarys = gpd.read_file(shapePaths + 'deoscounties.shp')
bigdeos = gpd.read_file(shapePaths + 'TRISTATE_OVERVIEW.shp')
inland_bays = gpd.read_file(shapePaths + 'InlandBays.shp')
state_outline = gpd.read_file(shapePaths + 'tristateMultiaddedPACo.shp')

# create cartopy instance of obscure projection
c = Proj('+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
oldproj = proj_to_cartopy(c)


# load in agwx dataset 

# load in agwx dataset 
bounds=(-76.2,38.3,-74.85, 40.3)

agwx_main = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSAG.nc")
agwx_main = agwx_main.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))

# create county data frames
chester_ds = xr.Dataset({"time": agwx_main['time'].values})
ncc_ds = xr.Dataset({"time": agwx_main['time'].values})
kent_ds = xr.Dataset({"time": agwx_main['time'].values})
sussex_ds = xr.Dataset({"time": agwx_main['time'].values})

# each geotiff would need to be a lone time slice...going to look into geopandas update
for co in range(0, len(deos_boundarys["NAME"])):
    county_outline = deos_boundarys.loc[[co], 'geometry']

    for var in agwx_main.data_vars:
        var_list = []
        for t in range(0,len(agwx_main.time.values)):
            da = xr.DataArray(agwx_main[var][t].values,dims=['latitude', 'longitude'],coords={'longitude': agwx_main.longitude.values, 'latitude' :agwx_main.latitude.values})
            da.rio.set_crs("epsg:4326")
            da.attrs['units'] = 'Fahrenheit'
            da.attrs['standard_name'] = 'Temperature'
            da.rio.set_spatial_dims('longitude', 'latitude')
            da.rio.to_raster('/Users/james/Downloads/' + var + '.tif', overwrite=True)
        
            xds = rioxarray.open_rasterio('/Users/james/Downloads/' + var +'.tif')
            # clip the interpolated data based on the shapefiles
            clipped = xds.rio.clip(county_outline.geometry.apply(mapping), xds.rio.crs, drop=True)
            cl = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, invert=True)
            ds_county = cl.mean("x")
            ds_county = ds_county.mean("y")
            var_list.append(ds_county.values[0])
        
        if co == 0:
            chester_ds[var] = (['time'], var_list)
        if co == 1:
            ncc_ds[var] = (['time'], var_list)
        if co == 2:
            kent_ds[var] = (['time'], var_list)
        if co == 3:
            sussex_ds[var] = (['time'], var_list)

    if co == 0:
        chester_ds.to_netcdf("")
        print('finished chester county')
    if co == 1:
        ncc_ds.to_netcdf("")
        print('finished new castle county')
    if co == 2:
        kent_ds.to_netcdf("")
        print('finished kent county')
    if co == 3:
        sussex_ds.to_netcdf("")
        print('finished sussex county')
