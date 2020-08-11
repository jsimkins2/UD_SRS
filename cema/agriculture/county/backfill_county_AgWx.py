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
from datetime import datetime, timedelta, date

# declare paths
shapePaths = "/home/james/mapLayers/"
colorPaths = "/home/james/colorramps/"
raster_path = "/home/sat_ops/deos/temp/"
my_dpi = 100
## define custom functions
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

# read in deos special shapefiles
deos_boundarys = gpd.read_file(shapePaths + 'deoscounties.shp')
bigdeos = gpd.read_file(shapePaths + 'TRISTATE_OVERVIEW.shp')
inland_bays = gpd.read_file(shapePaths + 'InlandBays.shp')
state_outline = gpd.read_file(shapePaths + 'tristateMultiaddedPACo.shp')

# create cartopy instance of obscure projection
c = Proj('+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
oldproj = proj_to_cartopy(c)


# load in agwx dataset 
bounds=(-76.2,38.3,-74.85, 40.3)

agwx_main = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSAG.nc")
agwx_main = agwx_main.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))

dsPrec = xr.open_dataset("http://thredds.demac.udel.edu/thredds/dodsC/NCEPIVQC.nc")
dsPrec = dsPrec.sel(lat=slice(bounds[1], bounds[3]), 
                    lon=slice(bounds[0],bounds[2]), 
                    time=slice(datetime.strptime("2010-01-01", "%Y-%m-%d"),
                              date.today()))
dsPrec = dsPrec.reindex(lat=list(reversed(dsPrec.lat)))
dsPrec = dsPrec.rename(name_dict= {'lat' : 'latitude'})
dsPrec = dsPrec.rename(name_dict= {'lon' : 'longitude'})
dsPrec = dsPrec.drop('crs')

# create county data frames
chester_ds = xr.Dataset({"time": agwx_main['time'].values})
ncc_ds = xr.Dataset({"time": agwx_main['time'].values})
kent_ds = xr.Dataset({"time": agwx_main['time'].values})
sussex_ds = xr.Dataset({"time": agwx_main['time'].values})

# each geotiff would need to be a lone time slice...going to look into geopandas update
for co in range(0, len(deos_boundarys["NAME"])):
    county_outline = deos_boundarys.loc[[co], 'geometry']
    print(co)
    for var in agwx_main.data_vars:
        var_list = []
        print(var)
        for t in range(0,len(agwx_main.time.values)):
            da = xr.DataArray(agwx_main[var][t].values,dims=['latitude', 'longitude'],coords={'longitude': agwx_main.longitude.values, 'latitude' :agwx_main.latitude.values})
            da.rio.set_crs("epsg:4326")
            da.attrs['units'] = 'Fahrenheit'
            da.attrs['standard_name'] = 'Temperature'
            da.rio.set_spatial_dims('longitude', 'latitude')
            da.rio.to_raster(raster_path + var + str(co) + str(t) + '.tif', overwrite=True)
            xds = rioxarray.open_rasterio(raster_path + var + str(co) + str(t) + '.tif')
            # clip the interpolated data based on the shapefiles
            clipped = xds.rio.clip(county_outline.geometry.apply(mapping), xds.rio.crs, drop=True)
            cl = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, invert=True)
            ds_county = cl.mean()
            var_list.append(round(ds_county.values.tolist(),2))
            # if we don't remove it, it won't overwrite properly
            os.system("/bin/rm " + raster_path + var + str(co) + str(t) + '.tif')
        if co == 0:
            chester_ds[var] = (['time'], var_list)
        if co == 1:
            ncc_ds[var] = (['time'], var_list)
        if co == 2:
            kent_ds[var] = (['time'], var_list)
        if co == 3:
            sussex_ds[var] = (['time'], var_list)
for co in range(0, len(deos_boundarys["NAME"])):
    county_outline = deos_boundarys.loc[[co], 'geometry']
    var = 'Precipitation_Flux'
    var_list = []
    for t in range(0,len(agwx_main.time.values)):
        da = xr.DataArray(dsPrec.Precipitation_Flux.sel(time = agwx_main.time.values[t], method = 'nearest').values,dims=['latitude', 'longitude'],coords={'longitude': agwx_main.longitude.values, 'latitude' :agwx_main.latitude.values})
        da.rio.set_crs("epsg:4326")
        da.attrs['units'] = 'Fahrenheit'
        da.attrs['standard_name'] = 'Temperature'
        da.rio.set_spatial_dims('longitude', 'latitude')
        da.rio.to_raster(raster_path + var + str(co) + str(t) + '.tif', overwrite=True)
        xds = rioxarray.open_rasterio(raster_path + var + str(co) + str(t) + '.tif')
        # clip the interpolated data based on the shapefiles
        clipped = xds.rio.clip(county_outline.geometry.apply(mapping), xds.rio.crs, drop=True)
        cl = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, invert=True)
        ds_county = cl.mean()
        var_list.append(round(ds_county.values.tolist(),2))
        # if we don't remove it, it won't overwrite properly
        os.system("/bin/rm " + raster_path + var + str(co) + str(t) + '.tif')
    if co == 0:
        chester_ds['NCEPstageIVPrecip'] = (['time'], var_list)
        chester_ds['NCEPstageIVPrecip'].attrs['units'] = dsPrec[var].attrs['units']
        chester_ds['NCEPstageIVPrecip'].attrs['long_name'] = dsPrec[var].attrs['long_name']
    if co == 1:
        ncc_ds['NCEPstageIVPrecip'] = (['time'], var_list)
        ncc_ds['NCEPstageIVPrecip'].attrs['units'] = dsPrec[var].attrs['units']
        ncc_ds['NCEPstageIVPrecip'].attrs['long_name'] = dsPrec[var].attrs['long_name']
    if co == 2:
        kent_ds['NCEPstageIVPrecip'] = (['time'], var_list)
        kent_ds['NCEPstageIVPrecip'].attrs['units'] = dsPrec[var].attrs['units']
        kent_ds['NCEPstageIVPrecip'].attrs['long_name'] = dsPrec[var].attrs['long_name']
    if co == 3:
        sussex_ds['NCEPstageIVPrecip'] = (['time'], var_list)
        sussex_ds['NCEPstageIVPrecip'].attrs['units'] = dsPrec[var].attrs['units']
        sussex_ds['NCEPstageIVPrecip'].attrs['long_name'] = dsPrec[var].attrs['long_name']


for var in chester_ds.variables:
    if var != 'time':
        if var != 'NCEPstageIVPrecip' :
            chester_ds[var].attrs['units'] = agwx_main[var].attrs['units']
            ncc_ds[var].attrs['units'] = agwx_main[var].attrs['units']
            kent_ds[var].attrs['units'] = agwx_main[var].attrs['units']
            sussex_ds[var].attrs['units'] = agwx_main[var].attrs['units']
            chester_ds[var].attrs['long_name'] = agwx_main[var].attrs['long_name']
            ncc_ds[var].attrs['long_name'] = agwx_main[var].attrs['long_name']
            kent_ds[var].attrs['long_name'] = agwx_main[var].attrs['long_name']
            sussex_ds[var].attrs['long_name'] = agwx_main[var].attrs['long_name']

chester_ds.to_netcdf("/data/DEOS/chester/chester_agwx_updated.nc", mode='w')
print('finished chester county')

ncc_ds.to_netcdf("/data/DEOS/ncc/ncc_agwx_updated.nc", mode='w')
print('finished new castle county')

kent_ds.to_netcdf("/data/DEOS/kent/kent_agwx_updated.nc", mode='w')
print('finished kent county')

sussex_ds.to_netcdf("/data/DEOS/sussex/sussex_agwx_updated.nc", mode='w')
print('finished sussex county')

chester_agwx.close()
ncc_agwx.close()
kent_agwx.close()
sussex_agwx.close()

os.system("/bin/mv /data/DEOS/chester/chester_agwx_updated.nc /data/DEOS/chester/chester_agwx.nc")
os.system("/bin/mv /data/DEOS/ncc/ncc_agwx_updated.nc /data/DEOS/ncc/ncc_agwx.nc")
os.system("/bin/mv /data/DEOS/kent/kent_agwx_updated.nc /data/DEOS/kent/kent_agwx.nc")
os.system("/bin/mv /data/DEOS/sussex/sussex_agwx_updated.nc /data/DEOS/sussex/sussex_agwx.nc")

