# deos geopandas to be run at home 
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import geopandas as gpd
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
from datetime import datetime, timedelta, date
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import rioxarray

# declare paths
#shapePaths = "/Users/James/Downloads/mapLayers/"
#colorPaths = "/Users/James/Downloads/colorramps/"
#tiffolder = '/Users/James/Downloads/'
# declare paths
shapePaths = "/home/james/mapLayers/"
colorPaths = "/home/james/colorramps/"
tiffolder = "/home/sat_ops/deos/static_tifs/"
my_dpi = 100
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

# Declare bounds of the data
bounds=(-76.2,38.3,-74.85, 40.3)
# open agwx main dataset
agwx_main = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSAG.nc")
agwx_main = agwx_main.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))

# open NCEP stage IV QC dataset
dsPrec = xr.open_dataset("http://thredds.demac.udel.edu/thredds/dodsC/NCEPIVQC.nc")
dsPrec = dsPrec.sel(lat=slice(bounds[1], bounds[3]), 
                    lon=slice(bounds[0],bounds[2]), 
                    time=slice(datetime.strptime("2014-01-01", "%Y-%m-%d"),
                              date.today()))

dsPrec = dsPrec.reindex(lat=list(reversed(dsPrec.lat)))
dsPrec = dsPrec.rename(name_dict= {'lat' : 'latitude'})
dsPrec = dsPrec.rename(name_dict= {'lon' : 'longitude'})
dsPrec = dsPrec.drop('crs')

sum_dict = dict(zip(['Heating Degree Days', 'Cooling Degree Days','DEOS Precip',  'Energy Density', 
                     'Reference Evapotranspiration', 'Growing Degree Days'],
                    ['HDD', 'CDD',  'dailyprecip', 'energyDens','refET', 'GDD']))
mean_dict = dict(zip(['Mean Temperature', 'Max Temperature', 'Min Temperature', 'Mean Wind Speed', 'Mean Dew Point',
                      'Mean Relative Humidity', 'Max Relative Humidity', 'Min Relative Humidity', 
                      'Mean Soil Temperature', 'Max Soil Temperature', 'Min Soil Temperature',
                      'Mean Volumetric Water Content','Max Volumetric Water Content', 'Min Volumetric Water Content',
                      'Mean Solar', 'Mean Wind Direction', 'Wind Gust', 'Min Wind Chill'],
                     ['meanTemp', 'maxTemp', 'minTemp', 'meanWS','meanDP','meanRH', 'maxRH', 'minRH',
                      'meanST', 'maxST', 'minST','meanVWC', 'maxVWC', 'minVWC','meanSolar', 'meanWD',
                      'dailyGust', 'dailyMinWC']))

max_dict = dict(zip(['Max Temperature','Max Soil Temperature', 'Wind Gust'],
                     ['maxTemp','maxST','dailyGust']))

min_dict = dict(zip(['Min Temperature', 'Min Soil Temperature','Min Wind Chill'],
                     ['minTemp', 'minST','dailyMinWC']))
                     
daysback_dict = dict(zip(['3 Months', '1 Month', '1 Week', '1 Day'], [90, 30, 7, 1]))
datasets = list(sum_dict.keys()) + list(mean_dict.keys()) 
nowdate=datetime.utcnow()
for var in datasets:
    for db in daysback_dict.keys():
        if any(var in s for s in mean_dict.keys()):
            df = agwx_main[mean_dict[var]]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Avg ' + df.units 
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.mean('time')
            dfvarname = "average_" + mean_dict[var]
        if any(var in s for s in sum_dict.keys()):
            df = agwx_main[sum_dict[var]]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Total ' + df.units
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.sum('time')
            dfvarname = "total_" + sum_dict[var]
        if any(var in s for s in max_dict.keys()):
            df = agwx_main[max_dict[var]]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Maximum ' + df.units 
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.max('time')
            dfvarname = "maximum_" + max_dict[var]
        if any(var in s for s in min_dict.keys()):
            df = agwx_main[min_dict[var]]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Minimum ' + df.units
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.min('time')
            dfvarname = "minimum_" + min_dict[var]
        # convert to geotiff so we can clip the extents
        df.rio.set_crs("epsg:4326")
        df.attrs['units'] = 'Fahrenheit'
        df.attrs['standard_name'] = 'Temperature'
        df.rio.set_spatial_dims('longitude', 'latitude')
        df.rio.to_raster(tiffolder + dfvarname  + str(daysback_dict[db]) + '.tif', overwrite=True)
        cl = rioxarray.open_rasterio(tiffolder + dfvarname + str(daysback_dict[db]) +'.tif')
        # clip the interpolated data based on the shapefiles
        #clipped = xds.rio.clip(deos_boundarys.geometry.apply(mapping), xds.rio.crs, all_touched=True,drop=False)
        #cl = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, all_touched=False,invert=True)
        if 'Temp' in dfvarname or 'ST' in dfvarname or 'DP' in dfvarname:
            cl.values[0] = ((cl.values[0] - 273.15)*(9/5)) + 32
            opLabel = dfvarname.split("_")[0] + ' Deg F'
        
        # create time label     
        timeLabel = datetime.strftime(time_recent, "%m-%d-%Y %H:%MZ")

        fig = plt.figure(figsize=(380/my_dpi, 772/my_dpi), dpi=my_dpi)
        ax = fig.add_subplot(111, projection=ccrs.Mercator())
        ax.set_extent([-76.15, -75.03, 38.44, 40.26], crs=ccrs.PlateCarree())
        for ind in range(0,len(bigdeos)):
                ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                              facecolor='silver', edgecolor='black')
        im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap='viridis',transform=ccrs.PlateCarree(),zorder=2)
        for ind in range(0,len(deos_boundarys)):
            ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                              facecolor='none', edgecolor='black', zorder=3, linewidth=1.5)
        for ind in range(0,len(inland_bays)):
            ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                              facecolor='white', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([state_outline['geometry'][74]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([bigdeos['geometry'][121]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([bigdeos['geometry'][120]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([bigdeos['geometry'][74]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([bigdeos['geometry'][75]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([state_outline['geometry'][117]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([state_outline['geometry'][103]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)




        #plt.text(-76.13, 38.523, var,horizontalalignment='left',color='black',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())
        plt.text(-76.11, 38.475, str(var + "\n  " + db + " back from\n " + 
                                     timeLabel),
                                     horizontalalignment='left',color='black',weight='bold',size=5.2,zorder=30,transform=ccrs.PlateCarree())
        #cbaxes = inset_axes(ax, width="3%", height="100%", pad='1.95%', loc=1) 
        cb = fig.colorbar(im, shrink=.7, pad=.02, label=opLabel)
        im1 = image.imread(shapePaths + "deos_logo.png")
        plt.figimage(im1, 18, 50 ,zorder=30, alpha=1)
        plt.savefig("/var/www/html/imagery/AgWx/weather/" + dfvarname + "_" + str(daysback_dict[db]) + ".png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)
        plt.close()

##################################################################
######################  Heat Index ###################### 
###################### ###################### ###################### 
# do the same for Heat Index which must be calculated
daysback_dict = dict(zip(['3 Months', '1 Month', '1 Week', '1 Day'], [90, 30, 7, 1]))
for db in daysback_dict.keys():
    TM = agwx_main['maxTemp']
    df = agwx_main['maxTemp']
    time_recent = pd.to_datetime(TM.time.values[-1])
    TM = TM.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
    df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
    TM.values = ((TM.values - 273.15)*(9/5)) + 32
    df.values = ((df.values - 273.15)*(9/5)) + 32
    TM.values = TM.where(TM.values < 112)
    df.values = df.where(df.values < 112)
    RH = agwx_main['maxRH']
    RH = RH.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
    RH.values = RH.where(RH.values < 101)
    # Heat Index formula grabbed from https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml
    HI = df
    for ts in range(0,len(TM.time.values)):
        HI.values[ts] = np.nan
        HI.values[ts] = 0.5 * (TM.values[ts] + 61.0 + [(TM.values[ts]-68.0)*1.2] + (RH.values[ts]*0.094))
        for r in range(0,HI.values[ts].shape[0]):
            for c in range(0,HI.values[ts].shape[1]):
                val = HI.values[ts][r,c]
                if val > 80:
                    val = -42.379 + 2.04901523*TM.values[ts][r,c] + 10.14333127*RH.values[ts][r,c] - .22475541*TM.values[ts][r,c]*RH.values[ts][r,c] - .00683783*TM.values[ts][r,c]*TM.values[ts][r,c] - .05481717*RH.values[ts][r,c]*RH.values[ts][r,c] + .00122874*TM.values[ts][r,c]*TM.values[ts][r,c]*RH.values[ts][r,c] + .00085282*TM.values[ts][r,c]*RH.values[ts][r,c]*RH.values[ts][r,c] - .00000199*TM.values[ts][r,c]*TM.values[ts][r,c]*RH.values[ts][r,c]*RH.values[ts][r,c]
                if RH.values[ts][r,c] < 13 and TM.values[ts][r,c] > 80 and TM.values[ts][r,c] < 112:
                    val = val - ((13-RH.values[ts][r,c])/4)*np.sqrt((17-np.abs(TM.values[ts][r,c]-95.))/17)
                if RH.values[ts][r,c] > 85 and TM.values[ts][r,c] > 80 and TM.values[ts][r,c] < 87:
                    val = val + ((RH.values[ts][r,c]-85)/10) * ((87-TM.values[ts][r,c])/5)
                HI.values[ts][r,c] = val
    
    df = HI.max('time')
    dfvarname = "maximum_HeatIndex"
    opLabel = 'Max ' + 'Heat Index'
    # convert to geotiff so we can clip the extents
    df.rio.set_crs("epsg:4326")
    df.attrs['units'] = 'Fahrenheit'
    df.attrs['standard_name'] = 'Temperature'
    df.rio.set_spatial_dims('longitude', 'latitude')
    df.rio.to_raster(tiffolder + dfvarname  + str(daysback_dict[db]) + '.tif', overwrite=True)
    cl = rioxarray.open_rasterio(tiffolder + dfvarname + str(daysback_dict[db]) +'.tif')
    
    # create time label     
    timeLabel = datetime.strftime(time_recent, "%m-%d-%Y %H:%MZ")

    fig = plt.figure(figsize=(380/my_dpi, 772/my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(111, projection=ccrs.Mercator())
    ax.set_extent([-76.15, -75.03, 38.44, 40.26], crs=ccrs.PlateCarree())
    for ind in range(0,len(bigdeos)):
            ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                          facecolor='silver', edgecolor='black')
    im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap='viridis',transform=ccrs.PlateCarree(),zorder=2)
    for ind in range(0,len(deos_boundarys)):
        ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='black', zorder=3, linewidth=1.5)
    for ind in range(0,len(inland_bays)):
        ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                          facecolor='white', edgecolor='black',zorder=3, linewidth=1.5)
    ax.add_geometries([state_outline['geometry'][74]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
    ax.add_geometries([bigdeos['geometry'][121]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=1.5)
    ax.add_geometries([bigdeos['geometry'][120]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
    ax.add_geometries([bigdeos['geometry'][74]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
    ax.add_geometries([bigdeos['geometry'][75]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
    ax.add_geometries([state_outline['geometry'][117]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
    ax.add_geometries([state_outline['geometry'][103]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)

    #plt.text(-76.13, 38.523, var,horizontalalignment='left',color='black',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())
    plt.text(-76.11, 38.475, str("Heat Index" + "\n  " + db + " back from\n " + 
                                 timeLabel),
                                 horizontalalignment='left',color='black',weight='bold',size=5.2,zorder=30,transform=ccrs.PlateCarree())
    #cbaxes = inset_axes(ax, width="3%", height="100%", pad='1.95%', loc=1) 
    cb = fig.colorbar(im, shrink=.7, pad=.02, label=opLabel)
    im1 = image.imread(shapePaths + "deos_logo.png")
    plt.figimage(im1, 18, 50 ,zorder=30, alpha=1)
    plt.savefig("/var/www/html/imagery/AgWx/weather/" + dfvarname + "_" + str(daysback_dict[db]) + ".png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)
    plt.close()




datasets = ['Reference Evapotranspiration', 'NCEP Stage IV Precip', 'NCEP Stage IV Precip - DEOS RefET']
daysback_dict = dict(zip(['18 Months', '12 Months', '6 Months', '3 Months', '1 Month', '1 Week', '1 Day'], [540, 360, 180, 90, 30, 7, 1]))
for var in datasets:
    for db in daysback_dict.keys():
        if var == 'Reference Evapotranspiration':
            df = agwx_main['refET']
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Total ' + df.units
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.sum('time')
            dfvarname = "total_refET"
            
        if var == 'NCEP Stage IV Precip':
            time_recent = pd.to_datetime(dsPrec.time.values[-1])
            df = dsPrec
            opLabel = 'Total mm'
            df = df.sel(time=slice(time_recent - timedelta(days=(daysback_dict[db] + 1)), time_recent - timedelta(days=1)))
            df = df.sum('time')
            df = df['Precipitation_Flux']
            #df.values = df.values * (1/0.03937007874) convert to inches when ready to 
            dfvarname = 'ncepIVprecip'

        if var == 'NCEP Stage IV Precip - DEOS RefET':
            time_recent = pd.to_datetime(dsPrec.time.values[-1])
            dfref = agwx_main['refET']
            df1 = dfref.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df1 = df1.sum('time')
            df2 = dsPrec.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df2 = df2.sum('time')
            df = df2 - df1.values
            df['Precip - ET'] = df.Precipitation_Flux
            df = df.drop('Precipitation_Flux')
            df = df['Precip - ET']
            dfvarname = 'ncepIVprecip_deosET'
            opLabel = 'Total mm'
        
        # convert to geotiff so we can clip the extents
        df.rio.set_crs("epsg:4326")
        df.attrs['units'] = 'Fahrenheit'
        df.attrs['standard_name'] = 'Temperature'
        df.rio.set_spatial_dims('longitude', 'latitude')
        df.rio.to_raster(tiffolder + dfvarname  + str(daysback_dict[db]) + '.tif', overwrite=True)
        cl = rioxarray.open_rasterio(tiffolder + dfvarname + str(daysback_dict[db]) +'.tif')

        # create time label     
        timeLabel = datetime.strftime(time_recent, "%m-%d-%Y %H:%MZ")

        fig = plt.figure(figsize=(380/my_dpi, 772/my_dpi), dpi=my_dpi)
        ax = fig.add_subplot(111, projection=ccrs.Mercator())
        ax.set_extent([-76.15, -75.03, 38.44, 40.26], crs=ccrs.PlateCarree())
        for ind in range(0,len(bigdeos)):
                ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                              facecolor='silver', edgecolor='black')
        im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap='viridis',transform=ccrs.PlateCarree(),zorder=2)
        for ind in range(0,len(deos_boundarys)):
            ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                              facecolor='none', edgecolor='black', zorder=3, linewidth=1.5)
        for ind in range(0,len(inland_bays)):
            ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                              facecolor='white', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([state_outline['geometry'][74]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([bigdeos['geometry'][121]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([bigdeos['geometry'][120]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([bigdeos['geometry'][74]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([bigdeos['geometry'][75]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([state_outline['geometry'][117]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)
        ax.add_geometries([state_outline['geometry'][103]], oldproj, facecolor='silver', edgecolor='black',zorder=3, linewidth=1.5)

        #plt.text(-76.13, 38.523, var,horizontalalignment='left',color='black',weight='bold',size=9,zorder=30,transform=ccrs.PlateCarree())
        plt.text(-76.11, 38.475, str(var + "\n  " + db + " back from\n " + 
                                     timeLabel),
                                     horizontalalignment='left',color='black',weight='bold',size=5.2,zorder=30,transform=ccrs.PlateCarree())
        #cbaxes = inset_axes(ax, width="3%", height="100%", pad='1.95%', loc=1) 
        cb = fig.colorbar(im, shrink=.7, pad=.02, label=opLabel)
        im1 = image.imread(shapePaths + "deos_logo.png")
        plt.figimage(im1, 18, 50 ,zorder=30, alpha=1)
        plt.savefig("/var/www/html/imagery/AgWx/water_quantity/" + dfvarname + "_" + str(daysback_dict[db]) + ".png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)
        plt.close()

