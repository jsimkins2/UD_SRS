# static ag weather maps hosted on /var/www/html/imagery/AgWx/
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
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
import pandas as pd
import matplotlib.patheffects as path_effects
import time
import matplotlib.image as image
from datetime import datetime, timedelta, date
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import rioxarray
from matplotlib.offsetbox import AnchoredText

# declare paths
#shapePaths = "/Users/James/Downloads/mapLayers/"
#colorPaths = "/Users/James/Downloads/colorramps/"
#tiffolder = '/Users/James/Downloads/'
# declare paths
shapePaths = "/home/james/mapLayers/"
colorPaths = "/home/james/colorramps/"
tiffolder = "/home/sat_ops/deos/static_tifs/"
my_dpi = 100

# define own colorbar
startcolor = '#8B4513'
midcolor = '#FFFFFF'
endcolor = '#008000'
own_cmap1 = mpl.colors.LinearSegmentedColormap.from_list( 'own2', [startcolor, midcolor, endcolor] )

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

sum_dict = dict(zip(['Heating Degree Days', 'Cooling Degree Days','DEOS Precip',  'Energy Density', 'Growing Degree Days'],
                    ['HDD', 'CDD',  'dailyprecip', 'energyDens', 'GDD']))
mean_dict = dict(zip(['Mean Temperature', 'Maxa Temperature', 'Mina Temperature', 'Mean Wind Speed', 'Mean Dew Point',
                      'Mean Relative Humidity', 'Max Relative Humidity', 'Min Relative Humidity', 
                      'Mean Soil Temperature', 'Maxa Soil Temperature', 'Mina Soil Temperature',
                      'Mean Volumetric Water Content','Max Volumetric Water Content', 'Min Volumetric Water Content',
                      'Mean Solar', 'Mean Wind Direction', 'Winda Gust', 'Mina Wind Chill', 'Dailya Max HI'],
                     ['meanTemp', 'maxTemp', 'minTemp', 'meanWS','meanDP','meanRH', 'maxRH', 'minRH',
                      'meanST', 'maxST', 'minST','meanVWC', 'maxVWC', 'minVWC','meanSolar', 'meanWD',
                      'dailyGust', 'dailyMinWC', 'maxHI']))

max_dict = dict(zip(['Max Temperature','Max Soil Temperature', 'Wind Gust', 'Daily Max HI'],
                     ['maxTemp','maxST','dailyGust', 'maxHI']))

min_dict = dict(zip(['Min Temperature', 'Min Soil Temperature','Min Wind Chill'],
                     ['minTemp', 'minST','dailyMinWC']))

nowtime = datetime.utcnow()
ytd = pd.to_datetime(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d")) -  pd.to_datetime(nowtime)
daysback_dict = dict(zip(['YTD', '3 Months', '1 Month', '1 Week', '1 Day'], [np.int(np.abs(ytd.days)), 90, 30, 7, 1]))
datasets = list(sum_dict.keys()) + list(mean_dict.keys()) + list(max_dict.keys()) + list(min_dict.keys())

nowdate=datetime.utcnow()
for var in datasets:
    for db in daysback_dict.keys():
        if any(var in s for s in mean_dict.keys()):
            df = agwx_main[mean_dict[var]]
            dfvarname = "average_" + mean_dict[var]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Avg ' + df.units 
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.mean('time')
            cmap = 'coolwarm'
            if var[3] == 'a':
                var_plot_name = var[0:3] + var[4:]
            elif var[4] == 'a':
                var_plot_name = var[0:4] + var[5:]
            elif var[5] == 'a':
                var_plot_name = var[0:5] + var[6:]
            else:
                var_plot_name = var
            
        if any(var in s for s in sum_dict.keys()):
            df = agwx_main[sum_dict[var]]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Total ' + df.units
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.sum('time')
            dfvarname = "total_" + sum_dict[var]
            var_plot_name = var
            if var == 'Cooling Degree Days':
                cmap = 'Spectral'
            else:
                cmap = 'Spectral_r'
            
        if any(var in s for s in max_dict.keys()):
            df = agwx_main[max_dict[var]]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Maximum ' + df.units 
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.max('time')
            dfvarname = "maximum_" + max_dict[var]
            var_plot_name = var
            if max_dict[var] == 'maxHI':
                dfvarname = "maximum_HeatIndex"
            cmap = 'coolwarm'
            
        if any(var in s for s in min_dict.keys()):
            df = agwx_main[min_dict[var]]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Minimum ' + df.units
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.min('time')
            dfvarname = "minimum_" + min_dict[var]
            var_plot_name = var
            cmap = 'coolwarm'
            # convert to geotiff so we can clip the extents
        df.rio.set_crs("epsg:4326")
        df.attrs['units'] = 'Fahrenheit'
        df.attrs['standard_name'] = 'Temperature'
        df.rio.set_spatial_dims('longitude', 'latitude')
        df.rio.to_raster(tiffolder + dfvarname  + str(daysback_dict[db]) + '.tif', overwrite=True)
        cl = rioxarray.open_rasterio(tiffolder + dfvarname + str(daysback_dict[db]) +'.tif')

        if 'Temp' in dfvarname or 'ST' in dfvarname or 'HeatIndex' in dfvarname or 'DP' in dfvarname or 'HI' in dfvarname:
            cl.values[0] = ((cl.values[0] - 273.15)*(9/5)) + 32
            opLabel = dfvarname.split("_")[0] + ' (Deg F)'
        
        if 'Gust' in dfvarname or 'WS' in dfvarname:
            cl.values[0] = cl.values[0] * 2.23694
            opLabel = dfvarname.split("_")[0] + ' (mph)'
        
        if 'precip' in dfvarname:
            cl.values[0] = cl.values[0] * 0.0393701
            opLabel = 'Total (inches)'
            
        if 'HeatIndex' in dfvarname or 'HI' in dfvarname:
            cl.values[0][cl.values[0] < 80] = np.nan
        # create time label     
        timeLabel = datetime.strftime(time_recent, "%m-%d-%Y %H:%MZ")
        
        fig = plt.figure(figsize=(380/my_dpi, 772/my_dpi), dpi=my_dpi)
        ax = fig.add_subplot(111, projection=ccrs.Mercator())
        ax.set_extent([-76.15, -75.03, 38.44, 40.26], crs=ccrs.PlateCarree())
        for ind in range(0,len(bigdeos)):
                ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                              facecolor='silver', edgecolor='black')
        im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap=cmap,transform=ccrs.PlateCarree(),zorder=2)
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

        if db == 'YTD':
            plt.text(-76.11, 38.475, str(var_plot_name + "\n" + 'Year To Date'),horizontalalignment='left',color='black',weight='bold',size=5.2,zorder=30,transform=ccrs.PlateCarree())
            cb = fig.colorbar(im, shrink=.7, pad=.02, label=opLabel)
            im1 = image.imread(shapePaths + "deos_logo.png")
            plt.figimage(im1, 18, 50 ,zorder=30, alpha=1)
            plt.savefig("/var/www/html/imagery/AgWx/weather/" + dfvarname + "_YTD.png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)
        else:
            plt.text(-76.11, 38.475, str(var_plot_name + "\n" + db + " back from\n" + 
                             timeLabel),
                             horizontalalignment='left',color='black',weight='bold',size=5.2,zorder=30,transform=ccrs.PlateCarree())
            cb = fig.colorbar(im, shrink=.7, pad=.02, label=opLabel)
            im1 = image.imread(shapePaths + "deos_logo.png")
            plt.figimage(im1, 18, 50 ,zorder=30, alpha=1)
            plt.savefig("/var/www/html/imagery/AgWx/weather/" + dfvarname + "_" + str(daysback_dict[db]) + ".png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)
        plt.close()



datasets = ['Reference Evapotranspiration', 'dailyprecip', 'NCEP Stage IV Precip', 'NCEP Stage IV Precip - DEOS RefET']
daysback_dict = dict(zip(['18 Months', '12 Months', 'YTD', '6 Months', '3 Months', '1 Month', '1 Week', '1 Day'], [540, 360,np.int(np.abs(ytd.days)), 180, 90, 30, 7, 1]))
cmap = 'BrBG'
for var in datasets:
    for db in daysback_dict.keys():
        if var == 'Reference Evapotranspiration':
            df = agwx_main['refET']
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Total ' + df.units
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.sum('time')
            dfvarname = "total_refET"
        
        if var =='dailyprecip':
            df = agwx_main['dailyprecip']
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Total ' + df.units
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.sum('time')
            dfvarname = "DEOSprecip"
            
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
        
        # convert to inches
        cl.values[0] = cl.values[0] * 0.0393701
        opLabel = 'Total (inches)'
        
        fig = plt.figure(figsize=(380/my_dpi, 772/my_dpi), dpi=my_dpi)
        ax = fig.add_subplot(111, projection=ccrs.Mercator())
        ax.set_extent([-76.15, -75.03, 38.44, 40.26], crs=ccrs.PlateCarree())
        for ind in range(0,len(bigdeos)):
                ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                              facecolor='silver', edgecolor='black')
        im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap=cmap,transform=ccrs.PlateCarree(),zorder=2)
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

        if db == 'YTD':
            plt.text(-76.11, 38.475, str(var + "\n" + 'Year To Date'),horizontalalignment='left',color='black',weight='bold',size=5.2,zorder=30,transform=ccrs.PlateCarree())
            cb = fig.colorbar(im, shrink=.7, pad=.02, label=opLabel)
            im1 = image.imread(shapePaths + "deos_logo.png")
            plt.figimage(im1, 18, 50 ,zorder=30, alpha=1)
            plt.savefig("/var/www/html/imagery/AgWx/water_quantity/" + dfvarname + "_YTD.png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)
        else:
            plt.text(-76.11, 38.475, str(var + "\n" + db + " back from\n" + 
                                         timeLabel),
                                         horizontalalignment='left',color='black',weight='bold',size=5.2,zorder=30,transform=ccrs.PlateCarree())
            cb = fig.colorbar(im, shrink=.7, pad=.02, label=opLabel)
            im1 = image.imread(shapePaths + "deos_logo.png")
            plt.figimage(im1, 18, 50 ,zorder=30, alpha=1)
            plt.savefig("/var/www/html/imagery/AgWx/water_quantity/" + dfvarname + "_" + str(daysback_dict[db]) + ".png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)

        plt.close()
        

# now create the departure maps 
# load in the datasets
climo = xr.open_dataset("/data/DEOS/doy_climatology/deos_doy_climatology.nc")
nowtime = datetime.utcnow()
ytd = pd.to_datetime(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d")) -  pd.to_datetime(nowtime)
daysback_dict = dict(zip(['YTD', '3 Months', '1 Month', '1 Week', '1 Day'], [np.int(np.abs(ytd.days)), 90, 30, 7, 1]))
datasets = list(sum_dict.keys()) + list(mean_dict.keys()) + ['Reference Evapotranspiration', 'NCEP Stage IV Precip']

nowdate=datetime.utcnow()
for var in datasets:
    for db in daysback_dict.keys():
        if any(var in s for s in mean_dict.keys()):
            df = agwx_main[mean_dict[var]]
            dfvarname = "Average_Difference_" + mean_dict[var]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Difference from Normal (' + df.units + ')'
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.mean('time')
            if 'HeatIndex' in dfvarname or 'HI' in dfvarname:
                df = df.where(df > 299.817)
            cf = climo[mean_dict[var]]
            cf_min = (time_recent.timetuple().tm_yday - daysback_dict[db]) if (time_recent.timetuple().tm_yday - daysback_dict[db]) > 0 else 0
            cf = cf.isel(dayofyear=slice(cf_min, time_recent.timetuple().tm_yday))
            cf = cf.mean('dayofyear')
            if 'HeatIndex' in dfvarname or 'HI' in dfvarname:
                cf = cf.where(cf > 299.817)
            df = df - cf
            cmap = 'coolwarm'
            if var[3] == 'a':
                var_plot_name = var[0:3] + var[4:]
            elif var[4] == 'a':
                var_plot_name = var[0:4] + var[5:]
            elif var[5] == 'a':
                var_plot_name = var[0:5] + var[6:]
            else:
                var_plot_name = var
            

        if any(var in s for s in sum_dict.keys()):
            df = agwx_main[sum_dict[var]]
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Difference from Normal (' + df.units + ')'
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.sum('time')
            dfvarname = "total_difference_" + sum_dict[var]
            cf = climo[sum_dict[var]]
            cf_min = (time_recent.timetuple().tm_yday - daysback_dict[db]) if (time_recent.timetuple().tm_yday - daysback_dict[db]) > 0 else 0
            cf = cf.isel(dayofyear=slice(cf_min, time_recent.timetuple().tm_yday))
            cf = cf.sum('dayofyear')
            df = df - cf
            var_plot_name = var
            if var == 'Cooling Degree Days':
                cmap = 'Spectral'
            if var == 'DEOS Precip':
                cmap = 'BrBG'
                df.values = df.values * (0.03937007874)
                dfvarname = "total_difference_DEOSPrecip"
                opLabel = 'Difference from Normal (inches)'
            else:
                cmap = 'Spectral_r'

        if var == 'Reference Evapotranspiration':
            df = agwx_main['refET']
            time_recent = pd.to_datetime(df.time.values[-1])
            opLabel = 'Difference from Normal (' + df.units + ')'
            df = df.sel(time=slice(time_recent - timedelta(days=daysback_dict[db]), time_recent))
            df = df.sum('time')
            cf = climo['refET']
            cf_min = (time_recent.timetuple().tm_yday - daysback_dict[db]) if (time_recent.timetuple().tm_yday - daysback_dict[db]) > 0 else 0
            cf = cf.isel(dayofyear=slice(cf_min, time_recent.timetuple().tm_yday))
            cf = cf.sum('dayofyear')
            df = df - cf
            dfvarname = "total_difference_refET"
            var_plot_name = var
            
        if var == 'NCEP Stage IV Precip':
            time_recent = pd.to_datetime(dsPrec.time.values[-1])
            df = dsPrec
            opLabel = 'Difference from Normal (inches)'
            df = df.sel(time=slice(time_recent - timedelta(days=(daysback_dict[db] + 1)), time_recent - timedelta(days=1)))
            df = df.sum('time')
            df = df['Precipitation_Flux']
            #df.values = df.values * (1/0.03937007874) convert to inches when ready to
            cf = climo['NCEPstageIVPrecip']
            cf_min = (time_recent.timetuple().tm_yday - daysback_dict[db]) if (time_recent.timetuple().tm_yday - daysback_dict[db]) > 0 else 0
            cf = cf.isel(dayofyear=slice(cf_min, time_recent.timetuple().tm_yday))
            cf = cf.sum('dayofyear')
            df = df - cf.values
            df.values = df.values * (0.03937007874)
            dfvarname = 'total_difference_ncepIVprecip'
            var_plot_name = var
            cmap = 'BrBG'
            
        # convert to geotiff so we can clip the extents
        df.rio.set_crs("epsg:4326")
        df.attrs['units'] = 'Fahrenheit'
        df.attrs['standard_name'] = 'Temperature'
        df.rio.set_spatial_dims('longitude', 'latitude')
        df.rio.to_raster(tiffolder + dfvarname  + str(daysback_dict[db]) + '.tif', overwrite=True)
        cl = rioxarray.open_rasterio(tiffolder + dfvarname + str(daysback_dict[db]) +'.tif')

        if 'Temp' in dfvarname or 'ST' in dfvarname or 'HeatIndex' in dfvarname or 'DP' in dfvarname or 'HI' in dfvarname:
            cl.values[0] = ((cl.values[0] - 273.15)*(9/5))
            opLabel = 'Difference (Deg F)'
        
        if 'Gust' in dfvarname or 'Speed' in dfvarname:
            cl.values[0] = cl.values[0] * 2.23694
            opLabel = 'Difference (mph)'
    
        # create time label     
        timeLabel = datetime.strftime(time_recent, "%m-%d-%Y %H:%MZ")
        tem_vmin = cl.min()
        tem_vmax = cl.max()
        if abs(tem_vmin) >= abs(tem_vmax):
            vmin = tem_vmin
            vmax = (-1)*tem_vmin
        if abs(tem_vmax) > abs(tem_vmin):
            vmax = tem_vmax
            vmin = (-1)*tem_vmax

        fig = plt.figure(figsize=(380/my_dpi, 772/my_dpi), dpi=my_dpi)
        ax = fig.add_subplot(111, projection=ccrs.Mercator())
        ax.set_extent([-76.15, -75.03, 38.44, 40.26], crs=ccrs.PlateCarree())
        for ind in range(0,len(bigdeos)):
                ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                              facecolor='silver', edgecolor='black')
        im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap=cmap,transform=ccrs.PlateCarree(),zorder=2,
                         vmin=vmin, vmax=vmax)
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

        if db == 'YTD':
            plt.text(-76.11, 38.475, str(var_plot_name + "\n" + 'Year To Date\n' + 'Difference'),horizontalalignment='left',color='black',weight='bold',size=5.2,zorder=30,transform=ccrs.PlateCarree())
            cb = fig.colorbar(im, shrink=.7, pad=.02, label=opLabel)
            im1 = image.imread(shapePaths + "deos_logo.png")
            plt.figimage(im1, 18, 50 ,zorder=30, alpha=1)
            plt.savefig("/var/www/html/imagery/AgWx/departures/" + dfvarname + "_YTD.png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)
        else:
            plt.text(-76.11, 38.475, str(var_plot_name + "\n" + db + " Difference\n" + 
                             timeLabel),
                             horizontalalignment='left',color='black',weight='bold',size=5.2,zorder=30,transform=ccrs.PlateCarree())
            cb = fig.colorbar(im, shrink=.7, pad=.02, label=opLabel)
            im1 = image.imread(shapePaths + "deos_logo.png")
            plt.figimage(im1, 18, 50 ,zorder=30, alpha=1)
            plt.savefig("/var/www/html/imagery/AgWx/departures/" + dfvarname + "_" + str(daysback_dict[db]) + ".png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)
        plt.close()

# cumulative county maps
# open up the county agwx datasets
# list of dataframes
countydf = ['chester_agwx.nc', 'ncc_agwx.nc', 'kentc_agwx.nc', 'sussex_agwx.nc']
doydf = ['chester_agwx_climatology.nc', 'ncc_agwx_climatology.nc', 'kent_agwx_climatology.nc', 'sussex_agwx_climatology.nc']
co_names = ['Chester', 'New Castle', 'Kent', 'Sussex']
nowtime = datetime.utcnow()
ytd = pd.to_datetime(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d")) -  pd.to_datetime(nowtime)
daysback_dict = dict(zip(['YTD'], [np.int(np.abs(ytd.days))]))#dict(zip(['YTD', '3 Months', '1 Month', '1 Week', '1 Day'], [np.int(np.abs(ytd.days)), 90, 30, 7, 1]))

# Precipitation Map
nowdate=datetime.utcnow()
for db in daysback_dict.keys():
    for co in range(0,len(countydf)):
        df = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['dailyprecip']
        time_recent = pd.to_datetime(df.time.values[-1])
        df = df.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        df.values = df.values * (0.0393701)
        cf = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['dailyprecip']
        cf = cf.isel(dayofyear=slice(0, 366))
        cf.values = cf.values * (0.0393701)
        ncep = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['NCEPstageIVPrecip']
        ncep = ncep.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        ncep.values = ncep.values * (0.0393701)
        ncepClim = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['NCEPstageIVPrecip']
        ncepClim.values = ncepClim.values * (0.0393701)
        ncepClim = ncepClim.isel(dayofyear=slice(0, 366))
        #Set X range here:
        left = date(nowdate.year, 1, 1)  #Makes it easy to quickly change the range
        right = date(nowdate.year, 12, 31)
        datelist = pd.date_range(str(str(nowdate.year) + "-01-01"), str(str(nowdate.year) + "-12-31")).tolist()

        fig = plt.figure(figsize=(12,8))
        # Create subplot of 
        plt.subplot(2, 1, 1)
        plt.plot(datelist, np.nancumsum(cf.values), linestyle='dashed', c="lightseagreen",  label="Climatology Cum. Sum")
        plt.plot(df.time.values, np.nancumsum(df.values), c="darkgreen",  label="Observed Cum. Sum")
        plt.plot(df.time.values, df.values, c="darkgreen", label="Observed (in/day)", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values) > 0 else "-"
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + " (in)")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Precipitation')
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Precipitation Totals (in)', labelpad = 2)
        plt.legend()
        plt.grid()
        fig.tight_layout()
        #Create subplot of ncep
        plt.subplot(2,1,2)
        plt.plot(datelist, np.nancumsum(ncepClim.values), linestyle='dashed', c="deepskyblue",  label="Climatology")
        plt.plot(ncep.time.values, np.nancumsum(ncep.values), c="blue",  label="Observed Cum. Sum")
        plt.plot(ncep.time.values, ncep.values, c="blue", label="Observed (in/day)", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + " (in)")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - NCEP Stage IV Precipitation')
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Precipitaiton Totals (in)', labelpad = 2)
        plt.legend()
        plt.grid()
        im1 = image.imread(shapePaths + "deos_logo.png")
        plt.figimage(im1, 1450, 75 ,zorder=30, alpha=1)
        fig.tight_layout()
        plt.savefig("/var/www/html/imagery/AgWx/county/" + co_names[co] + "_YTD_precipitation.png",bbox_inches='tight',pad_inches = 0.1,dpi=my_dpi*1.3)


# HDD and CDD 
nowdate=datetime.utcnow()
for db in daysback_dict.keys():
    for co in range(0,len(countydf)):
        df = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['HDD']
        time_recent = pd.to_datetime(df.time.values[-1])
        df = df.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        cf = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['HDD']
        cf = cf.isel(dayofyear=slice(0, 366))
        ncep = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['CDD']
        ncep = ncep.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        ncepClim = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['CDD']
        ncepClim = ncepClim.isel(dayofyear=slice(0, 366))
        #Set X range here:
        left = date(nowdate.year, 1, 1)  #Makes it easy to quickly change the range
        right = date(nowdate.year, 12, 31)
        datelist = pd.date_range(str(str(nowdate.year) + "-01-01"), str(str(nowdate.year) + "-12-31")).tolist()


        fig = plt.figure(figsize=(12,8))
        # Create subplot of 
        plt.subplot(2, 1, 1)
        plt.plot(datelist, np.nancumsum(cf.values), linestyle='dashed', c="coral",  label="Climatology")
        plt.plot(df.time.values, np.nancumsum(df.values), c="maroon",  label="Observed Cum. Sum")
        plt.plot(df.time.values, df.values, c="maroon", label="Observed", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + "")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Heating Degree Days')
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Total Days', labelpad = 2)
        plt.legend()
        plt.grid()
        fig.tight_layout()
        #Create subplot of ncep
        plt.subplot(2,1,2)
        plt.plot(datelist, np.nancumsum(ncepClim.values), linestyle='dashed', c="deepskyblue",  label="Climatology")
        plt.plot(ncep.time.values, np.nancumsum(ncep.values), c="blue",  label="Observed Cum. Sum")
        plt.plot(ncep.time.values, ncep.values, c="blue", label="Observed", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.round(np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + "")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Cooling Degree Days')
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Total Days', labelpad = 1)
        plt.legend()
        plt.grid()
        im1 = image.imread(shapePaths + "deos_logo.png")
        plt.figimage(im1, 1450, 75 ,zorder=30, alpha=1)
        fig.tight_layout()
        plt.savefig("/var/www/html/imagery/AgWx/county/" + co_names[co] + "_YTD_HDD_CDD.png",bbox_inches='tight',pad_inches = 0.1,dpi=my_dpi*1.3)

# GDD and Energy Density
nowdate=datetime.utcnow()
for db in daysback_dict.keys():
    for co in range(0,len(countydf)):
        df = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['GDD']
        time_recent = pd.to_datetime(df.time.values[-1])
        df = df.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        cf = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['GDD']
        cf = cf.isel(dayofyear=slice(0, 366))
        ncep = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['energyDens']
        ncep = ncep.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        ncepClim = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['energyDens']
        ncepClim = ncepClim.isel(dayofyear=slice(0, 366))
        #Set X range here:
        left = date(nowdate.year, 1, 1)  #Makes it easy to quickly change the range
        right = date(nowdate.year, 12, 31)
        datelist = pd.date_range(str(str(nowdate.year) + "-01-01"), str(str(nowdate.year) + "-12-31")).tolist()


        fig = plt.figure(figsize=(12,8))
        # Create subplot of 
        plt.subplot(2, 1, 1)
        plt.plot(datelist, np.nancumsum(cf.values), linestyle='dashed', c="lightseagreen",  label="Climatology")
        plt.plot(df.time.values, np.nancumsum(df.values), c="darkgreen",  label="Observed Cum. Sum")
        plt.plot(df.time.values, df.values, c="darkgreen", label="Observed", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + "")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Growing Degree Days', pad=5)
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Total Days', labelpad = 4)
        plt.legend()
        plt.grid()
        fig.tight_layout()
        #Create subplot of ncep
        plt.subplot(2,1,2)
        plt.plot(datelist, np.nancumsum(ncepClim.values), linestyle='dashed', c="deepskyblue",  label="Climatology")
        plt.plot(ncep.time.values, np.nancumsum(ncep.values), c="blue",  label="Observed Cum. Sum")
        plt.plot(ncep.time.values, ncep.values, c="blue", label="Observed", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.round(np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round((np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values))/(1e9),2)) + " (J^9)")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Energy Density', pad=5)
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Total Energy (J)', labelpad = 2)
        plt.legend()
        plt.grid()
        im1 = image.imread(shapePaths + "deos_logo.png")
        plt.figimage(im1, 1450, 75 ,zorder=30, alpha=1)
        fig.tight_layout()
        plt.savefig("/var/www/html/imagery/AgWx/county/" + co_names[co] + "_YTD_GDD_Energy.png",bbox_inches='tight',pad_inches = 0.1,dpi=my_dpi*1.3)

