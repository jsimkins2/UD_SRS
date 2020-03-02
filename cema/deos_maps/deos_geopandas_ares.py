# deos geopandas to be run at home
import os
from datetime import datetime, timedelta
stat = os.stat("/home/deos/web_push/map_data.json")
mod_time = datetime.utcfromtimestamp(stat.st_mtime)
nowtime = datetime.utcnow()

diff_time = nowtime - mod_time

# if the time difference between now and the map_data.json is less than 70 seconds than update the maps
# the cronjob will be set to every minute
if diff_time.seconds < 70:
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

    # declare paths
    shapePaths = "/home/map_maker/deos_maps/mapLayers/"
    colorPaths = "/home/map_maker/deos_maps/colorramps/"
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

    # grab data from json files
    deos_data = pd.read_json("http://128.175.28.202/deos_json/map_data.json")
    loc_deos = pd.read_json("http://128.175.28.202/deos_json/station_metadata.json")
    # index is the station nunumbers
    station_id = list(deos_data.index)
    # need to find time value
    date_deos = deos_data.columns[0]
    dst = "EST" if time.localtime().tm_isdst==0 else "EDT"
    zuluDIFF = 5 if dst=='EST' else 4
    date_deos_est = date_deos - timedelta(hours=zuluDIFF)
    deos_dateSTR = str("{0:0=2d}".format(date_deos_est.month) + '/' + "{0:0=2d}".format(date_deos_est.day) + '/' + str(date_deos_est.year) + ' ' + "{0:0=2d}".format(date_deos_est.hour) + ':' + "{0:0=2d}".format(date_deos_est.minute) + ' ' + dst)

    # create a dict of station IDs with the name of the station
    station_dict = {}
    for s in list(loc_deos.columns):
        station_dict[s] = loc_deos[s]['station_id']

    rev_station_dict = dict(zip(station_dict.values(),station_dict.keys()))

    nameDict = dict(zip(['Gage Precipitation (60)','Air Temperature','Dew Point Temperature','Wind Speed','Wind Direction','Barometric Pressure','Relative humidity', '24 Hour Precipitation', 'Peak Wind Gust Speed (60)'], ['precip', 'airT', 'dewP', 'wspeed', 'wdir', 'pres', 'rh', '24hr_precip', 'gust']))
    fancyDict = dict(zip(list(nameDict.keys()), ['1-hr Rain (in)', 'Air Temperature (F)', 'Dew Point (F)', '5-min Wind Speed', 'Wind', 'Pressure (mb)', 'Relative Humidity (%)', '24-hr Rain (in)', '1-hr Peak Wind Gust (mph)']))
    for var in list(nameDict.keys()):
        lats=list()
        lons=list()
        for key in station_id:
            lats.append(loc_deos[rev_station_dict[str(key)]]['latitude'])
            lons.append(loc_deos[rev_station_dict[str(key)]]['longitude'])

        # add in four corners to expand the interpolated grid
        lons = lons + list([-76.35,-76.35, -74.68,  -74.68])
        lats = lats + list([38.0, 40.6, 38.0, 40.6])
        if var == 'Wind Direction' or var == 'Wind Chill':
            pass
        else:
            # target grid to interpolate to
            temp = list()
            # Begin the for loop to place a location on a map
            for s in station_id:
                if str(s) in list(station_dict.values()):
                    ID = str(s)
                try:
                    temp.append(float(deos_data[date_deos][s][var]))
                except:
                    temp.append(np.nan)

            # just in case one of the gauges is out
            try:
                t1 = float(deos_data[date_deos][2321][var])
            except:
                pass
                try:
                    t1 = float(deos_data[date_deos][2932][var])
                except:
                    t1 = np.nanmean(temp)

            try:
                t2 = float(deos_data[date_deos][2999][var])
            except:
                pass
                try:
                    t2 = float(deos_data[date_deos][2982][var])
                except:
                    t2 = np.nanmean(temp)

            try:
                t3 = float(deos_data[date_deos][2304][var])
            except:
                pass
                try:
                    t3 = float(deos_data[date_deos][2747][var])
                except:
                    t3 = np.nanmean(temp)

            try:
                t4 = float(deos_data[date_deos][2983][var])
            except:
                pass
                try:
                    t4 = float(deos_data[date_deos][2984][var])
                except:
                    t4 = np.nanmean(temp)

            temp = temp + list([t1,t2,t3,t4])

            lons=np.array(lons)
            lats=np.array(lats)
            temp = np.array(temp)

            if var == 'Gage Precipitation (60)':
                temp = temp/25.4
                temp[temp < 0] = np.nan
                temp[temp > 15] = np.nan
                vmin=0
                vmax=12
                rounder = 2
                txt_cmap =  pd.read_csv(colorPaths + 'rain_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')

            if var == '24 Hour Precipitation':
                temp = temp/25.4
                temp[temp < 0] = np.nan
                temp[temp > 15] = np.nan
                vmin=0
                vmax=12
                rounder = 2
                txt_cmap =  pd.read_csv(colorPaths + 'rain_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')

            if var == 'Air Temperature':
                temp = temp - 273.15
                temp = (temp*(9/5)) + 32
                temp[temp < -30] = np.nan
                temp[temp > 120] = np.nan
                airT_array = temp
                vmin=-30
                vmax=120
                rounder = 1
                txt_cmap =  pd.read_csv(colorPaths + 'at_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')

            if var == 'Dew Point Temperature':
                temp = temp - 273.15
                temp = (temp*(9/5)) + 32
                temp[temp < 0] = np.nan
                temp[temp > 90] = np.nan
                dewP_array = temp
                vmin=-40
                vmax=120
                rounder = 1
                txt_cmap =  pd.read_csv(colorPaths + 'at_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')

            if var == 'Barometric Pressure':
                vmin=940
                vmax=1050
                temp[temp < 930] = np.nan
                temp[temp > 1100] = np.nan
                rounder = 0
                txt_cmap =  pd.read_csv(colorPaths + 'bp_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')

            if var == 'Relative humidity':
                vmin=0
                vmax=100
                temp[temp < 0] = np.nan
                temp[temp > 100] = np.nan
                rounder = 0
                txt_cmap =  pd.read_csv(colorPaths + 'rh_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')

            if var == 'Wind Speed':
                temp = temp * 2.237
                temp[temp < 0] = np.nan
                temp[temp > 90] = np.nan
                vmin=0
                vmax=52
                rounder = 1
                txt_cmap =  pd.read_csv(colorPaths + 'ws_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')

            if var == 'Peak Wind Gust Speed (60)':
                temp = temp * 2.237
                temp[temp < 0] = np.nan
                temp[temp > 90] = np.nan
                vmin=0
                vmax=120
                rounder = 1
                txt_cmap =  pd.read_csv(colorPaths + 'ws_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')

            lons,lats, temp = remove_nan_observations(lons,lats, temp)
            x = np.linspace(min(lons), max(lons), 750)
            y = np.linspace(min(lats), max(lats), 750)
            xi,yi = np.meshgrid(x,y)
            # interpolate
            #zi = griddata((lons,lats),temp,(xi,yi),method='cubic')
            # try the idw interpolation scheme
            xi, yi = xi.flatten(), yi.flatten()

            # Calculate IDW
            zi = linear_rbf(lons,lats,temp,xi,yi)
            zi=zi.reshape((len(x), len(y)))

            da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y})
            da.rio.set_crs("epsg:4326")
            da.attrs['units'] = 'Fahrenheit'
            da.attrs['standard_name'] = 'Temperature'
            da.rio.set_spatial_dims('lon', 'lat')
            da.rio.to_raster('/home/map_maker/deos_maps/temp/' + nameDict[var] + '.tif', overwrite=True)

            xds = rioxarray.open_rasterio('/home/map_maker/deos_maps/temp/' + nameDict[var] +'.tif')
            # clip the interpolated data based on the shapefiles
            clipped = xds.rio.clip(deos_boundarys.geometry.apply(mapping), xds.rio.crs, drop=True)
            cl = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, invert=True)



            ### Call the function make_cmap which returns your colormap
            raw_rgb = []
            for i in range(0,len(txt_cmap)):
                raw_rgb.append(tuple([txt_cmap['r'][i], txt_cmap['g'][i], txt_cmap['b'][i]]))
            cmap= make_cmap(raw_rgb, bit=True)
            bounds = []

            for i in range(1,len(list(txt_cmap['bound']))):
                lin=np.linspace(list(txt_cmap['bound'])[i-1],list(txt_cmap['bound'])[i], 40)
                bounds.extend(list(lin))

            norm = BoundaryNorm(bounds,ncolors=cmap.N)
            #or lons[l] != -75.913585

            fig = plt.figure(figsize=(380/my_dpi, 772/my_dpi), dpi=my_dpi)
            ax = fig.add_subplot(111, projection=ccrs.Mercator())
            ax.set_extent([-76.35, -75, 38.44, 40.26], crs=ccrs.PlateCarree())
            for ind in range(0,len(bigdeos)):
                    ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                                  facecolor='silver', edgecolor='black')
            im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap=cmap,norm=norm,transform=ccrs.PlateCarree(),zorder=2)
            for l in range(0,len(lons)):
                if var == 'Relative humidity':
                    if lons[l] != -76.35 and lons[l] != -74.68 and lons[l] != -75.062685 and lons[l] != -75.118033 and lons[l] != -75.247235 and lons[l] != -75.640685 and lons[l] != -75.727202:
                        text = plt.text(lons[l],lats[l],str(int(round(temp[l], rounder))), size=6.5,weight='bold',transform=ccrs.PlateCarree(),zorder=7)
                        text.set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='white'),path_effects.Normal()])
                if var == 'Barometric Pressure':
                    if lons[l] != -75.7311 and lons[l] != 75.6108 and lons[l] != -75.2472 and lons[l] != -75.118033 and lons[l] != -76.35 and lons[l] != -74.68 and lons[l] != -75.076665 and lons[l] != -75.118033 and lons[l] != -75.247235 and lons[l] != -75.640685 and lons[l] != -75.527755 and lons[l] != -75.682511 and lons[l] != -75.727202:
                        text = plt.text(lons[l],lats[l],str(int(round(temp[l], rounder))), size=5.8,weight='bold',verticalalignment='center',
                        horizontalalignment='center',transform=ccrs.PlateCarree(),zorder=5)
                        text.set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='white'),path_effects.Normal()])
                if var == 'Gage Precipitation (60)' or var == '24 Hour Precipitation':
                    if lons[l] != -76.35 and lons[l] != -74.68 and lons[l] != -75.062685 and lons[l] != -75.118033 and lons[l] != -75.247235 and lons[l] != -75.640685 and lons[l] != -75.527755 and lons[l] != -75.118033 and lons[l] != -75.148629 and lons[l] != -75.727202:
                        text = plt.text(lons[l],lats[l],str('{:.2f}'.format(round(temp[l], rounder))), size=6.5,weight='bold',verticalalignment='center',
                        horizontalalignment='center',transform=ccrs.PlateCarree(),zorder=5)
                        text.set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='white'),path_effects.Normal()])
                if var == 'Air Temperature' or var == 'Dew Point Temperature' or var == 'Wind Speed' or var == 'Peak Wind Gust Speed (60)':
                    if lons[l] != -76.35 and lons[l] != -74.68 and lons[l] != -75.062685 and lons[l] != -75.118033 and lons[l] != -75.247235 and lons[l] != -75.640685 and lons[l] != -75.527755 and lons[l] != -75.118033 and lons[l] != -75.148629 and lons[l] != -75.727202:
                        text = plt.text(lons[l],lats[l],str('{:.1f}'.format(round(temp[l], rounder))), size=6.5,weight='bold',verticalalignment='center',
                        horizontalalignment='center',transform=ccrs.PlateCarree(),zorder=5)
                        text.set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='white'),path_effects.Normal()])
            for ind in range(0,len(deos_boundarys)):
                ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                                  facecolor='none', edgecolor='gray', zorder=3, linewidth=0.5)
            for ind in range(0,len(inland_bays)):
                ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                                  facecolor='white', edgecolor='black',zorder=3, linewidth=0.5)
            ax.add_geometries([state_outline['geometry'][74]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=0.5)
            ax.add_geometries([bigdeos['geometry'][121]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=0.5)
            #plt.title(nameDict[var])
            if var == 'Peak Wind Gust Speed (60)':
                plt.text(-76.13, 38.503, fancyDict[var],horizontalalignment='left',weight='bold',color='white',size=5,zorder=30,transform=ccrs.PlateCarree())
                plt.text(-76.13, 38.473, deos_dateSTR,horizontalalignment='left',weight='bold',color='white',size=6,zorder=30,transform=ccrs.PlateCarree())
            else:
                plt.text(-76.13, 38.503, fancyDict[var],horizontalalignment='left',weight='bold',color='white',size=6,zorder=30,transform=ccrs.PlateCarree())
                plt.text(-76.13, 38.473, deos_dateSTR,horizontalalignment='left',weight='bold',color='white',size=6,zorder=30,transform=ccrs.PlateCarree())

            im1 = image.imread(shapePaths + "deos_logo.png")
            plt.figimage(im1, 24, 40 ,zorder=30, alpha=1)
            plt.savefig("/home/map_maker/deos_maps/imagery/deos_" + nameDict[var] + ".png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)



    var = 'Wind Direction'
    lats=list()
    lons=list()
    for key in station_id:
        lats.append(loc_deos[rev_station_dict[str(key)]]['latitude'])
        lons.append(loc_deos[rev_station_dict[str(key)]]['longitude'])

    # add in four corners to expand the interpolated grid
    lons = lons + list([-76.35,-76.35, -74.68,  -74.68])
    lats = lats + list([38.3, 40.3, 38.3, 40.3])

    temp = list()
    # Begin the for loop to place a location on a map
    for s in station_id:
        if str(s) in list(station_dict.values()):
            ID = str(s)
        try:
            temp.append(float(deos_data[date_deos][s][var]))
        except:
            temp.append(np.nan)

    # just in case one of the gauges is out
    try:
        t1 = float(deos_data[date_deos][2321][var])
    except:
        t1 = np.nanmean(temp)

    try:
        t2 = float(deos_data[date_deos][2980][var])
    except:
        t2 = np.nanmean(temp)

    try:
        t3 = float(deos_data[date_deos][2304][var])
    except:
        t3 = np.nanmean(temp)

    try:
        t4 = float(deos_data[date_deos][2983][var])
    except:
        t4 = np.nanmean(temp)

    temp = temp + list([t1,t2,t3,t4])

    lons=np.array(lons)
    lats=np.array(lats)
    temp = np.array(temp)
    temp[temp < 0] = np.nan
    temp[temp > 360] = np.nan
    # alright, now inteprolate wind speed values to put underneath wind directions
    ws = list()
    # Begin the for loop to place a location on a map
    for s in station_id:
        if str(s) in list(station_dict.values()):
            ID = str(s)
        try:
            ws.append(float(deos_data[date_deos][s]['Wind Speed']))
        except:
            ws.append(np.nan)

    # just in case one of the gauges is out
    try:
        t1 = float(deos_data[date_deos][2321]['Wind Speed'])
    except:
        t1 = np.nanmean(ws)

    try:
        t2 = float(deos_data[date_deos][2980]['Wind Speed'])
    except:
        t2 = np.nanmean(ws)

    try:
        t3 = float(deos_data[date_deos][2304]['Wind Speed'])
    except:
        t3 = np.nanmean(ws)

    try:
        t4 = float(deos_data[date_deos][2983]['Wind Speed'])
    except:
        t4 = np.nanmean(ws)


    ws = ws + list([t1,t2,t3,t4])

    lons=np.array(lons)
    lats=np.array(lats)
    ws = np.array(ws)
    ws=ws * 2.237
    # define other stuff
    txt_cmap =  pd.read_csv(colorPaths + 'ws_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')
    ws[ws < 0] = np.nan
    ws[ws > 90] = np.nan

    ws_with_nans = ws
    ws_w_nans = ws
    lons_with_nans=lons
    lats_with_nans=lats
    lons,lats, ws = remove_nan_observations(lons,lats, ws)
    x = np.linspace(min(lons), max(lons), 750)
    y = np.linspace(min(lats), max(lats), 750)
    xi,yi = np.meshgrid(x,y)
    # interpolate
    #zi = griddata((lons,lats),temp,(xi,yi),method='cubic')
    # try the idw interpolation scheme
    xi, yi = xi.flatten(), yi.flatten()

    # Calculate IDW
    zi = linear_rbf(lons,lats,ws,xi,yi)
    zi=zi.reshape((len(x), len(y)))

    da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y})
    da.rio.set_crs("epsg:4326")
    da.attrs['units'] = 'Fahrenheit'
    da.attrs['standard_name'] = 'Temperature'
    da.rio.set_spatial_dims('lon', 'lat')
    da.rio.to_raster('/home/map_maker/deos_maps/temp/' + nameDict[var] + '.tif', overwrite=True)

    xds = rioxarray.open_rasterio('/home/map_maker/deos_maps/temp/' + nameDict[var] +'.tif')
    # clip the interpolated data based on the shapefiles
    clipped = xds.rio.clip(deos_boundarys.geometry.apply(mapping), xds.rio.crs, drop=True)
    cl = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, invert=True)




    ### Call the function make_cmap which returns your colormap
    raw_rgb = []
    for i in range(0,len(txt_cmap)):
        raw_rgb.append(tuple([txt_cmap['r'][i], txt_cmap['g'][i], txt_cmap['b'][i]]))
    cmap= make_cmap(raw_rgb, bit=True)
    bounds = []

    for i in range(1,len(list(txt_cmap['bound']))):
        lin=np.linspace(list(txt_cmap['bound'])[i-1],list(txt_cmap['bound'])[i], 40)
        bounds.extend(list(lin))

    norm = BoundaryNorm(bounds,ncolors=cmap.N)
    def wind_components(wspd, wdir):
        rad = 4.0*np.arctan(1)/180.
        u = -wspd*np.sin(rad*wdir)
        v = -wspd*np.cos(rad*wdir)
        return u, v

    # now remove the lon values that we don't want
    ws_with_nans[(lons_with_nans == -76.35) | (lons_with_nans == -74.68) | (lons_with_nans == -75.062685)  | (lons_with_nans == -75.640685)] = np.nan
    temp[(lons_with_nans == -76.35) | (lons_with_nans == -74.68) | (lons_with_nans == -75.062685)   | (lons_with_nans == -75.640685)] = np.nan
    lats_with_nans[(lons_with_nans == -76.35) | (lons_with_nans == -74.68) | (lons_with_nans == -75.062685) |  (lons_with_nans == -75.640685)] = np.nan
    lons_with_nans[(lons_with_nans == -76.35) | (lons_with_nans == -74.68) | (lons_with_nans == -75.062685) |  (lons_with_nans == -75.640685)] = np.nan

    ws_with_nans = ws_with_nans[~np.isnan(temp)]
    lats_with_nans = lats_with_nans[~np.isnan(temp)]
    lons_with_nans = lons_with_nans[~np.isnan(temp)]
    temp = temp[~np.isnan(temp)]

    mws = np.ma.masked_where(ws_with_nans < 1, ws_with_nans)
    mtemp = np.ma.masked_where(ws_with_nans < 1, temp)
    mlats = np.ma.masked_where(ws_with_nans < 1, lats_with_nans)
    mlons = np.ma.masked_where(ws_with_nans < 1, lons_with_nans)

    # calculate the u v components of the wind
    u,v = wind_components(mws.compressed(), mtemp.compressed())
    normu = u / np.sqrt(u**2 + v**2)
    normv = v / np.sqrt(u**2 + v**2)


    fig = plt.figure(figsize=(380/my_dpi, 772/my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(111, projection=ccrs.Mercator())
    ax.set_extent([-76.35, -75, 38.44, 40.26], crs=ccrs.PlateCarree())

    for ind in range(0,len(bigdeos)):
            ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                          facecolor='silver', edgecolor='black')
    im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap=cmap,norm=norm,transform=ccrs.PlateCarree(),zorder=2)
    ax.quiver(mlons.compressed(),mlats.compressed(),normu,normv,transform=ccrs.PlateCarree(),zorder=5,scale_units='xy', alpha=0.8, headwidth=4.5, pivot='middle')
    plt.plot(np.ma.masked_array(lons_with_nans, mask=~mlons.mask),np.ma.masked_array(lats_with_nans, ~mlats.mask), 'ko',fillstyle='none',
             mew=1.3, markersize=5.5,transform=ccrs.PlateCarree(),zorder=5,alpha=0.8)

    for ind in range(0,len(deos_boundarys)):
        ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='gray', zorder=3, linewidth=0.5)
    for ind in range(0,len(inland_bays)):
        ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                          facecolor='white', edgecolor='black',zorder=3, linewidth=0.5)
    ax.add_geometries([state_outline['geometry'][74]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=0.5)
    ax.add_geometries([bigdeos['geometry'][121]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=0.5)
    plt.text(-76.13, 38.503, 'Wind',horizontalalignment='left',color='white',weight='bold',size=6,zorder=30,transform=ccrs.PlateCarree())
    plt.text(-76.13, 38.473, deos_dateSTR,horizontalalignment='left',weight='bold',color='white',size=6,zorder=30,transform=ccrs.PlateCarree())

    im1 = image.imread(shapePaths + "deos_logo.png")
    plt.figimage(im1, 24, 40 ,zorder=30, alpha=1)
    plt.savefig("/home/map_maker/deos_maps/imagery/deos_wdir.png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)


    ############################## feels like mapping ######################################
    ############################## feels like mapping ######################################
    ############################## feels like mapping ######################################
    ############################## feels like mapping ######################################

    var = 'Wind Chill'
    temp = list()
    # Begin the for loop to place a location on a map
    for s in station_id:
        if str(s) in list(station_dict.values()):
            ID = str(s)
        try:
            temp.append(float(deos_data[date_deos][s][var]))
        except:
            temp.append(np.nan)

    # just in case one of the gauges is out
    try:
        t1 = float(deos_data[date_deos][2321][var])
    except:
        t1 = np.nanmean(temp)

    try:
        t2 = float(deos_data[date_deos][2980][var])
    except:
        t2 = np.nanmean(temp)

    try:
        t3 = float(deos_data[date_deos][2304][var])
    except:
        t3 = np.nanmean(temp)

    try:
        t4 = float(deos_data[date_deos][2983][var])
    except:
        t4 = np.nanmean(temp)

    temp = temp + list([t1,t2,t3,t4])
    lons=np.array(lons)
    lats=np.array(lats)
    temp = np.array(temp)
    temp = temp - 273.15
    temp = (temp*(9/5)) + 32
    temp[temp < -50] == np.nan
    temp[temp > 40] == np.nan

    ### now read in heat index if it exists
    var = 'Heat Index'
    lats=list()
    lons=list()
    for key in station_id:
        lats.append(loc_deos[rev_station_dict[str(key)]]['latitude'])
        lons.append(loc_deos[rev_station_dict[str(key)]]['longitude'])

    # add in four corners to expand the interpolated grid
    lons = lons + list([-76.35,-76.35, -74.68,  -74.68])
    lats = lats + list([38.3, 40.3, 38.3, 40.3])

    temp2 = list()
    # Begin the for loop to place a location on a map
    for s in station_id:
        if str(s) in list(station_dict.values()):
            ID = str(s)
        try:
            temp2.append(float(deos_data[date_deos][s][var]))
        except:
            temp2.append(np.nan)

    # just in case one of the gauges is out
    try:
        t1 = float(deos_data[date_deos][2321][var])
    except:
        t1 = np.nanmean(temp2)

    try:
        t2 = float(deos_data[date_deos][2980][var])
    except:
        t2 = np.nanmean(temp2)

    try:
        t3 = float(deos_data[date_deos][2304][var])
    except:
        t3 = np.nanmean(temp2)

    try:
        t4 = float(deos_data[date_deos][2983][var])
    except:
        t4 = np.nanmean(temp2)

    temp2 = temp2 + list([t1,t2,t3,t4])
    lons=np.array(lons)
    lats=np.array(lats)
    temp2 = np.array(temp2)

    temp2 = temp2 - 273.15
    temp2 = (temp2*(9/5)) + 32
    temp2[temp2 < 80] == np.nan
    temp2[temp2 > 130] == np.nan
    vmin=-30
    vmax=120
    rounder = 1
    txt_cmap =  pd.read_csv(colorPaths + 'at_ramp.txt', header=None,names=['bound', 'r', 'g', 'b', 'a'],delimiter=' ')

    # combine heat index and wind chill arrays into one masked array
    m1 = np.ma.masked_invalid(temp)
    m2 = np.ma.masked_invalid(temp2)
    feelsLike = np.ma.filled(m1,m2)

    # replace certain feels like values with air temperature given conditions
    for i in range(0,len(feelsLike)-4):
        if np.isnan(feelsLike[i]) == True:
            if np.isnan(dewP_array[i]) == False:
                if np.isnan(airT_array[i]) == False:
                    if ws_w_nans[i] >= 0:
                        feelsLike[i] = airT_array[i]

    for i in range(len(feelsLike)-4,len(feelsLike)):
        if np.isnan(feelsLike[i]) == True:
            feelsLike[i] = airT_array[i]


    # remove the nan values of feelsLike from the applicable arrays
    lats = lats[~np.isnan(feelsLike)]
    lons = lons[~np.isnan(feelsLike)]
    feelsLike = feelsLike[~np.isnan(feelsLike)]

    x = np.linspace(min(lons), max(lons), 750)
    y = np.linspace(min(lats), max(lats), 750)
    xi,yi = np.meshgrid(x,y)
    # interpolate
    #zi = griddata((lons,lats),temp,(xi,yi),method='cubic')
    # try the idw interpolation scheme
    xi, yi = xi.flatten(), yi.flatten()

    # Calculate IDW
    zi = linear_rbf(lons,lats,feelsLike,xi,yi)
    zi=zi.reshape((len(x), len(y)))

    da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y})
    da.rio.set_crs("epsg:4326")
    da.attrs['units'] = 'Fahrenheit'
    da.attrs['standard_name'] = 'Temperature'
    da.rio.set_spatial_dims('lon', 'lat')
    da.rio.to_raster('/home/map_maker/deos_maps/temp/FeelsLike.tif', overwrite=True)

    xds = rioxarray.open_rasterio('/home/map_maker/deos_maps/temp/FeelsLike.tif')
    # clip the interpolated data based on the shapefiles
    clipped = xds.rio.clip(deos_boundarys.geometry.apply(mapping), xds.rio.crs, drop=True)
    cl = clipped.rio.clip(inland_bays.geometry.apply(mapping), oldproj.proj4_init, drop=False, invert=True)



    ### Call the function make_cmap which returns your colormap
    raw_rgb = []
    for i in range(0,len(txt_cmap)):
        raw_rgb.append(tuple([txt_cmap['r'][i], txt_cmap['g'][i], txt_cmap['b'][i]]))
    cmap= make_cmap(raw_rgb, bit=True)
    bounds = []

    for i in range(1,len(list(txt_cmap['bound']))):
        lin=np.linspace(list(txt_cmap['bound'])[i-1],list(txt_cmap['bound'])[i], 40)
        bounds.extend(list(lin))

    norm = BoundaryNorm(bounds,ncolors=cmap.N)
    #or lons[l] != -75.913585

    fig = plt.figure(figsize=(380/my_dpi, 772/my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(111, projection=ccrs.Mercator())
    ax.set_extent([-76.35, -75, 38.44, 40.26], crs=ccrs.PlateCarree())
    for ind in range(0,len(bigdeos)):
            ax.add_geometries([bigdeos['geometry'][ind]], oldproj,
                          facecolor='silver', edgecolor='black')
    im=ax.pcolormesh(cl['x'].values,cl['y'].values,cl.values[0],cmap=cmap,norm=norm,transform=ccrs.PlateCarree(),zorder=2)
    for l in range(0,len(lons)):
        if lons[l] != -76.35 and lons[l] != -74.68 and lons[l] != -75.062685 and lons[l] != -75.118033 and lons[l] != -75.247235 and lons[l] != -75.640685 and lons[l] != -75.527755 and lons[l] != -75.118033 and lons[l] != -75.148629 and lons[l] != -75.727202:
            text = plt.text(lons[l],lats[l],str('{:.1f}'.format(round(feelsLike[l], rounder))), size=6.5,weight='bold',verticalalignment='center',
            horizontalalignment='center',transform=ccrs.PlateCarree(),zorder=5)
            text.set_path_effects([path_effects.Stroke(linewidth=2.5, foreground='white'),path_effects.Normal()])
    for ind in range(0,len(deos_boundarys)):
        ax.add_geometries([deos_boundarys['geometry'][ind]], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='gray', zorder=3, linewidth=0.5)
    for ind in range(0,len(inland_bays)):
        ax.add_geometries([inland_bays['geometry'][ind]], oldproj,
                          facecolor='white', edgecolor='black',zorder=3, linewidth=0.5)
    ax.add_geometries([state_outline['geometry'][74]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=0.5)
    ax.add_geometries([bigdeos['geometry'][121]], oldproj, facecolor='none', edgecolor='black',zorder=3, linewidth=0.5)
    #plt.title(nameDict[var])
    plt.text(-76.13, 38.503, 'Feels Like (F)',horizontalalignment='left',weight='bold',color='white',size=6,zorder=30,transform=ccrs.PlateCarree())
    plt.text(-76.13, 38.473, deos_dateSTR,horizontalalignment='left',weight='bold',color='white',size=6,zorder=30,transform=ccrs.PlateCarree())

    im1 = image.imread(shapePaths + "deos_logo.png")
    plt.figimage(im1, 24, 40 ,zorder=30, alpha=1)
    plt.savefig("/home/map_maker/deos_maps/imagery/deos_feelslike.png",bbox_inches='tight',pad_inches = 0,dpi=my_dpi*1.3)
else:
    print("no new data yet..waiting...")
