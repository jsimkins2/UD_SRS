import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from metpy.plots import USCOUNTIES
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import pyproj
from pyproj import Proj
from matplotlib.colors import Normalize
from cartopy.feature import ShapelyFeature
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
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

    if proj.is_latlong():
        return ccrs.PlateCarree()

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

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

min_lon = -76.2 #lons['data'].min() + 2.5
min_lat = 38.3 #lats['data'].min() + 2
max_lat = 39.9 #lats['data'].max() - 2
max_lon = -74.6

data = pd.read_csv("/usr/local/src/cfms/cfms/water_levels.csv")
data.index = data.station

# shapefile
shp = shpreader.Reader('/home/james/cfms/proj_info/cfms_watersheds.shp')
shp2 = shpreader.Reader('/home/james/cfms/inland_bays/cfms_inland_bays.shp')
c = Proj('+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +datum=NAD83 +units=m +no_defs ')
oldproj = proj_to_cartopy(c)
subplot_kw = dict(projection=ccrs.Mercator())

cmap_cfms = ['#F7FBFF', '#E9F2FA', '#DDEAF6', '#D1E2F2', '#C3DAEE', '#AFD1E7', '#9BC7E0','#7FB9DA', '#66ABD4', '#509CCC', '#3D8DC3', '#2B7CBB', '#1D6BB0','#0C5AA2', '#07468B', '#08306B']
colrng = [2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.5,5,5.5,6,6.01]

fig, ax = plt.subplots(figsize=(7, 11),subplot_kw=subplot_kw)
ax.set_extent((min_lon, max_lon, min_lat, max_lat))
ax.background_patch.set_visible(False)
ax.outline_patch.set_visible(False)
for record, state in zip(shp.records(), shp.geometries()):
    tem = data.maxpred[str(record.attributes['station']) + '       ']
    col_ = cmap_cfms[find_nearest(colrng, tem)]
    ax.add_geometries([state], oldproj,
                      facecolor=col_, edgecolor='black', linewidth=0.5)

plt.savefig('/var/www/html/applications/cfms/home_kmz/cfms_plot.png', transparent=True, bbox_inches='tight', dpi=150)


fig, ax = plt.subplots(figsize=(7, 11),subplot_kw=subplot_kw)
ax.set_extent((min_lon, max_lon, min_lat, max_lat))
ax.background_patch.set_visible(False)
ax.outline_patch.set_visible(False)
for record, state in zip(shp2.records(), shp2.geometries()):
    tem = data.maxpred[str(record.attributes['station']) + '       ']
    col_ = cmap_cfms[find_nearest(colrng, tem)]
    ax.add_geometries([state], oldproj,
                      facecolor=col_, edgecolor='black', linewidth=0.5)

plt.savefig('/var/www/html/applications/cfms/home_kmz/cfms_inland_bays.png', transparent=True, bbox_inches='tight', dpi=150)