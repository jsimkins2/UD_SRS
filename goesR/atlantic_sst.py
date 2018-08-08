
# Velocity Script
# By James Simkins

# The velocity is very shakey, so I need to try and fix tha that issue
import matplotlib
import matplotlib.pyplot as plt
from siphon.catalog import TDSCatalog, get_latest_access_url
import numpy as np
from datetime import datetime, timedelta
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# query data
nowdate = datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/SeaSurfaceTemperature/FullDisk/' + \
                  str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')
dataset_name = sorted(cat.datasets.keys())[-1]
dataset = cat.datasets[dataset_name]

# load netcdf and read variables
nc = dataset.remote_access()
sst = np.array(nc.variables['SST'][:,:])
sst[np.isnan(sst)] = -1
sst = np.ma.array(sst)
sst[sst < 0] = np.ma.masked

X = nc.variables['x'][:]
Y = nc.variables['y'][:]

# define projections
proj_var = nc.variables['goes_imager_projection']
globe = ccrs.Globe(ellipse='sphere', semimajor_axis=proj_var.semi_major_axis,
                   semiminor_axis=proj_var.semi_minor_axis)

# define reprojection target
proj = ccrs.Mercator(central_longitude=proj_var.longitude_of_projection_origin, 
                     latitude_true_scale=proj_var.latitude_of_projection_origin,
                     globe=globe)

# Plot
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.coastlines(resolution='50m', color='black')
ax.add_feature(cfeature.STATES, linestyle=':', edgecolor='black')
ax.add_feature(cfeature.BORDERS, linewidth=2, edgecolor='black')
im = ax.imshow(sst, extent=(X.min(), X.max(), Y.min(), Y.max()), origin='upper')


# try again, this time with pyproj
from pyproj import Proj     
p = Proj(proj='geos', h=proj_var.perspective_point_height, lon_0=proj_var.longitude_of_projection_origin, sweep=proj_var.sweep_angle_axis)

X = nc.variables['x'][:] * proj_var.perspective_point_height
Y = nc.variables['y'][:] * proj_var.perspective_point_height
XX, YY = np.meshgrid(X,Y)
lons, lats = p(XX, YY, inverse=True)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.coastlines(resolution='50m', color='black')
ax.add_feature(cfeature.STATES, linestyle=':', edgecolor='black')
ax.add_feature(cfeature.BORDERS, linewidth=2, edgecolor='black')
im = ax.imshow(sst, extent=(lons.min(), lons.max(), lats.min(), lats.max()), origin='upper')



