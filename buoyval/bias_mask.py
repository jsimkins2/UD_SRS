# now use the bias data to create a bias mask that Sarah can use to interpolate
from scipy.interpolate import griddata
import rioxarray
import xarray as xr
import pyproj
from pyproj import Proj
import numpy as np

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
    
stations = {'44065' : 'NY Harbor',
            '44009' : 'Delaware Bay',
            '44091' : 'Barnegat',
            '44025' : 'URC', 
            '44089' : 'LLC', 
            '44402' : 'LRC',
            '44069' : 'ULC'
    }

NYharbias = 0.2090 #40.369 -73.703 91, 183
DBayBias =  0.3461 #38.457 -74.702 196, 128
BarnBias = -.0198 #39.778 -73.769 123, 180
LLC = .1502
URC = .2884

lons = np.array([-73.703, -74.702, -73.769, -76.98615, -72.00014, -72.00014, -76.98615])
lats = np.array([40.369, 38.457, 39.778, 37.01693, 41.99443, 37.01693, 41.99443])
temp = np.array([NYharbias, DBayBias, BarnBias, LLC, URC, LLC, URC])
bounds = [-76.98615, -72.00014, 37.01693, 41.99443]

x = np.linspace(bounds[1], bounds[0], 278)
y = np.linspace(bounds[2], bounds[3], 276)
y = list(y)
y.append((y[-1] + (y[-1] - y[-2])))
y.append((y[-1] + (y[-1] - y[-2])*2))
y = np.array(y)
xi,yi = np.meshgrid(x,y)
# interpolate
#zi = griddata((lons,lats),temp,(xi,yi),method='nearest')
# try the idw interpolation scheme
xi, yi = xi.flatten(), yi.flatten()
#zi = griddata((lons,lats),temp,(xi,yi),method='cubic')
# Calculate IDW
zi = linear_rbf(lons,lats,temp,xi,yi)
zi=zi.reshape((len(x), len(y)))
#zi = np.transpose(zi)
zi = np.delete(zi, 277, 0)
zi = np.delete(zi, 276, 0)
da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y[0:-2]})
da.rio.set_crs("epsg:4326")
da.attrs['units'] = 'Kelvin'
da.attrs['standard_name'] = 'Temperature'
da.rio.set_spatial_dims('lon', 'lat')
da.rio.to_raster('Downloads/idw_sarah.tif', overwrite=True)
da.to_netcdf('Downloads/idw_sarah.nc')

x = np.linspace(bounds[1], bounds[0], 278)
y = np.linspace(bounds[2], bounds[3], 276)
y = list(y)
y.append((y[-1] + (y[-1] - y[-2])))
y.append((y[-1] + (y[-1] - y[-2])*2))
y = np.array(y)
xi,yi = np.meshgrid(x,y)
zi = griddata((lons,lats),temp,(xi,yi),method='nearest')
# try the idw interpolation scheme
#zi = np.transpose(zi)
zi = np.delete(zi, 277, 0)
zi = np.delete(zi, 276, 0)
da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y[0:-2]})
da.rio.set_crs("epsg:4326")
da.attrs['units'] = 'Kelvin'
da.attrs['standard_name'] = 'Temperature'
da.rio.set_spatial_dims('lon', 'lat')
da.rio.to_raster('Downloads/nearest_sarah.tif', overwrite=True)
da.to_netcdf('Downloads/nearest_sarah.nc')


x = np.linspace(bounds[1], bounds[0], 278)
y = np.linspace(bounds[2], bounds[3], 276)
y = list(y)
y.append((y[-1] + (y[-1] - y[-2])))
y.append((y[-1] + (y[-1] - y[-2])*2))
y = np.array(y)
xi,yi = np.meshgrid(x,y)
zi = griddata((lons,lats),temp,(xi,yi),method='linear')
# try the idw interpolation scheme
#zi = np.transpose(zi)
zi = np.delete(zi, 277, 0)
zi = np.delete(zi, 276, 0)
da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y[0:-2]})
da.rio.set_crs("epsg:4326")
da.attrs['units'] = 'Kelvin'
da.attrs['standard_name'] = 'Temperature'
da.rio.set_spatial_dims('lon', 'lat')
da.rio.to_raster('Downloads/linear_sarah.tif', overwrite=True)
da.to_netcdf('Downloads/linear_sarah.nc')


x = np.linspace(bounds[1], bounds[0], 278)
y = np.linspace(bounds[2], bounds[3], 276)
y = list(y)
y.append((y[-1] + (y[-1] - y[-2])))
y.append((y[-1] + (y[-1] - y[-2])*2))
y = np.array(y)
xi,yi = np.meshgrid(x,y)
zi = griddata((lons,lats),temp,(xi,yi),method='cubic')
# try the idw interpolation scheme
#zi = np.transpose(zi)
zi = np.delete(zi, 277, 0)
zi = np.delete(zi, 276, 0)
da = xr.DataArray(zi,dims=['lat', 'lon'],coords={'lon': x, 'lat' :y[0:-2]})
da.rio.set_crs("epsg:4326")
da.attrs['units'] = 'Kelvin'
da.attrs['standard_name'] = 'Temperature'
da.rio.set_spatial_dims('lon', 'lat')
da.rio.to_raster('Downloads/cubic_sarah.tif', overwrite=True)
da.to_netcdf('Downloads/cubic_sarah.nc')





