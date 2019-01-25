import xarray as xr
import metpy
from datetime import datetime
from xarray.backends import NetCDF4DataStore
from netCDF4 import Dataset, num2date
from siphon.catalog import TDSCatalog
from pyproj import Proj
import cartopy.crs as ccrs
import os.path
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import os.path
import os
import sys

datadir = '/home/sat_ops/goesR/radar/prectype/hrrr_temp/'
nowdate = datetime.utcnow()
cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/latest.xml')
dataset_name = sorted(cat.datasets.keys())[-1]
dataset = cat.datasets[dataset_name]

def trim_data(lats, lons, ref, boundinglat, boundinglon):
    for i in range(lats.shape[0]):
        for j in range(lats.shape[1]):
            if (lats[i][j]>boundinglat[0] and lats[i][j]<boundinglat[1] and \
                lons[i][j]>boundinglon[0] and lons[i][j]<boundinglon[1]):
                pass
            else:
                # lats[i] = np.nan; lons[i] = np.nan
                ref[i][j] = np.nan
    return ref

if os.path.isfile(datadir + dataset_name[:-5] + 'nc') == False:
    os.system("rm " + datadir + "*")
    print('downloading ' + dataset_name)
    ds = dataset.remote_access(service='OPENDAP')
    ds = NetCDF4DataStore(ds)
    ds = xr.open_dataset(ds)
    tempiso = ds.metpy.parse_cf('Temperature_isobaric')
    hproj = tempiso.metpy.cartopy_crs
    hproj = Proj(hproj.proj4_init)
    
    
    # Write the netcdf file
    
    f = Dataset(datadir + dataset_name[:-5] + 'nc','w', format='NETCDF4') #'w' stands for write
    # dimensions
    f.createDimension('x', tempiso['x'].size)
    f.createDimension('y', tempiso['y'].size)
    f.createDimension('time', 1)
    
    # variables
    
    x = f.createVariable('x', 'f4', 'x')
    x.standard_name = tempiso['x'].standard_name
    x.units = tempiso['x'].units
    
    y = f.createVariable('y', 'f4', 'y')
    y.standard_name = tempiso['y'].standard_name
    y.units = tempiso['y'].units
    
    time = f.createVariable('time', 'f8', 'time')
    time.standard_name = tempiso['time'][1].standard_name
    time.long_name = tempiso['time'][1].long_name
    
    proj = f.createVariable('projection', 'f4')
    proj.crs = 'lambert_conformal_conic'
    proj.srs = hproj.srs
    
    tsurf = f.createVariable('temperature_surface', 'f4', ('time', 'y', 'x'))
    tsurf.long_name = 'Temperature @ Surface'
    tsurf.standard_name = 'temperature_surface'
    tsurf.units = 'Kelvin'
    
    t925 = f.createVariable('temperature_925', 'f4', ('time', 'y', 'x'))
    t925.long_name = 'Temperature @ 925mb'
    t925.standard_name = 'temperature_925'
    t925.units = 'Kelvin'
    
    t850 = f.createVariable('temperature_850', 'f4', ('time', 'y', 'x'))
    t850.long_name = 'Temperature @ 850mb'
    t850.standard_name = 'temperature_850'
    t850.units = 'Kelvin'
    
    ht500 = f.createVariable('height_500', 'f4', ('time', 'y', 'x'))
    ht500.long_name = 'Geopotential_height_isobaric @ 500mb'
    ht500.standard_name = 'height_500'
    ht500.units = 'gpm'
    
    ht1000 = f.createVariable('height_1000', 'f4', ('time', 'y', 'x'))
    ht1000.long_name = 'Geopotential_height_isobaric @ 1000mb'
    ht1000.standard_name = 'height_1000'
    ht1000.units = 'gpm'
    
    # fill the data
    x[:] = tempiso['x'].values
    y[:] = tempiso['y'].values
    time[:] = tempiso['time'][1].values
    tsurf[:]= ds.metpy.parse_cf('Temperature_height_above_ground')[1][0].values
    t925[:]=tempiso[1][3].values
    t850[:]=tempiso[1][2].values
    ht500[:]=ds.metpy.parse_cf("Geopotential_height_isobaric")[1][0].values
    ht1000[:]=ds.metpy.parse_cf("Geopotential_height_isobaric")[1][3].values
    f.close()
    
    # now run the reprojecting script 
    os.system("Rscript /home/sat_ops/goesR/github/UD_SRS/hrrr/reproject_hrrr.R " + dataset_name[:-5] + 'nc')
    
    datadir = "/home/sat_ops/goesR/radar/prectype/hrrr_temp/"
    filenames = [f for f in listdir(datadir) if isfile(join(datadir, f))]
    hrrrdata = datadir + ''.join(([f for f in filenames if f[0:3]=='rep']))

    ds = xr.open_dataset(hrrrdata)
    # parse the temperature at various heights
    tsurf = ds.metpy.parse_cf('temperature_surface')
    t850 = ds.metpy.parse_cf('temperature_850')
    t925 = ds.metpy.parse_cf('temperature_925')
    ht1000 = ds.metpy.parse_cf("height_1000")
    ht500 = ds.metpy.parse_cf("height_500")
    thick = ht500 - ht1000
    
    hproj = ccrs.Geodetic()
    hproj = Proj(hproj.proj4_init)
    
    # trim the data to save space
    lons, lats = np.meshgrid(tsurf['longitude'], tsurf['latitude'])
    hrrr_t850 = trim_data(lats, lons, ma.getdata(t850), boundinglat, boundinglon)
    hrrr_t925 = trim_data(lats, lons, ma.getdata(t925), boundinglat, boundinglon)
    hrrr_tsurf = trim_data(lats, lons, ma.getdata(tsurf), boundinglat, boundinglon)
    thick = trim_data(lats, lons, ma.getdata(thick), boundinglat, boundinglon)
    
    # we have to ravel these for scipy interpolate
    rav_lats = lats.ravel()
    rav_lons = lons.ravel()
    rav_t850 = hrrr_t850.ravel()
    rav_t925 = hrrr_t925.ravel()
    rav_tsurf = hrrr_tsurf.ravel()
    rav_thick = thick.ravel()
    
    #Grid Data using scipy interpolate
    grid_lons = np.linspace(boundinglon[0],boundinglon[1],nlon)
    grid_lats = np.linspace(boundinglat[0],boundinglat[1],nlat)
    glon,glat = np.meshgrid(grid_lons,grid_lats)
    grid850= griddata((rav_lons,rav_lats),rav_t850,(glon,glat),method='linear')
    grid925 = griddata((rav_lons,rav_lats),rav_t925,(glon,glat),method='linear')
    gridsurf = griddata((rav_lons,rav_lats),rav_tsurf,(glon,glat),method='linear')
    gridthick = griddata((rav_lons,rav_lats),rav_thick,(glon,glat),method='linear')
    
    # create a masked array for each precipitation type
    rain = (gridsurf > 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
    rain = np.ma.masked_array(gref, ~rain)
    ice = (grid850 > 273.15) & (grid925 > 273.15) & (gridsurf < 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
    ice = np.ma.masked_array(gref, ~ice)
    sleet = (grid850 > 273.15) & (grid925 < 273.15) & (gridsurf < 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
    sleet = np.ma.masked_array(gref, ~sleet)
    snow = (grid850 < 273.15) & (grid925 < 273.15) & (gridsurf < 273.15) &  (np.isfinite(gref)) & (np.isfinite(grid850)) & (np.isfinite(grid925)) & (np.isfinite(gridsurf))
    snow = np.ma.masked_array(gref, ~snow) 

    np.save(datadir + 'rain.npy', rain)
    np.save(datadir + 'ice.npy', ice)
    np.save(datadir + 'sleet.npy', sleet)
    np.save(datadir + 'snow.npy', snow)
    