# this is the primary script responsible for downloading and extracting JPL GOES-16 SST data 
import os
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
from datetime import datetime, timedelta
from netCDF4 import Dataset

# define helper function for removing problem attributes - I don't know if we use this at all but I'm afraid to take it out... 4/15/2019
def remove_problematic_attrs(ds):
    for variable in ds.variables.values():
        if 'coordinates' in variable.attrs:
            del variable.attrs['coordinates']
            
nowdate = datetime.utcnow()
nowdate = nowdate.replace(second=0, microsecond=0)
temdate = nowdate - timedelta(days=3)
mainURL = 'https://thredds.jpl.nasa.gov/thredds/dodsC/OceanTemperature/ABI_G16-STAR-L3C-v2.70.nc'

# Begin loop
# JPL Thredds Server call, use try here in case there isn't any data for the whole day

jpl = xr.open_dataset(mainURL,decode_coords=False)
jpl = jpl.sel(time=slice(temdate,nowdate))
jpl = jpl.sel(lat=slice(52,16), lon=slice(-100,-50))
datetimes_jpl = []
    
for t in jpl.time.values:
    t_ = t.timetuple()
    datetimes_jpl.append(str("{0:0=4d}".format(t_.tm_year) + "_" + "{0:0=2d}".format(t_.tm_mon) + "{0:0=2d}".format(t_.tm_mday) + "_" + "{0:0=2d}".format(t_.tm_hour) + "{0:0=2d}".format(t_.tm_min)))

for sstDate in datetimes_jpl:
    sstTime = datetime.strptime(sstDate, "%Y_%m%d_%H%M")
    print("processing" + str(sstTime))
    if os.path.isfile("/data/GOES/GOES-R/jpl_sst/" + str(sstTime.year) + "/jpl_goesSST_" + str(sstDate) + ".nc") == False:
        jpl_tem = jpl.sel(time=sstTime, method='nearest')
        
        f = Dataset("/data/GOES/GOES-R/jpl_sst/" + str(sstTime.year) + "/jpl_goesSST_" + str(sstDate) + ".nc", 
                    'w', format='NETCDF4')  # 'w' stands for write
        # dimensions
        f.createDimension('longitude', jpl_tem['lon'].size)
        f.createDimension('latitude', jpl_tem['lat'].size)
        f.createDimension('time', 1)
        
        # variables
        longitude = f.createVariable('longitude', 'f4', 'longitude')
        longitude.standard_name = jpl_tem['lon'].attrs['long_name']
        longitude.units = jpl_tem['lon'].attrs['units']
        
        latitude = f.createVariable('latitude', 'f4', 'latitude')
        latitude.standard_name = jpl_tem['lat'].attrs['long_name']
        latitude.units = jpl_tem['lat'].attrs['units']
        
        time = f.createVariable('time', 'f8', 'time')
        time.standard_name = "time"
        time.long_name = "EPOCH Time"
        time.units = "seconds since 1970-01-01T00:00:00Z"
        
        proj = f.createVariable('projection', 'f4')
        proj.proj4_string = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        proj.epsg_code = 4326
        
        sst = f.createVariable(
            'sea_surface_temperature', 'f4', ('time', 'latitude', 'longitude'))
        sst.long_name = jpl_tem.variables['sea_surface_temperature'].attrs['long_name']
        sst.standard_name = jpl_tem.variables['sea_surface_temperature'].attrs['standard_name']
        sst.units = jpl_tem.variables['sea_surface_temperature'].attrs['units']
        sst.comment = jpl_tem.variables['sea_surface_temperature'].attrs['comment']
        
        dqf = f.createVariable(
            'quality_level', 'f4', ('time', 'latitude', 'longitude'))
        dqf.long_name = jpl_tem.variables['quality_level'].attrs['long_name']
        dqf.flag_values = jpl_tem.variables['quality_level'].attrs['flag_values']
        dqf.flag_meanings = jpl_tem.variables['quality_level'].attrs['flag_meanings']
        
        dt = f.createVariable(
            'dt_analysis', 'f4', ('time', 'latitude', 'longitude'))
        dt.long_name = jpl_tem.variables['dt_analysis'].attrs['long_name']
        dt.units = jpl_tem.variables['dt_analysis'].attrs['units']
        dt.comment = jpl_tem.variables['dt_analysis'].attrs['comment']
        
        zen = f.createVariable(
            'satellite_zenith_angle', 'f4', ('time', 'latitude', 'longitude'))
        zen.long_name = jpl_tem.variables['satellite_zenith_angle'].attrs['long_name']
        zen.units = jpl_tem.variables['satellite_zenith_angle'].attrs['units']
        zen.comment = jpl_tem.variables['satellite_zenith_angle'].attrs['comment']
        
        sses = f.createVariable(
            'sses_bias', 'f4', ('time', 'latitude', 'longitude'))
        sses.long_name = jpl_tem.variables['sses_bias'].attrs['long_name']
        sses.units = jpl_tem.variables['sses_bias'].attrs['units']
        sses.comment = jpl_tem.variables['sses_bias'].attrs['comment']
        
        sses_dev = f.createVariable(
            'sses_standard_deviation', 'f4', ('time', 'latitude', 'longitude'))
        sses_dev.long_name = jpl_tem.variables['sses_standard_deviation'].attrs['long_name']
        sses_dev.units = jpl_tem.variables['sses_standard_deviation'].attrs['units']
        sses_dev.comment = jpl_tem.variables['sses_standard_deviation'].attrs['comment']
        
        wind = f.createVariable(
            'wind_speed', 'f4', ('time', 'latitude', 'longitude'))
        wind.long_name = jpl_tem.variables['wind_speed'].attrs['long_name']
        wind.units = jpl_tem.variables['wind_speed'].attrs['units']
        wind.comment = jpl_tem.variables['wind_speed'].attrs['comment']
        
    
        # data
        latitude[:] = jpl_tem['lat'].values
        longitude[:] = jpl_tem['lon'].values
        sst[:] = jpl_tem['sea_surface_temperature'].values
        dqf[:] = jpl_tem['quality_level'].values
        dt[:] = jpl_tem['dt_analysis'].values
        zen[:] = jpl_tem['satellite_zenith_angle'].values
        sses[:] = jpl_tem['sses_bias'].values
        sses_dev[:] = jpl_tem['sses_standard_deviation'].values
        wind[:] = jpl_tem['wind_speed'].values
        time[:] = sstTime.timestamp()
        
        # metadata
        f.creator_name = "James Simkins"
        f.creator_email = "simkins@udel.edu"
        f.institution = "University of Delaware Ocean Exploration, Remote Sensing and Biogeography Group (ORB)"
        f.url = "http://orb.ceoe.udel.edu/"
        f.source = "NOAA/NESDIS/STAR GHRSST- http://www.star.nesdis.noaa.gov - SST "
        f.groundstation = "University of Delaware, Newark, Center for Remote Sensing"
        f.summary = "NOAA/NESDIS/STAR GHRSST GOES16 SST product, fit to UDel ORB lab NW Atlantic extent, reprojected to EPSG:4326."
        f.acknowledgement = "These data were provided by Group for High Resolution Sea Surface Temperature (GHRSST) and the National Oceanic and Atmospheric Administration (NOAA)"
        f.close()


# place the most recent days data (all 24 hours) in a file without resampling to daily res.
today = pd.date_range(pd.datetime.today(), pd.datetime.today()).tolist()
goes_nc = xr.open_dataset(
    "http://basin.ceoe.udel.edu/thredds/dodsC/GOESJPLSST.nc")
goes_nc = goes_nc.sel(time=datetime.strftime(today[0].date(), '%Y-%m-%d'))
goes_nc = goes_nc.drop('projection')
goes_nc = goes_nc.drop('dt_analysis')
goes_nc = goes_nc.drop('satellite_zenith_angle')
goes_nc = goes_nc.drop('sses_standard_deviation')
goes_nc = goes_nc.drop('wind_speed')
goes_nc['sea_surface_temperature'] = goes_nc['sea_surface_temperature'] - goes_nc['sses_bias']
goes_nc = goes_nc.drop('sses_bias')
newtimestamp = goes_nc.time.values[-1]
for t in range(len(goes_nc.time.values)):
    x = goes_nc['sea_surface_temperature'][t]
    goes_nc['sea_surface_temperature'][t] = x.where(goes_nc['quality_level'][t] == 5)
    #x = goes_nc['DQF'][t]
    #goes_nc['DQF'][t] = x.where(goes_nc['DQF'][t] == 3)


goes_nc['sst'] = goes_nc['sea_surface_temperature']
goes_nc = goes_nc.drop(['sea_surface_temperature'])

outpath = "/data/GOES/GOES-R/1day/"
goes_nc.to_netcdf(path=outpath + '/' + str(today[0].year) + '/GOES16_SST_1day_' + str(today[0].year) + str("{0:0=3d}".format(
    today[0].dayofyear)) + '_' + str("{0:0=2d}".format(today[0].month)) + str("{0:0=2d}".format(today[0].day)) + '.nc', format='NETCDF4')

# resample to a daily composite
goes_nc = goes_nc.drop(['quality_level'])
goes_nc = goes_nc.resample(time='1D').mean('time')

newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
goes_nc.time.values = np.array([newtimestamp], dtype='float64')
goes_nc.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'

outpath = "/data/GOES/GOES-R/daily_composite/"
goes_nc.to_netcdf(path=outpath + '/' + str(today[0].year) + '/GOES16_SST_dailycomposite_' + str(today[0].year) + str("{0:0=3d}".format(
    today[0].dayofyear)) + '_' + str("{0:0=2d}".format(today[0].month)) + str("{0:0=2d}".format(today[0].day)) + '.nc', format='NETCDF4')


# now make a rolling 1 day aka last 24 hours IN CELSIUS
goes_nc = xr.open_dataset(
    "http://basin.ceoe.udel.edu/thredds/dodsC/GOESJPLSST.nc")
# grab the last 24 hours of sst dataset
goes_nc = goes_nc.isel(time=range(-24, 0))
goes_nc = goes_nc.drop('projection')
goes_nc = goes_nc.drop('dt_analysis')
goes_nc = goes_nc.drop('satellite_zenith_angle')
goes_nc = goes_nc.drop('sses_standard_deviation')
goes_nc = goes_nc.drop('wind_speed')
goes_nc['sea_surface_temperature'] = goes_nc['sea_surface_temperature'] - goes_nc['sses_bias']
goes_nc = goes_nc.drop('sses_bias')
newtimestamp = goes_nc.time.values[-1]
for t in range(len(goes_nc.time.values)):
    x = goes_nc['sea_surface_temperature'][t]
    goes_nc['sea_surface_temperature'][t] = x.where(goes_nc['quality_level'][t] == 5)

goes_nc = goes_nc.drop(['quality_level'])
goes_nc['sea_surface_temperature'] = goes_nc['sea_surface_temperature'] - 273.15
dat = goes_nc
dat['sst'] = dat['sea_surface_temperature']
dat = dat.drop(['sea_surface_temperature'])
dat = dat.where(dat.sst > -1)
dat.attrs['units'] = 'Celsius'
dat = dat.mean('time')
newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
x = dat.assign_coords(time=newtimestamp)
dat = x.expand_dims('time')
dat.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
outpath = "/data/GOES/GOES-R/rolling_1day/"
dat.to_netcdf(path=outpath + 'GOES16_SST_rolling_1day.nc',
                  format='NETCDF4')


# now make a rolling 1 day aka last 24 hours IN FAHRENHEIT
goes_nc = xr.open_dataset(
    "http://basin.ceoe.udel.edu/thredds/dodsC/GOESJPLSST.nc")
# grab the last 24 hours of sst dataset
goes_nc = goes_nc.isel(time=range(-24, 0))
goes_nc = goes_nc.drop('projection')
goes_nc = goes_nc.drop('dt_analysis')
goes_nc = goes_nc.drop('satellite_zenith_angle')
goes_nc = goes_nc.drop('sses_standard_deviation')
goes_nc = goes_nc.drop('wind_speed')
goes_nc['sea_surface_temperature'] = goes_nc['sea_surface_temperature'] - goes_nc['sses_bias']
goes_nc = goes_nc.drop('sses_bias')
newtimestamp = goes_nc.time.values[-1]
for t in range(len(goes_nc.time.values)):
    x = goes_nc['sea_surface_temperature'][t]
    goes_nc['sea_surface_temperature'][t] = x.where(goes_nc['quality_level'][t] == 5)

goes_nc = goes_nc.drop(['quality_level'])
goes_nc['sea_surface_temperature'] = (goes_nc['sea_surface_temperature'] - 273.15)*(9/5) + 32
dat = goes_nc
dat['sst'] = dat['sea_surface_temperature']
dat = dat.drop(['sea_surface_temperature'])
dat = dat.where(dat.sst > 31)
dat.attrs['units'] = 'Fahrenheit'
dat = dat.mean('time')
newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
x = dat.assign_coords(time=newtimestamp)
dat = x.expand_dims('time')
dat.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
outpath = "/data/GOES/GOES-R/rolling_1day_fahrenheit/"
dat.to_netcdf(path=outpath + 'GOES16_SST_rolling_1day_fahrenheit.nc',
                  format='NETCDF4')
