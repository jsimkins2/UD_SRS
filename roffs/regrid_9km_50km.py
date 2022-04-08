#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 15:45:02 2022

@author: james
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe


ds_50k = xr.open_dataset("/Users/james/Downloads/50km_test.nc")
ds_50k.sst.plot()
ds_9k = xr.open_dataset("/Users/james/Downloads/avhrr_8day.1997.concat.nc4")
ds_9k.sst.plot()

def add_corners(xrds, lonvar="lon", latvar="lat"):
    lon_centers = xrds[lonvar].values
    lat_centers = xrds[latvar].values
    # should rename lat/lon variables here
    #slice down xrds because corner length must be greater than center length
    xrds = xrds.isel(lon=slice(1,len(lon_centers)-1), lat=slice(1,len(lat_centers)-1))
    
    #expanded lon corners
    lon_corners = 0.5 * (lon_centers[1:] + lon_centers[:-1])
    lat_corners = 0.5 * (lat_centers[1:]+ lat_centers[:-1])
    
    xrds['lon_b'] = lon_corners
    xrds['lat_b'] = lat_corners
    xrds
    
    return xrds

ds_9k = add_corners(xrds=ds_9k)
ds_50k = add_corners(xrds=ds_50k)
# both grids are regular, so we can just keep lat/lon coordinate names
# create the regridder: this takes in 9k and regrids to 50k
regridder = xe.Regridder(ds_9k, ds_50k, "conservative")

# regrid 9k to 50k ()
dr_out = regridder(ds_9k.sst)
dr_out.to_netcdf("bilinear_fill_50km.nc")
dr_out[0].shape == ds_50k.sst.shape
dr_out.name ="sst"

dr_out.to_netcdf("/Users/james/Downloads/conservative_fill_50km.nc")


ds_50k = xr.open_dataset("/Users/james/Downloads/50km_test.nc")
ds_50k.sst.plot()
ds_9k = xr.open_dataset("/Users/james/Downloads/avhrr_8day.1997.concat.nc4")
ds_9k.sst.plot()
ds_9k.sst[0].plot(vmin=0,vmax=35)

regridder = xe.Regridder(ds_9k, ds_50k, "bilinear")

# regrid 9k to 50k ()
dr_out = regridder(ds_9k.sst)

dr_out[0].plot(vmin=0,vmax=35)
dr_out.name ="sst"

dr_out.to_netcdf("/Users/james/Downloads/bilinear_fill_50km.nc")

