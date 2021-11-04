#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 21:29:07 2021

@author: james
"""
import xarray as xr
import time

# original at 9km
orig = xr.open_dataset("/Users/james/Downloads/S2008349.nc4")

# coarsen to 45km
start_time = time.time()
coarse = orig.coarsen(lat=5, lon=5, boundary='pad').mean()
elapsed = time.time() - start_time

# elapsed ~ 5 seconds

x = xr.open_dataset('/Users/james/Downloads/seawifs_resample.nc')