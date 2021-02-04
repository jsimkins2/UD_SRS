# experiment with blended viirs (both suomi npp and noaa20) chlorophyll products
# from JPL
# ERDDAP:https://coastwatch.pfeg.noaa.gov/erddap/griddap/nesdisVHNnoaaSNPPnoaa20chlaDaily.graph?chl_oci[(2020-08-04T12:00:00Z)][(0.0)][(89.75625):(-89.75626)][(-179.9812):(179.9813)]&.draw=surface&.vars=longitude%7Clatitude%7Cchl_oci&.colorBar=%7C%7C%7C%7C%7C&.bgColor=0xffccccff
# THREDDS: 
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import pandas as pd

area2 = [36.5, 40, -75.25, -72.0]
blended_main = xr.open_dataset("https://www.star.nesdis.noaa.gov/thredds/dodsC/chlociVIIRSnpp-n20GlobalDailyWW00/V2020218_D1_NPP-N20_WW00_chloci.nc")
viirs_blend = blended_main.sel(lat=slice(area2[1],area2[0]), lon=slice(area2[2], area2[3]))


# paths
outpath = "/home/james/roffs/"
#landmask.to_netcdf(path=outpath + 'landmask_roffs_area' + str(a) + '.nc', format='NETCDF3_CLASSIC')
goes_nc.to_netcdf(path=outpath + 'roffs_' +  'area' + str(a) + '_' + str(daysback[d]) + "day"  + '.nc', format='NETCDF3_CLASSIC')

