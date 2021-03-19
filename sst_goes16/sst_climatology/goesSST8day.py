import pandas as pd
import os
import os.path
import time
import xarray as xr
from datetime import datetime as datetime
from datetime import timedelta
import numpy as np


start_date = datetime.utcnow() - timedelta(days=10)
end_date = datetime.utcnow()
daterange = pd.date_range(start_date, end_date)

workdir = '/data/aquaGoesSST/C'
# now make a rolling 1 day aka last 24 hours
goes_main = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/GOESJPLSST.nc")

for d in daterange:
    # grab the last 24 hours of sst dataset
    goes_nc = goes_main.sel(time=slice(d.to_pydatetime() - timedelta(days=8), d.to_pydatetime()))
    newtimestamp = goes_nc.time.values[-1]
    newtimestamp = (newtimestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    
    if os.path.isfile("/data/aquaGoesSST/C/" + str(time.strftime('%Y', time.localtime(newtimestamp))) + "/SSTanomalyGoesAqua8dayCelsius" + str(time.strftime('%Y%m%d', time.localtime(newtimestamp))) + ".nc") == False:
        goes_nc = goes_nc.mean('time')
        print(str(time.strftime('%Y%m%d', time.localtime(newtimestamp))))
        x = goes_nc.assign_coords(time=newtimestamp)
        goes_nc = x.expand_dims('time')
        goes_nc.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00'
        outpath = "/home/sat_ops/goesR/jpl_sst/sstClimatology/"
        goes_nc.to_netcdf(path=outpath + 'GOES16_SST_8day.nc')
        os.system("Rscript /home/sat_ops/goesR/github/UD_SRS/sst_goes16/sst_climatology/goesVsAquaClimatologyCelsius.R")
        os.system("Rscript /home/sat_ops/goesR/github/UD_SRS/sst_goes16/sst_climatology/goesVsAquaClimatologyFahrenheit.R")
        os.system("rm /home/sat_ops/goesR/jpl_sst/sstClimatology/GOES16_SST_8day.nc")

