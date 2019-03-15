# goes16 sst quality control algorithm
import xarray as xr
import metpy
from datetime import datetime
from os import listdir
from os.path import isfile, join
import os
import sys
import numpy as np 
# test - seems that the lower right hand corner is a solid region as far as data goes.
# the algorithm processes from North Pole to South Pole and many times the failed datasets
# are due to lack of time for algorithm to complete. Thus, we check the lower right hand corner.
# we will check this region for gaps in the data. Instead of deleting the passes, let's just 
# add a suspect folder and place them in there so they aren't just deleted


nowdate = datetime.utcnow()
datadir = "/data/GOES/GOES-R/sst/" + str(2018) + "/"
filenames = [f for f in listdir(datadir) if isfile(join(datadir, f))]

for f in range(0,len(filenames)):
    ds = xr.open_dataset(datadir + filenames[f])
    sst = ds.metpy.parse_cf("SST")
    sst = sst.sel(longitude=slice(-88, -50), latitude=slice(16,18))
    sst = sst.where(sst.values > 270, np.nan)
    sst = sst.fillna(-999)
    if np.percentile(sst.values, 60) == -999:
        print(sst.time.values)
        move_files = "mv " + datadir + filenames[f] + " /data/GOES/GOES-R/suspect/"
        os.system(move_files)
    if np.percentile(sst.values, 60) != -999:
        dqf = ds.metpy.parse_cf("DQF")
        dqf = dqf.sel(longitude=slice(-65, -50), latitude=slice(37,45))
        dqf = dqf.where(dqf.values > 1, np.nan)
        dqf = dqf.fillna(-999)
        if np.percentile(dqf.values, 60) == -999:
            move_files = "mv " + datadir + filenames[f] + " /data/GOES/GOES-R/suspect/"
            os.system(move_files)


'''
badfiles2 = []
ds = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
for t in range(0,len(ds.time.values)):
    sst = ds.metpy.parse_cf("SST")[t]
    #sst = sst.sel(time=datetime.datetime.strptime('2019-03-08 16:28:18', '%Y-%m-%d %H:%M:%S'))
    sst = sst.sel(longitude=slice(-88, -50), latitude=slice(16,18))
    sst = sst.where(sst.values > 270, np.nan)
    sst = sst.fillna(-999)
    if np.percentile(sst.values, 60) == -999:
        badfiles.append(sst.time.values)



ds = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/goes_r_sst.nc")
for t in range(0,len(ds.time.values)):
    dqf = ds.metpy.parse_cf("DQF")[t]
    dqf = dqf.sel(longitude=slice(-65, -50), latitude=slice(37,45))
    dqf = dqf.where(dqf.values > 1, np.nan)
    dqf = dqf.fillna(-999)
    if np.percentile(dqf.values, 60) == -999:
        badfiles2.append(dqf.time.values)

#from itertools import groupby
#[len(list(v)) for k,v in groupby(sst.values) if k==-999]
'''

