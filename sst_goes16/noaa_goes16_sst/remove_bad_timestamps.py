import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import warnings
import datetime as datetime
x = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/GOESNOAASST.nc")




bad_values = sorted(x.time.values)[0:2]

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

idx_bad = list()
for b in bad_values:
    idx_bad.append(find_nearest(array = x.time.values, value = b))


x.time.values[12196]


x.time.values[2638]
dt64 = np.datetime64(x.time.values[2638])
ts = (dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
datetime.datetime.utcfromtimestamp(ts).timetuple().tm_yday
