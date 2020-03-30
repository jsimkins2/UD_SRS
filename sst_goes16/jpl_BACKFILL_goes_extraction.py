
# Download JPL version of GOES-16 SST data
# this runs on /home/sat_ops/goesR/sstClimatology/mkanomaly.sh
import os
import numpy as np
import xarray as xr
import netCDF4
from datetime import datetime, timedelta

def remove_problematic_attrs(ds):
    for variable in ds.variables.values():
        if 'coordinates' in variable.attrs:
            del variable.attrs['coordinates']
            
nowdate = datetime.utcnow()
nowdate = nowdate.replace(second=0, microsecond=0)
temdate = datetime.strptime("2018-01-01", "%Y-%m-%d")
# Use both servers to double check
mainURL = 'https://thredds.jpl.nasa.gov/thredds/dodsC/OceanTemperature/ABI_G16-STAR-L3C-v2.70.nc'

# Begin loop
# JPL Thredds Server call, use try here in case there isn't any data for the whole day
try:
    jpl = xr.open_dataset(mainURL,decode_coords=False)
    jpl = jpl.sel(time=slice(temdate,nowdate))
    jpl = jpl.sel(lat=slice(52,16), lon=slice(-100,-50))
    datetimes_jpl = []
    
    for t in jpl.time.values:
        t_ = t.timetuple()
        datetimes_jpl.append(str("{0:0=4d}".format(t_.tm_year) + "_" + "{0:0=2d}".format(t_.tm_mon) + "{0:0=2d}".format(t_.tm_mday) + "_" + "{0:0=2d}".format(t_.tm_hour) + "{0:0=2d}".format(t_.tm_min)))
    
    for sstDate in datetimes_jpl:
        sstTime = datetime.strptime(sstDate, "%Y_%m%d_%H%M")
        if os.path.isfile("/data/GOES/GOES-R/jpl_sst/" + str(sstTime.year) + "/jpl_goesSST_" + str(sstDate) + ".nc") == False:
            jpl_tem = jpl.sel(time=sstTime, method='nearest')
            remove_problematic_attrs(jpl_tem)
            #add in a time dimensinone
            outfile = "/data/GOES/GOES-R/jpl_sst/" + str(sstTime.year) + "/jpl_goesSST_" + str(sstDate) + ".nc"
            jpl_tem.to_netcdf(path=outfile, format='NETCDF3_CLASSIC')
            #add in a time dimensinon
            ncap2 -s \"defdim("time",1);time[time]=74875.0;time@long_name="Time";time@units="seconds since 1970-01-01"' -O jpl_goesSST_2020_0330_1200.nc out.nc
            addTime = "ncap2 -s " + "'" + "defdim(\"time\",1);time[time]=" + str(sstTime.timestamp()) + ";time@long_name=\"Time\";time@units=\"seconds since 1970-01-01\"" + "'" + " -O " + outfile + " " + outfile
            os.system(addTime)
            print("downloaded and saved " + str(sstDate) + ".nc")
        else:
            print("already downloaded this file")
except:
    pass

        
