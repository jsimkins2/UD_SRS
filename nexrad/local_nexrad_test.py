import matplotlib
import matplotlib.pyplot as plt
from siphon.catalog import TDSCatalog, get_latest_access_url
import urllib
from netCDF4 import Dataset, num2date
from matplotlib import ticker
import matplotlib as mpl
from dateutil import tz
import time
from time import mktime
import os.path
from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.image as image
from datetime import datetime, timedelta
from matplotlib.patches import Rectangle
import pyart
from siphon.radarserver import RadarServer, get_radarserver_datasets
import numpy.ma as ma
import netCDF4
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import shapely.geometry as sgeom
import xarray
import cartopy.crs as ccrs
site = 'KDOX'
nowtime = datetime.utcnow().replace(second=0, microsecond=0)
#https://thredds-jumbo.unidata.ucar.edu/thredds/catalog/nexrad/level2/catalog.html
cat = TDSCatalog('https://thredds-test.unidata.ucar.edu/thredds/catalog/nexrad/level2/' + site + '/' + str(nowtime.year) + str(nowtime.month).zfill(2) + str(nowtime.day).zfill(2) + '/catalog.xml')
#https://thredds-test.unidata.ucar.edu/thredds/dodsC/nexrad/level2/KDOX/20220830/Level2_KDOX_20220830_2042.ar2v
dataset = list(cat.datasets.values())[-1]


os.system('wget ' + dataset.access_urls['HTTPServer'] + ' /home/sat_ops/goesR/radar/data/radar.ar2v')

radar = pyart.io.read_nexrad_archive('/home/sat_ops/goesR/radar/data/radar.ar2v')

pyart.io.read_nexrad_archive('https://thredds-test.unidata.ucar.edu/thredds/fileServer/nexrad/level2/KDOX/20220830/Level2_KDOX_20220830_0005.ar2v')

my_gf = pyart.filters.GateFilter(radar)
my_gf.exclude_above('differential_reflectivity', 6, exclude_masked=False)
my_gf.exclude_below('reflectivity', 10, exclude_masked=False)
my_gf.exclude_below('cross_correlation_ratio', 1, exclude_masked=False)
my_ds_gf = pyart.correct.despeckle_field(radar, 'reflectivity', gatefilter=my_gf)
# Here we see reflectivity values below zero masked.
fig = plt.figure(figsize=[8, 8])
display = pyart.graph.RadarMapDisplay(radar)
display.plot_ppi_map('reflectivity', sweep=0, resolution='50m',
                     vmin=10,min_lon=-79, max_lon=-76,
                     min_lat=38, max_lat=40,vmax=60,
                     projection=ccrs.PlateCarree())#, gatefilter=my_ds_gf)
plt.show()




radar.fields['reflectivity']







#### MMIDATLANTIC TESTING ####


nowdate = datetime.utcnow()
# https://thredds.ucar.edu/thredds/catalog/grib/NCEP/MRMS/BaseRef/latest.html
cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/MRMS/CONUS/BaseRef/catalog.html')

nexrad_name = cat.datasets['Full Collection Dataset']
nexrad = nexrad_name.remote_access(use_xarray=True)
proj_var = nexrad.variables['LatLon_Projection']
created_plot = False

fileind = [-1]
for i in fileind:
    refltime=i
    tdim = nexrad['MergedBaseReflectivityQC_altitude_above_msl'].dims[0]
    if tdim == 'time1':
        refl = nexrad['MergedBaseReflectivityQC_altitude_above_msl'].isel(time1=refltime, altitude_above_msl=0)
        timestamp = pd.Timestamp(nexrad.refvalidtime1.values[i]).to_pydatetime()
    if tdim == 'time2':
        refl = nexrad['MergedBaseReflectivityQC_altitude_above_msl'].isel(time2=refltime, altitude_above_msl=0)
        timestamp = pd.Timestamp(nexrad.refvalidtime2.values[i]).to_pydatetime()
    if tdim == 'time3':
        refl = nexrad['MergedBaseReflectivityQC_altitude_above_msl'].isel(time3=refltime, altitude_above_msl=0)
        timestamp = pd.Timestamp(nexrad.refvalidtime3.values[i]).to_pydatetime()
    if tdim == 'time4':
        refl = nexrad['MergedBaseReflectivityQC_altitude_above_msl'].isel(time4=refltime, altitude_above_msl=0)
        timestamp = pd.Timestamp(nexrad.refvalidtime4.values[i]).to_pydatetime()








import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates

import cartopy.crs as ccrs
import pyart
import pandas as pd

import nexradaws
import tempfile
import pytz

templocation = tempfile.mkdtemp()

import cartopy.feature as cfeature
from metpy.plots import USCOUNTIES



### Define the radar, start time and end time
radar_id = 'KDOX'
start = pd.Timestamp(2022,8,30,18,0).tz_localize(tz='UTC')
end = pd.Timestamp(2022,8,30,20,15).tz_localize(tz='UTC')

### Bounds of map we want to plot
min_lon = -93.25
max_lon = -88.
min_lat = 40.35
max_lat = 43.35



### Define the radar, start time and end time
radar_id = 'KDVN'
start = pd.Timestamp(2020,8,10,16,30).tz_localize(tz='UTC')
end = pd.Timestamp(2020,8,10,21,0).tz_localize(tz='UTC')

### Bounds of map we want to plot
min_lon = -93.25
max_lon = -88.
min_lat = 40.35
max_lat = 43.35
## download these files
#results = conn.download(scans[0:2], templocation)
#### and get the data
conn = nexradaws.NexradAwsInterface()
scans = conn.get_avail_scans_in_range(start, end, radar_id)
print("There are {} scans available between {} and {}\n".format(len(scans), start, end))
print(scans[0:4])

## download these files
#results = conn.download(scans[0:2], templocation)
results = conn.download(scans, templocation)





radar = pyart.io.read_nexrad_archive('/Users/james/Downloads/Level2_KDOX_20220830_2042.ar2v')




