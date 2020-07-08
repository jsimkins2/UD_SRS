import xarray as xr
import pyart
from siphon.catalog import TDSCatalog, get_latest_access_url
from datetime import datetime, timedelta

x = xr.open_dataset("http://thredds-jetstream.unidata.ucar.edu/thredds/dodsC/nexrad/level3/N0H/DOX/20200625/Level3_DOX_N0H_20200625_0007.nids")
radar = pyart.io.read_nexrad_cdm(ds.access_urls['OPENDAP'])

radar = pyart.io.read_nexrad_cdm(dataset.access_urls['OPENDAP'])

date = datetime.utcnow() - timedelta(days=1)
print(date)

site = 'DOX'
nowtime = datetime.utcnow().replace(second=0, microsecond=0)
cat = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/nexrad/level3/'
                 f'N0H/DOX/{date:%Y%m%d}/catalog.xml')
dataset = list(cat.datasets.values())[0]

radar = pyart.io.read_nexrad_level3(dataset.access_urls['OPENDAP'])

radar = pyart.io.read_nexrad_level3("Downloads/Level3_DOX_N0H_20200625_0007.nids")

add_metpy_logo(fig, 190, 85, size='large')

import matplotlib.pyplot as plt
import numpy as np

from metpy.cbook import get_test_data
from metpy.io import Level3File
from metpy.plots import add_metpy_logo, add_timestamp, colortables
    # Open the file
name = 'Downloads/Level3_DOX_N0H_20200625_0007.nids'
f = Level3File(name)

# Pull the data out of the file object
datadict = f.sym_block[0][0]

# Turn into an array, then mask
data = np.ma.array(datadict['data'])
data[data == 0] = np.ma.masked

# Grab azimuths and calculate a range based on number of gates
az = np.array(datadict['start_az'] + [datadict['end_az'][-1]])
rng = np.linspace(0, f.max_range, data.shape[-1] + 1)

# Convert az,range to x,y
xlocs = rng * np.sin(np.deg2rad(az[:, np.newaxis]))
ylocs = rng * np.cos(np.deg2rad(az[:, np.newaxis]))

# Plot the data
norm, cmap = colortables.get_with_steps(ctable, 16, 16)
ax.pcolormesh(xlocs, ylocs, data, norm=norm, cmap=cmap)
ax.set_aspect('equal', 'datalim')
ax.set_xlim(-40, 20)
ax.set_ylim(-30, 30)
add_timestamp(ax, f.metadata['prod_time'], y=0.02, high_contrast=True)

plt.show()