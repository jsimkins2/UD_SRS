import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
import numpy as np

from metpy.cbook import get_test_data
from metpy.interpolate import (interpolate_to_grid, remove_nan_observations,
                               remove_repeat_coordinates)

from metpy.plots import ctables

levels = list(range(-30, 120, 1))
cmap = ctables.registry.get_colortable('rainbow')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
from scipy.interpolate import griddata

# target grid to interpolate to
temp = list()
# Begin the for loop to place a location on a map
for s in station_id:
    if str(s) in list(station_dict.values()):
        ID = str(s)
        try:
            temp.append(float(deos_data[date_deos][s]['Air Temperature']))
        except:
            temp.append(np.nan)


lats=list()
lons=list()
for key in station_id:
    # grag latitude and longitudes location of station
    lats.append(loc_deos[rev_station_dict[str(key)]]['latitude'])
    lons.append(loc_deos[rev_station_dict[str(key)]]['longitude'])

lons=np.array(lons)
lats=np.array(lats)
temp = np.array(temp)

temp = temp - 273.15
temp = (temp*(9/5)) + 32

lons,lats, temp = remove_nan_observations(lons,lats, temp)

x = np.linspace(min(lons), max(lons), 5000)
y = np.linspace(min(lats), max(lats), 5000)
xi,yi = np.meshgrid(x,y)
# interpolate
zi = griddata((lons,lats),temp,(xi,yi),method='cubic')

# mask out the field

# plot
fig = plt.figure()
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
plt.pcolormesh(xi,yi,zi, cmap=cmap, norm=norm,transform=ccrs.PlateCarree())
plt.plot(lons,lats,'k.')
plt.xlabel('xi',fontsize=16)
plt.ylabel('yi',fontsize=16)

