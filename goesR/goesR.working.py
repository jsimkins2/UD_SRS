%matplotlib inline
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
from pyproj import Proj                  # Note to self: pyproj not installed on meteo19.
                             



C_file = "http://thredds.demac.udel.edu/thredds/dodsC/GOESR_FD.nc?band01[1][1:962][1:1011],band02[1][1:962][1:1011],band03[1][1:962][1:1011]"
C = Dataset(C_file, 'r')

R = np.sqrt(C.variables['band02'][0, 0:961, 0:1010]) # Band 2 is red (0.64 um)
G = np.sqrt(C.variables['band03'][0, 0:961, 0:1010]) # Band 3 is "green" (0.865 um)
B = np.sqrt(C.variables['band01'][0, 0:961, 0:1010]) # Band 1 is blue (0.47 um)

# "True Green" is some linear interpolation between the three channels
G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G

# The final RGB array :)
RGB = np.dstack([R, G_true, B])

plt.figure(figsize=[10, 8])
plt.imshow(RGB)
plt.title(DATE)