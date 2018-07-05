import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from siphon.catalog import TDSCatalog
import urllib
from netCDF4 import Dataset, num2date
from matplotlib import ticker
import matplotlib as mpl
import pyart
from dateutil import tz
import time
from time import mktime
import os.path
import numpy as np
import matplotlib.image as image
import datetime
############# Initial Set Up ##################
workdir = "/home/sat_ops/goesR/radar/"
imgdir = "imgconus/"

datadir = "/home/sat_ops/goesR/data/"
# Go to the Unidata Thredds Server for the Current Day
nowdate = datetime.utcnow()
cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/TotalPrecipitableWater/CONUS/current/catalog.xml')
dataset_name = sorted(cat.datasets.keys())[-1]
dataset = cat.datasets[dataset_name]
pwnc = dataset.remote_access()
